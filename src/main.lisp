;;; Copyright 2025 Steve Losh and contributors
;;; SPDX-License-Identifier: GPL-3.0-or-later

(in-package :minimera)

;;;; Configuration ------------------------------------------------------------
(defparameter *k* 8)
(defparameter *w* 16)

(defparameter *foldback-position-epsilon* 50)
(defparameter *intercept-epsilon* 30)
(defparameter *gap-epsilon* 300)
(defparameter *minimum-cluster-length* 40)
(defparameter *minimum-foldback-length-absolute* 50)
(defparameter *minimum-foldback-length-relative* 0.05)

(defparameter *monotony-threshold* 0.50)
(defparameter *minimum-qscore* 9.0)
(defparameter *dorado-mean-qscore* nil)

(defparameter *plot-foldbacks* nil)
(defparameter *plot-normal* nil)
(defparameter *low-quality-threshold* 3)
(defparameter *low-quality-window-size* 10)

(defparameter *worker-threads* 6)
(defparameter *output-directory* "results")
(defparameter *annotation-path* nil)
(defparameter *annotation-tracks* (list))
(defparameter *report-progress* t)
(defparameter *interactive* t)
(defparameter *input-format* :fastq)


;;;; Types --------------------------------------------------------------------
(deftype u8 ()
  '(unsigned-byte 8))

(deftype a8 ()
  '(simple-array u8 (*)))

(defun str->a8 (str)
  (map 'a8 #'char-code str))


;;;; Minimizers ---------------------------------------------------------------
(declaim (inline make-minmimizer))

(defstruct (minimizer (:constructor make-minimizer (hash pos))
                      (:conc-name nil))
  (hash 0 :type u64)
  (pos  0 :type u64))

(defmethod print-object ((o minimizer) s)
  (print-unreadable-object (o s :type t)
    (format s "~16,'0X ~D" (hash o) (pos o))))


(defun nreverse-complement-bytes (seq)
  "Mutate and return the reverse complement of `seq`."
  (declare (optimize (speed 3) (safety 1) (debug 1)))
  (check-type seq a8)
  (map-into seq (lambda (base)
                  (ecase base
                    (#.(char-code #\A) (char-code #\T))
                    (#.(char-code #\C) (char-code #\G))
                    (#.(char-code #\G) (char-code #\C))
                    (#.(char-code #\T) (char-code #\A))))
            seq)
  (nreverse seq))

(defun reverse-complement-bytes (seq)
  "Return a fresh a8 of the reverse complement of `seq`."
  (nreverse-complement-bytes (copy-seq seq)))


(defun-inline base-bits (base)
  (ecase base
    (#.(char-code #\A) #b00)
    (#.(char-code #\C) #b01)
    (#.(char-code #\G) #b10)
    (#.(char-code #\T) #b11)))


(defun compute-initial-chunk (w sequence)
  (declare (optimize (speed 3) (safety 1) (debug 1)))
  (check-type sequence a8)
  (iterate (with-result result = 0)
           (for base :in-vector sequence :below w)
           (zapf result (wrap64 (logior (ash % 2) (base-bits base))))))

(defun minimizer (k w chunk window-start)
  (declare (optimize (speed 3) (safety 1) (debug 1)))
  (check-type k (integer 1 30))
  (check-type w (integer 1 31))
  (check-type window-start (and fixnum (integer 0)))
  (iterate
    (with k2 = (* 2 k))
    (with w2 = (* 2 w))
    (declare (type fixnum k2 w2))
    (for window/pos :from window-start)
    (declare (type fixnum window/pos))
    (for chunk/pos :from (- w2 k2) :downto 0 :by 2)
    (for kmer = (ldb (byte k2 chunk/pos) chunk))
    (for hash = (phihash-chunk kmer))
    (finding (make-minimizer hash window/pos) :minimizing hash)))

(defun minimizers (k w sequence)
  "Return a list of the (`k`, `w`, φ)-mimimizer sketches of `sequence`.

  `sequence` must be an array of FASTQ bytes.

  "
  (declare (optimize (speed 3) (safety 1) (debug 1)))
  (check-type k (integer 1 30))
  (check-type w (integer 1 31))
  (check-type sequence a8)
  (iterate
    (with k2 = (* 2 k))
    (with w2 = (* 2 w))
    (with w-k = (- w k))
    (declare (type fixnum k2 w2 w-k))
    (with prev = nil)
    (with candidate = (make-minimizer 0 0))
    (for start :from 0 :to (- (length sequence) w))
    (for end :from w)
    (declare (type fixnum start end))
    (when (and prev (< (pos prev) start))
      ;; If we've moved the window off the previous best minimizer, discard it.
      (setf prev nil))
    (for chunk
         :first (compute-initial-chunk w sequence)
         :then (_ chunk
                 (ash _ 2) ; move
                 (ldb (byte w2 0) _) ; mask off old start base
                 (logior _ (base-bits (aref sequence (1- end)))) ; add new base
                 wrap64)) ; promise sbcl it's small enough
    (if prev
      ;; If we still have a previous minimizer, we just compare the new one and
      ;; check if it's better than that one.
      (progn
        (setf (hash candidate) (phihash-chunk (ldb (byte k2 0) chunk))
              (pos candidate) (+ start w-k))
        (when (< (hash candidate) (hash prev))
          (let ((next (copy-minimizer candidate)))
            (collect next)
            (setf prev next))))
      ;; Otherwise we don't have a previous winning minimizer, either because
      ;; we're at the beginning of the sequence or it fell off the moving
      ;; window.  Just recompute the window from scratch here.
      (let ((next (minimizer k w chunk start)))
        (collect next)
        (setf prev next)))))


;;;; Monotony -----------------------------------------------------------------
(defun estimate-monotony (minimizers)
  (_ minimizers
    (proportions _ :key #'hash)
    ;; TODO make a version of this supporting :result-type 'vector
    (alexandria:hash-table-values _)
    (coerce _ 'vector)
    (sort _ #'>)
    ;; Sum the proportions of the top 0.1% of minimizers.
    (take (max (truncate (length minimizers) 1000) 2) _)
    summation))


;;;; Minimizer Hits/Matches ---------------------------------------------------
(defstruct (hit (:constructor make-hit (c l1 l2)))
  (c nil :type fixnum)
  (l1 nil :type fixnum)
  (l2 nil :type fixnum))

(defmethod print-object ((o hit) s)
  (print-unreadable-object (o s :type t)
    (format s "~D ~D ~D" (hit-c o) (hit-l1 o) (hit-l2 o))))

(defun hit< (a b)
  (declare (optimize (speed 3) (safety 1) (debug 1)))
  (check-type a hit)
  (check-type b hit)
  ;; ordering hierarchy is (< c) (> l2) (< l1)
  (let ((x (hit-c a))
        (y (hit-c b)))
    (cond ((< x y) t)
          ((< y x) nil)
          (t (let ((x (hit-l2 a))
                   (y (hit-l2 b)))
               (cond ((> x y) t)
                     ((> y x) nil)
                     (t (let ((x (hit-l1 a))
                              (y (hit-l1 b)))
                          (cond ((< x y) t)
                                ((< y x) nil)
                                (t nil))))))))))

(defun hits (ms1 ms2)
  "Return a vector of matching minimizer hits from sketches `ms1` and `ms2`.

  The vector will be sorted in parameter-space-intercept-order.

  "
  (let ((index-1 (group-by #'hash ms1 :test #'equal :map #'pos))
        (index-2 (group-by #'hash ms2 :test #'equal :map #'pos)))
    (_ (iterate
         (for (minimizer locs-1) :in-hashtable index-1)
         (for locs-2 = (gethash minimizer index-2))
         (appending (alexandria:map-product (lambda (l1 l2)
                                              (make-hit (- l2 l1) l1 l2))
                                            locs-1 locs-2)
                    :into result)
         (returning (sort (coerce result 'vector) #'hit<))))))


;;;; Clustering ---------------------------------------------------------------
(defun break-gaps (cluster)
  "Return a list of new clusters from `cluster`, breaking gaps larger than `epsilon`.

  The list may contain just a single element (the original `cluster`) if no gaps
  are present.

  "
  (iterate
    (with epsilon = *gap-epsilon*)
    (with last-i = (1- (length cluster)))
    (with start = 0)
    (for hit :in-vector cluster :with-index i)
    (for next-i = (1+ i))
    (for next-hit = (if (= i last-i)
                      nil
                      (aref cluster next-i)))
    (when (or (null next-hit)
              (> (- (hit-l1 next-hit) (hit-l1 hit)) epsilon)
              (> (- (hit-l2 next-hit) (hit-l2 hit)) epsilon))
      (collect (subseq cluster start next-i))
      (setf start next-i))))


(defun longest-increasing-subsequence (cluster)
  "Return the longest increasing subsequence of hits in `cluster`."
  (declare (optimize (speed 3) (safety 1) (debug 1))
           (type (simple-array hit (*)) cluster))
  ;; Used https://github.com/privet-kitty/cl-competitive/blob/master/module/lis.lisp
  ;; as a starting point, which is public domain and implements the LIS algorithm
  ;; described at https://en.wikipedia.org/wiki/Longest_increasing_subsequence#Efficient_algorithms
  ;; for speed, because profiling showed that a naive implementation using
  ;; Patience Sort from Rosetta code was the performance bottleneck.
  (let* ((n (length cluster))
         (dp (make-array (+ 1 n)
               :element-type 'fixnum
               :initial-element most-positive-fixnum))
         (end 0))
    (declare ((integer 0 #.array-dimension-limit) end))
    (labels ((bisect (ng ok value)
               (declare (type (integer -1 #.array-dimension-limit) ng ok)
                        (type fixnum value))
               (if (<= (- ok ng) 1)
                 ok
                 (let ((mid (ash (+ ng ok) -1)))
                   (if (< (aref dp mid) value)
                     (bisect mid ok value)
                     (bisect ng mid value))))))
      (let ((prevs (make-array n :element-type 'fixnum :initial-element -1))
            (indices (make-array n :element-type '(integer 0 #.most-positive-fixnum))))
        (loop :for i :below n
              :for val = (hit-l2 (aref cluster i))
              :for pos = (bisect -1 end val)
              :do (progn (setf (aref dp pos) val
                               (aref indices pos) i)
                         (unless (zerop pos)
                           (setf (aref prevs i)
                                 (aref indices (- pos 1))))
                         (when (= end pos)
                           (incf end))))
        (let ((result (make-array end)))
          (unless (zerop end)
            (loop :with index = (aref indices (- end 1))
                  :for i :from (- end 1) :downto 0
                  :do (setf (aref result i) (aref cluster index)
                            index (aref prevs index))))
          result)))))

(defun filter-clusters (clusters min-length)
  "Return a list of clusters, with gaps broken and small clusters removed."
  (let* ((deduped (mapcar (lambda (cluster)
                            (longest-increasing-subsequence
                              (sort cluster #'< :key #'hit-l1)))
                          clusters))
         (degapped (alexandria:mappend #'break-gaps deduped))
         (results (remove-if (lambda (cluster)
                               (< (length cluster) min-length))
                             degapped)))
    results))

(defun cluster (hits &key (min-length *minimum-cluster-length*))
  "Cluster the sorted `hits`, returning a list of vectors of `hit`s."
  (unless (zerop (length hits))
    (iterate
      (with last-i = (1- (length hits)))
      (with start = 0)
      (with epsilon = *intercept-epsilon*)
      (for hit :in-vector hits :with-index i)
      (for next-i = (1+ i))
      (for next-hit = (if (= i last-i)
                        nil
                        (aref hits next-i)))
      ;; (for stop = (or (null next-hit)
      ;;                 (> (hit-c next-hit) 500)))
      (when (or (null next-hit)
                ;; stop
                (> (- (hit-c next-hit) (hit-c hit)) epsilon))
        (collect (subseq hits start next-i) :into results)
        (setf start next-i))
      ;; (when stop
      ;;   (finish))
      (returning (filter-clusters results min-length)))))


;;;; Foldback Detection -------------------------------------------------------
(defun cluster-bounds (cluster)
  "Return (as multiple values) the bounds of the cluster's region on both reads."
  (assert (plusp (length cluster)))
  (let ((lo (aref cluster 0))
        (hi (aref cluster (1- (length cluster)))))
    (values (hit-l1 lo) (hit-l1 hi)
            (hit-l2 lo) (hit-l2 hi))))

(defun cluster-lengths (cluster)
  "Return (as multiple values) the length of the cluster's region on both reads."
  (assert (plusp (length cluster)))
  (let ((lo (aref cluster 0))
        (hi (aref cluster (1- (length cluster)))))
    (values (- (hit-l1 hi) (hit-l1 lo))
            (- (hit-l2 hi) (hit-l2 lo)))))

(defun region-in-foldback-position-p (read-length cluster)
  "Return whether `cluster`'s region is in the expected position for a foldback region."
  (multiple-value-bind (r1s r1e r2s r2e) (cluster-bounds cluster)
    (declare (ignore r1s r2e))
    (and (<= (abs (- r2s 0)) *foldback-position-epsilon*)
         (<= (abs (- r1e read-length)) *foldback-position-epsilon*))))

(defun region-long-enough-p (read-length cluster)
  "Return whether `cluster`'s region is long enough for a foldback region."
  (let ((min-length (max *minimum-foldback-length-absolute*
                         (truncate (* read-length *minimum-foldback-length-relative*)))))
    (multiple-value-bind (l1 l2) (cluster-lengths cluster)
      (and (>= l1 min-length)
           (>= l2 min-length)))))

(defun find-foldback-clusters (read-length clusters)
  "Return candidates for a foldback region from the given clusters.

  May return an empty list if no clusters look like possible foldback regions.

  "
  (remove-if-not (alexandria:conjoin
                   (curry #'region-in-foldback-position-p read-length)
                   (curry #'region-long-enough-p read-length))
                 clusters))

(defun compute-foldback-position (cluster)
  (multiple-value-bind (start end) (cluster-bounds cluster)
    (truncate (+ start end) 2)))


;;;; Annotations --------------------------------------------------------------
(defun make-smith-waterman-grid (rows cols)
  (make-array (list rows cols)
    :element-type 'fixnum
    :initial-element 0))

(defun find-best-score (grid)
  (iterate outer
    (with (rows cols) = (array-dimensions grid))
    (for r :from 0 :below rows)
    (iterate (for c :from 0 :below cols)
             (in outer (finding (cons r c) :maximizing (aref grid r c))))))

(defun extract-smith-waterman-path (grid prev)
  (iterate
    (for coord
         :first (find-best-score grid)
         :then (aref prev (car coord) (cdr coord)))
    (until (zerop (aref grid (car coord) (cdr coord))))
    (collect (cons (1- (car coord))
                   (1- (cdr coord)))
             :at start)))

(defun smith-waterman (seq1 seq2)
  (check-type seq1 a8)
  (check-type seq2 a8)
  (iterate
    (with rows = (1+ (length seq1)))
    (with cols = (1+ (length seq2)))
    (with grid = (make-smith-waterman-grid rows cols))
    (with prev = (make-array (list rows cols) :initial-element nil))
    (for row :from 1 :below rows)
    (for ch1 :in-array seq1)
    (iterate
      (for col :from 1 :below cols)
      (for ch2 :in-array seq2)
      (for row-1 = (1- row))
      (for col-1 = (1- col))
      (for match = (+ (aref grid row-1 col-1) (if (= ch1 ch2) 1 -1)))
      (for gap-1 = (+ (aref grid row col-1) -2))
      (for gap-2 = (+ (aref grid row-1 col) -2))
      (for best = (max match gap-1 gap-2))
      (setf (aref grid row col) best
            (aref prev row col) (cond ((= best match) (cons row-1 col-1))
                                      ((= best gap-1) (cons row col-1))
                                      ((= best gap-2) (cons row-1 col)))))
    (returning (extract-smith-waterman-path grid prev))))

(defun nucleotide-char-p (char)
  (member char '(#\A #\C #\T #\G)))

(defclass* (annotation-track :conc-name track-) ()
  ((name)
   (sequence)
   (id-sequence-bytes)
   (rc-sequence-bytes)
   (id-minimizers)
   (rc-minimizers)))

(defmethod print-object ((o annotation-track) s)
  (print-unreadable-object (o s :type t)
    (format s "~A" (track-name o))))

(defun read-fasta-entry (stream &optional (eof-error-p t) eof-value)
  (with-eof-handled (stream eof-error-p eof-value)
    (let* ((header-line (read-line stream))
           (name (progn (assert (char= #\> (char header-line 0)) ()
                          "Bad annotation header line ~S, must start with >." header-line)
                        (subseq header-line 1)))
           (chunks (iterate (when (char= #\> (peek-char nil stream nil #\>))
                              (finish))
                            (for line :in-stream stream :using #'read-line)
                            (collect line)))
           (sequence (apply #'concatenate 'string chunks))
           (id-sequence-bytes (progn
                                (assert (every #'nucleotide-char-p sequence) ()
                                  "Bad annotation sequence ~A, invalid sequence content (must be ACTG only)." name)
                                (str->a8 sequence)))
           (rc-sequence-bytes (reverse-complement-bytes id-sequence-bytes))
           (id-minimizers (minimizers *k* *w* id-sequence-bytes))
           (rc-minimizers (minimizers *k* *w* rc-sequence-bytes)))
      (make-instance 'annotation-track
        :name name
        :sequence sequence
        :id-sequence-bytes id-sequence-bytes
        :rc-sequence-bytes rc-sequence-bytes
        :id-minimizers id-minimizers
        :rc-minimizers rc-minimizers))))

(defun parse-annotation-file (path)
  (with-open-file (stream path)
    (iterate
      (for track :in-stream stream :using #'read-fasta-entry)
      (collect track))))

(defun match-short-annotation-track (track seq)
  (flet ((bounds (alignment)
           (cons (reduce #'min alignment :key #'car)
                 (reduce #'max alignment :key #'car))))
    (let* ((id-align (smith-waterman seq (track-id-sequence-bytes track)))
           (rc-align (smith-waterman seq (track-rc-sequence-bytes track))))
      (list (when (> (length id-align) 10)
              (list (bounds id-align)))
            (when (> (length rc-align) 10)
              (list (bounds rc-align)))))))

(defun match-long-annotation-track (track ms)
  (flet ((bounds (cluster)
           (cons (reduce #'min cluster :key #'hit-l1)
                 (reduce #'max cluster :key #'hit-l1))))
    (let* ((id-hits (hits ms (track-id-minimizers track)))
           (rc-hits (hits ms (track-rc-minimizers track)))
           (id-clusters (cluster id-hits :min-length 12))
           (rc-clusters (cluster rc-hits :min-length 12)))
      (list (map 'list #'bounds id-clusters)
            (map 'list #'bounds rc-clusters)))))


(defun match-annotation-track (track ms seq)
  (if (> (length (track-sequence track)) 50)
    (match-long-annotation-track track ms)
    (match-short-annotation-track track seq)))


;;;; Hallucination ------------------------------------------------------------
(eval-when (:compile-toplevel :load-toplevel :execute)
  (defun q->p% (q)
    ;; P = 10 ^ (-Q/10)
    (expt 10.0d0 (/ (- q) 10.0d0)))

  (defun make-q-table% ()
    (let ((result (make-array (list 61) :element-type 'double-float)))
      (loop :for q :from 0 :to 60
            :do (setf (aref result q) (q->p% q)))
      result)))

(defun-inline q->p (q)
  "Return the error probability (as a double-float) corresponding to Q score `q`.

   `q` must be an `(integer 0 60)`, and a lookup table is used to make this
   computation fast.  If you need to convert arbitrary Q scores (e.g. a mean
   Q score) use `p->q%`.

  "
  (aref #.(make-q-table%) q))

(defun-inline p->q (p)
  "Return the Q score (as a double-float) corresponding to error probability `p`."
  ;; Q = -10 * log10(P)
  (* -10.0d0 (log p 10.0d0)))

(defun compute-mean-qscore (qs &key start end)
  (if (and (null start) (null end))
    (let ((truncate-read (and *dorado-mean-qscore* (> (length qs) 60))))
      (compute-mean-qscore qs :start (if truncate-read 60 0)))
    (let* ((start (or start 0))
           (end (or end (length qs)))
           (length (- end start)))
      (if (zerop length)
        0.0d0
        ;; Compute this in probability space.  Q scores are logarithmic so averaging
        ;; them doesn't make much sense.
        (coerce (p->q (/ (reduce #'+ qs :start start :end end :key #'q->p) length))
                'double-float)))))

(defun compute-foldback-qscores (qs cluster)
  (multiple-value-bind (foldback-start foldback-end) (cluster-bounds cluster)
    (let* ((foldback-mid (compute-foldback-position cluster))
           (prefix-qscore (unless (zerop foldback-start)
                            (compute-mean-qscore qs :start 0 :end foldback-start)))
           (foldback-prefix-qscore (compute-mean-qscore qs :start foldback-start :end foldback-mid))
           (foldback-suffix-qscore (compute-mean-qscore qs :start foldback-mid   :end foldback-end)))
      (list prefix-qscore foldback-prefix-qscore foldback-suffix-qscore))))

(defun compute-llqr (quality-scores &key
                     (window-size *low-quality-window-size*)
                     (low-quality-threshold *low-quality-threshold*))
  "Return the longest low-quality region in `quality-scores`."
  (declare (optimize (speed 3) (space 1) (safety 1) (debug 1)))
  (check-type quality-scores a8)
  (check-type low-quality-threshold (integer 1 60))
  (check-type window-size (integer 1 4096))
  (flet ((q->p% (q) (q->p q))) ; only inline the lookup table once
    (iterate
      (declare (type fixnum       q run i w w-1 best best-end)
               (type double-float p sum sum-threshold))

      (with w = window-size)
      (with w-1 = (1- window-size))
      (with sum = 0.0d0)
      (with sum-threshold = (* w (q->p% low-quality-threshold)))
      (with best = 0)
      (with best-end = -1)
      (with run = 0)

      (for q :in-vector quality-scores :with-index i)
      (for p = (q->p% q))

      (incf sum p)

      (cond ((< i w-1) (next-iteration))
            ((>= i w) (decf sum (q->p% (aref quality-scores (- i w))))))

      (if (> sum sum-threshold) ; in p(error) space, higher is *worse*
        (incf run)
        (setf run 0))

      (when (> run best)
        (setf best run
              best-end (1+ i)))

      (finally (return
                 (if (plusp best)
                   ;; Convert length in windows to length/coords in bases.
                   (let* ((region-length (+ best w-1))
                          (region-end best-end)
                          (region-start (- region-end region-length)))
                     (declare (type fixnum region-length region-start region-end))
                     (values region-length region-start region-end))
                   (values 0 nil nil)))))))


;;;; Plotting -----------------------------------------------------------------
(defun index-hits (clusters)
  "Index the minimizers hits inside each cluster of `clusters`.

  Return a hash table of `{hit: n}` where `n` is a cluster number (starting at 1)
  for all hits.

  This will be a large hash table, but it only used to color hits by cluster
  when plotting.

  "
  (iterate
    (for i :from 1)
    (for hits :in clusters)
    (doseq (hit hits)
      (collect-hash (hit i)))))

(defun choose-tics (sequence &aux (length (length sequence)))
  "Return an appropriate value for axis tics based on the length of `sequence`."
  (cond
    ((< length    500)    50)
    ((< length   5000)   100)
    ((< length  10000)   500)
    ((< length  50000)  1000)
    ((< length 100000)  5000)
    (t 10000)))


(defun parse-move-table (move-table)
  (nreverse
    (loop
      :with stride = (aref move-table 0)
      :with result = (make-array (reduce #'+ move-table :start 1)
                                 :element-type 'fixnum
                                 :initial-element 0)
      :with ri = (1- (length result))
      :for mi :from (1- (length move-table)) :downto 1
      :do (progn
            (when (minusp ri)
              (assert (null (find-if-not #'zerop move-table :start 1 :end (1+ mi))))
              (loop-finish))
            (incf (aref result ri) stride)
            (ecase (aref move-table mi)
              (0 (progn))
              (1 (decf ri))))
      :finally (return result))))

(defun bin-move-table-signal (move-table sequence-length bins)
  (check-type move-table (simple-array fixnum (*)))
  (iterate
    (with signal = (parse-move-table move-table))
    (with bin-width = (truncate sequence-length bins))
    (with d = (coerce bin-width 'double-float))
    (repeat bins)
    (for bin-start :from 0 :by bin-width)
    (for bin-end = (+ bin-start bin-width))
    (for bin-sum = (reduce #'+ signal :start bin-start :end bin-end))
    (collect (/ bin-sum d) :result-type 'vector)))

(defun-inline cluster-color (cluster)
  (case cluster
    ((nil) (values #x00 #x00 #x00))
    ;; https://brand.umich.edu/design-resources/colors/
    (1 (values #x9A #x33 #x24))
    (2 (values #x75 #x98 #x8d))
    (3 (values #xA5 #xA5 #x08))
    (4 (values #x00 #xB2 #xA9))
    (5 (values #x70 #x20 #x82))
    (6 (values #x9B #x9A #x6D))
    (7 (values #xD8 #x60 #x18))
    (8 (values #x57 #x52 #x94))
    (t (values #x2F #x65 #xA7))))

;; The plot should look like:
;;
;; +------------------------------------------------------+
;; | Title...                pad                          |
;; |   +----------------------------------------------+   |
;; |   |                                              |   |
;; |   |                                              |   |
;; |   |                                              |   |
;; |   |                                              |   |
;; |   |                                              |   |
;; |   |                                              |   |
;; |   |                                              |   |
;; |   |                                              |   |
;; |pad|                Minimizer Plot                |pad|
;; |   |                                              |   |
;; |   |                                              |   |
;; |   |                                              |   |
;; |   |                                              |   |
;; |   |                                              |   |
;; |   |                                              |   |
;; |   |                                              |   |
;; |   |                                              |   |
;; |   +----------------------------------------------+   |
;; |                         pad                          |
;; |   +----------------------------------------------+   |
;; |   |             Annotation Plots                 |   |
;; |   +----------------------------------------------+   |
;; |                         pad                          |
;; |   +----------------------------------------------+   |
;; |   |                 Quality Plot                 |   |
;; |   +----------------------------------------------+   |
;; |                         pad                          |
;; |   +----------------------------------------------+   |
;; |   |                 Signal Speed                 |   |
;; |   +----------------------------------------------+   |
;; |                       1/2 pad                        |
;; | Legend1               1/2 pad                        |
;; | Legend2               1/2 pad                        |
;; +------------------------------------------------------+
;;
;; mp = minimizer plot
;; ap = annotation plot
;; qp = quality plot
;; sp = signal plot
;; lg = legend
;; ti = title
;; s = start
;; e = end
;; w = width
;; h = height

(defun render-minimizers (sequence quality-scores move-table result)
  (let* ((has-minimizers (not (slot-boundp result 'skip-reasons)))
         (has-annotations (not (null *annotation-tracks*)))
         (foldback-point (foldback-point result))
         (id (read-id result))
         (minimizers (read-minimizers result))
         (llqr (llqr result))
         (llqr-start (llqr-start result))
         (llqr-end (llqr-end result))
         (annotation-track-height 14)
         (annotation-colors '#1=(((#x1F #x78 #xB4) (#x33 #xA0 #x2C))
                                 ((#xE3 #x1A #x1C) (#x6A #x3D #x9A))
                                 ((#xB1 #x59 #x28) (#xFF #x7F #x00)) . #1#))
         (size 512)
         (pad 26)
         (pad/2 (truncate pad 2))
         (pad/4 (truncate pad 4))
         ;; Always using the llqma window size results in too many pixels drawn
         ;; for larger reads.  todo fix this by being smarter about drawing?
         (avg-window (max *low-quality-window-size*
                          (truncate (length sequence) 200)))
         (w size)
         (tix 5) ; title
         (tiy 5)
         (mpw (- size (* 2 pad))) ; minimizer plot
         (mph (- size (* 2 pad)))
         (mpx pad)
         (mpxs (+ mpx 0))
         (mpxe (+ mpx mpw))
         (mpy pad)
         (mpys (+ mpy 0))
         (mpye (+ mpy mph))
         (msx (truncate w 2)) ; minimizers skipped message
         (msy (+ mpys (truncate (- mpye mpys) 2)))
         (apw (- size (* 2 pad))) ; annotation plots
         (aph (* annotation-track-height (length *annotation-tracks*)))
         (apx pad)
         (apxs (+ apx 0))
         (apxe (+ apx apw))
         (apy (+ mph (* 2 pad) -6))
         (apys (+ apy 0))
         ;; (apye (+ apy aph))
         (qpw (- size (* 2 pad))) ; quality plot
         (qph (truncate size 10))
         (qpx pad)
         (qpxs (+ qpx 0))
         (qpxe (+ qpx qpw))
         (qpy (if has-annotations
                (+ apy aph pad/4)
                (+ mpye pad/2)))
         (qpys (+ qpy 0))
         (qpye (+ qpy qph))
         (spw (- size (* 2 pad))) ; signal plot
         (sph (if move-table (truncate size 10) 0))
         (spx pad)
         (spxs (+ spx 0))
         (spxe (+ spx spw))
         (spy (+ qpye pad))
         ;; (spys (+ spy 0))
         (spye (+ spy sph))
         (signal (when move-table
                   (bin-move-table-signal move-table
                                          (length sequence)
                                          (min (length sequence) spw))))
         (lgh pad) ; legend
         (h (+ pad mph pad
               (if has-annotations (+ aph pad/2) 0)
               qph
               (if move-table (+ sph pad) 0)
               lgh))
         (lg1x 5)
         (lg1y (- h pad))
         (lg2x 5)
         (lg2y (- h pad/2))
         (png (make-instance 'zpng:png
                :color-type :truecolor
                :width w
                :height h))
         (img (zpng:data-array png))
         (len (length sequence))
         (idx (when has-minimizers
                (index-hits (foldback-clusters result))))
         (pp 1)
         (pixels-per-base (/ mpw len))
         (tic-bases (choose-tics sequence))
         (tic-width (truncate (* tic-bases pixels-per-base)))
         (tic-length (truncate pad 3)))
    (dotimes (i (array-total-size img)) ; fill white
      (setf (row-major-aref img i) #xFF))
    (labels ((base->x (base) (values (truncate (map-range 0 len mpx (+ mpx mpw) base))))
             (base->y (base) (values (truncate (map-range 0 len (+ mpy mph) mpy base))))
             (quality->y (q) (values (truncate (map-range 0 50 qpye qpys q))))
             (hit-x (hit) (base->x (hit-l1 hit)))
             (hit-y (hit) (base->y (hit-l2 hit)))
             (ink (x y &optional (a #x30))
               "Darken a pixel."
               (dotimes (ch 3)
                 (setf (aref img y x ch) (max 0 (- (aref img y x ch) a)))))
             (draw (x y r &optional (g r) (b g))
               "Set a pixel to a color."
               (setf (aref img y x 0) r
                     (aref img y x 1) g
                     (aref img y x 2) b))
             (draw-unclustered-point (x y)
               (when (< len 20000)
                 (ink (- x pp) y)
                 (ink (+ x pp) y)
                 (ink x (- y pp))
                 (ink x (+ y pp)))
               (ink x y #x50))
             (draw-point (x y r &optional (g r) (b g))
               (draw (- x pp) y r g b)
               (draw (+ x pp) y r g b)
               (draw x (- y pp) r g b)
               (draw x (+ y pp) r g b)
               (draw x y r g b))
             (char-width (char)
               (let ((pixels (gethash char *m5x7*)))
                 (assert pixels () "No glyph available for ~S" char)
                 (array-dimension pixels 0)))
             (string-width (string)
               (+ (1- (length string))
                  (summation string :key #'char-width)))
             (draw-char (x y char)
               (let ((pixels (gethash char *m5x7*)))
                 (assert pixels () "No glyph available for ~S" char)
                 (destructuring-bind (width height) (array-dimensions pixels)
                   (dotimes (dx width)
                     (dotimes (dy height)
                       (when (plusp (aref pixels dx dy))
                         (draw (+ x dx) (+ y dy) #x00))))
                   width)))
             (draw-string (x y string &key centered)
               (when centered
                 (decf x (truncate (string-width string) 2)))
               (iterate
                 (for char :in-string string)
                 (for w = (draw-char x y char))
                 (incf x (1+ w)))))
      ;; Title
      (draw-string tix tiy (format nil "Read ~A (~A)"
                                   id (string-downcase (classification result))))
      ;; Legend
      (draw-string lg1x lg1y (format nil "tics = ~:Dbp, read length = ~:D, ~A"
                                     tic-bases len
                                     (if foldback-point
                                       (format nil "foldback point = ~:D" foldback-point)
                                       "no foldback detected")))
      (draw-string lg2x lg2y (format nil "mean Q = ~,1F, llqr = ~D, monotony = ~,2F"
                                     (mean-qscore result)
                                     llqr
                                     (monotony result)))
      (if (not has-minimizers)
        (let ((y msy))
          (draw-string msx y "Minimizer computation skipped due to:" :centered t)
          (incf y 13)
          (when (member :monotony (skip-reasons result))
            (draw-string msx y (format nil "High monotony (~,2F >= ~,2F)"
                                       (monotony result) *monotony-threshold*)
                         :centered t)
            (incf y 13))
          (when (member :qscore (skip-reasons result))
            (draw-string msx y (format nil "Low mean Q score (~,1F < ~,1F)"
                                       (mean-qscore result) *minimum-qscore*)
                         :centered t)
            (incf y 13)))
        (progn
          ;; X/Y axes and tic marks
          (do-range ((x mpxs mpxe)) (draw x mpye #x66))
          (do-range ((y mpys mpye)) (draw mpxs y #x66))
          (iterate (for tx :from (+ mpxs tic-width) :to mpxe :by tic-width) ; xticx
                   (loop :for ty :from mpye
                         :repeat tic-length
                         :do (draw tx ty #x66)))
          (iterate (for ty :from (- mpye tic-width) :downto mpxs :by tic-width) ; ytics
                   (loop :for tx :downfrom mpxs
                         :repeat tic-length
                         :do (draw tx ty #x66)))
          ;; Draw dashed line at foldback point
          (when foldback-point
            (iterate (with x = (base->x foldback-point))
                     (for y :from mpys :to qpye :by 4)
                     (draw x y      #xBF #xB0 #x86)
                     (draw x (1+ y) #xBF #xB0 #x86)))

          ;; Unclustered hits
          (doseq (hit (foldback-hits result))
            (let ((cluster (gethash hit idx)))
              (unless cluster
                (draw-unclustered-point (hit-y hit) (hit-x hit)))))
          ;; Clustered hits
          (doseq (hit (foldback-hits result))
            (let ((cluster (gethash hit idx)))
              (when cluster
                (multiple-value-call #'draw-point
                  (hit-y hit) (hit-x hit) (cluster-color cluster)))))))
      ;; Annotation tracks
      (loop :for track in *annotation-tracks*
            :for y :from apys :by annotation-track-height
            :for (id-clusters rc-clusters) = (match-annotation-track track minimizers sequence)
            :for name = (track-name track)
            :for ((idR idG idB) (rcR rcG rcB)) :in annotation-colors
            :do (progn
                  (draw-string 2 (- y 4)
                               (subseq name 0 (min 4 (length name))))
                  (do-range ((x apxs apxe))
                    (draw x y #xCC))
                  (do-range ((x (base->x 0)
                                (min (base->x (length (track-sequence track)))
                                     (1- size))))
                    (draw x y #x00))
                  (loop :for (pos1 . pos2) :in id-clusters
                        :for x1 = (base->x pos1)
                        :for x2 = (base->x pos2)
                        :do (do-irange ((x x1 x2)
                                        (dy -3 -1))
                              (draw x (+ y dy) idR idG idB)))
                  (loop :for (pos1 . pos2) :in rc-clusters
                        :for x1 = (base->x pos1)
                        :for x2 = (base->x pos2)
                        :do (do-irange ((x x1 x2)
                                        (dy 1 3))
                              (draw x (+ y dy) rcR rcG rcB)))))
      ;; Quality scores labels
      (let ((q20y (quality->y 20))
            (lqy (quality->y *low-quality-threshold*)))
        (do-range ((x qpxs qpxe))
          (draw x (1- qpys) #xCC)
          (draw x (1+ qpye) #xCC))
        (do-range ((x (1- qpxs) (1+ qpxe)))
          (draw x q20y #x00 #xAA #x00)
          (draw x lqy #xAA #x00 #x00))
        (draw-string (- qpxs (min pad 20)) (- qpys 4) "Q50")
        (draw-string (- qpxs (min pad 20)) (- q20y 3) "Q20")
        (draw-string (- qpxs (min pad 15)) (- qpye 2) "Q0"))
      ;; Highlight low-quality region by marking beneath the x-axis
      (when (plusp llqr)
        (do-range ((dy 0 4))
          (do-range ((x (base->x llqr-start) (base->x llqr-end)))
            (draw x (+ 1 dy qpye) #xCC #x00 #x00))))
      ;; Quality scores
      (iterate (with qs = (make-array avg-window))
               (for b :from 0)
               (for q :in-vector quality-scores)
               (setf (aref qs (mod b avg-window)) q)
               (for x = (base->x b))
               (for y = (quality->y q))
               (ink x y #x66)
               (when (> b avg-window)
                 (for xa = (base->x b))
                 (for ya = (quality->y (/ (summation qs) avg-window)))
                 (draw xa ya #xFF 0 0)))
      ;; Signal speed
      (flet ((signal->dy (s)
               (min (truncate s) 200)))
        (when move-table
          (do-range ((x spxs spxe)
                     (s 0 50 5))
            (draw x (- spye (signal->dy s)) #xCC))
          (if (>= (length sequence) spw)
            (loop :for x :from spxs :below spxe
                  :for i :from 0
                  :for s = (aref signal i)
                  :do (loop
                        :for dy :from (signal->dy s) :downto 0
                        :do (draw x (- spye dy) #x00 #x00 #XFF)))
            nil))))
    (ensure-directories-exist (format nil "~A/plots/" *output-directory*))
    (zpng:write-png png (format nil "~A/plots/~A.png" *output-directory* id)
                    :if-exists :supersede)))


;;;; Main Entry Point ---------------------------------------------------------
(defclass* minimera-result ()
  (read-id
   read-length
   read-minimizers
   classification
   monotony
   mean-qscore
   (prefix-mean-qscore :initform nil)
   (foldback-prefix-mean-qscore :initform nil)
   (foldback-suffix-mean-qscore :initform nil)
   skip-reasons
   (foldback-point :initform nil)
   foldback-hits
   foldback-clusters
   llqr
   (llqr-start :initform nil)
   (llqr-end :initform nil)
   processing-time))

(defun foldbackp (result)
  (eql (classification result) :foldback))


(defun process-read% (id sequence quality-scores)
  (check-type sequence a8)
  (multiple-value-bind (llqr llqr-start llqr-end)
      (compute-llqr quality-scores)
    (let* ((read-length (length sequence))
           (m1 (minimizers *k* *w* sequence))
           (monotony (estimate-monotony m1))
           (mean-qscore (compute-mean-qscore quality-scores))
           (skip-reasons (list)))
      (when (>= monotony *monotony-threshold*)
        (push :monotony skip-reasons))
      (when (< mean-qscore *minimum-qscore*)
        (push :qscore skip-reasons))
      (if skip-reasons
        (make-instance 'minimera-result
          :read-id id
          :read-length read-length
          :read-minimizers m1
          :classification :failed
          :monotony monotony
          :mean-qscore mean-qscore
          :skip-reasons skip-reasons
          :llqr llqr
          :llqr-start llqr-start
          :llqr-end llqr-end)
        ;; CAREFUL HERE, this nreverse breaks using this function interactively
        (let* ((m2 (minimizers *k* *w* (nreverse-complement-bytes sequence)))
               (hits (hits m1 m2))
               (clusters (cluster hits))
               (results (find-foldback-clusters read-length clusters))
               (result (alexandria:extremum results #'> :key #'cluster-lengths))
               (classification (if result :foldback :normal))
               (foldback-position (when result
                                    (compute-foldback-position result)))
               (foldback-qscores (when result
                                   (compute-foldback-qscores quality-scores result))))
          (make-instance 'minimera-result
            :read-id id
            :read-length read-length
            :read-minimizers m1
            :classification classification
            :monotony monotony
            :mean-qscore mean-qscore
            :prefix-mean-qscore (when foldback-qscores (first foldback-qscores))
            :foldback-prefix-mean-qscore (when foldback-qscores (second foldback-qscores))
            :foldback-suffix-mean-qscore (when foldback-qscores (third foldback-qscores))
            :llqr llqr
            :llqr-start llqr-start
            :llqr-end llqr-end
            :foldback-point foldback-position
            :foldback-hits hits
            :foldback-clusters clusters))))))

(defun process-read (id sequence quality-scores move-table)
  (let ((result (process-read% id sequence quality-scores)))
    (when (if (foldbackp result)
            *plot-foldbacks*
            *plot-normal*)
      (render-minimizers sequence quality-scores move-table result))
    result))


;;;; Toplevel -----------------------------------------------------------------
(defparameter *csv-headers*
  (list "read-id"
        "read-length"
        "classification"
        "monotony"
        "foldback-point"
        "mean-qscore"
        "prefix-mean-qscore"
        "foldback-prefix-mean-qscore"
        "foldback-suffix-mean-qscore"
        "llqr"
        "llqr-start"
        "llqr-end"
        "processing-time-microsec"))


(defun output% (csv-stream result)
  (conserve:write-row (list (read-id result)
                            (princ-to-string (read-length result))
                            (string-downcase (classification result))
                            (format nil "~,4F" (monotony result))
                            (if (foldback-point result)
                              (princ-to-string (foldback-point result))
                              "")
                            (princ-to-string (mean-qscore result))
                            (princ-to-string (or (prefix-mean-qscore result) ""))
                            (princ-to-string (or (foldback-prefix-mean-qscore result) ""))
                            (princ-to-string (or (foldback-suffix-mean-qscore result) ""))
                            (princ-to-string (llqr result))
                            (princ-to-string (or (llqr-start result) ""))
                            (princ-to-string (or (llqr-end result) ""))
                            (princ-to-string (processing-time result)))
                      csv-stream))


(defun gettime ()
  (multiple-value-bind (sec microsec)
      (sb-ext:get-time-of-day)
    (+ (* 1000000 sec) microsec)))

(defun work-fastq% (fastq-read)
  (let* ((start (gettime))
         (id (faster:parse-sequence-id (faster:id fastq-read)))
         (seq (faster:seq fastq-read))
         (qs (faster:n-bytes-to-quality-scores (faster:qs fastq-read)))
         (result (process-read id seq qs nil))
         (end (gettime))
         (processing-time (- end start)))
    (setf (processing-time result) processing-time)
    result))

(defun work-bam% (alignment)
  (let* ((start (gettime))
         (id (faster/bam::alignment-read-name alignment))
         (seq (let ((seq (faster/bam::alignment-sequence/fastq alignment)))
                (when (faster/bam::revcompp alignment)
                  (nreverse-complement-bytes seq))
                seq))
         (qs (faster/bam::alignment-phred alignment))
         (move-table (let ((mv (faster/bam::alignment-tag alignment "mv")))
                       (when mv
                         (cdr (faster/bam::tag-value mv)))))
         (result (process-read id seq qs move-table))
         (end (gettime))
         (processing-time (- end start)))
    (setf (processing-time result) processing-time)
    result))


(defun run (filename)
  (when *annotation-path*
    (setf *annotation-tracks* (parse-annotation-file *annotation-path*)))
  (ensure-directories-exist (uiop:ensure-directory-pathname *output-directory*))
  (with-open-file (csv-stream (format nil "~A/foldbacks.csv" *output-directory*)
                              :direction :output
                              :if-exists :supersede)
    (conserve:write-row *csv-headers* csv-stream)
    (faster:run filename
                :work-function (ecase *input-format*
                                 (:fastq #'work-fastq%)
                                 (:bam #'work-bam%))
                :format *input-format*
                :output-function (curry #'output% csv-stream)
                :interactive *interactive*
                :report-progress *report-progress*
                :worker-threads *worker-threads*
                :queue-size (* 10 *worker-threads*))))

;; (setf *output-directory* "results/multiqs")
;; (setf *input-format* :bam)

;; (run "data/bams/20251216_1450_P2S-01935-A_PBE89432_059690b8.sup.unmapped.bam")
