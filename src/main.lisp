(in-package :minimera)

;;;; Configuration ------------------------------------------------------------
(defparameter *k* 8
  "Size of kmers (in base pairs) to use for minimizers.")

(defparameter *w* 16
  "Size (total) of windows (in base pairs) to use for minimizer sketches.")

(defparameter *foldback-position-epsilon* 50
  "How close (in base pairs) a cluster must be to the beginning of read 2 and end of read 1 to be considered as a foldback cluster.")

(defparameter *intercept-epsilon* 30
  "Epsilon (in y-intercept space) used to cluster colinear points during the initial clustering step.")

(defparameter *gap-epsilon* 300
  "Maximum width (in y-intercept space) of an allowable gap in a cluster.  Gaps wider than this will result in splitting the cluster.")

(defparameter *minimum-cluster-length* 40
  "Minimum number of allowable points in a cluster.  Clusters with fewer than this many hits will be removed.")

(defparameter *minimum-foldback-length-absolute* 50
  "Minimum length (in base pairs) of a foldback region.  Regions shorter than this will be excluded.")

(defparameter *minimum-foldback-length-relative* 0.05
  "Minimum length (as a fraction of total read length) of a foldback region.  Regions shorter than this will be excluded.")

(defparameter *plot-foldbacks* nil)

(defparameter *plot-normal* nil)

(defparameter *worker-threads* 6)
(defparameter *input-queue*  (lparallel.queue:make-queue :fixed-capacity 12))
(defparameter *output-queue* (lparallel.queue:make-queue :fixed-capacity 12))

;;;; Data Generation ----------------------------------------------------------
(defparameter *foldback-chance* 0.2)

(defun random-dna (n)
  (values (iterate (repeat n)
                   (collect (random-elt "ACTG") :result-type 'string))
          (format nil "~A real" (uuid:make-v4-uuid))))

(defparameter *snp-chance*    5/100)
(defparameter *indel-chance* 10/100)

(defun mutate-dna (seq)
  (iterate
    (for ch :in-string seq)
    (unless (randomp *indel-chance*)
      (collect (if (randomp *snp-chance*)
                 (random-elt "ATCG")
                 ch)
               :result-type 'string))
    (when (randomp *indel-chance*)
      (do-repeat (random-range 1 5)
        (collect (random-elt "ATCG") :result-type 'string)))))

(defun random-foldback (n)
  (let* ((foldback-fraction (random-range 0.3 1.0))
         (real-length (truncate (/ n (1+ foldback-fraction))))
         (foldback-length (- n real-length))
         (prefix (iterate (repeat real-length)
                          (collect (random-elt "ACTG") :result-type 'string)))
         (suffix (reverse-complement (subseq prefix (- (length prefix) foldback-length))))
         (adapter (random-dna (random-range 80 100))))
    ;; (format t "~v@A~%~v@A" real-length prefix real-length (reverse suffix))
    (values (concatenate 'string prefix adapter suffix)
            (format nil "~A foldback ~D" (uuid:make-v4-uuid) foldback-length)
            foldback-length)))

(defun random-quality (n)
  (iterate (repeat n)
           (collect (random-elt "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQ") :result-type 'string)))

(defun random-fastq-entry (n)
  (multiple-value-bind (seq id)
      (if (randomp *foldback-chance*)
        (random-foldback n)
        (random-dna n))
    (format nil ">~A~%~A~%+~%~A~%"
            (string-downcase id)
            seq
            (random-quality n))))

(defun write-random-fastq (filename entries min-length max-length)
  (with-open-file (f filename :direction :output :if-exists :supersede)
    (do-repeat entries
      (write-string (random-fastq-entry
                      (random-range-inclusive min-length max-length))
                    f))))

(defun print-foldback (seq n)
  (let ((width (- (length seq) n)))
    (format t "~v@A~%~v@A"
            width (subseq seq 0 width)
            width (reverse (subseq seq width))))
  (values))


;;;; Longest Increasing Subsequence -------------------------------------------

;;; Based on Rosetta code, cleaned up and extended to allow :key function:
;;;
;;; https://rosettacode.org/wiki/Longest_increasing_subsequence#Common_Lisp:_Using_the_Patience_Sort_approach

(defun patience-insert (item piles &key key)
  ;; piles is an array that will be mutated.  each pile is ((item . backpointer)
  ;; …) where backpointer points at the previous element (and this conveniently
  ;; happens to resemble vanilla list structure).
  (loop :for pile :across piles
        :and prev = nil :then pile
        :for i from 0
        :for top = (car pile)
        :do (when (<= (funcall key item)
                      (funcall key (car top)))
              (push (cons item (car prev)) (aref piles i))
              (return))
        :finally (vector-push-extend (list (cons item top)) piles)))

(defun longest-increasing-subsequence (sequence &key (key #'identity) (result-type 'list))
  "Return the longest increasing subsequence of elements from `sequence`.

  Elements will be compared with `<`.  `key` will be called on each first.

  The result will be a sequence of type `result-type`.

  "
  (let ((piles (make-array 256 :adjustable t :fill-pointer 0)))
    (map nil (lambda (item) (patience-insert item piles :key key)) sequence)
    (nreverse (coerce (car (aref piles (1- (length piles)))) result-type))))


;;;; Hashing ------------------------------------------------------------------

;;; https://naml.us/post/inverse-of-a-hash-function/

(deftype u64 ()
  `(integer 0 ,(1- (expt 2 64))))

(defun-inline wrap64 (n)
  (mod n (expt 2 64)))

(defun-inline shl (x bits)
  (logand (ash x bits)
          (1- (ash 1 64))))

(defun-inline shr (x bits)
  (logand (ash x (- bits))
          (1- (ash 1 64))))

(defun-inline hash64 (n)
  (declare (optimize (speed 3) (safety 1) (debug 1)))
  (etypecase n
    (fixnum
      (_ n
        (wrap64 (+ (lognot _) (shl _ 21)))
        (wrap64 (logxor _ (shr _ 24)))
        (wrap64 (+ (+ _ (shl _ 3))
                   (shl _ 8)))
        (wrap64 (logxor _ (shr _ 14)))
        (wrap64 (+ (+ _ (shl _ 2))
                   (shl _ 4)))
        (wrap64 (logxor _ (shr _ 28)))
        (wrap64 (+ _ (shl _ 31)))))))

(defun ihash64 (n)
  (declare (optimize (speed 3) (safety 1) (debug 1)))
  ;;; https://naml.us/post/inverse-of-a-hash-function/
  (etypecase n
    (fixnum
      (let* ((key n)

             ;;   // Invert key = key + (key << 31)
             ;;   tmp = key - (key <<31);
             ;;   key = key - (tmp <<31);
             (tmp (wrap64 (- key (shl key 31))))
             (key (wrap64 (- key (shl tmp 31))))

             ;;   // Invert key = key ^ (key >> 28)
             ;;   tmp = key ^ key >>28;
             ;;   key = key ^ tmp >>28;
             (tmp (wrap64 (logxor key (shr key 28))))
             (key (wrap64 (logxor key (shr tmp 28))))

             ;;   // Invert key *= 21
             ;;   key *= 14933078535860113213u;
             (key (wrap64 (* key 14933078535860113213)))

             ;;   // Invert key = key ^ (key >> 14)
             ;;   tmp = key ^ key >> 14;
             ;;   tmp = key ^ tmp >> 14;
             ;;   tmp = key ^ tmp >> 14;
             ;;   key = key ^ tmp >> 14;
             (tmp (wrap64 (logxor key (shr key 14))))
             (tmp (wrap64 (logxor key (shr tmp 14))))
             (tmp (wrap64 (logxor key (shr tmp 14))))
             (key (wrap64 (logxor key (shr tmp 14))))

             ;;   // Invert key *= 265
             ;;   key *= 15244667743933553977u;
             (key (wrap64 (* key 15244667743933553977)))

             ;;   // Invert key = key ^ (key >> 24)
             ;;   tmp = key ^ key >> 24;
             ;;   key = key ^ tmp >> 24;
             (tmp (wrap64 (logxor key (shr key 24))))
             (key (wrap64 (logxor key (shr tmp 24))))

             ;;   // Invert key = (~key) + (key << 21)
             ;;   tmp = ~key;
             ;;   tmp = ~( key -( tmp <<21));
             ;;   tmp = ~( key -( tmp <<21));
             ;;   key = ~( key -( tmp <<21));
             (tmp (wrap64 (lognot key)))
             (tmp (wrap64 (lognot (- key (shl tmp 21)))))
             (tmp (wrap64 (lognot (- key (shl tmp 21)))))
             (key (wrap64 (lognot (- key (shl tmp 21))))))
        key))))

(defun phihash-string (string &key (start 0) (end (length string)))
  "Return the φ-hash of the subsequence of `string` bounded by `start` and `end`."
  (iterate
    (with hash = 0)
    (for ch :in-string string :from start :below end)
    (zapf hash (_ %
                 (ash _ 2)
                 (logior _ (ecase ch
                             (#\A #b00)
                             (#\C #b01)
                             (#\G #b10)
                             (#\T #b11)))))
    (returning (hash64 hash))))

(defun-inline phihash-chunk (chunk)
  (hash64 chunk))


;;;; FASTQ (Reference) --------------------------------------------------------
(defun do-fastq (function stream)
  "Read a FASTQ from `stream` and call `(function id sequence)` on each entry."
  (iterate
    (until (eql :eof (peek-char nil stream nil :eof)))
    (assert (char= #\> (read-char stream)))
    (for id = (read-line stream))
    (for seq = (read-line stream))
    (assert (char= #\+ (read-char stream)))
    (peek-char #\newline stream) (read-char stream) ; +
    (peek-char #\newline stream) (read-char stream) ; quality scores
    (funcall function id seq)))


;;;; FASTQ (Optimized) --------------------------------------------------------
(deftype buffer-length ()
  `(integer 1 ,array-dimension-limit))

(deftype u8 ()
  '(unsigned-byte 8))

(deftype a8 ()
  '(simple-array u8 (*)))

(defstruct (fastq-read (:conc-name nil))
  (id nil :type (simple-array u8 (*)))
  (seq nil :type (simple-array u8 (*))))

(defun-inline copy-buffer (buffer start end)
  (let ((result (make-array (- end start) :element-type 'u8)))
    (replace result buffer :start2 start :end2 end)
    result))

(defun bs (buffer)
  (map 'string #'code-char buffer))

(declaim (ftype (function (list) a8) revcat-buffers))

(defun-inline revcat-buffers (buffers)
  (cond ((null buffers) (make-array 0 :element-type 'u8))
        ((null (rest buffers)) (first buffers))
        (t (iterate
             (with total = (reduce #'+ buffers :key #'length))
             (with-result result = (make-array total :element-type 'u8))
             (with end = total)
             (for buffer in buffers)
             (for n = (length buffer))
             (replace result buffer :start1 (- end n))
             (decf end n)))))

(defun map-fastq (function stream &key (buffer-length (* 10 (expt 2 20))))
  (declare (optimize (speed 3) (safety 1) (debug 1)))
  (check-type buffer-length buffer-length)
  (let ((buffer (make-array buffer-length :element-type 'u8))
        (i 0)
        (end 0)
        (refills 0))
    (declare (type a8 buffer)
             (type fixnum refills)
             (type (integer 0 #.array-total-size-limit) i end))
    (labels ((ensure-refill ()
               "Refill buffer if necessary and return whether we have data left."
               (if (>= i end)
                 (progn (setf i 0
                              end (read-sequence buffer stream))
                        (incf refills)
                        (not (zerop end))) ; eof if we didn't read anything
                 t)) ; not at end of buffer, definitely not eof yet
             (chomp-byte (expected-char)
               (unless (ensure-refill)
                 (error "Unexpected EOF, expected ~S." expected-char))
               (unless (= (char-code expected-char) (aref buffer i))
                 (error "Bad byte ~X, expected ~X." (char-code expected-char) (aref buffer i)))
               (incf i))
             (skip-line% (&key (allow-eof nil))
               "Skip the rest of the line, return number of bytes skipped (not including the NL)."
               (iterate
                 (with-result bytes-read = 0)
                 (unless (ensure-refill)
                   (if allow-eof
                     (finish)
                     (error "Unexpected EOF skipping line.")))
                 (for n = (position (char-code #\newline) buffer :start i))
                 (if n
                   (progn (incf bytes-read (- n i))
                          (setf i (1+ n))
                          (finish))
                   (progn (incf bytes-read (- end i))
                          (setf i end)))))
             (read-line% ()
               (iterate
                 (with result = nil)
                 (unless (ensure-refill)
                   (error "Unexpected EOF reading line."))
                 (for n = (position (char-code #\newline) buffer :start i))
                 (if n
                   (progn (push (copy-buffer buffer i n) result)
                          (setf i (1+ n))
                          (finish))
                   (progn (push (copy-buffer buffer i end) result)
                          (setf i end)))
                 (returning (revcat-buffers result)))))
      (iterate
        (while (ensure-refill))
        (chomp-byte #\>)
        (for id = (read-line%))
        (for seq = (read-line%))
        (chomp-byte #\+)
        (skip-line%) ; + line
        (for quality-length = (skip-line% :allow-eof t)) ; quality scores
        (unless (= quality-length (length seq))
          (error "Mismatched number of quality scores (~D) for read of length ~D."
                 quality-length (length seq)))
        (funcall function (make-fastq-read :id id :seq seq)))
      ;; (print refills)
      (values))))


;;;; Minimizers ---------------------------------------------------------------
(declaim (inline make-minmimizer))

(defstruct (minimizer (:constructor make-minimizer (hash pos))
                      (:conc-name nil))
  (hash 0 :type u64)
  (pos  0 :type u64))

(defmethod print-object ((o minimizer) s)
  (print-unreadable-object (o s :type t)
    (format s "~16,'0X ~D" (hash o) (pos o))))


;;;; Minimizers (Reference) ---------------------------------------------------
(defun reverse-complement (seq)
  "Return a fresh string of the reverse complement of `seq`."
  (nreverse (map 'string (lambda (ch)
                           (ecase ch
                             (#\A #\T)
                             (#\C #\G)
                             (#\G #\C)
                             (#\T #\A)))
                 seq)))

(defun minimizer (k window window-start)
  "Return the `k`-minimizer of `window`.

  Will return a cons of the (hashed) minimizer and its position in the sequence.

  `window-start` must the the start position of the window in the overall
  sequence, and is used to compute the position.

  "
  (iterate
    (for (kmer ks ke) :window k :on window)
    (for hash = (phihash-string kmer))
    (for pos = (+ window-start ks))
    (finding (make-minimizer hash pos) :minimizing hash)))

(defun minimizers (k w sequence)
  "Return a list of the (`k`, `w`, φ)-mimimizer sketch of `sequence`."
  (iterate
    (with prev)
    (for (window ws we) :window w :on sequence)
    (for minimizer = (minimizer k window ws))
    (unless (and prev (= (pos prev) (pos minimizer)))
      (collect minimizer)
      (setf prev minimizer))))


;;;; Minimizers (Optimized) ---------------------------------------------------
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
  (loop :with result = 0
        :for i :from 0 :below w
        :do (setf result (logior (ash result 2)
                                 (base-bits/id (aref sequence i))))
        :finally (return result)))

(defun minimizer/fast (k w chunk window-start)
  (declare (optimize (speed 3) (safety 1) (debug 1)))
  (check-type k (integer 1 31))
  (check-type w (integer 1 64))
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

(defun minimizers/fast (k w sequence)
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
                 (logior _ (base-bits (aref sequence (1- end)))))) ; add new end base
    (if prev
      ;; If we still have a previous minimizer, we just compare the new one and
      ;; check if it's better than that one.
      (progn
        (setf
          (hash candidate) (ldb (byte k2 0) chunk)
          (hash candidate) (phihash-chunk (hash candidate))
          (pos candidate) (+ start w-k))
        (when (< (hash candidate) (hash prev))
          (let ((next (copy-minimizer candidate)))
            (collect next)
            (setf prev next))))
      ;; Otherwise we don't have a previous winning minimizer, either because
      ;; we're at the beginning of the sequence or it fell off the moving
      ;; window.  Just recompute the window from scratch here.
      (let ((next (minimizer/fast k w chunk start)))
        (collect next)
        (setf prev next)))))


(defun check-implementation (k w data)
  (let ((bytes (map 'a8 #'char-code data)))
    (equalp (loop :for m in (minimizers k w data)
                  :collect (cons (hash m) (pos m)))
            (loop :for m in (minimizers/fast k w bytes)
                  :collect (cons (hash m) (pos m))))))


;;;; Minimizer Hits/Matches ---------------------------------------------------
(defstruct (hit (:constructor make-hit (c l1 l2)))
  (c nil :type fixnum)
  (l1 nil :type fixnum)
  (l2 nil :type fixnum))

(defmethod print-object ((o hit) s)
  (print-unreadable-object (o s :type t)
    (format s "~D ~D ~D" (hit-c o) (hit-l1 o) (hit-l2 o))))

(define-sorting-predicate hit<
  (#'< :key #'hit-c)
  (#'> :key #'hit-l2)
  (#'< :key #'hit-l1))

(defun hits (ms1 ms2)
  "Return a vector of matching minimizer hits from sketches `ms1` and `ms2`.

  The vector will be sorted in parameter-space-intercept-order.

  "
  (let ((index-1 (group-by #'hash ms1 :test #'equal :map #'pos))
        (index-2 (group-by #'hash ms2 :test #'equal :map #'pos)))
    (sort (coerce (iterate
                    (for (minimizer locs-1) :in-hashtable index-1)
                    (for locs-2 = (gethash minimizer index-2))
                    (appending (alexandria:map-product (lambda (l1 l2)
                                                         (make-hit (- l2 l1) l1 l2))
                                                       locs-1 locs-2)))
                  'vector)
          #'hit<)))



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

(defun filter-clusters (clusters)
  "Return a list of clusters, with gaps broken and small clusters removed."
  (let* ((min-length *minimum-cluster-length*)
         (deduped (mapcar (lambda (cluster)
                            (longest-increasing-subsequence
                              (sort cluster #'< :key #'hit-l1)
                              :key #'hit-l2
                              :result-type 'vector))
                          clusters))
         (degapped (alexandria:mappend #'break-gaps deduped))
         (results (remove-if (lambda (cluster)
                               (< (length cluster) min-length))
                             degapped)))
    ;; (prl (map 'list #'length clusters)
    ;;      (map 'list #'length deduped)
    ;;      (map 'list #'length degapped)
    ;;      (map 'list #'length results))
    results))

(defun cluster (hits)
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
      (when (or (null next-hit)
                (> (- (hit-c next-hit) (hit-c hit)) epsilon))
        (collect (subseq hits start next-i) :into results)
        (setf start next-i))
      (finally (return (filter-clusters results))))))


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

  May return an empty list of no clusters look like possible foldback regions.

  "
  (remove-if-not (alexandria:conjoin
                   (curry #'region-in-foldback-position-p read-length)
                   (curry #'region-long-enough-p read-length))
                 clusters))


;;;; Plotting -----------------------------------------------------------------
(defun index-hits (clusters)
  (iterate
    (for i :from 1)
    (for hits :in clusters)
    (format t "cluster ~D: ~D hits~%" i (length hits))
    (doseq (hit hits)
      (collect-hash (hit i)))))

(defun choose-tics (sequence &aux (length (length sequence)))
  (cond
    ((< length   500)    50)
    ((< length  1000)   100)
    ((< length  5000)   500)
    ((< length 10000)  1000)
    ((< length 50000)  5000)
    (t 10000)))

(defun plot-minimizers (info sequence hits clusters)
  (with-open-file (f "hits.csv" :direction :output :if-exists :supersede)
    (conserve:write-row (list "c" "l1" "l2" "cluster") f)
    (let ((cluster-index (index-hits clusters)))
      (doseq (hit hits)
        (conserve:write-row (list (princ-to-string (hit-c hit))
                                  (princ-to-string (hit-l1 hit))
                                  (princ-to-string (hit-l2 hit))
                                  (princ-to-string (gethash hit cluster-index "NA")))
                            f))))
  (write-line "Plotting…")
  (rscript-file "plot.R" (length sequence) (choose-tics sequence))
  (destructuring-bind (id class pos rev)
      (str:split "," info)
    (declare (ignore rev))
    (sh (list "cp" "matches.png" (format nil "plots/~A-~A-~A.png" class id pos)))))


;;;; Scoring ------------------------------------------------------------------
(defun foldbacks (sequence &key info (optimized nil))
  (if optimized
    (check-type sequence a8)
    (check-type sequence string))
  (let* ((read-length (length sequence))
         (hits (if optimized
                 (hits (minimizers/fast *k* *w* sequence)
                       (minimizers/fast *k* *w* (nreverse-complement-bytes sequence)))
                 (hits (minimizers *k* *w* sequence)
                       (minimizers *k* *w* (reverse-complement sequence)))))
         (clusters (cluster hits))
         (results (find-foldback-clusters read-length clusters)))
    ;; (format *debug-io* "~D minimizer hit~:P found (~F/kbp)~%"
    ;;         (length hits)
    ;;         (/ (length hits)
    ;;            (/ (length sequence) 1000)))
    ;; (format *debug-io* "~D cluster~:P found~%" (length clusters))
    ;; (format *debug-io* "Found ~D candidate foldback cluster~:P.~%" (length results))
    (if results
      (when *plot-foldbacks*
        (plot-minimizers info sequence hits clusters))
      (when *plot-normal*
        (plot-minimizers info sequence hits clusters)))
    results))


;;;; Toplevel -----------------------------------------------------------------
(defun compute-foldback-coordinates (cluster alignment-pos reverse?)
  (multiple-value-bind (start end) (cluster-bounds cluster)
    (let ((pos (truncate (+ start end) 2)))
      (values pos ; coord in read
              (if reverse? ; coord in genome
                (- alignment-pos pos)
                (+ alignment-pos pos))))))

(defun parse-read-info (read-info)
  ;; 44c97e05-69a2-401a-a722-ac5c3c9080b1,normal,chr21:19180531,True
  (destructuring-bind (read-name class alignment-location reversed) (str:split "," read-info)
    (destructuring-bind (contig pos) (str:split ":" alignment-location)
      (values read-name
              class
              contig
              (parse-integer pos)
              (alexandria:eswitch (reversed :test #'string=)
                ("True" t)
                ("False" nil))))))

(defun run-read (read-info sequence &key optimized)
  (let* ((candidates (foldbacks sequence :info read-info :optimized optimized))
         (candidate (alexandria:extremum candidates #'> :key #'cluster-lengths))
         (read-info (if optimized
                      (map 'string #'code-char read-info)
                      read-info)))
    (multiple-value-bind (name class contig pos rev) (parse-read-info read-info)
      (declare (ignore class))
      (if candidate
        (multiple-value-bind (read-coord genome-coord) (compute-foldback-coordinates candidate pos rev)
          (cons (conserve:write-row (list name "true" (princ-to-string read-coord)) nil)
                (let ((conserve:*delimiter* #\tab))
                  (conserve:write-row (list contig
                                            (princ-to-string genome-coord)
                                            (princ-to-string (1+ genome-coord))
                                            (format nil "~A" read-info))
                                      nil))))
        (cons (conserve:write-row
                (list name "false" "")
                nil)
              nil)))))

(defun run/slow (filename &key (output-basename "results"))
  (with-open-file (input-stream filename :direction :input)
    (with-open-file (output-stream (format nil "~A.csv" output-basename) :direction :output :if-exists :supersede)
      (conserve:write-row (list "read-id" "is-foldback" "foldback-point")
                          output-stream)
      (with-open-file (bed-stream (format nil "~A.bed" output-basename) :direction :output :if-exists :supersede)
        (do-fastq (lambda (info seq)
                    ;; (format *debug-io* "~%Checking record ~A~%" info)
                    (destructuring-bind (csv-out . bed-out)
                        (run-read info seq :optimized nil)
                      (write-string csv-out output-stream)
                      (when bed-out
                        (write-string bed-out bed-stream))))
                  input-stream)))))

(defparameter *input-done* nil)
(defparameter *work-done* nil)

(defun run/fast (filename &key (output-basename "results"))
  (setf *input-done* nil
        *work-done* nil)
  (bt2:make-thread
    (lambda ()
      (with-open-file (input-stream filename :direction :input :element-type 'u8)
        (map-fastq (lambda (fastq-read)
                     (lparallel.queue:push-queue fastq-read *input-queue*))
                   input-stream))
      (setf *input-done* t))
    :name "FASTQ Reader")
  (dotimes (i *worker-threads*)
    (bt2:make-thread
      (lambda ()
        (loop (multiple-value-bind (fastq-read found)
                  (lparallel.queue:try-pop-queue *input-queue* :timeout 1)
                (cond ((and (not found) *input-done*) (return))
                      (found (lparallel.queue:push-queue
                               (run-read (id fastq-read) (seq fastq-read) :optimized t)
                               *output-queue*))
                      (t (progn)))))
        (setf *work-done* t))
      :name (format nil "Worker Thread ~D" (1+ i))))
  (bt2:join-thread
    (bt2:make-thread
      (lambda ()
        (with-open-file (output-stream (format nil "~A.csv" output-basename) :direction :output :if-exists :supersede)
          (conserve:write-row (list "read-id" "is-foldback" "foldback-point") output-stream)
          (with-open-file (bed-stream (format nil "~A.bed" output-basename) :direction :output :if-exists :supersede)
            (loop (multiple-value-bind (out found)
                      (lparallel.queue:try-pop-queue *output-queue* :timeout 1)
                    (cond ((and (not found) *work-done*) (return))
                          (found (destructuring-bind (csv-out . bed-out) out
                                   (write-string csv-out output-stream)
                                   (when bed-out
                                     (write-string bed-out bed-stream))))
                          (t (progn))))))))
      :name "Output Writer")))



#; Scratch --------------------------------------------------------------------

(defparameter *k* 8)
(defparameter *w* 16)
(defparameter *foldback-position-epsilon* 50)
(defparameter *intercept-epsilon* 30)
(defparameter *gap-epsilon* 300)
(defparameter *minimum-cluster-length* 40)
(defparameter *minimum-foldback-length-absolute* 50)
(defparameter *minimum-foldback-length-relative* 0.05)


(defparameter *snp-chance*    3/100)
(defparameter *indel-chance*  5/100)

(defparameter *seq* (mutate-dna (random-foldback 20000)))
(defparameter *seq* (random-foldback 20000))
(defparameter *seq* (random-dna 20000))

(plot-minimizers *seq*
               :k 8 :w 16
               :intercept-epsilon 30
               :cluster-min 10
               :gap-epsilon 300
               )

(plot-minimizers (alexandria:read-file-into-string "data/okay2")
               :k 8 :w 16
               :intercept-epsilon 30
               :cluster-min 10
               :gap-epsilon 300)

(plot-minimizers (pbpaste)
               :k 8 :w 16
               :intercept-epsilon 30
               :cluster-min 10
               :gap-epsilon 300)

(run "data/9b58328c3b631816942cbe400a807241935897e6.manual-test.2.fastq" :output-basename "results.2")

(run "out.fastq" :output-basename "results.unclassified")

(run "out2.fastq" :output-basename "results.unclassified.2")

(time (run/slow "big.fastq" :output-basename "bench.out"))

(time (run/fast "big.fastq" :output-basename "bench.fast.out"))


;; (defparameter *test-seq* (random-dna 1000))

;; (check-implementation 6 10 *test-seq*)

;; (minimizers 4 12 *test-seq*)
;; (minimizers/fast 4 12 (map 'a8 #'char-code *test-seq*))

;; (defparameter *data/string* (random-dna 10000000))
;; (defparameter *data/bytes* (map 'a8 #'char-code *data/string*))

;; (time (loop :repeat 1 :do (minimizers 15 25 *data/string*)))

;; (profile (loop :repeat 10 :do (minimizers/fast 15 25 *data/bytes*)))

;; (check-implementation 15 25 *data/string*)
