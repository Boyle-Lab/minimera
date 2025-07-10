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


(defparameter *plot-foldbacks* nil
  "Whether to generate plots (with the accompanying R script) of reads determined to be foldbacks.")

(defparameter *plot-normal* nil
  "Whether to generate plots (with the accompanying R script) of reads determined *not* to be foldbacks.")


(defparameter *worker-threads* 6)

(defparameter *output-directory* "results")


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
                 (error "Bad byte ~2,'0X (~S), expected ~2,'0X (~S)."
                        (aref buffer i) (code-char (aref buffer i))
                        (char-code expected-char) expected-char))
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
        (chomp-byte #\@)
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
  (declare (optimize (speed 3) (safety 1) (debug 1)))
  (check-type sequence a8)
  (iterate (with-result result = 0)
           (for base :in-vector sequence :below w)
           (zapf result (wrap64 (logior (ash % 2) (base-bits base))))))

(defun minimizer/fast (k w chunk window-start)
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
      (returning (filter-clusters results)))))


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


;;;; Plotting -----------------------------------------------------------------
(defparameter *plotting-code* (alexandria:read-file-into-string "src/plot.R"))

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
    ((< length   500)    50)
    ((< length  1000)   100)
    ((< length  5000)   500)
    ((< length 10000)  1000)
    ((< length 50000)  5000)
    (t 10000)))

(defun plot-minimizers (id sequence hits clusters)
  (uiop:with-temporary-file (:stream f
                             :pathname p
                             :direction :output
                             :prefix "hits."
                             :type "csv")
    (conserve:write-row (list "c" "l1" "l2" "cluster") f)
    (let ((cluster-index (index-hits clusters)))
      (doseq (hit hits)
        (conserve:write-row (list (princ-to-string (hit-c hit))
                                  (princ-to-string (hit-l1 hit))
                                  (princ-to-string (hit-l2 hit))
                                  (princ-to-string (gethash hit cluster-index "NA")))
                            f)))
    (finish-output f)
    (ensure-directories-exist (format nil "~A/plots/" *output-directory*))
    (rscript *plotting-code*
             (length sequence)
             (choose-tics sequence)
             (princ-to-string p)
             (format nil "~A/plots/~A.png" *output-directory* id))))



;;;; Main Entry Point ---------------------------------------------------------
(defun foldbacks (id sequence &key (optimized nil))
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
    (if results
      (when *plot-foldbacks* (plot-minimizers id sequence hits clusters))
      (when *plot-normal*    (plot-minimizers id sequence hits clusters)))
    results))


;;;; Toplevel -----------------------------------------------------------------
(defun compute-foldback-position (cluster)
  (multiple-value-bind (start end) (cluster-bounds cluster)
    (truncate (+ start end) 2)))

(defun run-read (id sequence &key optimized)
  (let* ((candidates (foldbacks id sequence :optimized optimized))
         (candidate (alexandria:extremum candidates #'> :key #'cluster-lengths)))
    (if candidate
      (list id "true" (princ-to-string (compute-foldback-position candidate)))
      (list id "false" ""))))

(defun run/slow (filename)
  (with-open-file (input-stream filename :direction :input)
    (with-open-file (output-stream (format nil "~A.csv" *output-directory*) :direction :output :if-exists :supersede)
      (conserve:write-row (list "read-id" "is-foldback" "foldback-point") output-stream)
      (do-fastq (lambda (id seq)
                  (conserve:write-row (run-read id seq :optimized nil) output-stream))
                input-stream))))


(defparameter *interactive* t)
(defparameter *input-done* nil)
(defparameter *work-done*  nil)

(defparameter *input-queue*  (lparallel.queue:make-queue :fixed-capacity (* 4 *worker-threads*)))
(defparameter *output-queue* (lparallel.queue:make-queue :fixed-capacity (* 4 *worker-threads*)))

(defmacro do-thread (name &body body)
  `(bt2:make-thread (lambda () ,@body) :name ,name))

(defun parse-read-id (bytes)
  ;; TODO: Remove or document the comma splitting here.
  (first (str:split "," (map 'string #'code-char bytes) :limit 2)))


(defun die (error)
  (adopt:print-error-and-exit
    error
    :exit-function (lambda (code)
                     (sb-ext:exit :code code
                                  :abort t
                                  :timeout 4))))

(defmacro with-exit-during-noninteractive-mode (&body body)
  "Run `body`, killing everything upon an error when running noninteractively."
  (with-gensyms (thunk err)
    `(let ((,thunk (lambda () ,@body)))
       (if *interactive*
           (funcall ,thunk)
           (handler-case (funcall ,thunk)
             (error (,err) (die ,err)))))))

(defun run/fast (filename)
  ;; Set up status variables.
  (setf *input-done* nil *work-done* nil)

  ;; Make sure the output directory exists.
  (ensure-directories-exist (uiop:ensure-directory-pathname *output-directory*))

  ;; Spawn single reader thread.
  (do-thread "Minimera FASTQ Reader"
    (with-exit-during-noninteractive-mode
      (flet ((run% (stream)
               (map-fastq (lambda (fastq-read)
                            (lparallel.queue:push-queue fastq-read *input-queue*))
                          stream)))
        (if (string= filename "-")
          (run% *standard-input*)
          (with-open-file (input-stream filename :direction :input :element-type 'u8)
            (run% input-stream))))
      (setf *input-done* t)))

  ;; Spawn N worker threads.
  (dotimes (i *worker-threads*)
    (do-thread (format nil "Minimera Worker ~D" (1+ i))
      (with-exit-during-noninteractive-mode
        (loop (multiple-value-bind (fastq-read found)
                  (lparallel.queue:try-pop-queue *input-queue* :timeout 1)
                (cond ((and (not found) *input-done*) (return))
                      (found (lparallel.queue:push-queue
                               (run-read (parse-read-id (id fastq-read))
                                         (seq fastq-read) :optimized t)
                               *output-queue*))
                      (t (progn)))))
        (setf *work-done* t))))

  ;; Spawn a single writer thread, and join it to wait for things to finish.
  (bt2:join-thread
    (do-thread "Minimera Output Writer"
      (with-open-file (output-stream (format nil "~A/foldbacks.csv" *output-directory*) :direction :output :if-exists :supersede)
        (conserve:write-row (list "read-id" "is-foldback" "foldback-point") output-stream)
        (loop (multiple-value-bind (out found)
                  (lparallel.queue:try-pop-queue *output-queue* :timeout 1)
                (cond ((and (not found) *work-done*) (return))
                      (found (conserve:write-row out output-stream))
                      (t (progn)))))))))


;;;; UI -----------------------------------------------------------------------
(adopt:define-string *documentation*
  "Minimera is a tool for detecting foldback chimeric reads in Oxford Nanopore ~
   data using minimizers.~@
   ~@
   Reads will be read from the specified fastq file (or stdin if - is given) and ~
   the results written to foldbacks.csv inside the specified output directory.  ~
   If plots are generated they will be inside a plots/ subdirectory of the ~
   output directory.")

(defparameter *o/help*
  (adopt:make-option 'help
    :help "Display help and exit."
    :long "help"
    :short #\h
    :reduce (constantly t)))

(defparameter *o/threads*
  (adopt:make-option 'threads
    :help "Number of worker threads to spawn. Does not include the reader and writer threads (default: 1)."
    :long "threads"
    :short #\j
    :parameter "N"
    :reduce #'adopt:last
    :initial-value 1
    :key #'parse-integer))

(defparameter *o/output-directory*
  (adopt:make-option 'output-directory
    :help "Output directory (required)."
    :long "output"
    :short #\o
    :parameter "PATH"
    :reduce #'adopt:last))

(defparameter *o/k*
  (adopt:make-option 'k
    :help "Size of kmers (in base pairs) to use for minimizers (default: 8)."
    :long "kmer-size"
    :short #\K
    :parameter "N"
    :reduce #'adopt:last
    :initial-value 8
    :key #'parse-integer))

(defparameter *o/w*
  (adopt:make-option 'w
    :help "Size (total) of windows (in base pairs) to use for minimizer sketches (default: 16)."
    :long "window-size"
    :short #\W
    :parameter "N"
    :reduce #'adopt:last
    :initial-value 16
    :key #'parse-integer))

(defparameter *o/foldback-position-epsilon*
  (adopt:make-option 'foldback-position-epsilon
    :help "How close (in base pairs) a cluster must be to the beginning of read 2 and end of read 1 to be considered as a foldback cluster (default: 50)."
    :long "foldback-position-epsilon"
    :short #\F
    :parameter "N"
    :reduce #'adopt:last
    :initial-value 50
    :key #'parse-integer))

(defparameter *o/intercept-epsilon*
  (adopt:make-option 'intercept-epsilon
    :help "Epsilon (in y-intercept space) used to cluster colinear points during the initial clustering step (default: 30)."
    :long "intercept-epsilon"
    :short #\I
    :parameter "N"
    :reduce #'adopt:last
    :initial-value 30
    :key #'parse-integer))

(defparameter *o/gap-epsilon*
  (adopt:make-option 'gap-epsilon
    :help "Maximum width (in y-intercept space) of an allowable gap in a cluster.  Gaps wider than this will result in splitting the cluster (default: 300)."
    :long "gap-epsilon"
    :short #\G
    :parameter "N"
    :reduce #'adopt:last
    :initial-value 300
    :key #'parse-integer))

(defparameter *o/minimum-cluster-length*
  (adopt:make-option 'minimum-cluster-length
    :help "Minimum number of allowable points in a cluster.  Clusters with fewer than this many hits will be removed (default: 40)."
    :long "minimum-cluster-length"
    :short #\C
    :parameter "N"
    :reduce #'adopt:last
    :initial-value 40
    :key #'parse-integer))

(defparameter *o/minimum-foldback-length-absolute*
  (adopt:make-option 'minimum-foldback-length-absolute
    :help "Minimum length (in base pairs) of a foldback region.  Regions shorter than this will be excluded (default: 50)."
    :long "minimum-foldback-length-absolute"
    :short #\A
    :parameter "N"
    :reduce #'adopt:last
    :initial-value 50
    :key #'parse-integer))

(defparameter *o/minimum-foldback-length-relative*
  (adopt:make-option 'minimum-foldback-length-relative
    :help "Minimum length (as a fraction of total read length) of a foldback region.  Regions shorter than this will be excluded (default: 0.05)."
    :long "minimum-foldback-length-relative"
    :short #\R
    :parameter "X"
    :reduce #'adopt:last
    :initial-value 0.05
    :key #'parse-float:parse-float))

(adopt:defparameters (*o/plot/foldbacks* *o/plot/no-foldbacks*)
  (adopt:make-boolean-options 'plot-foldbacks
    :long "plot-foldbacks"
    :help "Generate plots for all reads determined to be foldbacks."
    :help-no "Do not generate plots for reads determined to be foldbacks (the default)."))

(adopt:defparameters (*o/plot/normal* *o/plot/no-normal*)
  (adopt:make-boolean-options 'plot-normal
    :long "plot-normal"
    :help "Generate plots for all reads determined to be normal."
    :help-no "Do not generate plots for reads determined to be normal (the default)."))

(defparameter *ui*
  (adopt:make-interface
    :name "minimera"
    :usage "[OPTIONS] --output PATH FASTQ"
    :summary "detect foldback chimeric reads using minimizers"
    :help *documentation*
    :contents
    (list *o/help*
          *o/threads*
          *o/output-directory*
          (adopt:make-group 'algorithm
            :title "Algorithm Options"
            :options (list *o/k*
                           *o/w*
                           *o/foldback-position-epsilon*
                           *o/intercept-epsilon*
                           *o/gap-epsilon*
                           *o/minimum-cluster-length*
                           *o/minimum-foldback-length-absolute*
                           *o/minimum-foldback-length-relative*))
          (adopt:make-group 'plotting
            :title "Plotting Options"
            :help "Minimera can optionally generate plots of the minimizers for reads.  This requires Rscript and Tidyverse libraries to be installed.  Generating these plots is very slow, but they can be useful to debug edge cases."
            :options (list *o/plot/foldbacks*
                           *o/plot/no-foldbacks*
                           *o/plot/normal*
                           *o/plot/no-normal*)))))


(defun toplevel ()
  (adopt::quit-on-ctrl-c ()
    (multiple-value-bind (arguments options) (adopt:parse-options-or-exit *ui*)
      (when (gethash 'help options)
        (adopt:print-help-and-exit *ui*))
      (setf
        *interactive* nil
        *worker-threads* (gethash 'threads options)
        *k* (gethash 'k options)
        *w* (gethash 'w options)
        *foldback-position-epsilon* (gethash 'foldback-position-epsilon options)
        *intercept-epsilon* (gethash 'intercept-epsilon options)
        *gap-epsilon* (gethash 'gap-epsilon options)
        *minimum-cluster-length* (gethash 'minimum-cluster-length options)
        *minimum-foldback-length-absolute* (gethash 'minimum-foldback-length-absolute options)
        *minimum-foldback-length-relative* (gethash 'minimum-foldback-length-relative options)
        *plot-foldbacks* (gethash 'plot-foldbacks options)
        *plot-normal* (gethash 'plot-normal options)
        *output-directory* (gethash 'output-directory options))
      (handler-case
          (progn
            (assert *output-directory* () "Output directory must be specified.")
            (assert (= 1 (length arguments)) () "Exactly one .fastq file must be specified.")
            (assert (<= 1 *w* 31) () "Invalid window size ~A, must be in the range [1, 31]." *w*)
            (assert (<= 1 *k* 30) () "Invalid kmer size ~A, must be in the range [1, 30]." *k*)
            (assert (> *w* *k*) () "Window size (~A) must be larger than kmer size (~A)." *w* *k*)
            (assert (plusp *foldback-position-epsilon*) ()
              "Foldback position epsilon (~A) must be positive." *foldback-position-epsilon*)
            (assert (plusp *intercept-epsilon*) ()
              "Intercept epsilon (~A) must be positive." *intercept-epsilon*)
            (assert (plusp *gap-epsilon*) ()
              "Gap epsilon (~A) must be positive." *gap-epsilon*)
            (assert (plusp *minimum-cluster-length*) ()
              "Minimum cluster length (~A) must be positive." *minimum-cluster-length*)
            (assert (plusp *minimum-foldback-length-absolute*) ()
              "Minimum foldback length (absolute) (~A) must be positive." *minimum-foldback-length-absolute*)
            (assert (<= 0.0 *minimum-foldback-length-relative* 1.0) ()
              "Minimum foldback length (relative) (~A) must be in the range [0, 1]." *minimum-foldback-length-relative*)
            (run/fast (first arguments)))
        (error (e)
               (adopt:print-error-and-exit e))))))
