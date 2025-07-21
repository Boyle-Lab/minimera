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

(defparameter *monotony-threshold* 0.80)

(defparameter *plot-foldbacks* nil
  "Whether to generate plots (with the accompanying R script) of reads determined to be foldbacks.")

(defparameter *plot-normal* nil
  "Whether to generate plots (with the accompanying R script) of reads determined *not* to be foldbacks.")


(defparameter *worker-threads* 6)

(defparameter *output-directory* "results")


(defparameter *report-progress* nil)

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


;;;; Alignments ---------------------------------------------------------------
(defparameter *alignments* nil)

(defun parse-alignment-file (path)
  (with-open-file (stream path)
    (conserve:read-row stream) ; header
    (iterate
      (for (id contig pos rev flag cigar) :in-stream stream :using #'conserve:read-row)
      (declare (ignorable flag))
      (collect-hash (id (list contig
                              (parse-integer pos)
                              (string= rev "1")
                              cigar))
                    :test #'equal))))

(defun parse-cigar (cigar-string)
  (gathering
    (ppcre:do-register-groups ((#'parse-integer n) op)
        ("([0-9]+)?([^0-9])" cigar-string)
      (gather (cons (char op 0) (or n 1))))))

(defun alignment-length (cigar)
  (loop :for (op . n) :in cigar
        :summing (ecase op
                   ((#\S #\H #\I) 0)
                   ((#\M #\= #\X #\D) n))))

(defun compute-reference-offset (cigar foldback-point)
  (let ((offset 0))
    (destructuring-bind (op . n) (first cigar)
      (case op
        ;; Account for the weird start position in the case of soft clip.
        (#\S (decf offset n))
        ;; If it's hard clipped we wouldn't have included that bit in the read
        ;; sequence we searched for minimizers, I think?
        (#\H (progn (pop cigar)))))
    (iterate
      (with f = foldback-point)
      (for (op . n) :in cigar)
      (while (plusp f))
      (ecase op
        ((#\D) (incf offset n))
        ((#\M #\= #\X #\S) (progn (incf offset (min f n))
                                  (setf f (max 0 (- f n)))))
        ((#\I) (setf f (max 0 (- f n))))))
    offset))

(defun compute-reference-position (pos cigar-string foldback-point is-rev)
  (let ((cigar (parse-cigar cigar-string)))
    (if is-rev
        (let* ((pos (+ pos (alignment-length cigar)))
               (cigar (reverse cigar)))
          (- pos (compute-reference-offset cigar foldback-point)))
      (+ pos (compute-reference-offset cigar foldback-point)))))


;;;; Main Entry Point ---------------------------------------------------------
(defun compute-foldback-position (cluster)
  (multiple-value-bind (start end) (cluster-bounds cluster)
    (truncate (+ start end) 2)))


(defun find-foldback (id sequence &key (optimized nil))
  (if optimized
    (check-type sequence a8)
    (check-type sequence string))
  (let* ((read-length (length sequence))
         (m1 (if optimized
               (minimizers/fast *k* *w* sequence)
               (minimizers *k* *w* sequence)))
         (monotony (estimate-monotony m1)))
    (if (>= monotony *monotony-threshold*)
      (list id read-length :monotonous monotony nil)
      (let* ((m2 (if optimized
                   (minimizers/fast *k* *w* (nreverse-complement-bytes sequence))
                   (minimizers *k* *w* (reverse-complement sequence))))
             (hits (hits m1 m2))
             (clusters (cluster hits))
             (results (find-foldback-clusters read-length clusters))
             (result (alexandria:extremum results #'> :key #'cluster-lengths)))
        (if results
          (when *plot-foldbacks* (plot-minimizers id sequence hits clusters))
          (when *plot-normal*    (plot-minimizers id sequence hits clusters)))
        (if result
          (list id read-length :foldback monotony (compute-foldback-position result))
          (list id read-length :normal monotony nil))))))


;;;; Toplevel -----------------------------------------------------------------
(defparameter *csv-headers* (list "read-id" "read-length" "classification" "monotony" "foldback-point"))

(defun write-csv-result (result stream)
  (destructuring-bind (id read-length class monotony foldback-point) result
    (conserve:write-row (list id
                              (princ-to-string read-length)
                              (string-downcase class)
                              (format nil "~,4F" monotony)
                              (if foldback-point
                                (princ-to-string foldback-point)
                                ""))
                        stream)))

(defun write-bed-result (result stream)
  (destructuring-bind (id read-length class monotony foldback-point) result
    (declare (ignore read-length class monotony))
    (when foldback-point
      (when-let ((alignment-info (gethash id *alignments*)))
        (destructuring-bind (contig pos is-rev cigar) alignment-info
          (let ((conserve:*delimiter* #\tab)
                (reference-position (compute-reference-position pos cigar foldback-point is-rev)))
            (conserve:write-row (mapcar 'princ-to-string
                                        (list contig
                                              (max 0 (1- reference-position))
                                              (1+ reference-position)
                                              id))
                                stream)))))))

(defun run/slow (filename)
  (with-open-file (input-stream filename :direction :input)
    (with-open-file (output-stream (format nil "~A.csv" *output-directory*) :direction :output :if-exists :supersede)
      (conserve:write-row *csv-headers* output-stream)
      (do-fastq (lambda (id seq)
                  (write-csv-result (find-foldback id seq :optimized nil) output-stream))
                input-stream))))


(defparameter *interactive* t)
(defparameter *input-done* nil)
(defparameter *work-done*  nil)

(defparameter *input-queue*  (lparallel.queue:make-queue :fixed-capacity (* 2 *worker-threads*)))
(defparameter *output-queue* (lparallel.queue:make-queue :fixed-capacity (* 2 *worker-threads*)))
(defparameter *progress-queue* (lparallel.queue:make-queue))


(defun parse-read-id (bytes)
  ;; TODO: Remove or document the comma splitting here.
  (first (str:split #\tab (first (str:split "," (map 'string #'code-char bytes)
                                          :limit 2))
                    :limit 2)))


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


(defun run/reader% (fastq-filename)
  (with-exit-during-noninteractive-mode
    (flet ((run% (stream)
             (map-fastq (lambda (fastq-read)
                          (lparallel.queue:push-queue fastq-read *input-queue*))
                        stream)))
      (if (string= fastq-filename "-")
        (run% *standard-input*)
        (with-open-file (input-stream fastq-filename :direction :input :element-type 'u8)
          (run% input-stream))))
    (setf *input-done* t)))

(defun run/writer% ()
  (let ((bed-stream (when *alignments*
                      (open (format nil "~A/foldbacks.bed" *output-directory*)
                            :direction :output
                            :if-exists :supersede))))
    (unwind-protect
        (with-open-file (csv-stream (format nil "~A/foldbacks.csv" *output-directory*)
                                    :direction :output
                                    :if-exists :supersede)
          (conserve:write-row *csv-headers* csv-stream)
          (loop (multiple-value-bind (result found)
                    (lparallel.queue:try-pop-queue *output-queue* :timeout 1)
                  (cond ((and (not found) *work-done*) (return))
                        (found (progn
                                 (write-csv-result result csv-stream)
                                 (when bed-stream
                                   (write-bed-result result bed-stream))))
                        (t (progn))))))
      (when bed-stream
        (close bed-stream)))))

(defun run/worker% ()
  (with-exit-during-noninteractive-mode
    (loop (multiple-value-bind (fastq-read found)
              (lparallel.queue:try-pop-queue *input-queue* :timeout 1)
            (cond ((and (not found) *input-done*) (return))
                  (found (let ((id (parse-read-id (id fastq-read))))
                           (lparallel.queue:push-queue id *progress-queue*)
                           (lparallel.queue:push-queue
                             (find-foldback id (seq fastq-read) :optimized t)
                             *output-queue*)))
                  (t (progn)))))
    (setf *work-done* t)))

(defun run/progress% ()
  (with-exit-during-noninteractive-mode
    (loop (multiple-value-bind (update found)
              (lparallel.queue:try-pop-queue *progress-queue* :timeout 1)
            (cond ((and (not found) *work-done*) (return))
                  (found (when *report-progress*
                           (format *debug-io* "Processing: ~A~%" update)))
                  (t (progn)))))))

(defun run/fast (filename &key alignments)
  (setf *input-done* nil *work-done* nil)
  (ensure-directories-exist (uiop:ensure-directory-pathname *output-directory*))
  (when alignments
    (setf *alignments* (parse-alignment-file alignments)))
  ;; Spawn single reader thread.
  (bt2:make-thread (lambda () (run/reader% filename)) :name "Minimera FASTQ Reader")
  ;; Spawn N worker threads.
  (dotimes (i *worker-threads*)
    (bt2:make-thread (lambda () (run/worker%)) :name (format nil "Minimera Worker ~D" (1+ i))))
  ;; Spawn single progress thread.
  (bt2:make-thread (lambda () (run/progress%)) :name "Minimera Progress Reporter")
  ;; Spawn a single writer thread, and join it to wait for things to finish.
  (bt2:join-thread
    (bt2:make-thread (lambda () (run/writer%)) :name "Minimera Output Writer")))

#; Scratch --------------------------------------------------------------------

(defparameter *r* nil)

(with-open-file (f "data/hallucination.fastq" :element-type 'u8)
  (map-fastq (lambda (x) (setf *r* x)) f))

(doseq (m (minimizers/fast 8 15 (seq *r*)))
  (print m))

(_ (seq *r*)
  (minimizers/fast 8 15 _)
  (mapcar #'hash _)
  proportions
  alexandria:hash-table-alist
  (sort _ #'> :key #'cdr)
  (take 10 _)
  )

(_ (seq *r*)
  (minimizers/fast 8 15 _)
  (mapcar #'hash _)
  frequencies
  alexandria:hash-table-alist
  (sort _ #'> :key #'cdr)
  (take 10 _)
  )

(_ (seq *r*)
  (minimizers/fast 8 15 _)
  length
  (truncate _ 1000)
  )

