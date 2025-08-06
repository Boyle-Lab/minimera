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


(defparameter *report-progress* t)

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
    ((< length    500)    50)
    ((< length   5000)   100)
    ((< length  10000)   500)
    ((< length  50000)  1000)
    ((< length 100000)  5000)
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

(defparameter *png-size*    512)

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

(defun render-minimizers (id sequence classification foldback-point hits clusters)
  (assert (>= *png-size* 64))
  (let* ((size *png-size*)
         (png (make-instance 'zpng:png
                :color-type :truecolor
                :width size
                :height size))
         (img (zpng:data-array png))
         (len (length sequence))
         (idx (index-hits clusters))
         (pad (truncate size 20))
         (pp 1)
         (s pad)
         (e (- size pad))
         (pixels-per-base (/ (- e s) len))
         (tic-bases (choose-tics sequence))
         (tic-width (truncate (* tic-bases pixels-per-base)))
         (tic-length (truncate pad 3)))
    (dotimes (i (array-total-size img))
      (setf (row-major-aref img i) #xFF))
    (labels ((base->x (base) (values (truncate (map-range 0 len s e base))))
             (base->y (base) (values (truncate (map-range 0 len e s base))))
             (hit-x (hit) (base->x (hit-l1 hit)))
             (hit-y (hit) (base->y (hit-l2 hit)))
             (ink (x y &optional (a #x30))
               (dotimes (ch 3)
                 (setf (aref img y x ch) (max 0 (- (aref img y x ch) a)))))
             (draw (x y r &optional (g r) (b g))
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
             (draw-char (x y char)
               (let ((pixels (gethash char *m5x7*)))
                 (assert pixels () "No glyph available for ~S" char)
                 (destructuring-bind (width height) (array-dimensions pixels)
                   (dotimes (dx width)
                     (dotimes (dy height)
                       (when (plusp (aref pixels dx dy))
                         (draw (+ x dx) (+ y dy) #x00))))
                   width)))
             (draw-string (x y string)
               (iterate
                 (for char :in-string string)
                 (for w = (draw-char x y char))
                 (incf x (1+ w)))))
      (draw-string 5 5 (format nil "Read ~A (~A)" id (string-downcase classification)))
      (draw-string 5 (- size 10)
                   (format nil "tics = ~:Dbp, read length = ~:D, ~A"
                           tic-bases len
                           (if foldback-point
                             (format nil "foldback point = ~:D" foldback-point)
                             "no foldback detected")))
      (do-range ((x (1- s) e)) (draw x e #x66)) ; x axis
      (do-range ((y s e)) (draw (1- s) y #x66)) ; y axis
      (iterate (for tx :from (+ s tic-width) :to e :by tic-width) ; xticx
               (loop :for ty :from (1+ e)
                     :repeat tic-length
                     :do (draw tx ty #x66)))
      (iterate (for ty :from (- e tic-width) :downto s :by tic-width) ; ytics
               (loop :for tx :downfrom (- s 2)
                     :repeat tic-length
                     :do (draw tx ty #x66)))
      (when foldback-point ; draw dashed line at foldback point
        (iterate (with x = (base->x foldback-point))
                 (for y :from s :to e :by 4)
                 (draw x y      #xBF #xB0 #x86)
                 (draw x (1+ y) #xBF #xB0 #x86)))
      (doseq (hit hits) ; unclustered hits
        (let ((cluster (gethash hit idx)))
          (unless cluster
            (multiple-value-call #'draw-unclustered-point
              (hit-y hit) (hit-x hit)))))
      (doseq (hit hits) ; clustered hits
        (let ((cluster (gethash hit idx)))
          (when cluster
            (multiple-value-call #'draw-point
              (hit-y hit) (hit-x hit) (cluster-color cluster))))))
    (ensure-directories-exist (format nil "~A/plots/" *output-directory*))
    (zpng:write-png png (format nil "~A/plots/~A.png" *output-directory* id)
                    :if-exists :supersede)))


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
             (result (alexandria:extremum results #'> :key #'cluster-lengths))
             (classification (if result :foldback :normal))
             (foldback-position (when result
                                  (compute-foldback-position result))))
        (when (or (and results *plot-foldbacks*)
                  (and (null results) *plot-normal*))
          (render-minimizers id sequence classification foldback-position hits clusters))
        (list id read-length classification monotony foldback-position)))))


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

(defparameter *rs* nil)

(setf *output-directory* "results/repl")

(with-open-file (f "data/bench/bulk.fastq" :element-type 'u8)
  (setf *rs* (gathering-vector () (map-fastq #'gather f))))

(mapcar #'car (last (take 200 (sort (enumerate *rs*) #'> :key (lambda (x) (length (seq (cdr x)))))) 10))

(defparameter *r* (aref *rs* 709)) ; the longboi
(defparameter *r* (aref *rs* 124))
(pr (length (seq *r*)))

;; (map 'string 'code-char (seq *r*))

(progn
  (defparameter *min/id* (minimizers/fast 8 15 (seq *r*)))
  (defparameter *min/rc* (minimizers/fast 8 15 (reverse-complement-bytes (seq *r*))))
  (defparameter *hits* (hits *min/id* *min/rc*))
  (defparameter *clusters* (cluster *hits*)))


(_ (index-hits *clusters*)
  alexandria:hash-table-alist
  (mapcar #'cdr _)
  (frequencies _)
  )

(time (plot-minimizers "repltest-r" (seq *r*) *hits* *clusters*))
(time (render-minimizers "repltest-lisp" (seq *r*) :unknown 5000 *hits* *clusters*))

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

