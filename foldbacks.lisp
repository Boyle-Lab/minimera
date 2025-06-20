(eval-when (:compile-toplevel :load-toplevel :execute)
  (ql:quickload '(:losh :uuid :str :conserve)))

(defpackage :foldbacks
  (:use :cl :losh :iterate))

(in-package :foldbacks)

;;;; Data Generation ----------------------------------------------------------
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

(defun hash64 (n)
  (etypecase n
    (u64
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
  ;;; https://naml.us/post/inverse-of-a-hash-function/
  (etypecase n
    (u64
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

(defun phihash (string &key (start 0) (end (length string)))
  "Return the φ-hash of the subsequence of `string` bounded by `start` and `end`."
  (iterate
    (with hash = 0)
    (for ch :in-string string :from start :below end)
    (zapf hash (_ %
                 (ash _ 2)
                 (logior _ (ecase ch
                             (#\A #b00)
                             (#\C #b01)
                             (#\T #b10)
                             (#\G #b11)))))
    (returning (hash64 hash))))


;;;; FASTQ --------------------------------------------------------------------
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


;;;; Minimizers ---------------------------------------------------------------
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
  (alexandria:extremum
    (iterate
      (for (kmer ks ke) :window k :on window)
      (collect (cons (phihash kmer) (+ window-start ks))))
    #'< :key #'car))

(defun minimizers (k w sequence)
  "Return a list of the (`k`, `w`, φ)-mimimizer sketch of `sequence`."
  (iterate
    (with prev)
    (for (window ws we) :window w :on sequence)
    (for minimizer = (minimizer k window ws))
    (unless (and prev (= (cdr prev) (cdr minimizer)))
      (collect minimizer)
      (setf prev minimizer))))


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
  (let ((index-1 (group-by #'car ms1 :test #'equal :map #'cdr))
        (index-2 (group-by #'car ms2 :test #'equal :map #'cdr)))
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
(defun foldback-score (sequence &key
                       k w
                       intercept-epsilon
                       minimum-cluster-length
                       gap-epsilon
                       foldback-position-epsilon
                       foldback-length-absolute
                       foldback-length-relative
                       info
                       (plot nil))
  (let ((*intercept-epsilon* intercept-epsilon)
        (*minimum-cluster-length* minimum-cluster-length)
        (*gap-epsilon* gap-epsilon)
        (*foldback-position-epsilon* foldback-position-epsilon)
        (*minimum-foldback-length-absolute* foldback-length-absolute)
        (*minimum-foldback-length-relative* foldback-length-relative))
    (let* ((read-length (length sequence))
           (hits (hits (minimizers k w sequence)
                       (minimizers k w (reverse-complement sequence))))
           (clusters (cluster hits))
           (results (find-foldback-clusters read-length clusters)))
      (format *debug-io* "~D minimizer hit~:P found (~F/kbp)~%"
              (length hits)
              (/ (length hits)
                 (/ (length sequence) 1000)))
      (format *debug-io* "~D cluster~:P found~%"
              (length clusters))
      (format *debug-io* "Found ~D candidate foldback cluster~:P: ~A~%"
              (length results)
              results)
      (when plot
        (plot-minimizers info sequence hits clusters)
        ;; (break)
        )
      results)))


;;;; Toplevel -----------------------------------------------------------------
(defun run (filename &rest args)
  (with-open-file (f filename :direction :input)
    (gathering
      (do-fastq (lambda (id seq)
                  (format *debug-io* "~4%Checking record ~A~%" id)
                  (let ((foldback-regions (apply #'foldback-score
                                                 seq :info id
                                                 args)))
                    (gather (list id (length foldback-regions) foldback-regions))))
                f))))



#; Scratch --------------------------------------------------------------------

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

(foldback-score
  (alexandria:read-file-into-string "data/foldback2")
  :k 8 :w 16
  :intercept-epsilon 30
  :minimum-cluster-length 10
  :gap-epsilon 300
  :foldback-position-epsilon 100
  :foldback-length-absolute 50
  :foldback-length-relative 0.05
  :plot t)

;; (defun plot-minimizers (sequence &key
;;                         k w
;;                         (intercept-epsilon 10)
;;                         (cluster-min 10)
;;                         (gap-epsilon 10))
;;   (with-open-file (f "hits.csv" :direction :output :if-exists :supersede)
;;     (conserve:write-row (list "c" "l1" "l2" "cluster") f)
;;     (let* ((hits (hits (minimizers k w sequence)
;;                        (minimizers k w (reverse-complement sequence))))
;;            (clusters (index-hits (cluster hits
;;                                           :intercept-epsilon intercept-epsilon
;;                                           :min-length cluster-min
;;                                           :gap-epsilon gap-epsilon))))
;;       (doseq (hit hits)
;;         (conserve:write-row (list (princ-to-string (hit-c hit))
;;                                   (princ-to-string (hit-l1 hit))
;;                                   (princ-to-string (hit-l2 hit))
;;                                   (princ-to-string (gethash hit clusters "NA")))
;;                             f))))
;;   (write-line "Plotting…")
;;   (rscript-file "plot.R" (length sequence) (choose-tics sequence)))


(time
  (let ((*debug-io*
          ;; (make-broadcast-stream)
          *debug-io*
          ))
    (remove-if (lambda (result)
                 (null (third result)))
               (run "data/9b58328c3b631816942cbe400a807241935897e6.manual-test.2.fastq"
                    :k 8 :w 16
                    :intercept-epsilon 30
                    :minimum-cluster-length 10
                    :gap-epsilon 300
                    :foldback-position-epsilon 100
                    :foldback-length-absolute 50
                    :foldback-length-relative 0.05
                    :plot t))))

