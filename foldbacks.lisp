(eval-when (:compile-toplevel :load-toplevel :execute)
  (ql:quickload '(:losh :uuid :str :conserve)))

(defpackage :foldbacks
  (:use :cl :losh :iterate))

(in-package :foldbacks)

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

(defun longest-increasing-subsequence (sequence &key (key #'identity))
  (let ((piles (make-array 256 :adjustable t :fill-pointer 0)))
    (map nil (lambda (item) (patience-insert item piles :key key)) sequence)
    (nreverse (car (aref piles (1- (length piles)))))))


;;;; Reference Implementation -------------------------------------------------
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

(defun invertible-hash (n)
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

;; (format t "~16,'0X" (invertible-hash #x7ffffbffffdfffff))


(defun minhash (string &key (start 0) (end (length string)))
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
    (returning (invertible-hash hash))))


(defun do-fastq (function stream)
  (iterate
    (until (eql :eof (peek-char nil stream nil :eof)))
    (assert (char= #\> (read-char stream)))
    (for id = (read-line stream))
    (for seq = (read-line stream))
    (assert (char= #\+ (read-char stream)))
    (peek-char #\newline stream) (read-char stream) ; +
    (peek-char #\newline stream) (read-char stream) ; quality scores
    (funcall function id seq)))

(defun reverse-complement (seq)
  (nreverse (map 'string (lambda (ch)
                           (ecase ch
                             (#\A #\T)
                             (#\C #\G)
                             (#\G #\C)
                             (#\T #\A)))
                 seq)))

(define-sorting-predicate string-minimizer-<
  (#'string< :key #'car)
  (#'< :key #'cdr))

(defun minimizers (k w seq)
  (iterate
    (with prev)
    (for (window ws we) :window w :on seq)
    (for minimizer = (alexandria:extremum
                       (iterate
                         (for (kmer ks ke) :window k :on window)
                         (collect
                           (cons (minhash kmer) (+ ws ks))
                           ;; (cons kmer (+ ws ks))

                           ))
                       #'< :key #'car
                       ;; #'string-minimizer-<
                       ))
    (unless (and prev (= (cdr prev) (cdr minimizer)))
      (collect minimizer)
      (setf prev minimizer))))

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

(define-sorting-predicate cluster-hit-<
  (#'< :key #'hit-l1)
  (#'< :key #'hit-l2))

(define-sorting-predicate cluster-hit->
  (#'> :key (lambda (hit) (+ (* (hit-l1 hit) (hit-l1 hit))
                             (* (hit-l2 hit) (hit-l2 hit)))))
  (#'> :key #'hit-l1)
  (#'> :key #'hit-l2))

(defun cluster (hits &key (epsilon 10) (min-length 10))
  (unless (zerop (length hits))
    (iterate
      (with last-i = (1- (length hits)))
      (with start = 0)
      (for hit :in-vector hits :with-index i)
      (for next-i = (1+ i))
      (for next-hit = (if (= i last-i)
                          nil
                          (aref hits next-i)))
      (when (or (null next-hit)
                (> (- (hit-c next-hit) (hit-c hit))
                   epsilon))
        (when (>= (- next-i start) min-length)
          (collect (subseq hits start next-i) :into results))
        (setf start next-i))
      (finally
        (prl (mapcar #'length results))
        (let ((deduped (mapcar (lambda (cluster)
                                 (longest-increasing-subsequence (sort cluster #'< :key #'hit-l1) :key #'hit-l2))
                               results)))
          (prl (mapcar #'length deduped))
          (return deduped))))))

(defun foldback-score (seq &key (k 4) (w 12))
  (let* ((minimizers-id (minimizers k w seq))
         (minimizers-rc (minimizers k w (reverse-complement seq)))
         (hits (hits minimizers-id minimizers-rc)))
    (prl hits))
  0.0d0)

;; (defun run (filename)
;;   (with-open-file (f filename :direction :input)
;;     (gathering
;;       (do-fastq (lambda (id seq)
;;                       (identity (cons id (foldback-score seq))))
;;                     f))))



(defparameter *data/okay-1*     (alexandria:read-file-into-string "data/okay1"))
(defparameter *data/foldback-1* (alexandria:read-file-into-string "data/foldback1"))
(defparameter *data/foldback-2* (alexandria:read-file-into-string "data/foldback2"))
(defparameter *data/foldback-3* (alexandria:read-file-into-string "data/foldback2"))

;; (write-random-fastq "test.fastq" 100 200 400)

(defun index-hits (clusters)
  (iterate
    (for i :from 1)
    (for hits :in clusters)
    (format t "cluster ~D: ~D hits~%" i (length hits))
    (doseq (hit hits)
      (collect-hash (hit i)))))

(defun run-test (sequence &key k w (mut? t) (cluster-epsilon 10) (cluster-min 10))
  (let* ((seq (if mut? (mutate-dna sequence) sequence))
         (rlen (length seq)))
    (with-open-file (f "hits.csv" :direction :output :if-exists :supersede)
      (conserve:write-row (list "c" "l1" "l2" "cluster") f)
      (let* ((hits (hits (minimizers k w seq)
                         (minimizers k w (reverse-complement seq))))
             (clusters (index-hits (cluster hits :epsilon cluster-epsilon :min-length cluster-min))))
        (doseq (hit hits)
          (conserve:write-row (list (princ-to-string (hit-c hit))
                                    (princ-to-string (hit-l1 hit))
                                    (princ-to-string (hit-l2 hit))
                                    (princ-to-string (gethash hit clusters "NA")))
                              f))))
    (write-line "Plotting…")
    (rscript-file "plot.R" rlen)))


#; Scratch --------------------------------------------------------------------

(defparameter *seq* (mutate-dna (random-foldback 20000)))

(defparameter *seq* (random-foldback 20000))

(defparameter *seq* (random-dna 20000))

(run-test *seq* :k 12 :w 15 :mut? nil :cluster-epsilon 30 :cluster-min 10)


(run-test (alexandria:read-file-into-string "data/foldback4")
          :k 8 :w 20
          :mut? nil
          :cluster-epsilon 30 :cluster-min 20)

(run-test (alexandria:read-file-into-string "data/ligation-with-rep") :k 12 :w 25 :mut? nil)

(run-test (pbpaste) :k 12 :w 20 :mut? nil)
