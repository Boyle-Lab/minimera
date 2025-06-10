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
         (suffix (reverse-complement (subseq prefix (- (length prefix) foldback-length)))))
    ;; (format t "~v@A~%~v@A" real-length prefix real-length (reverse suffix))
    (values (concatenate 'string prefix suffix)
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


;;;; Reference Implementation -------------------------------------------------



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
                         (collect (cons kmer (+ ws ks))))
                       #'string-minimizer-<))
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
    (sort (iterate
            (for (minimizer locs-1) :in-hashtable index-1)
            (for locs-2 = (gethash minimizer index-2))
            (appending (alexandria:map-product (lambda (l1 l2)
                                                 (make-hit (- l2 l1) l1 l2))
                                               locs-1 locs-2)))
          #'hit<)))

(defun foldback-score (seq &key (k 4) (w 12))
  (let* ((minimizers-id (minimizers k w seq))
         (minimizers-rc (minimizers k w (reverse-complement seq)))
         (hits (hits minimizers-id minimizers-rc)))
    (prl hits))
  0.0d0)

(defun run (filename)
  (with-open-file (f filename :direction :input)
    (gathering
      (do-fastq (lambda (id seq)
                      (identity (cons id (foldback-score seq))))
                    f))))


#; Scratch --------------------------------------------------------------------

(write-random-fastq "test.fastq" 100 200 400)

(multiple-value-bind (seq id l) (random-foldback 100000)
  (declare (ignore l))
  (let* ((seq "AGTCCCCAGTTAAACACCAGGTTTGGGAATGCTGCGGCAGGAGGATCGCTTGAGCCCATGAGTTTGAGACCAGCCAGAACAATATAGTGAGACCCTGTCTCTAAAAAAAAAAAAAAAAAAAAAAAAATTTAAACACTTAGCTGAGGCATGGTGGGGCATGCCTGTAGTCCCAGCTACACTGGGAGGCTGTGGTAGGAGGGTCGTTTGAGCTTGGAATATTGAGGCTGTAGTGAATGGTGATCAAGCCACTGCACTCCAGCCTGGGTAACAGAGGGAGACTCTGTCTCATAAATCAAACGTTTTGTATAGATTCCCATAGAAGTGAGTTAGGCATCAATCATAGAGTTATTAGCCACTTTGGTGTCTACCTTGGGAGTAAAACATATAATAAGGGGCAGCTTTAAACCATCTCAATCAATAACCTCCAACTTCTCGAGAAGGTTCTTATTTCATGAATTTCTAAACAAGAGACTACCTAGATTAAGACATTTGGTGGACACCATTTTGAGATGAAGAATCTGGAGTGGGAAGAAGGGAGATCTCTACTTACTGAAGCTTCCCAATGACATAGTTAAGTGTTCCCCAAAAGAAACTTTAGAACAAGAATTTCATCGTGCCATATCTCTATGAAAAGGAATTTCTCTAAAAGAAAACAAAGGCAAACAAGTGATAATCTGATTCTCATGGGAAAGTTTTCATCATAAAAGAAAAAGAGTGCTGGTTGCCGTGGCTCACGTCTGTAATCCCAACAGTTTGGGAGGCTGAGGTGGGTGGATTACCTGAGGTCAGAAGTTCAAAAACAGCCTGGCCAACATGGTGAAACCCTGTCTCTACTGAAAATACAAAAATTAGCCAGGTGTGGTGGTGTGCACCTGTAGTCCCAGCTACTTGGGAGGCCGAGGCAGGAGAATCACTTGAACCCAGGAGGTGGGAGTTGCAGTAAGCCGAGATGGTGCCACTGCACTCCAGGCTAGATGACACAGTGTGACTCCATCTCAAAAAAAAGAAAAAAAGAAAAACAAAAAAGGGACAAAGTATACTGGTCCAAAAAAGAAGAAGGCAAGAAAAAAAGGACAAAGTATATTGTTTACTATCATCACAGTGAGATAGTCCCCCTTTGAGATTAGAAAATAACAGTATACTCAAAGTAACATTAATGAGAACCAACATAAAATAGACAACATTCACTATCTACAAAAGTAATCTGCACCAATTAGCAATGTATGAGCATGTGGTTGAGAATATTGTCTATAATATGTATGTTAGAAGGAAGAGACCTCAAGAAAATGGTCAGAGCTGGAAATGTAGATTAGGGAATCTAGGTCAAAGTTTTGAGATTTTAGGAGTCCTGAGAGAATTTAAAAAGAGAAATAGCCACCAGGCGTGGTGGCCACACCTGTAATCCCAGCACTTTGGGAGGCCAAGGCAGGAAGATCATGAGGTCAGGAGTTCAAGACCAGTCTGGCCAACATAGTGAAACCCCGTCTCTACTAAAAATACAAAAAATTAGCTGGGCTTGGTAACACATACCTGTAATCCTAGCTACTTGGGAGGCTGAGGCAGGAGAATCGCTTGAACACACGAGGTGGAGGTTGCATTGAGCTGAGATCATGCCACTGCACTCCAGACTGGGTGACAGTGGGAGACTCCATCTCAAAACAAAAAACGAAAACAAAACAAAACAAAAAAACAGAAAAGGATAGGACTAAAGAACAGAGGTTGCTAAATTTAGAAATGAAGCAGGATCAGAGGAATAGGAAGGGATAGGGTAAAGAACAGATGTTGCTGCATTTAGAAAGGAAGGGGGTCAGAGGAATAGAAAGGGATAGGGCTGAAGAATAGAGGTCGCTGCATTTATAAATGAAGCCGAGTCAGGAATAGAAAGAGATAGGGCTGAAGAACAGAGGTTGCTGCATTTAGAAATGAAGTAGGGTCAGAGGAATAGAAAGGGATAGGGCTGAAGAACAGAGGTCACTGCATTTAGAAAGGAAGCGGGGTCAGAGGAGCAGAGGGAGCATTTGGTCACTGCTCTGCTGAGCAAAGCAGGATAAAGTCCTCCATGACCCTTGGACTTCTATATTGGAATTATTAAAAATCAGATTTCAGAATAAAAAACACAATAAGTGATGAAAAATGGATTTCTGAATGAGATCATGTGTCATAGAGTCCAATGGAAGGGGAGAAACAGGATATTAGAAAAGCCACAAAAAGTAGACAAAAGTTGTTTTTGTTGATTATAGA")
         (k 4)
         (w 12)
         (mut (mutate-dna seq)))
    (prl seq id)
    ;; (print-foldback seq l)
    (flet ((run (seq filename)
             (with-open-file (f filename :direction :output :if-exists :supersede)
               (conserve:write-row (list "c" "l1" "l2") f)
               (dolist (hit (hits (minimizers k w seq)
                                  (minimizers k w (reverse-complement seq))))
                 (conserve:write-row (list (princ-to-string (hit-c hit))
                                           (princ-to-string (hit-l1 hit))
                                           (princ-to-string (hit-l2 hit)))
                                     f)))

             ))
      (run seq "hits-seq.csv")
      (run mut "hits-mut.csv"))))



(run "test.fastq")
