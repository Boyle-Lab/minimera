(in-package :minimera)

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

