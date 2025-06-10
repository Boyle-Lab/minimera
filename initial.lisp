(eval-when (:compile-toplevel :load-toplevel :execute)
  (ql:quickload '(:losh :uuid :str)))

(defpackage :foldbacks
  (:use :cl :losh :iterate))

(in-package :foldbacks)


;;;; Data Generation ----------------------------------------------------------
(defun random-dna (n)
  (iterate (repeat n)
    (collect (random-elt "ACTG") :result-type 'string)))

(defun random-quality (n)
  (iterate (repeat n)
           (collect (random-elt "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQ") :result-type 'string)))

(defun random-fastq-entry (n)
  (format nil ">~A~%~A~%+~%~A~%"
          (string-downcase (princ-to-string (uuid:make-v4-uuid)))
          (random-dna n)
          (random-quality n)))

(defun write-random-fastq (filename entries min-length max-length)
  (with-open-file (f filename :direction :output :if-exists :supersede)
    (do-repeat entries
      (write-string (random-fastq-entry
                      (random-range-inclusive min-length max-length))
                    f))))


;;;; Reference Implementation -------------------------------------------------
(defun ref/map-fastq (function stream)
  (iterate
    (until (eql :eof (peek-char nil stream nil :eof)))
    (assert (char= #\> (read-char stream)))
    (for id = (read-line stream))
    (for seq = (read-line stream))
    (assert (char= #\+ (read-char stream)))
    (peek-char #\newline stream) (read-char stream) ; +
    (peek-char #\newline stream) (read-char stream) ; quality scores
    (funcall function id seq)))

(defun ref/gc-content (id seq)
  (_ seq
    (count-if (lambda (ch) (or (char= #\G ch) (char= #\C ch))) _)
    (coerce _ 'double-float)
    (let ((n (length seq)))
      (if (zerop n)
        0.0d0
        (/ _ n)))
    (cons id _)))

(defun ref/run (filename)
  (with-open-file (f filename :direction :input)
    (gathering
      (ref/map-fastq (lambda (id seq)
                       (gather (ref/gc-content id seq)))
                     f))))

(define-sorting-predicate string-minimizer-<
  (#'string< :key #'car)
  (#'< :key #'cdr))

(defun ref/minimizers (k w id seq)
  (declare (ignore id))
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

(ref/minimizers 4 10 "foo" "CCTAGCGGAACTGAGTCCTCTGGCTGTAC")

(defun ref)


;;;; Reference Implementation -------------------------------------------------

;; uint64_t hash(uint64_t key) {
;;   key = (~key) + (key << 21); // key = (key << 21) - key - 1;
;;   key = key ^ (key >> 24);
;;   key = (key + (key << 3)) + (key << 8); // key * 265
;;   key = key ^ (key >> 14);
;;   key = (key + (key << 2)) + (key << 4); // key * 21
;;   key = key ^ (key >> 28);
;;   key = key + (key << 31);
;;   return key;
;; }

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


(hex (minhash "GACTTATAG" :start 1 :end 8) 64)


;;;; FASTQ Implementation -----------------------------------------------------

(defconstant +buffer-size+ (* 10 (expt 2 20)))

;; (defconstant +buffer-size+ 4096)

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

(defun map-fastq (function stream)
  (declare (optimize (speed 3) (safety 1) (debug 1)))
  (let ((buffer (make-array +buffer-size+ :element-type 'u8))
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


;;;; Optimized Implementation -------------------------------------------------
(defun gc-content (read)
  (declare (optimize (speed 3) (safety 1) (debug 1)))
  (_ (loop :for b :across (seq read)
           :counting (or (= (char-code #\G) b)
                         (= (char-code #\C) b)))
    (coerce _ 'double-float)
    (let ((n (length (seq read))))
      (if (zerop n)
        0.0d0
        (/ _ n)))
    (cons (bs (id read)) _)))

(defun run (filename)
  (with-open-file (f filename :direction :input :element-type 'u8)
    (gathering
      (map-fastq (lambda (read)
                   (gather (gc-content read)))
                 f))))


#; Scratch --------------------------------------------------------------------

(defun check (filename)
  (let ((result-ref (time (ref/run filename)))
        (result-opt (time (run filename))))
    (assert (equalp result-ref
                    result-opt))))

(check "scratch.fastq")

(check "stress.fastq")

(check "small.fastq")

(write-random-fastq "small.fastq" 10 150 150)

(write-random-fastq "big.fastq" 15000 200 80000)

(profile (progn (ref/run "big.fastq") (values)))

(profile (progn (run "big.fastq") (values)))


(check "big.fastq")
