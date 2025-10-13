(in-package :minimera)

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
    (for qs = (read-line stream))
    (funcall function id seq qs)))


;;;; FASTQ (Optimized) --------------------------------------------------------
(deftype buffer-length ()
  `(integer 1 ,array-dimension-limit))

(deftype u8 ()
  '(unsigned-byte 8))

(deftype a8 ()
  '(simple-array u8 (*)))

(defstruct (fastq-read (:conc-name nil))
  (id nil :type (simple-array u8 (*)))
  (seq nil :type (simple-array u8 (*)))
  (qs nil :type (simple-array u8 (*))))

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
        (for quality-scores = (read-line%))
        (map-into quality-scores (lambda (q) (- q (char-code #\!))) quality-scores)
        (unless (= (length quality-scores) (length seq))
          (error "Mismatched number of quality scores (~D) for read of length ~D."
                 (length quality-scores) (length seq)))
        (funcall function (make-fastq-read :id id :seq seq :qs quality-scores)))
      ;; (print refills)
      (values))))

