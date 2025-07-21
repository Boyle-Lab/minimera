(in-package :minimera)

;;;; SAM ----------------------------------------------------------------------

;; https://samtools.github.io/hts-specs/SAMv1.pdf

(defun assert-char (expected-char stream)
  (let ((ch (peek-char nil stream)))
    (assert (char= expected-char ) () "Expected ~S but got ~S" expected-char ch)))

(defun parse-tag (tag)
  (destructuring-bind (tag data) (str:split #\: tag :limit 2)
    (cons tag data)))

(defun read-sam-header-line (stream)
  (assert-char #\@ stream)
  (let ((line (read-line stream)))
    (destructuring-bind (record-type data) (str:split #\tab line :limit 2)
      (if (string= record-type "@CO")
        (cons :comment (subseq line 4))
        (let ((data (str:split #\tab data)))
          (cons (alexandria:switch (record-type :test #'string=)
                  ("@HD" :file-metadata)
                  ("@SQ" :reference-sequence)
                  ("@RG" :read-group)
                  ("@PG" :program)
                  (t (error "Malformed header line ~S" line)))
                (mapcar #'parse-tag data)))))))

(defun skip-sam-headers (stream)
  (loop :while (char= #\@ (peek-char nil stream nil :eof))
        :do (read-line stream)))

(defstruct (alignment)
  qname
  flag
  rname
  pos
  mapq
  cigar
  rnext
  pnext
  tlen
  seq
  qual
  optional-fields)

(defmethod print-object ((o alignment) s)
  (print-unreadable-object (o s :type t)
    (format s "~A ~12,'0B ~5A Q~2,'0D ~A ~A:~D +~D field~:P"
            (alignment-qname o)
            (alignment-flag o)
            (cond ((unmappedp o) "UNMAP")
                  ((supplementaryp o) "SUPPL")
                  ((secondaryp o) "SEC")
                  (t "PRIM"))
            (alignment-mapq o)
            (str:shorten 12 (alignment-seq o))
            (alignment-rname o)
            (alignment-pos o)
            (length (alignment-optional-fields o)))))

(defun multi-segment-p             (alignment) (logbitp  0 (alignment-flag alignment)))
(defun all-aligned-p               (alignment) (logbitp  1 (alignment-flag alignment)))
(defun unmappedp                   (alignment) (logbitp  2 (alignment-flag alignment)))
(defun next-unmapped-p             (alignment) (logbitp  3 (alignment-flag alignment)))
(defun reverse-complemented-p      (alignment) (logbitp  4 (alignment-flag alignment)))
(defun next-reverse-complemented-p (alignment) (logbitp  5 (alignment-flag alignment)))
(defun first-segment-p             (alignment) (logbitp  6 (alignment-flag alignment)))
(defun last-segment-p              (alignment) (logbitp  7 (alignment-flag alignment)))
(defun secondaryp                  (alignment) (logbitp  8 (alignment-flag alignment)))
(defun not-passing-p               (alignment) (logbitp  9 (alignment-flag alignment)))
(defun duplicatep                  (alignment) (logbitp 10 (alignment-flag alignment)))
(defun supplementaryp              (alignment) (logbitp 11 (alignment-flag alignment)))

(defun or-* (form)
  (if (string= "*" form)
    nil
    form))

(defmacro parse-and-check (value type &optional (parser 'parse-integer))
  (with-gensyms (result)
    `(let ((,result (,parser ,value)))
       (assert (typep ,result ',type))
       ,result)))

(defun parse-u8  (string) (parse-and-check string (unsigned-byte  8)))
(defun parse-i8  (string) (parse-and-check string (signed-byte    8)))
(defun parse-u16 (string) (parse-and-check string (unsigned-byte 16)))
(defun parse-i16 (string) (parse-and-check string (signed-byte   16)))
(defun parse-u32 (string) (parse-and-check string (unsigned-byte 32)))
(defun parse-i32 (string) (parse-and-check string (signed-byte   32)))
(defun parse-fl  (string) (parse-and-check string float parse-float:parse-float))

(defun parse-numeric-array (payload)
  (destructuring-bind (type &rest data) (str:split #\, payload)
    (map 'vector
         (alexandria:eswitch (type :test #'string=)
           ("c" #'parse-i8)
           ("C" #'parse-u8)
           ("s" #'parse-i16)
           ("S" #'parse-u16)
           ("i" #'parse-i32)
           ("I" #'parse-u32)
           ("f" #'parse-fl))
         data)))

(defun parse-field-value (type value)
  (assert (= 1 (length type)) () "Malformed optional field type: ~S" type)
  (ecase (char type 0)
    (#\A (progn
           (assert (= 1 (length type)) () "Malformed character field: ~S" value)
           (cons :char (char value 0))))
    (#\i (cons :int (parse-integer value)))
    (#\f (cons :float (parse-float:parse-float value)))
    (#\Z (cons :str value))
    (#\H (progn
           (assert (evenp (length value)) () "Malformed hex array: ~S" value)
           (cons :hex (iterate
                        (for i :from 0 :below (length value) :by 2)
                        (collect (parse-integer (subseq value i (+ i 2)) :radix 16)
                                 :result-type vector)))))
    (#\B (cons :numeric (parse-numeric-array value)))))

(defun parse-optional-field (field)
  (destructuring-bind (tag type value)
      (str:split #\: field :limit 3)
    (cons tag (parse-field-value type value))))

(defun read-sam-record (stream &optional eof-error-p eof-value)
  (with-eof-handled (stream eof-error-p eof-value)
    (let* ((line (read-line stream))
           (fields (str:split #\tab line)))
      (assert (>= (length fields) 11) () "Malformed SAM line: ~S" line)
      (destructuring-bind
          (qname flag rname pos mapq cigar rnext pnext tlen seq qual &rest optional-fields)
          fields
        (make-alignment
          :qname (or-* qname)
          :flag (parse-integer flag)
          :rname (or-* rname)
          :pos (parse-integer pos)
          :mapq (parse-integer mapq)
          :cigar (or-* cigar)
          :rnext (or-* rnext)
          :pnext (parse-integer pnext)
          :tlen (parse-integer tlen)
          :seq (or-* seq)
          :qual (or-* qual)
          :optional-fields (mapcar #'parse-optional-field optional-fields))))))

#; Scratch --------------------------------------------------------------------




(time (with-open-file
          (f "data/small.sam" :direction :input :element-type 'character)
        ;; (skip-sam-headers f)
        (loop :while (char= #\@ (peek-char nil f nil :eof))
              :do (read-sam-header-line f))
        (loop :until (eql :eof (peek-char nil f nil :eof))
              :collect (read-sam-record f))))
