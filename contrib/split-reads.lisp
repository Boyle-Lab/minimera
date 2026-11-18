(eval-when (:compile-toplevel :load-toplevel :execute)
  (ql:quickload '(:conserve :faster :adopt :local-time)))

(defpackage :split-reads
  (:use :cl)
  (:export :build :toplevel))

(in-package :split-reads)


;;;; Configuration ------------------------------------------------------------
(defparameter *input-format* :fastq)
(defparameter *worker-threads* 1)
(defparameter *interactive* nil)


;;;; Types --------------------------------------------------------------------
(deftype u8 () '(unsigned-byte 8))
(deftype a8 () '(simple-array u8 (*)))


;;;; Utilities ----------------------------------------------------------------
(declaim (inline make-a8))

(defun make-a8 (size &rest make-array-args)
  (apply #'make-array size :element-type 'u8 make-array-args))


;;;; Split Point File ---------------------------------------------------------
(defun parse-split-point (string)
  (parse-integer string))

(defun parse-split-points (strings)
  (sort (mapcar #'parse-split-point strings) #'>))

(defun parse-splits (splits-path)
  (let ((result (make-hash-table :test #'equal)))
    (with-open-file (csv-stream splits-path)
      (loop :for line = (conserve:read-row csv-stream nil :eof)
            :until (eql line :eof)
            :do (setf (gethash (first line) result)
                      (parse-split-points (rest line)))))
    result))


;;;; FASTQ Splitting ----------------------------------------------------------
(defun split/fastq (record point)
  (let* ((read-id (faster:parse-sequence-id (faster:id record)))
         (prefix-id (map 'a8 #'char-code (format nil "~A_split-0" read-id)))
         (suffix-id (map 'a8 #'char-code (format nil "~A_split-1" read-id)))
         (prefix-seq (subseq (faster:seq record) 0 point))
         (suffix-seq (subseq (faster:seq record) point nil))
         (prefix-qs (subseq (faster:qs record) 0 point))
         (suffix-qs (subseq (faster:qs record) point nil)))
    (list (list prefix-id prefix-seq prefix-qs)
          (list suffix-id suffix-seq suffix-qs))))


;;;; BAM Splitting ------------------------------------------------------------
(defun parse-skips (string)
  (let ((result (make-array (count #\, string) :element-type 'fixnum))
        (start 1))
    (dotimes (i (length result))
      (multiple-value-bind (n pos)
          (parse-integer string :start start :junk-allowed t)
        (setf start (1+ pos)
              (aref result i) n)))
    result))

(defun parse-mods (string)
  (if (digit-char-p (char string 0))
    (list (parse-integer string))
    (coerce string 'list)))

(defstruct mm
  (base #\nul :type base-char)
  (strand #\+ :type base-char)
  (mods nil :type list)
  (skip-meaning nil :type (or null base-char))
  (skips #() :type (simple-array fixnum (*))))

(defun parse-methylation-tag (data)
  (let ((result nil))
    (ppcre:do-register-groups (base strand mods skip-meaning skips)
        ("([ACGTUN])([-+])([a-zA-Z]+|[0-9]+)([.?])?((?:,[0-9]+)*);" data)
      (push (make-mm :base (char base 0)
                     :strand (char strand 0)
                     :mods (parse-mods mods)
                     :skip-meaning (when skip-meaning
                                     (char skip-meaning 0))
                     :skips (parse-skips skips))
            result))
    (nreverse result)))

(defun render-methylation-tags (mms)
  (with-output-to-string (s)
    (dolist (mm mms)
      (format s "~A~A~A~A"
              (mm-base mm)
              (mm-strand mm)
              (let ((mods (mm-mods mm)))
                (if (integerp (first mods))
                  (first mods)
                  (concatenate 'string mods)))
              (or (mm-skip-meaning mm) ""))
      (loop :for n :across (mm-skips mm)
            :do (format s ",~D" n))
      (write-char #\; s))))

(defun count-bases-for-modification-splitting (seq split-point)
  "Return the count of each base up to the split point.

  The result is an alist of `((#\A . …) (#\C . …) (#\G . …) (#\T . …) (#\N . …))`.

  Note that `:n` is the number of *all* bases, not specifically Ns produced by
  the basecaller.  See the first paragraph of the MM tag's documentation in the
  SAM specfication section 1.7 more details.

  "
  (check-type seq a8)
  (let ((na 0) (nc 0) (ng 0) (nt 0) (nn split-point))
    (loop :for b :across seq
          :repeat split-point
          :do (ecase b
                (1 (incf na))
                (2 (incf nc))
                (4 (incf ng))
                (8 (incf nt))))
    (list (cons #\A na)
          (cons #\C nc)
          (cons #\G ng)
          (cons #\T nt)
          (cons #\N nn))))

(defun split-skips (skips target-base-count)
  ; split a skips array into two skips arrays P(refix) and S(uffix).
  ; target-base-count is the total count of the target bases before the split point

  ; first chomp off all the obviously-in-P skips
  (destructuring-bind (end . chomped)
      (loop :with n = 0
            :with i = 0
            :for skip :across skips
            :while (<= (+ n (1+ skip)) target-base-count)
            :do (progn (incf n (1+ skip))
                       (incf i))
            :finally (return (cons i n)))

    (if (= end (length skips))
      ;; if we have no more entries, all of the methylation happened before the
      ;; split point and the second array will be empty
      (cons (copy-seq skips) (make-array 0 :element-type 'fixnum))
      ;; otherwise we need to reduce the first skip in the second array by
      ;; however many bases are left in the prefix after the final modded base
      (let ((P (subseq skips 0 end))
            (S (subseq skips end)))
        (decf (aref S 0) (- target-base-count chomped))
        (cons P S)))))

(defun split-mm (mm counts)
  (destructuring-bind (skips-prefix . skips-suffix)
      (split-skips (mm-skips mm) (alexandria:assoc-value counts (mm-base mm)))
    (let ((P (copy-mm mm))
          (S (copy-mm mm)))
      (setf
        (mm-mods P) (copy-seq (mm-mods P))
        (mm-mods S) (copy-seq (mm-mods S))
        (mm-skips P) skips-prefix
        (mm-skips S) skips-suffix)
      (cons P S))))

(defun split-mms (mms seq split-point)
  (let ((counts (count-bases-for-modification-splitting seq split-point)))
    (mapcar (lambda (mm) (split-mm mm counts)) mms)))

(defun lengths-to-start-ends (lengths)
  (let ((start 0)
        (result '()))
    (map 'list (lambda (length &aux (end (+ start length)))
                 (push (cons start end) result)
                 (setf start end))
         lengths)
    (nreverse result)))

(defun mm-quality-length (mm)
  (* (length (mm-mods mm))
     (length (mm-skips mm))))

(defun split-ml (ml split-mms)
  (loop
    :with p-ml = (make-a8 (reduce #'+ split-mms :key (alexandria:compose #'mm-quality-length #'car)))
    :with s-ml = (make-a8 (reduce #'+ split-mms :key (alexandria:compose #'mm-quality-length #'cdr)))
    :with ml-start = 0
    :with p-start = 0
    :with s-start = 0
    :for (p-mm . s-mm) :in split-mms
    :for p-length = (mm-quality-length p-mm)
    :for s-length = (mm-quality-length s-mm)
    :for p-end = (+ p-start p-length)
    :for s-end = (+ s-start s-length)
    :do (progn
          (replace p-ml ml :start1 p-start :end1 p-end :start2 ml-start)
          (incf ml-start p-length)
          (incf p-start p-length)
          (replace s-ml ml :start1 s-start :end1 s-end :start2 ml-start)
          (incf ml-start s-length)
          (incf s-start s-length))
    :finally (return (cons p-ml s-ml))))

(defun split-modifications (mms ml mn seq split-point)
  "Split the base modification information, returning prefix and suffix."
  (let* ((split-mms (split-mms mms seq split-point))
         (split-ml (split-ml ml split-mms))
         (seq-length (length seq))
         (p-mn split-point)
         (s-mn (- seq-length split-point))
         (p (list (mapcar #'car split-mms) (car split-ml) p-mn))
         (s (list (mapcar #'cdr split-mms) (cdr split-ml) s-mn)))
    (assert (= mn seq-length (+ p-mn s-mn)) ()
            "Sanity check of MN:i:~D tag failed, seq was length ~D, prefix ~D suffix ~D."
            mn seq-length p-mn s-mn)
    (cons p s)))

(defun unmap-alignment (alignment)
  (setf (faster/bam::alignment-reference-id      alignment)        -1
        (faster/bam::alignment-pos               alignment)        -1
        (faster/bam::alignment-mapq              alignment)         0
        (faster/bam::alignment-bin               alignment)      4680
        (faster/bam::alignment-flag              alignment)         4 ; todo do this more properly
        (faster/bam::alignment-next-pos          alignment)        -1
        (faster/bam::alignment-next-reference-id alignment)        -1
        (faster/bam::alignment-cigar             alignment)  (vector)
        (faster/bam::alignment-tlen              alignment)         0))


(defconstant +max-qscore+ 93)

(eval-when (:compile-toplevel :load-toplevel :execute)
  (defun q->p% (q)
    ;; P = 10 ^ (-Q/10)
    (expt 10.0d0 (/ (- q) 10.0d0)))

  (defun make-q-table% ()
    (let ((result (make-array (list (1+ +max-qscore+)) :element-type 'double-float)))
      (loop :for q :from 0 :to +max-qscore+
            :do (setf (aref result q) (q->p% q)))
      result)))


(declaim (inline q->p p->q))

(defun q->p (q)
  "Return the error probability (as a double-float) corresponding to Q score `q`.

   `q` must be an `(integer 0 60)`, and a lookup table is used to make this
   computation fast.  If you need to convert arbitrary Q scores (e.g. a mean
   Q score) use `p->q%`.

  "
  (aref #.(make-q-table%) q))

(defun p->q (p)
  "Return the Q score (as a double-float) corresponding to error probability `p`."
  ;; Q = -10 * log10(P)
  (* -10.0d0 (log p 10.0d0)))

(defun mean-qs (phred &aux (length (length phred)))
  ;; TODO: Port this in from Minimera instead of copy/pasting?
  ;; TODO: Should we attempt to doradoize the computation?  I don't think so,
  ;; because this might be a split of a split of a split, so it's hard to know
  ;; where the *original* start was.  And Dorado's dropping of the original
  ;; bases is… somewhat janky anyway.
  (if (zerop length)
    0.0f0
    (coerce (p->q (/ (reduce #'+ phred :key #'q->p) length))
            'double-float)))

(defun estimate-read-duration (original-duration original-length new-length)
  (* (/ new-length original-length) original-duration))

(defun format-dorado-timestring (timestamp)
  (local-time:format-rfc3339-timestring nil timestamp
                                        :use-zulu nil
                                        :timezone local-time:+utc-zone+))

(defun estimate-suffix-start-time (original-start-time suffix-duration)
  (local-time:timestamp+ original-start-time (truncate (* 1000000000 suffix-duration)) :nsec))

(defun split-alignment (alignment split-point)
  (flet ((set-tag (alignment tag-id tag-type value)
           (let ((tag (faster/bam::alignment-tag alignment tag-id)))
             (when tag
               (assert (eql (faster/bam::tag-type tag) tag-type) ()
                 "Mismatched tag type ~S for tag ~A when splitting alignment."
                 tag-type tag)
               (setf (faster/bam::tag-value tag) value))))
         (remove-tag (alignment tag-id)
           (alexandria:removef (faster/bam::alignment-tags alignment) tag-id
                               :key #'faster/bam::tag-id :test #'string=)))
    (let* ((p (faster/bam::copy-alignment alignment))
           (s (faster/bam::copy-alignment alignment))
           (duration (faster/bam::tag-value (faster/bam::alignment-tag alignment "du")))
           (seq (faster/bam::alignment-seq alignment))
           (phred (faster/bam::alignment-phred alignment))
           (p-len split-point)
           (s-len (- (length seq) split-point))
           (p-seq (subseq seq 0 split-point))
           (s-seq (subseq seq split-point nil))
           (p-phred (subseq phred 0 split-point))
           (s-phred (subseq phred split-point nil))
           (p-start-time (local-time:parse-rfc3339-timestring (faster/bam::tag-value (faster/bam::alignment-tag alignment "st"))))
           (p-dur (estimate-read-duration duration (length seq) p-len))
           (s-dur (estimate-read-duration duration (length seq) s-len))
           (s-start-time (estimate-suffix-start-time p-start-time s-dur)))
      (unmap-alignment p)
      (unmap-alignment s)
      ;;x . ts  - i  - signal index start  - 10
      ;;x . ns  - i  - signal index end  - 2027
      ;;x . sm  - f  - scaling midpoint  - 769.8418
      ;;x . sd  - f  - scaling dispersion  - 125.65758
      ;;x . sv  - Z  - scaling version  - "pa"
      ;;x . mv  - B  - move table - (#\c . #(6 1 0 0  ...))
      ;;  . mx  - i  - read mux  - 4
      ;;  . ch  - i  - read channel - 2161
      ;;  . rn  - i  - read number - 33
      ;;  . fn  - Z  - filename  - "PBE89432_skip_059690b8_3033dae4_0.pod5"
      ;;  . RG  - Z  - run/model/barcoding - "3033dae4f3d039e5bba20f3215a8d405ca61e3c7_dna_r10.4.1_e8.2_400bps_sup@v5.0.0"
      ;;! . qs  - f  - mean basecall q-score - 6.9754105
      ;;! . du  - f  - duration of read in seconds   - 0.4054
      ;;! . st  - Z  - read start time  - "2025-12-16T20:52:21.914+00:00"
      ;;!   MN  - i  - length of seq at the time methylation was written  - 71
      ;;!   MM  - Z  -   - "C+h?,9;C+m?,9;"
      ;;!   ML  - B  - methylation quality scores  - (#\C . #(57 26))
      ;;?   dx  - i  - is duplex? (boolean)  - 0
      (dolist (id '("ts" "ns" "sm" "sd" "sv" "mv"))
        (remove-tag p id)
        (remove-tag s id))
      (setf (faster/bam::alignment-read-name p) (str:concat (faster/bam::alignment-read-name p) "_split-0")
            (faster/bam::alignment-read-name s) (str:concat (faster/bam::alignment-read-name s) "_split-1")
            (faster/bam::alignment-tags p) (faster/bam::copy-alignment-tags (faster/bam::alignment-tags p))
            (faster/bam::alignment-tags s) (faster/bam::copy-alignment-tags (faster/bam::alignment-tags s))
            (faster/bam::alignment-seq p) p-seq
            (faster/bam::alignment-seq s) s-seq
            (faster/bam::alignment-phred p) p-phred
            (faster/bam::alignment-phred s) s-phred)
      (set-tag p "qs" #\f (mean-qs p-phred))
      (set-tag s "qs" #\f (mean-qs s-phred))
      (set-tag p "du" #\f p-dur)
      (set-tag s "du" #\f s-dur)
      (set-tag p "st" #\Z (format-dorado-timestring p-start-time))
      (set-tag s "st" #\Z (format-dorado-timestring s-start-time))
      (cons p s))))


(defun split/bam (alignment point)
  (destructuring-bind (prefix . suffix) (split-alignment alignment point)
    (list prefix suffix)))


;;;; Work ---------------------------------------------------------------------
(defun work/fastq (splits record)
  (let* ((read-id (faster:parse-sequence-id (faster:id record)))
         (split-points (gethash read-id splits)))
    (if split-points
      (progn
        (assert (= (length split-points) 1) () "TODO: Support multiple split points")
        (split/fastq record (first split-points)))
      (list))))

(defun work/bam (splits alignment)
  (let* ((read-id (faster/bam::alignment-read-name alignment))
         (split-points (gethash read-id splits)))
    (if split-points
      (progn
        (assert (= (length split-points) 1) () "TODO: Support multiple split points")
        (split/bam alignment (first split-points)))
      (list))))


;;;; Output -------------------------------------------------------------------
(defun write-char-code (character &optional (stream *standard-output*))
  (write-byte (char-code character) stream))

(defun output/fastq (output-stream records)
  (dolist (record records)
    (destructuring-bind (id seq qs) record
      (write-char-code #\@ output-stream)
      (write-sequence id output-stream)
      (write-char-code #\newline output-stream)
      (write-sequence seq output-stream)
      (write-char-code #\newline output-stream)
      (write-char-code #\+ output-stream)
      (write-char-code #\newline output-stream)
      (write-sequence qs output-stream)
      (write-char-code #\newline output-stream))))

(defun output/bam (bgzf alignments)
  (dolist (alignment alignments)
    (faster/bgzf::write-bgzf bgzf (faster/bam::alignment-bytes alignment))))


;;;; Running ------------------------------------------------------------------
(defun run/fastq (fastq-path splits output-path)
  (with-open-file (output-stream output-path :element-type 'u8 :direction :output :if-exists :supersede)
    (faster:run fastq-path
                :format :fastq
                :work-function (alexandria:curry #'work/fastq splits)
                :output-function (alexandria:curry #'output/fastq output-stream)
                :interactive *interactive*
                :report-progress nil
                :worker-threads *worker-threads*)))

(defun run/bam (bam-path splits output-path)
  (with-open-file (output-stream output-path :element-type 'u8 :direction :output :if-exists :supersede)
    (let ((bgzf (faster/bgzf::make-bgzf-compressor output-stream)))
      (faster/bgzf::write-bgzf bgzf (faster/bam::bam-magic-bytes))
      (faster:run bam-path
                  :format :bam
                  :process-header (lambda (header)
                                    (faster/bgzf::write-bgzf bgzf (faster/bam::bam-header-bytes header)))
                  :process-references (lambda (references)
                                        (faster/bgzf::write-bgzf bgzf (faster/bam::bam-references-bytes references)))
                  :work-function (alexandria:curry #'work/bam splits)
                  :output-function (alexandria:curry #'output/bam bgzf)
                  :interactive *interactive*
                  :report-progress t
                  :worker-threads *worker-threads*)
      (faster/bgzf::finish-bgzf bgzf))))

(defun run (reads-path splits-path output-path &key (interactive t))
  (let ((*interactive* interactive)
        (splits (parse-splits splits-path)))
    (ecase *input-format*
      (:fastq (run/fastq reads-path splits output-path))
      (:bam (run/bam reads-path splits output-path)))))


;;;; Toplevel -----------------------------------------------------------------
(defun toplevel ()
  (run (rest (adopt:argv)) :interactive nil))

(defun build ()
  (sb-ext:save-lisp-and-die "build/split-reads"
    :toplevel #'toplevel
    :executable t
    :compression t))


#; Scratch --------------------------------------------------------------------
(run "../data/giab-tiny/tiny.fastq"
     "../results/tinytest-fastq/split-points.csv"
     "../results/tinytest-fastq/split-reads.fastq")

(let ((*input-format* :bam))
  (run "../data/giab-tiny/tiny.bam"
       "../results/tinytest-fastq/split-points.csv"
       "../results/tinytest-fastq/split-reads.bam"))
