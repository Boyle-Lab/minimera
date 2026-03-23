(ql:quickload '(:alexandria :losh :iterate))

(defpackage :quick-fastq
  (:use :cl :iterate :losh)
  (:export :toplevel :build))

(in-package :quick-fastq)

;; data is represented as a conses of (bases . quality-scores)

(defun phred-char (q)
  (code-char (+ (char-code #\!) q)))

(defun random-qscore ()
  (phred-char (+ 18 (random 4))))

(defun random-base ()
  (random-elt "ACTG"))

(defun random-dna (n)
  (cons (coerce (gimme n (random-base)) 'string)
        (coerce (gimme n (random-qscore)) 'string)))

(defun rev (dna)
  (cons (reverse (car dna))
        (reverse (cdr dna))))

(defun complement-base (base)
  (ecase base
    (#\A #\T)
    (#\T #\A)
    (#\C #\G)
    (#\G #\C)))

(defun comp (dna)
  (cons (map 'string #'complement-base (car dna))
        (cdr dna)))

(defun copy (dna)
  (cons (copy-seq (car dna))
        (copy-seq (cdr dna))))

(defun concat (a b)
  (cons (concatenate 'string (car a) (car b))
        (concatenate 'string (cdr a) (cdr b))))

(defun revcomp (dna)
  (rev (comp dna)))

(defun qual (q dna)
  (cons (car dna)
        (make-string (length (cdr dna)) :initial-element (phred-char q))))

(defun take-first (n dna)
  (cons (subseq (car dna) 0 n)
        (subseq (cdr dna) 0 n)))

(defun take-last (n dna &aux (len (length (car dna))))
  (cons (subseq (car dna) (- len n) len)
        (subseq (cdr dna) (- len n) len)))

(defun mutate (base)
  (random-elt (ecase base
                (#\A "TCG")
                (#\T "ACG")
                (#\C "ATG")
                (#\G "ATC"))))

(defun add-snp (freq dna)
  (iterate (with (seq . qs) = (copy dna))
           (for b :in-string seq :with-index i)
           (when (randomp freq)
             (setf (aref seq i) (mutate b)))
           (returning (cons seq qs))))

(defun add-ins (freq dna)
  (iterate (for b :in-string (car dna))
           (for q :in-string (cdr dna))
           (collect b :into seq)
           (collect q :into qs)
           (when (randomp freq)
             (collect (random-base) :into seq)
             (collect (random-qscore) :into qs))
           (returning (cons (coerce seq 'string) (coerce qs 'string)))))

(defun add-del (freq dna)
  (iterate (for b :in-string (car dna))
           (for q :in-string (cdr dna))
           (unless (randomp freq)
             (collect b :into seq)
             (collect q :into qs))
           (returning (cons (coerce seq 'string) (coerce qs 'string)))))

(defun add-err (freq dna)
  (add-ins freq (add-del freq (add-snp freq dna))))

(defun literal (seq)
  (cons (copy-seq seq)
        (coerce (gimme (length seq) (random-qscore)) 'string)))

(defun repeat (n seq)
  (reduce #'concat (gimme n seq)))

(defun run (binding-forms fastq-form)
  (let ((bindings (make-hash-table)))
    (labels ((r (form) (rev (eval-form form)))
             (c (form) (comp (eval-form form)))
             (rc (form) (revcomp (eval-form form)))
             (f (n form) (take-first n (eval-form form)))
             (l (n form) (take-last n (eval-form form)))
             (q (n form) (qual n (eval-form form)))
             (snp (freq form) (add-snp freq (eval-form form)))
             (ins (freq form) (add-ins freq (eval-form form)))
             (del (freq form) (add-del freq (eval-form form)))
             (err (freq form) (add-err freq (eval-form form)))
             (rep (n form) (repeat n (eval-form form)))
             (eval-form (form)
               (typecase form
                 (integer (random-dna form))
                 (string (literal form))
                 (vector (reduce #'concat form :key #'eval-form))
                 (symbol (gethash form bindings))
                 (list (destructuring-bind (op . args) form
                         (if (char= #\Q (char (symbol-name op) 0))
                           (apply #'q
                                  (parse-integer (subseq (symbol-name op) 1))
                                  args)
                           (apply (ecase op
                                    ((rep tr) #'rep)
                                    ((rev r) #'r)
                                    ((comp c) #'c)
                                    ((revcomp rc) #'rc)
                                    ((first f) #'f)
                                    ((last l) #'l)
                                    ((snp) #'snp)
                                    ((ins) #'ins)
                                    ((del) #'del)
                                    ((err) #'err))
                                  args))))
                 (t form))))

      ;; Process the bindings
      (dolist (binding-form binding-forms)
        (destructuring-bind (symbol form) binding-form
          (case symbol
            (:seed (setf *random-state* (sb-ext:seed-random-state form)))
            (t (setf (gethash symbol bindings) (eval-form form))))))

      ;; Process the FASTQ form using the bindings and print it
      (destructuring-bind (seq . qs)
          (eval-form fastq-form)
        (format nil "@quickfastq~%~A~%+~%~A" seq qs)))))

(defun read-form ()
  (let ((*package* (find-package :quick-fastq)))
    (read)))

(defun toplevel ()
  (setf *random-state* (make-random-state t))
  ; todo safe-read this
  (let ((binding-forms (read-form))
        (fastq-form (read-form)))
    (write-string (run binding-forms fastq-form))))

(defun build ()
  (sb-ext:save-lisp-and-die "build/qfq"
    :toplevel #'toplevel
    :executable t
    :compression t))
