(ql:quickload :minimera)

(defpackage :fastq-stats
  (:use :cl)
  (:export :build :toplevel))

(in-package :fastq-stats)

(defun work (fs)
  (let ((minimera::*dorado-mean-qscore* nil))
    (list
      (faster:parse-sequence-id (faster:id fs))
      (princ-to-string (length (faster:seq fs)))
      (princ-to-string (minimera::compute-mean-qscore (faster:n-bytes-to-quality-scores (faster:qs fs)))))))

(defun output (result)
  (conserve:write-row result))

(defun run (path &key (interactive t))
  (conserve:write-row (list "read id" "read length" "mean qscore"))
  (faster:run path
              :work-function #'work
              :output-function #'output
              :interactive interactive
              :report-progress nil
              :worker-threads 4))

(defun toplevel ()
  (run (rest (adopt:argv)) :interactive nil))

(defun build ()
  (sb-ext:save-lisp-and-die "build/fastq-stats"
    :toplevel #'toplevel
    :executable t
    :compression t))

