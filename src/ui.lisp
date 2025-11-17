;;; Copyright 2025 Steve Losh and contributors
;;; SPDX-License-Identifier: GPL-3.0-or-later

(in-package :minimera)

;;;; UI -----------------------------------------------------------------------
(defparameter *version* "unknown")

(adopt:define-string *documentation*
  "Minimera is a tool for detecting foldback chimeric and other problematic ~
   reads in Oxford Nanopore data using minimizers.~@
   ~@
   Reads will be read from the specified fastq file (or stdin if - is given) and ~
   the results written to foldbacks.csv inside the specified output directory.  ~
   If plots are generated they will be inside a plots/ subdirectory of the ~
   output directory.")

(adopt:define-string *manual*
  "~A~@
   ~@
   The resulting foldbacks.csv file will contain the following columns:~@
   ~@
   1. read-id: the ID of the read~@
   ~@
   2. read-length: the length (in base pairs) of the read~@
   ~@
   3. classification: one of foldback/monotonous/normal~@
   ~@
   4. monotony: monotony score for the read~@
   ~@
   5. foldback-point: estimated position of the foldback point for foldbacks, empty otherwise~@
   ~@
   6. llqma: length of the longest low-quality moving average run in the read~@
   ~@
   7. llqma-start: start of the longest low-quality moving average run, blank if no run exists~@
   ~@
   8. llqma-end: end of the longest low-quality moving average run, blank if no run exists~@
   ~@
   9. processing-time-microsec: how long minimera took to process this read~@
   ~@
   New columns may be added after these in future versions of Minimera, so when ~
   processing foldbacks.csv make sure to allow extra columns if you want to ~
   remain forwards compatible."
  *documentation*)


(defparameter *o/help*
  (adopt:make-option 'help
    :help "Display help and exit."
    :long "help"
    :short #\h
    :reduce (constantly t)))

(defparameter *o/version*
  (adopt:make-option 'version
    :help "Display version and exit."
    :long "version"
    :short #\v
    :reduce (constantly t)))


(defparameter *o/threads*
  (adopt:make-option 'threads
    :help "Number of worker threads to spawn. Does not include the reader and writer threads (default: 1)."
    :long "threads"
    :short #\j
    :parameter "N"
    :reduce #'adopt:last
    :initial-value 1
    :key #'parse-integer))

(defparameter *o/output-directory*
  (adopt:make-option 'output-directory
    :help "Output directory (required)."
    :long "output"
    :short #\o
    :parameter "PATH"
    :reduce #'adopt:last))

(defparameter *o/k*
  (adopt:make-option 'k
    :help "Size of kmers (in base pairs) to use for minimizers (default: 8)."
    :long "kmer-size"
    :short #\K
    :parameter "N"
    :reduce #'adopt:last
    :initial-value 8
    :key #'parse-integer))

(defparameter *o/w*
  (adopt:make-option 'w
    :help "Size (total) of windows (in base pairs) to use for minimizer sketches (default: 16)."
    :long "window-size"
    :short #\W
    :parameter "N"
    :reduce #'adopt:last
    :initial-value 16
    :key #'parse-integer))

(defparameter *o/foldback-position-epsilon*
  (adopt:make-option 'foldback-position-epsilon
    :help "How close (in base pairs) a cluster must be to the beginning of read 2 and end of read 1 to be considered as a foldback cluster (default: 50)."
    :long "foldback-position-epsilon"
    :short #\F
    :parameter "N"
    :reduce #'adopt:last
    :initial-value 50
    :key #'parse-integer))

(defparameter *o/intercept-epsilon*
  (adopt:make-option 'intercept-epsilon
    :help "Epsilon (in y-intercept space) used to cluster colinear points during the initial clustering step (default: 30)."
    :long "intercept-epsilon"
    :short #\I
    :parameter "N"
    :reduce #'adopt:last
    :initial-value 30
    :key #'parse-integer))

(defparameter *o/gap-epsilon*
  (adopt:make-option 'gap-epsilon
    :help "Maximum width (in y-intercept space) of an allowable gap in a cluster.  Gaps wider than this will result in splitting the cluster (default: 300)."
    :long "gap-epsilon"
    :short #\G
    :parameter "N"
    :reduce #'adopt:last
    :initial-value 300
    :key #'parse-integer))

(defparameter *o/minimum-cluster-length*
  (adopt:make-option 'minimum-cluster-length
    :help "Minimum number of allowable points in a cluster.  Clusters with fewer than this many hits will be removed (default: 40)."
    :long "minimum-cluster-length"
    :short #\C
    :parameter "N"
    :reduce #'adopt:last
    :initial-value 40
    :key #'parse-integer))

(defparameter *o/minimum-foldback-length-absolute*
  (adopt:make-option 'minimum-foldback-length-absolute
    :help "Minimum length (in base pairs) of a foldback region.  Regions shorter than this will be excluded (default: 50)."
    :long "minimum-foldback-length-absolute"
    :short #\A
    :parameter "N"
    :reduce #'adopt:last
    :initial-value 50
    :key #'parse-integer))

(defparameter *o/minimum-foldback-length-relative*
  (adopt:make-option 'minimum-foldback-length-relative
    :help "Minimum length (as a fraction of total read length) of a foldback region.  Regions shorter than this will be excluded (default: 0.05)."
    :long "minimum-foldback-length-relative"
    :short #\R
    :parameter "X"
    :reduce #'adopt:last
    :initial-value 0.05
    :key #'parse-float:parse-float))

(defparameter *o/monotony-threshold*
  (adopt:make-option 'monotony-threshold
    :help "Monotony score above which reads will be classified as monotonous and not analyzed for foldback finding (default: 0.80)."
    :long "monotony-threshold"
    :short #\M
    :parameter "X"
    :reduce #'adopt:last
    :initial-value 0.80
    :key #'parse-float:parse-float))

(defparameter *o/low-quality-threshold*
  (adopt:make-option 'low-quality-threshold
    :help "Threshold (inclusive) for low-quality Phred scores when computing LLQMA (default: 3)."
    :long "low-quality-threshold"
    :short #\Q
    :parameter "Q"
    :reduce #'adopt:last
    :initial-value 3
    :key #'parse-integer))

(defparameter *o/low-quality-window-size*
  (adopt:make-option 'low-quality-window-size
    :help "Size of moving average window when computing LLQMA (default: 10)."
    :long "low-quality-window-size"
    :short #\L
    :parameter "N"
    :reduce #'adopt:last
    :initial-value 10
    :key #'parse-integer))

(adopt:defparameters (*o/plot/foldbacks* *o/plot/no-foldbacks*)
  (adopt:make-boolean-options 'plot-foldbacks
    :long "plot-foldbacks"
    :help "Generate plots for all reads determined to be foldbacks."
    :help-no "Do not generate plots for reads determined to be foldbacks (the default)."))

(adopt:defparameters (*o/plot/normal* *o/plot/no-normal*)
  (adopt:make-boolean-options 'plot-normal
    :long "plot-normal"
    :help "Generate plots for all reads determined to be normal."
    :help-no "Do not generate plots for reads determined to be normal (the default)."))

(adopt:defparameters (*o/progress* *o/no-progress*)
  (adopt:make-boolean-options 'progress
    :long "progress"
    :help "Report progress as reads are classified (the default)."
    :help-no "Do not report progress."
    :initial-value t))


(defparameter *examples*
  '(("Run minimera on a FASTQ:"
     . "minimera path/to/data.fastq --output ./results/")
    ("Run minimera with 8 worker threads, plotting any foldbacks:"
     . "minimera path/to/data.fastq --threads 8 --plot-foldbacks --output ./results/")))

(defparameter *ui*
  (adopt:make-interface
    :name "minimera"
    :usage "[OPTIONS] --output=PATH FASTQ"
    :summary "detect foldback chimeric reads using minimizers"
    :help *documentation*
    :manual *manual*
    :examples *examples*
    :contents
    (list *o/help*
          *o/version*
          *o/threads*
          *o/progress*
          *o/no-progress*
          *o/output-directory*
          (adopt:make-group 'foldback
            :title "Foldback Options"
            :options (list *o/k*
                           *o/w*
                           *o/foldback-position-epsilon*
                           *o/intercept-epsilon*
                           *o/gap-epsilon*
                           *o/minimum-cluster-length*
                           *o/minimum-foldback-length-absolute*
                           *o/minimum-foldback-length-relative*))
          (adopt:make-group 'monotony
            :title "Monotony Options"
            :options (list *o/monotony-threshold*))
          (adopt:make-group 'llqma
            :title "Low Quality Moving Average Options"
            :options (list *o/low-quality-window-size*
                           *o/low-quality-threshold*))
          (adopt:make-group 'plotting
            :title "Plotting Options"
            :help "Minimera can optionally generate plots of the minimizers for reads.  Generating these plots is slow, but they can be useful to debug edge cases."
            :options (list *o/plot/foldbacks*
                           *o/plot/no-foldbacks*
                           *o/plot/normal*
                           *o/plot/no-normal*)))))


(defun toplevel (&optional argv)
  (declare (ignore argv)) ; buildapp junk
  (adopt::quit-on-ctrl-c ()
    (multiple-value-bind (arguments options) (adopt:parse-options-or-exit *ui*)
      (cond ((gethash 'help options)
             (adopt:print-help-and-exit *ui*))
            ((gethash 'version options)
             (format t "minimera ~A~%" *version*)
             (adopt:exit)))
      (setf
        *interactive* nil
        *worker-threads* (gethash 'threads options)
        *k* (gethash 'k options)
        *w* (gethash 'w options)
        *foldback-position-epsilon* (gethash 'foldback-position-epsilon options)
        *intercept-epsilon* (gethash 'intercept-epsilon options)
        *gap-epsilon* (gethash 'gap-epsilon options)
        *minimum-cluster-length* (gethash 'minimum-cluster-length options)
        *minimum-foldback-length-absolute* (gethash 'minimum-foldback-length-absolute options)
        *minimum-foldback-length-relative* (gethash 'minimum-foldback-length-relative options)
        *monotony-threshold* (gethash 'monotony-threshold options)
        *low-quality-threshold* (gethash 'low-quality-threshold options)
        *low-quality-window-size* (gethash 'low-quality-window-size options)
        *plot-foldbacks* (gethash 'plot-foldbacks options)
        *plot-normal* (gethash 'plot-normal options)
        *output-directory* (gethash 'output-directory options)
        *report-progress* (gethash 'progress options))
      (handler-case
          (progn
            (assert *output-directory* () "Output directory must be specified.")
            (assert (= 1 (length arguments)) () "Exactly one .fastq file must be specified.")
            (assert (<= 1 *w* 31) () "Invalid window size ~A, must be in the range [1, 31]." *w*)
            (assert (<= 1 *k* 30) () "Invalid kmer size ~A, must be in the range [1, 30]." *k*)
            (assert (> *w* *k*) () "Window size (~A) must be larger than kmer size (~A)." *w* *k*)
            (assert (<= 1 *low-quality-threshold* 60) ()
              "Invalid low-quality threshold ~A, must be in the range [1, 60]." *low-quality-threshold*)
            (assert (<= 1 *low-quality-window-size* 4096) ()
              "Invalid low-quality window size ~A, must be in the range [1, 4096]." *low-quality-window-size*)
            (assert (plusp *foldback-position-epsilon*) ()
              "Foldback position epsilon (~A) must be positive." *foldback-position-epsilon*)
            (assert (plusp *intercept-epsilon*) ()
              "Intercept epsilon (~A) must be positive." *intercept-epsilon*)
            (assert (plusp *gap-epsilon*) ()
              "Gap epsilon (~A) must be positive." *gap-epsilon*)
            (assert (plusp *minimum-cluster-length*) ()
              "Minimum cluster length (~A) must be positive." *minimum-cluster-length*)
            (assert (plusp *minimum-foldback-length-absolute*) ()
              "Minimum foldback length (absolute) (~A) must be positive." *minimum-foldback-length-absolute*)
            (assert (<= 0.0 *minimum-foldback-length-relative* 1.0) ()
              "Minimum foldback length (relative) (~A) must be in the range [0, 1]." *minimum-foldback-length-relative*)
            (assert (<= 0.0 *monotony-threshold* 1.0) ()
              "Monotony threshold (~A) must be in the range [0, 1]." *monotony-threshold*)
            (run (first arguments)))
        (error (e)
               (adopt:print-error-and-exit e))))))
