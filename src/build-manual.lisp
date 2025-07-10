(ql:quickload :minimera)

(with-open-file (stream "build/minimera.1" :direction :output :if-exists :supersede)
  (adopt:print-manual minimera:*ui* :stream stream))
