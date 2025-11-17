;;; Copyright 2025 Steve Losh and contributors
;;; SPDX-License-Identifier: GPL-3.0-or-later

(ql:quickload :minimera)

(with-open-file (stream "build/minimera.1" :direction :output :if-exists :supersede)
  (adopt:print-manual minimera:*ui* :stream stream))

(with-open-file (stream "DOCUMENTATION.markdown" :direction :output :if-exists :supersede)
  (adopt::print-manual/md minimera:*ui* :stream stream))
