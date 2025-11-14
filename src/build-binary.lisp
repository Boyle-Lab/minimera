;;; Copyright 2025 Steve Losh and contributors
;;; SPDX-License-Identifier: GPL-3.0-or-later

(require :asdf)
(require :uiop)
(asdf:load-system :minimera)

#+sbcl
(progn
  (sb-ext:gc :full t)
  (sb-ext:save-lisp-and-die
    "build/minimera"
    :executable t
    :compression t
    :toplevel #'minimera:toplevel
    :save-runtime-options :accept-runtime-options))

#-sbcl
(error "Only SBCL is supported.")
