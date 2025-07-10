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
    :save-runtime-options t))

#-sbcl
(error "Only SBCL is supported.")
