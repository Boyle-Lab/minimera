(asdf:defsystem :minimera
  :description "A tool for detecting problematic reads in Oxford Nanopore data."
  :author "Steve Losh <steve@stevelosh.com>"
  :homepage "https://github.com/Boyle-Lab/minimera"

  :license "GPL-3.0-or-later"
  :version "0.4.0"

  :depends-on (:adopt
               :alexandria
               :cl-netpbm
               :conserve
               :faster
               :iterate
               :losh
               :parse-float
               :zpng)

  :serial t
  :components ((:module "src" :serial t
                :components ((:file "package")
                             (:file "hash")
                             (:file "patience")
                             (:file "m5x7")
                             (:file "main")
                             (:file "ui")))))


