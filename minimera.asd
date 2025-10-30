(asdf:defsystem :minimera
  :description "TODO"
  :author "Steve Losh <steve@stevelosh.com>"
  :homepage "TODO"

  :license "TODO"
  :version "0.0.1"

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


