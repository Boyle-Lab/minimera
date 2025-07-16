(asdf:defsystem :minimera
  :description "TODO"
  :author "Steve Losh <steve@stevelosh.com>"
  :homepage "TODO"

  :license "TODO"
  :version "0.0.1"

  :depends-on (:losh :iterate :alexandria :uuid :str :conserve :lparallel :adopt :parse-float)


  :serial t
  :components ((:module "src" :serial t
                :components ((:file "package")
                             (:file "fastq")
                             (:file "main")
                             (:file "ui")))))


