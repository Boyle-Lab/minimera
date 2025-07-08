(asdf:defsystem :foldbacks
  :description "TODO"
  :author "Steve Losh <steve@stevelosh.com>"
  :homepage "https://docs.stevelosh.com/adopt/"

  :license "TODO"
  :version "0.0.1"

  :depends-on (:losh :iterate :alexandria :uuid :str :conserve :lparallel)


  :serial t
  :components ((:module "src" :serial t
                :components ((:file "package")
                             (:file "foldbacks")))))


