(in-package :asdf-user)

(defsystem "ham-single-crack"
  :depends-on ("cl-mpm/examples/single-crack"
               "vgplot"
               "swank.live"
               "magicl")
  :serial t
  :components ((:file "single-crack")))
