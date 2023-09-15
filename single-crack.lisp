;; (restrict-compiler-policy 'speed 3 3)
;; (restrict-compiler-policy 'debug 0 0)
;; (restrict-compiler-policy 'safety 0 0)
;; (setf *block-compile-default* t)
;(push :magicl.use-mkl *features*)
;(pushnew :cl-mpm-fbar *features*)
;; (setf *features* (delete :cl-mpm-pic *features*))
(ql:quickload "cl-mpm/examples/single-crack" :silent t)
(defpackage :ham-single-crack
  (:use :cl
        :cl-mpm/examples/single-crack))
(sb-ext:restrict-compiler-policy 'speed  0 0)
(sb-ext:restrict-compiler-policy 'debug  3 3)
(sb-ext:restrict-compiler-policy 'safety 3 3)
(in-package :ham-single-crack)

;; (ql:quickload "magicl")
;; (ql:quickload "cl-mpm")
;; (ql:quickload "cl-mpm/examples/single-crack" :silent t)
;; (in-package :cl-mpm/examples/single-crack)
(defparameter *debug* t)

(defparameter *ice-height* 0d0)
(defparameter *ice-length* 0d0)
(defparameter *crack-depth* 0d0)
(defparameter *original-crack-height* 0d0)
(defparameter *last-min* 0d0)
(defparameter *water-height* 0d0)
(defparameter *crack-water-width* 0d0)
(defparameter *meltwater-fill* 0.20d0)
(defparameter *floor-bc* nil)
(defparameter *sim* nil)
(defparameter *t* 0)
(defparameter *crack-water-bc* nil)
(defparameter *crack-water-height* 0d0)
(defparameter *sim-step* 0d0)
(defparameter *velocity* '())
(defparameter *time* '())
(defparameter *t* 0)
(defparameter *x* 0d0)
(defparameter *x-pos* '())
(defparameter *cfl-max* '())
(defparameter *max-stress* '())
(defparameter *sim-step* 0)

(defun length-from-def (sim mp dim)
  (let* ((mp-scale 2)
         (h-initial (magicl:tref (cl-mpm/particle::mp-domain-size mp) dim 0))
                                        ;(h-initial  (/ (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)) mp-scale))
         )
    h-initial
                                        ;(* (magicl:tref (cl-mpm::mp-deformation-gradient mp) dim dim) h-initial)
    ))
(defun max-stress (mp)
  (declare (optimize (speed 0) (debug 3)))
  (multiple-value-bind (l v) (magicl:eig (cl-mpm::voight-to-matrix (cl-mpm/particle:mp-stress mp)))
    (apply #'max l)
    ;; (- (max 0d0 (apply #'max l))
    ;;    (max 0d0 (apply #'min l)))
    ;; (cl-mpm/fastmath::voigt-tensor-reduce-simd (cl-mpm/particle::mp-velocity-rate mp))
    ;; (magicl:tref (cl-mpm/particle::mp-velocity-rate mp) 2 0)
    ;; (cl-mpm/particle::mp-damage-ybar mp)
    ;; (cl-mpm/constitutive::effective-strain-rate (cl-mpm/particle::mp-eng-strain-rate mp))
    ;; (cl-mpm/particle::mp-time-averaged-visc mp)
    ;; (magicl:tref (cl-mpm/particle::mp-stress mp) 0 0)
    )
  )

(defun plot (sim &optional (plot :stress))
  (declare (optimize (speed 0) (debug 3)))
  (vgplot:format-plot t "set palette defined (0 'blue', 1 'red')")
  (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh sim)))
         (ms-x (first ms))
         (ms-y (second ms))
         )
    (vgplot:format-plot t "set object 1 rect from ~f,0 to ~f,~f fc rgb 'blue' fs transparent solid 0.5 noborder behind" *crack-water-width* ms-x *water-height*)
    (when *crack-water-bc*
      (vgplot:format-plot t "set object 2 rect from 0,0 to ~f,~f fc rgb 'blue' fs transparent solid 0.5 noborder behind" *crack-water-width* *crack-water-height*))
    )

  (multiple-value-bind (x y c stress-y lx ly e density temp vx ybar dist)
    (loop for mp across (cl-mpm:sim-mps sim)
          collect (magicl:tref (cl-mpm::mp-position mp) 0 0) into x
          collect (magicl:tref (cl-mpm::mp-position mp) 1 0) into y
          collect (magicl:tref (cl-mpm::mp-velocity mp) 0 0) into vx
          collect (length-from-def sim mp 0) into lx
          collect (length-from-def sim mp 1) into ly
          collect (if (slot-exists-p mp 'cl-mpm/particle::damage) (cl-mpm/particle:mp-damage mp) 0) into c
          collect (if (slot-exists-p mp 'cl-mpm/particle::damage-ybar) (cl-mpm/particle::mp-damage-ybar mp) 0) into ybar
          collect (if (slot-exists-p mp 'cl-mpm/particle::temperature) (cl-mpm/particle:mp-temperature mp) 0) into temp
          collect (if (slot-exists-p mp 'cl-mpm/particle::strain-energy-density) (cl-mpm/particle::mp-strain-energy-density mp) 0) into e
          collect (/ (cl-mpm/particle:mp-mass mp) (cl-mpm/particle:mp-volume mp)) into density
          ;; collect (cl-mpm/particle:mp-volume mp) into density
          collect (max-stress mp) into stress-y
          ;; collect (local-dist sim mp) into dist
          collect 0d0 into dist
          finally (return (values x y c stress-y lx ly e density temp vx ybar dist)))

    (let* ((node-x '())
           (node-y '())
           (node-c '())
           (mesh (cl-mpm:sim-mesh *sim*))
           (nodes (cl-mpm/mesh:mesh-nodes mesh))
           )
      (dotimes (i (array-total-size nodes))
        (let ((n (row-major-aref nodes i)))
          (with-accessors ((index cl-mpm/mesh:node-index)
                           (boundary cl-mpm/mesh::node-boundary-node)
                           (boundary-s cl-mpm/mesh::node-boundary-scalar)
                           (active cl-mpm/mesh::node-active)
                           )

              n
            (let ((n-ratio (/ (cl-mpm/mesh::node-volume n) (cl-mpm/mesh::node-volume-true n))))
              (when active
                (if boundary
                    ;; (print index)
                    (destructuring-bind (x y) (cl-mpm/mesh:index-to-position mesh index)
                      ;; (push (nth 0 index) node-x)
                      ;; (push (nth 1 index) node-y)
                      (push x node-x)
                      (push y node-y)
                      (push 2 node-c)
                      ;; (push n-ratio node-c)
                      ;; (push boundary-s node-c)
                      )
                    (destructuring-bind (x y) (cl-mpm/mesh:index-to-position mesh index)
                      (push x node-x)
                      (push y node-y)
                      ;; (push n-ratio node-c)
                      (push 0 node-c)
                      ;; (push boundary-s node-c)
                      ))
                ;; (push (/ (cl-mpm/mesh::node-volume n) (cl-mpm/mesh::node-volume-true n)) node-c)
                )))))
      (cond
        ((eq plot :point)
         ;; (vgplot:format-plot t "set cbrange [0:1]")
         ;; (vgplot:plot x y ";;with points pt 7")
                                        ;(vgplot:format-plot t "set cbrange [0:2]")
         (vgplot:format-plot t "set style fill empty")
         (if node-x
             ;; (multiple-value-bind (xp yp) (max-spacial (map 'list #'identity (cl-mpm::sim-mps *sim*)))
             (progn
               (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min node-c) (+ 0.01 (apply #'max node-c)))
               (vgplot:plot
                ;; x y ";;with points pt 7"
                x y lx ly ";;with ellipses"
                node-x node-y node-c ";;with points pt 7 lc palette"
                ;; (list xp) (list yp) ";;with points"
                ))
               ;)
             (vgplot:plot
              x y lx ly ";;with ellipses")
             )
         )
        ((eq plot :contact)
         (let* ((contact-points (cl-mpm/penalty::collect-contact-points-bc (cl-mpm:sim-mesh *sim*) (cl-mpm:sim-mps *sim*) *floor-bc*))
                (c-x (mapcar (lambda (p) (magicl:tref p 0 0)) contact-points))
                (c-y (mapcar (lambda (p) (magicl:tref p 1 0)) contact-points))
                (c-c (mapcar (lambda (p) 1d0) contact-points))
           )
           ;; (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min stress-y) (+ 1e-20 (apply #'max stress-y)))
           (vgplot:format-plot t "set cbrange [~f:~f]" 0d0 (+ 1e-40 (apply #'max c)))
         ;; (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min node-c) (+ 0.01 (apply #'max node-c)))
         (if c-x
             ;; (multiple-value-bind (xp yp) (max-spacial (map 'list #'identity (cl-mpm::sim-mps *sim*)))
             (vgplot:plot
              ;; x y ";;with points pt 7"
              x y c ";;with points pt 7 lc palette"
              c-x c-y ";;with points pt 7"
              ;; (list xp) (list yp) ";;with points"
              )
             (vgplot:plot
              ;; x y lx ly ";;with ellipses"
              x y c ";;with points pt 7 lc palette"
              )
             ))
         )
        ((eq plot :damage)
         (vgplot:format-plot t "set style fill solid")
         ;; (vgplot:format-plot t "set cbrange [0:1]")
         ;; (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min c) (+ 1e-6 (apply #'max c)))
         (vgplot:format-plot t "set cbrange [~f:~f]" 0d0 (+ 1e-6 (apply #'max c)))
         ;; (vgplot:format-plot t "set cbrange [~f:~f]" (+ 0d0 (apply #'min c)) (+ 1e-6 (apply #'max c)))
         ;; (vgplot:plot x y c ";;with points pt 7 lc palette")
         (vgplot:plot x y lx ly c ";;with ellipses lc palette")
         )
        ((eq plot :dist)
         (vgplot:format-plot t "set cbrange [~f:~f]" 0d0 (+ 1e-6 (apply #'max dist)))
         (vgplot:plot x y dist ";;with points pt 7 lc palette")
         )
        ((eq plot :velocity)
         (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min vx) (+ 0.01 (apply #'max vx)))
         (vgplot:plot x y vx ";;with points pt 7 lc palette")
         )
        ((eq plot :temperature)
         (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min temp) (+ 0.01 (apply #'max temp)))
         ;; (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min c) (+ 0.01 (apply #'max c)))
         (vgplot:plot x y temp ";;with points pt 7 lc palette")
         )
        ((eq plot :energy)
         (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min e) (+ 0.01 (apply #'max e)))
         (vgplot:plot x y e ";;with points pt 7 lc palette")
         )
        ((eq plot :stress)
         (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min stress-y) (+ 1e-40 (apply #'max stress-y)))
         (vgplot:plot x y stress-y ";;with points pt 7 lc palette"))
        ((eq plot :damage-ybar)
         (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min ybar) (+ 1e-20 (apply #'max ybar)))
         (vgplot:plot x y ybar ";;with points pt 7 lc palette"))
        ((eq plot :deformed)
         ;; (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min stress-y) (+ 0.01 (apply #'max stress-y)))
         (vgplot:plot x y lx ly ";;with ellipses"))
        ((eq plot :density)
         (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min density) (+ 0.01 (apply #'max density)))
         (vgplot:plot x y e ";;with points pt 7 lc palette")
         )))
    )
  ;; (vgplot:format-plot t "replot (~f*x + ~f)~%" *sliding-slope* *sliding-offset*)
  (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh sim)))
         (ms-x (first ms))
         (ms-y (second ms))
         )
    (vgplot:axis (list 0 ms-x
                       0 ms-y))
    (vgplot:format-plot t "set size ratio ~f" (/ ms-y ms-x)))
    (let ((h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*))))
      (vgplot:format-plot t "set ytics ~f" h)
      (vgplot:format-plot t "set xtics ~f" h))

  (vgplot:replot))

(defun rectangle-sdf (position size)
  (lambda (pos)
      (let* ((position (magicl:from-list position '(2 1) :type 'double-float))
             (dist-vec (magicl:.- (magicl:map! #'abs (magicl:.- pos position))
                                  (magicl:from-list size '(2 1) :type 'double-float))))

        (+ (sqrt (magicl::sum
                  (magicl:map! (lambda (x) (* x x))
                               (magicl:map! (lambda (x) (max 0d0 x)) dist-vec))))
           (min (max (magicl:tref dist-vec 0 0)
                     (magicl:tref dist-vec 1 0)
                     ) 0d0)))))

(defun setup-test-column (size block-size block-offset &optional (e-scale 1d0) (mp-scale 1d0)
                          &rest mp-args)
  (declare (optimize (speed 0)))
  (let* ((sim (cl-mpm/setup::make-block (/ 1 e-scale)
                                        (mapcar (lambda (s) (* s e-scale)) size)
                                        #'cl-mpm/shape-function:make-shape-function-bspline
                                        'cl-mpm/damage::mpm-sim-damage
                                        ))

         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         (angle 0d0)
         (density *ice-density*)
         (elements (mapcar (lambda (s) (* e-scale (/ s 2))) size)))
    (progn
      (let ((block-position
              (mapcar #'+ (list (* h-x (- (+ (/ 1 (* 2 mp-scale))) 0))
                                (* h-y (+ (/ 1d0 (* 2d0 mp-scale)))))
                      block-offset)))
        (setf (cl-mpm:sim-mps sim)
              (cl-mpm/setup::make-mps-from-list
               (cl-mpm/setup::make-block-mps-list
                block-position
                block-size
                (mapcar (lambda (e) (* e e-scale mp-scale)) block-size)
                density
                ;'cl-mpm/particle::particle-viscoplastic-damage
                'cl-mpm/particle::particle-elastic-damage
                :E 1.5d9
                :nu 0.3250d0
                ;; :nu 0.0000d0
                ;; :visc-factor 111d6
                ;; :visc-power 3d0

                :initiation-stress 0.1d6
                :damage-rate 0d0
                :critical-damage 0.50d0
                :local-length 20d0
                ;; :local-length-damaged 20d0
                ;; :local-length-damaged 20d0
                :local-length-damaged 0.1d0
                :damage 0.0d0

                :gravity -9.8d0
                ;; :gravity-axis (magicl:from-list '(0.5d0 0.5d0) '(2 1))
                :index 0
                ;; :damage-model
                ;; (cl-mpm/damage::make-damage-model-creep
                ;;  :damage-rate 1d-4
                ;;  :initiation-stress 0.01d6
                ;;  )
                ;:angle angle
                )))
        )
      (let ((mass-scale 1d1))
        (setf (cl-mpm::sim-mass-scale sim) mass-scale)
        (setf (cl-mpm:sim-damping-factor sim)
              (* 0.1d0 mass-scale)
              ;; 1d0
              ;; (* 0.00000001d0 mass-scale)
              ;; 0.1d0
              ;; 0.01d0
              ;; 0.1d0
              ;; 0.0d0
              ;; 100d0
              )
        )
      (setf (cl-mpm:sim-mass-filter sim) 1d-15)
      (setf (cl-mpm::sim-allow-mp-split sim) nil)
      (setf (cl-mpm::sim-allow-mp-damage-removal sim) t)
      (setf (cl-mpm::sim-nonlocal-damage sim) t)
      (setf (cl-mpm::sim-enable-damage sim) nil)
      (setf (cl-mpm::sim-mp-damage-removal-instant sim) nil)
      (setf (cl-mpm:sim-dt sim) 1d-4)
      (setf (cl-mpm:sim-bcs sim) (make-array 0))
      (setf (cl-mpm:sim-bcs sim)
            (cl-mpm/bc::make-outside-bc-var
             (cl-mpm:sim-mesh sim)
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(0 nil)))
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(0 nil)))
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(nil 0)))
             (lambda (i) (cl-mpm/bc:make-bc-fixed (mapcar #'+ i '(0 0)) '(nil 0)))
             ))
      (format t "Bottom level ~F~%" h-y)
      (defparameter *ocean-fill* 0.50d0)
      (let* ((terminus-size (+ (second block-size) (* 0d0 (first block-size))))
             (ocean-x 1000)
            ;; (ocean-y (+ h-y (* 0.90d0 0.0d0 terminus-size)))
             ;; (ocean-y (* (round (* 0.5d0 terminus-size) h-y) h-y))
             (ocean-y (+ 0d0 (* *ocean-fill* terminus-size)))
            ;(angle -1d0)
            )

        (format t "Ocean level ~a~%" ocean-y)
        (defparameter *water-height* ocean-y)
        (defparameter *meltwater-fill* 0.20d0)
        (defparameter *floor-bc*
          (cl-mpm/penalty::make-bc-penalty-point-normal
           sim
           (magicl:from-list (list (sin (- (* pi (/ angle 180d0))))
                                   (cos (+ (* pi (/ angle 180d0))))) '(2 1))
           (magicl:from-list (list 00d0 (+ 1d0 h-y)) '(2 1))
           (* *ice-density* 1d3)
           0.0d0
           ))
        (defparameter *crack-water-height* terminus-size)
        (setf *crack-water-height* 0d0)
        (defparameter *crack-water-width* (- (first block-size) h-x))
        (defparameter *crack-water-bc*
          (cl-mpm/buoyancy::make-bc-buoyancy-clip
           sim
           ocean-y
           ;*water-density*
           1000d0
           (lambda (pos datum)
             (and
              (< (magicl:tref pos 0 0) *crack-water-width*)))))
        (setf (cl-mpm::sim-bcs-force-list sim)
              (list
               (cl-mpm/bc:make-bcs-from-list
                (list
                 (cl-mpm/buoyancy::make-bc-buoyancy-clip
                  sim
                  ocean-y
                  *water-density*
                  (lambda (pos datum)
                    (and
                     (>= (magicl:tref pos 0 0) *crack-water-width*)
                     )))))
               (cl-mpm/bc:make-bcs-from-list
                (list
                 *crack-water-bc*))
               (cl-mpm/bc:make-bcs-from-list
                (list
                 (cl-mpm/bc::make-bc-closure
                  '(0 0)
                  (lambda ()
                    (cl-mpm/buoyancy::set-pressure-all sim *crack-water-bc*)))))
               ))
        )
      (let ((normal (magicl:from-list (list (sin (- (* pi (/ angle 180d0))))
                                            (cos (+ (* pi (/ angle 180d0))))) '(2 1))))
        (defparameter *sliding-slope* 0d0)
        (defparameter *sliding-offset* (- h-y (* 0d0 0d0))))
      sim)))

(defparameter *ice-density* 917d0)
(defparameter *water-density* 1020d0)
;; (defparameter *ice-density* 900)
;; (defparameter *water-density* 1000)
;Setup
;; (defparameter *run-sim* nil)
(defun setup ()
  (declare (optimize (speed 0)))
  (defparameter *run-sim* nil)
  (let* ((mesh-size 5)
         (mps-per-cell 2)
         (shelf-height 125)
         (shelf-length 500)
         (offset (list 0 0))
         )
    (defparameter *sim*
      (setup-test-column (list (+ shelf-length (* 1 shelf-height))
                                                 (+ shelf-height 100))
                         (list shelf-length shelf-height)
                         (mapcar #'+ offset (list 000 (* 0 mesh-size)))
                         (/ 1 mesh-size) mps-per-cell))

    ;;Delete all the plotted frames
    (ensure-directories-exist (uiop:merge-pathnames* "./output/"))
    (ensure-directories-exist (uiop:merge-pathnames* "./outframes/"))
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
    (let ((cut-depth 20))
      (cl-mpm/setup::remove-sdf
       *sim*
       (rectangle-sdf
        (list (* 0.5d0 shelf-length)
              shelf-height)
        (list
         5
         cut-depth
         )))
      (defparameter *ice-height* shelf-height)
      (defparameter *ice-length* shelf-length)
      (defparameter *crack-depth* 0d0)
      (defparameter *original-crack-height* (- shelf-height cut-depth))
      (defparameter *last-min* *original-crack-height*)
      )
    (loop for mp across (cl-mpm:sim-mps *sim*)
          when
          (or
           (<= (magicl:tref (cl-mpm/particle:mp-position mp) 0 0) (* 0.25d0 shelf-length))
           (>= (magicl:tref (cl-mpm/particle:mp-position mp) 0 0) (* 0.75d0 shelf-length))
           )
          do (setf  (cl-mpm/particle::mp-damage-rate mp) 0d0))
    )
  (print (cl-mpm:sim-dt *sim*))
  (format t "Sim MPs: ~a~%" (length (cl-mpm:sim-mps *sim*)))
  (defparameter *velocity* '())
  (defparameter *time* '())
  (defparameter *t* 0)
  (defparameter *x* 0d0)
  (defparameter *x-pos* '())
  (defparameter *cfl-max* '())
  (defparameter *sim-step* 0)
  (setf *max-stress* '())
  ;(defparameter *run-sim* t)
  (defparameter *load-mps*
    (let* ((mps (cl-mpm:sim-mps *sim*))
           (least-pos
              (apply #'max (loop for mp across mps
                                 collect (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)))))
      (loop for mp across mps when (>= (magicl:tref (cl-mpm/particle:mp-position mp) 0 0) (- least-pos 0.001))
              collect mp)))
  (with-accessors ((mps cl-mpm:sim-mps))
      *sim*
    (let ((x-min (loop for mp across mps
                       minimize (magicl:tref
                                 (cl-mpm/particle:mp-position mp)
                                 0 0))))
      (defparameter *far-field-mps*
        (loop for mp across mps
              when (= x-min (magicl:tref
                         (cl-mpm/particle:mp-position mp)
                         0 0))
                collect mp)))
    (let ((x-max (loop for mp across mps
                       maximize (magicl:tref
                                 (cl-mpm/particle:mp-position mp)
                                 0 0))))
      (defparameter *terminus-mps*
        (loop for mp across mps
              when (= x-max (magicl:tref
                             (cl-mpm/particle:mp-position mp)
                             0 0))
                collect mp)))
    ;; (increase-load *sim* *terminus-mps*
    ;;                (magicl:from-list (list (* 1d4) 0d0) '(2 1)))
    )
  ;; (defparameter *dist-mp* (nth 0 *terminus-mps*))
  ;; (increase-load *sim* *load-mps* 1)
  ;; (increase-load *sim* *load-mps* 100)
  )





(defparameter *run-sim* nil)

(defun get-crack-waterlevel (fill-percent)
  (let ((min-damage 0.4d0))
    (let ((min *original-crack-height*))
      ;; (cl-mpm::iterate-over-nodes-serial
      ;;  (cl-mpm:sim-mesh *sim*)
      ;;  (lambda (node)
      ;;    (let ((ypos (magicl:tref (cl-mpm/mesh::node-position node) 1 0)))
      ;;      (when (and
      ;;             (< (magicl:tref (cl-mpm/mesh::node-position node) 0 0) *crack-water-width*)
      ;;             ;; (> (cl-mpm/mesh::node-damage node) 0.4d0)
      ;;             (< (cl-mpm/mesh::node-volume node) (* 0.9d0 (cl-mpm/mesh::node-volume-true node)))
      ;;             ;; (cl-mpm/mesh::node-boundary-node node)
      ;;             )
      ;;        (setf min (min ypos min))
      ;;        ))
      ;;    ))
      (loop for mp across (cl-mpm:sim-mps *sim*)
            when
            (and
             (< (magicl:tref (cl-mpm/particle::mp-position mp) 0 0) *crack-water-width*)
             ;; (>= (magicl:tref (cl-mpm/particle:mp-position mp) 0 0) (* (- 0.5d0 crack-width) *ice-length*))
             ;; (<= (magicl:tref (cl-mpm/particle:mp-position mp) 0 0) (* (+ 0.5d0 crack-width) *ice-length*))
             (> (cl-mpm/particle::mp-damage mp) min-damage)
             )
            do (setf min (min (magicl:tref (cl-mpm/particle:mp-position mp) 1 0) min)))
      (setf min (min *last-min* min))
      (setf *last-min* min)

      (with-accessors ((mesh cl-mpm:sim-mesh))
          *sim*
        (let ((h (cl-mpm/mesh::mesh-resolution mesh)))
          (let* ((crack-depth (- *ice-height* min))
                 (v (+ min
                       (* fill-percent
                          crack-depth
                          )))
                 ;;Normalise to mesh
                 ;(v (* (round v h) h))
                 )
            ;; (print v)
            (setf *crack-water-height* v)
            (setf (cl-mpm/buoyancy::bc-buoyancy-datum *crack-water-bc*) v)
            ;; v
            (- *ice-height* crack-depth)
            ))))))
(defun run ()
  (cl-mpm/output:save-vtk-mesh (merge-pathnames "output/mesh.vtk")
                          *sim*)
  (defparameter *run-sim* t)
  (let* ((target-time 1d0)
         (dt-scale 1d0)
         (substeps (floor target-time (cl-mpm::sim-dt *sim*))))

    (cl-mpm::update-sim *sim*)
    (time (let* ((dt-e (* dt-scale (cl-mpm::calculate-min-dt *sim*)))
                 (substeps-e (floor target-time dt-e)))
            (format t "CFL dt estimate: ~f~%" dt-e)
            (format t "CFL step count estimate: ~D~%" substeps-e)
            (setf (cl-mpm:sim-dt *sim*) dt-e)
            (setf substeps substeps-e)))
    (format t "Substeps ~D~%" substeps)
    (time (loop for steps from 0 to 30
                while *run-sim*
                do
                   (progn
                     (let ((base-damping 1d1))
                       ;; (let ((damping-step 10)
                       ;;       (damping-final 1d-3)
                       ;;       (damping-speed 50))
                       ;;   (when (> steps damping-step)
                       ;;     (setf (cl-mpm::sim-damping-factor *sim*)
                       ;;           (max
                       ;;            (* (exp (/ (- damping-step steps) 5)) base-damping)
                       ;;            1d-3
                       ;;            )
                       ;;           ;; dt-scale 0.1d0
                       ;;           ;; target-time 1d-2
                       ;;           )
                       ;;     ))
                       (when (= steps 3)
                         (progn
                           (setf (cl-mpm::sim-enable-damage *sim*) t)
                           (setf (cl-mpm::sim-damping-factor *sim*)
                                 0d0)
                           ;;       ;; 1d0
                           ;;       ;; base-damping
                           ;;       ;; (* base-damping 0.1)
                           ;;       ;; dt-scale 0.1d0
                           ;;       ;; target-time 1d-2
                           ;;       )
                           (let* ((crack-width ;; 0.1d0
                                        (/ 20d0 *ice-length*)
                                               )
                                  (init-stress
                                   (loop for mp across (cl-mpm:sim-mps *sim*)
                                         when
                                         (and
                                          (>= (magicl:tref (cl-mpm/particle:mp-position mp) 0 0) (* (- 0.5d0 crack-width) *ice-length*))
                                          (<= (magicl:tref (cl-mpm/particle:mp-position mp) 0 0) (* (+ 0.5d0 crack-width) *ice-length*)))
                                         maximize (cl-mpm/particle::mp-damage-ybar mp)
                                          )
                                         )
                                  (init-stress-reduced (* 0.10d0 init-stress))
                                  (init-stress-reduced 0.1d6)
                                  )
                             (format t "Bounds ~F - ~F~%" (* (- 0.5d0 crack-width) *ice-length*) (* (+ 0.5d0 crack-width) *ice-length*))
                             (format t "Init stress found at :~F MPa~%" (* 1d-6 init-stress-reduced))
                           (loop for mp across (cl-mpm:sim-mps *sim*)
                                 when
                                 (and
                                  (>= (magicl:tref (cl-mpm/particle:mp-position mp) 0 0) (* (- 0.5d0 crack-width) *ice-length*))
                                  (<= (magicl:tref (cl-mpm/particle:mp-position mp) 0 0) (* (+ 0.5d0 crack-width) *ice-length*))
                                  ;; (<= (magicl:tref (cl-mpm/particle:mp-position mp) 1 0) (+ *original-crack-height* 0))
                                  )
                                 do (setf  (cl-mpm/particle::mp-damage-rate mp) 1d-2
                                           (cl-mpm/particle::mp-initiation-stress mp) init-stress-reduced
                                           ))
                           )))
                       (when (= steps 0)
                         (progn
                           (setf (cl-mpm::sim-enable-damage *sim*) nil
                                 (cl-mpm::sim-damping-factor *sim*)
                                 base-damping
                                 ;; (+ (* 1d0 (cl-mpm::sim-mass-scale *sim*) ;;       (exp (- steps))
                                 ;;       )
                                 ;;    base-damping)
                                 )
                           ))
                       )
                     (let* ((dt-e (* dt-scale (cl-mpm::calculate-min-dt *sim*)))
                            (substeps-e (floor target-time dt-e)))
                       (format t "CFL dt estimate: ~f~%" dt-e)
                       (format t "CFL step count estimate: ~D~%" substeps-e)
                       (setf (cl-mpm:sim-dt *sim*) dt-e)
                       (setf substeps substeps-e))
                     (format t "Damping ~F~%" (cl-mpm::sim-damping-factor *sim*))
                     ;; (when (> steps 2)
                     ;;   (setf dt-scale 1d0)
                     ;;     )
                     (format t "Step ~d ~%" steps)
                     (when nil;*debug*
                       (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
                       (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_nodes_~5,'0d.vtk" *sim-step*)) *sim*))
                     ;; (cl-mpm/output:save-csv (merge-pathnames (format nil "output/sim_~5,'0d.csv" *sim-step*)) *sim*)

                     (push *t* *time*)
                     (push
                      (loop for mp across (cl-mpm:sim-mps *sim*)
                            maximize (max-stress mp))
                      *max-stress*)
                     ;; (let ((cfl (find-max-cfl *sim*)))
                     ;;   (push cfl *cfl-max*))
                     (setf *x*
                           (loop for mp across (cl-mpm:sim-mps *sim*)
                                 maximize (magicl:tref (cl-mpm/particle::mp-displacement mp) 0 0)))
                     (push
                      *x*
                      *x-pos*)
                     (let ((cfl 0))
                       (time (dotimes (i substeps)
                               (cl-mpm::update-sim *sim*)
                               (when *crack-water-bc*
                                 (defparameter *crack-depth* (get-crack-waterlevel *meltwater-fill*)))
                               (setf *t* (+ *t* (cl-mpm::sim-dt *sim*)))))

                       (format t "Crack depth: ~F~%" *crack-depth*)
                       (format t "Crack depth %: ~F~%" (- 1 (/ *crack-depth* *ice-height*)))
                       (format t "CFL: ~f~%" cfl)
                       (push cfl *cfl-max*)
                       (let* ((dt-e (* dt-scale (cl-mpm::calculate-min-dt *sim*)))
                              (substeps-e (floor target-time dt-e)))
                         (format t "CFL dt estimate: ~f~%" dt-e)
                         (format t "CFL step count estimate: ~D~%" substeps-e)
                         (setf (cl-mpm:sim-dt *sim*) dt-e)
                         (setf substeps substeps-e))
                         )
                     (incf *sim-step*)
                     ;(when *debug*
                     ;  (plot *sim*)
                     ;  (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" *sim-step*))
                     ;                     :terminal "png size 1920,1080"
                     ;                     )
                     ;  (swank.live:update-swank)
                     ;  (sleep .01))
                     ))))
  (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
  (when *crack-water-bc*
    (defparameter *crack-depth* (get-crack-waterlevel *meltwater-fill*)))
  (format t "Crack depth: ~F~%" *crack-depth*)
  (format t "Crack depth %: ~F~%" (- 1 (/ *crack-depth* *ice-height*)))
  )


;(setf lparallel:*kernel* (lparallel:make-kernel 32 :name "custom-kernel"))
;; (setf lparallel:*kernel* (lparallel:make-kernel 8 :name "custom-kernel"))


(defun test-bounds ()
  (let ((mesh (cl-mpm:sim-mesh *sim*)))
    (time-form
     1000000
     ;; (dotimes (i 100000))
     (progn
       (cl-mpm/mesh:in-bounds mesh '(0 0))))))

(defun test-all-meltwater ()
  (loop for mw from 0.0d0 to 1d0 by 0.1d0
        for i from 0
        do (progn
             (setup)
             (defparameter *meltwater-fill* mw)
             (format t "MP count:~D~%" (length (cl-mpm:sim-mps *sim*)))
             (run)
             (format t "Ocean: ~F - Meltwater: ~F - Crack depth %: ~F~%" *ocean-fill* mw (- 1 (/ *crack-depth* *ice-height*)))
             (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" i)) *sim*)
             )))

(defparameter *debug* nil)
 (if *debug*
   (progn
     ;(setf lparallel:*kernel* (lparallel:make-kernel 8 :name "custom-kernel"))
     (setf lparallel:*kernel* (lparallel:make-kernel 32 :name "custom-kernel"))
     (defparameter *run-sim* nil)
     (setup)
     (format t "MP count:~D~%" (length (cl-mpm:sim-mps *sim*)))
     (run))
   (progn
     ;(setf lparallel:*kernel* (lparallel:make-kernel 32 :name "custom-kernel"))
     (setf lparallel:*kernel* (lparallel:make-kernel 32 :name "custom-kernel"))
     (test-all-meltwater)
     ;(setup)
     ;(format t "MP count:~D~%" (length (cl-mpm:sim-mps *sim*)))
     ;(run)
     ))
