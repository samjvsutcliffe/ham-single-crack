(restrict-compiler-policy 'speed 3 3)
(restrict-compiler-policy 'debug 0 0)
(restrict-compiler-policy 'safety 0 0)
(setf *block-compile-default* t)
;(push :magicl.use-mkl *features*)
;(pushnew :cl-mpm-fbar *features*)
(setf *features* (delete :cl-mpm-pic *features*))
(ql:quickload "magicl")
(ql:quickload "cl-mpm")
(asdf:compile-system :cl-mpm :force T)
(ql:quickload "cl-mpm/examples/single-crack")
(asdf:compile-system :cl-mpm/examples/single-crack :force T)
(in-package :cl-mpm/examples/single-crack)

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
                :E 1d9
                :nu 0.3250d0
                ;; :visc-factor 111d6
                ;; :visc-power 3d0

                :initiation-stress 0.1d6
                :damage-rate 0d0
                :critical-damage 0.50d0
                :local-length 50d0
                :local-length-damaged 50d0
                ;; :local-length-damaged 20d0
                ;; :local-length-damaged 0.1d0
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
      (let ((mass-scale 1d0))
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
      (setf (cl-mpm::sim-mp-damage-removal-instant sim) t)
      (setf (cl-mpm:sim-dt sim) 1d-4)
      (setf (cl-mpm:sim-bcs sim) (make-array 0))
      (setf (cl-mpm:sim-bcs sim)
            (cl-mpm/bc::make-outside-bc-var
             (cl-mpm:sim-mesh sim)
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(0 nil)))
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(0 nil)))
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(nil 0)))
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(nil 0)))
             ))
      (format t "Bottom level ~F~%" h-y)
      (let* ((terminus-size (+ (second block-size) (* 0d0 (first block-size))))
             (ocean-x 1000)
            ;; (ocean-y (+ h-y (* 0.90d0 0.0d0 terminus-size)))
             (ocean-y (* 0.0d0 terminus-size))
            ;(angle -1d0)
            )

        (format t "Ocean level ~a~%" ocean-y)
        (defparameter *water-height* ocean-y)
        (defparameter *meltwater-fill* 0.50d0)
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
           *water-density*
           (lambda (pos datum)
             (and
              (< (magicl:tref pos 0 0) *crack-water-width*)
              ;; (< (magicl:tref pos 1 0) datum)
              ))))

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
                     )))
                 ;; (cl-mpm/buoyancy::make-bc-buoyancy
                 ;;  sim
                 ;;  ocean-y
                 ;;  *water-density*
                 ;;  )
                 *crack-water-bc*
                 ))
               ;; (cl-mpm/bc:make-bcs-from-list
               ;;  (list *floor-bc*)
               ;;  )
               ))
        ;; (setf (cl-mpm::sim-bcs-force-list sim)
        ;;       (list
        ;;        (cl-mpm/bc:make-bcs-from-list
        ;;         (list
        ;;          (cl-mpm/buoyancy::make-bc-buoyancy
        ;;           sim
        ;;           ocean-y
        ;;           *water-density*
        ;;           )
        ;;          ))
        ;;        (cl-mpm/bc:make-bcs-from-list
        ;;         (list *floor-bc*)
        ;;         )
        ;;        ))
        )
      (let ((normal (magicl:from-list (list (sin (- (* pi (/ angle 180d0))))
                                            (cos (+ (* pi (/ angle 180d0))))) '(2 1))))
        (defparameter *sliding-slope* 0d0)
        (defparameter *sliding-offset* (- h-y (* 0d0 0d0))))
      sim)))

(defparameter *ice-density* 900)
(defparameter *water-density* 1000)
;; (defparameter *ice-density* 900)
;; (defparameter *water-density* 1000)
;Setup
(defun setup ()
  (declare (optimize (speed 0)))
  (defparameter *run-sim* nil)
  (let* ((mesh-size 50)
         (mps-per-cell 2)
         (shelf-height 120)
         (shelf-length 600)
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
         10
         cut-depth
         )))
      (defparameter *ice-height* shelf-height)
      (defparameter *ice-length* shelf-length)
      (defparameter *original-crack-height* (- shelf-height cut-depth))
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
  (defparameter *dist-mp* (nth 0 *terminus-mps*))
  ;; (increase-load *sim* *load-mps* 1)
  ;; (increase-load *sim* *load-mps* 100)
  )
(defparameter *water-height* 0d0)





(defparameter *run-sim* nil)

(defun get-crack-waterlevel (fill-percent)
  (let ((min-damage 0.5d0))
    (let ((min *original-crack-height*))
      (cl-mpm::iterate-over-nodes-serial
       (cl-mpm:sim-mesh *sim*)
       (lambda (node)
         (let ((ypos (magicl:tref (cl-mpm/mesh::node-position node) 1 0)))
           (when (and
                  (< (magicl:tref (cl-mpm/mesh::node-position node) 0 0) *crack-water-width*)
                  ;; (> (cl-mpm/mesh::node-damage node) 0.4d0)
                  (< (cl-mpm/mesh::node-volume node) (* 0.9d0 (cl-mpm/mesh::node-volume-true node)))
                  ;; (cl-mpm/mesh::node-boundary-node node)
                  )
             (setf min (min ypos min))
             ))
         ))
      (let ((v (+ min
                  (* fill-percent
                     (- *ice-height* min)))))
        ;; (print v)
        (setf *crack-water-height* v)
        (setf (cl-mpm/buoyancy::bc-buoyancy-datum *crack-water-bc*) v)
        v
        ))))
(defun run ()
  (cl-mpm/output:save-vtk-mesh (merge-pathnames "output/mesh.vtk")
                          *sim*)
  (defparameter *run-sim* t)
  (let* ((target-time 0.1d0)
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
    (time (loop for steps from 0 to 100
                while *run-sim*
                do
                   (progn
                     (let ((base-damping 1d0))
                       (when (= steps 30)
                         (progn
                           (setf (cl-mpm::sim-enable-damage *sim*) t)
                           (setf (cl-mpm::sim-damping-factor *sim*) base-damping
                                 ;; dt-scale 0.1d0
                                 ;; target-time 1d-2
                                 )
                           (let* ((crack-width (/ 50d0 *ice-length*))
                                  (init-stress
                                   (loop for mp across (cl-mpm:sim-mps *sim*)
                                         when
                                         (and
                                          (>= (magicl:tref (cl-mpm/particle:mp-position mp) 0 0) (* (- 0.5d0 crack-width) *ice-length*))
                                          (<= (magicl:tref (cl-mpm/particle:mp-position mp) 0 0) (* (+ 0.5d0 crack-width) *ice-length*)))
                                         maximize (cl-mpm/particle::mp-damage-ybar mp)
                                          )
                                         )
                                  (init-stress-reduced (* 0.1d0 init-stress))
                                  )
                             (format t "Bounds ~F - ~F~%" (* (- 0.5d0 crack-width) *ice-length*) (* (+ 0.5d0 crack-width) *ice-length*))
                             (format t "Init stress found at :~F MPa~%" (* 1d-6 init-stress-reduced))
                           (loop for mp across (cl-mpm:sim-mps *sim*)
                                 when
                                 (and
                                  (>= (magicl:tref (cl-mpm/particle:mp-position mp) 0 0) (* (- 0.5d0 crack-width) *ice-length*))
                                  (<= (magicl:tref (cl-mpm/particle:mp-position mp) 0 0) (* (+ 0.5d0 crack-width) *ice-length*))
                                  )
                                 do (setf  (cl-mpm/particle::mp-damage-rate mp) 1d-1
                                           (cl-mpm/particle::mp-initiation-stress mp) init-stress-reduced
                                           ))
                           )))
                       (when (= steps 0)
                         (progn
                           (setf (cl-mpm::sim-enable-damage *sim*) t
                                 (cl-mpm::sim-damping-factor *sim*)
                                 base-damping
                                 ;; (+ (* 1d0 (cl-mpm::sim-mass-scale *sim*)
                                 ;;       (exp (- steps))
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
                     (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
                     (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_nodes_~5,'0d.vtk" *sim-step*)) *sim*)
                     (cl-mpm/output:save-csv (merge-pathnames (format nil "output/sim_~5,'0d.csv" *sim-step*)) *sim*)

                     (push *t* *time*)
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
                               (defparameter *crack-depth* 0d0)
                               (when *crack-water-bc*
                                 (defparameter *crack-depth* (get-crack-waterlevel *meltwater-fill*)))
                               (setf *t* (+ *t* (cl-mpm::sim-dt *sim*)))))

                       (format t "Crack depth: ~F~%" *crack-depth*)
                       (format t "Crack depth %: ~F~%" (- 1 (/ *crack-depth* *ice-height*)))
                       (format t "CFL: ~f~%" cfl)
                       (push cfl *cfl-max*)
                       ;; (let* ((dt-e (* dt-scale (cl-mpm::calculate-min-dt *sim*)))
                       ;;        (substeps-e (floor target-time dt-e)))
                       ;;   (format t "CFL dt estimate: ~f~%" dt-e)
                       ;;   (format t "CFL step count estimate: ~D~%" substeps-e)
                       ;;   (setf (cl-mpm:sim-dt *sim*) dt-e)
                       ;;   (setf substeps substeps-e))
                         )
                     (incf *sim-step*)
                     (cl-mpm/examples/single-crack::plot *sim*)
                     (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" *sim-step*))
                                        :terminal "png size 1920,1080"
                                        )
                     (swank.live:update-swank)
                     (sleep .01)
                     ))))
  (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
  )


;(setf lparallel:*kernel* (lparallel:make-kernel 32 :name "custom-kernel"))
(setf lparallel:*kernel* (lparallel:make-kernel 4 :name "custom-kernel"))
(require :sb-sprof)
(defparameter *run-sim* nil)
(setup)
(format t "MP count:~D~%" (length (cl-mpm:sim-mps *sim*)))
(run)

