;; This file defines a guix package.
;; You can use it to launch an interactive development shell with
;; ```
;; guix shell
;; ```
;; or to build the package with
;; ```
;; guix build -f guix.scm
;; ```
(define-module (gnu packages dune-dpg)
  #:use-module (srfi srfi-1)
  #:use-module (guix gexp)
  #:use-module ((guix licenses) #:prefix license:)
  #:use-module (guix packages)
  #:use-module (guix git-download)
  #:use-module (guix build-system cmake)
  #:use-module (guix download)
  #:use-module (gnu packages)
  #:use-module (gnu packages algebra)
  #:use-module (gnu packages boost)
  #:use-module (gnu packages gcc)
  #:use-module (gnu packages maths)
  #:use-module (gnu packages multiprecision)
  #:use-module (gnu packages pkg-config)
  #:use-module (gnu packages python))

(define %source-dir (dirname (current-filename)))

(define-public dune-dpg
  (package
    (name "dune-dpg")
    (version "git")
    (source (local-file %source-dir
                        #:recursive? #t
                        #:select? (git-predicate %source-dir)))
    (build-system cmake-build-system)
    (arguments
     `(#:phases
       (modify-phases %standard-phases
         (add-after 'build 'build-tests
           (lambda* (#:key make-flags parallel-build? #:allow-other-keys)
             (apply invoke "make" "build_tests"
                    `(,@(if parallel-build?
                            `("-j" ,(number->string (parallel-job-count)))
                            '())
                      ,@make-flags)))))))
    (inputs
     (list dune-common
           dune-geometry
           dune-istl
           dune-localfunctions
           dune-functions
           dune-typetree
           dune-grid
           dune-subgrid
           boost
           eigen
           python
           suitesparse
           metis
           openblas
           gmp))
    (native-inputs
     (list gfortran pkg-config))
    (home-page "https://gitlab.dune-project.org/felix.gruber/dune-dpg")
    (synopsis "Discontinuous Petrov–Galerkin finite element solver")
    (description "The dune-dpg library allows to solve partial differential
equations using the Discontinuous Petrov–Galerkin finite element method.")
    ;; GPL version 2 with "runtime exception" to make it behave like LGPLv2.
    (license license:gpl2+)))

;; Return the dune-dpg package to make it easier to use with guix shell.
dune-dpg
