#!/bin/bash
# checked with shellcheck
# formatted with shfmtv3.36.0

action="run"

basedir=basedir=$(git rev-parse --show-toplevel)

export PLUMED_STACK_TRACE=1
export PLUMED_EXCEPTIONS_DEBUG=1

dirs=(
  regtest/symfunc/rt20
  regtest/symfunc/rt20-mpi
  regtest/symfunc/rt24
  regtest/symfunc/rt-atomic-smac
  regtest/symfunc/rt-averaged-q6
  regtest/symfunc/rt-averaged-q6-spAspB
  regtest/symfunc/rt-average-vec
  regtest/symfunc/rt-averaging-test
  # regtest/symfunc/rt-behler
  # regtest/symfunc/rt-behler2
  regtest/symfunc/rt-center
  regtest/symfunc/rt-center2
  regtest/symfunc/rt-center3
  regtest/symfunc/rt-cmap-times-stack
  regtest/symfunc/rt-combine
  regtest/symfunc/rt-compose-vector-biased
  regtest/symfunc/rt-coordination-powers
  regtest/symfunc/rt-coordination-powers-multi
  regtest/symfunc/rt-dumpmulti
  regtest/symfunc/rt-hexactic
  regtest/symfunc/rt-localav-test
  regtest/symfunc/rt-local-crystal
  regtest/symfunc/rt-matrix-tricks
  regtest/symfunc/rt-nbonds-one-q6
  regtest/symfunc/rt-nbonds-q6
  regtest/symfunc/rt-new-tetra
  regtest/symfunc/rt-ntwo-lq6
  regtest/symfunc/rt-numder
  regtest/symfunc/rt-parse-only
  regtest/symfunc/rt-product
  regtest/symfunc/rt-q6
  regtest/symfunc/rt-q6-nochain
  regtest/symfunc/rt-q6-subset
  regtest/symfunc/rt-simple-center
  regtest/symfunc/rt-simple-center2
  regtest/symfunc/rt-simplecubic
  regtest/symfunc/rt-single-groupb
  regtest/symfunc/rt-smac
  regtest/symfunc/rt-small-lq6
  regtest/symfunc/rt-small-q6
  regtest/symfunc/rt-splane
  # regtest/symfunc/rt-structure-factor
  regtest/symfunc/rt-two-average-q1
  regtest/symfunc/rt-two-coord-stack-mult-q1
  regtest/symfunc/rt-two-coord-stack-mult-q1-av
  regtest/symfunc/rt-two-coord-stack-mult-q1-withnorm
  regtest/symfunc/rt-two-coord-stack-mult-q1-withnorm-matsum
  regtest/symfunc/rt-two-filter
  regtest/symfunc/rt-two-lq1
  regtest/symfunc/rt-two-lq6
  regtest/symfunc/rt-two-mols
  regtest/symfunc/rt-two-nonorm
  regtest/symfunc/rt-two-norm
  regtest/symfunc/rt-two-q6
  regtest/symfunc/rt-urea

  regtest/targetGPU/rt-symfunc-averaged-q6
  regtest/targetGPU/rt-symfunc-averaged-q6-spAspB
  regtest/targetGPU/rt-symfunc-nbonds-one-q6
  regtest/targetGPU/rt-symfunc-nbonds-q6
  regtest/targetGPU/rt-symfunc-ntwo-lq6
  regtest/targetGPU/rt-symfunc-q6
  regtest/targetGPU/rt-symfunc-q6-nochain
  regtest/targetGPU/rt-symfunc-q6-subset
  regtest/targetGPU/rt-symfunc-small-lq6
  regtest/targetGPU/rt-symfunc-small-q6
  regtest/targetGPU/rt-symfunc-two-lq6
  regtest/targetGPU/rt-symfunc-two-q6
  regtest/targetGPU/rt-symfunc-simplecubic

  regtest/targetGPU/rt-symfunc-spherical-harmonics

  # regtest/adjmat/rt-eigvecs
  # regtest/adjmat/rt-join-eigvecs
  # regtest/adjmat/rt-waterbridge2
  # regtest/basic/rt-combine-periodic
  # regtest/basic/rt-drmsd
  # regtest/basic/rt-fix-229
  # regtest/basic/rt-fix-235
  # regtest/basic/rt-ghbfix
  # regtest/basic/rt-lepton-asmjit
  # regtest/basic/rt-lepton
  # regtest/basic/rt-maxent-3
  # regtest/basic/rt-molfile-1
  # regtest/basic/rt-molfile-2
  # regtest/basic/rt-multi-1
  # regtest/basic/rt-multirmsd
  # regtest/basic/rt0
  # regtest/basic/rt11b
  # regtest/basic/rt11c
  # regtest/basic/rt12
  # regtest/basic/rt16
  # regtest/basic/rt19-mpi
  # regtest/basic/rt19
  # regtest/basic/rt26
  # regtest/basic/rt4
  # regtest/basic/rt42
  # regtest/basic/rt42d
  # regtest/basic/rt63b
  # regtest/basic/rt71
  # regtest/basic/rt71b
  # regtest/gridtools/rt-grad-2
  # regtest/gridtools/rt-kldiv
  # regtest/gridtools/rt-spherical-integral-2
  # regtest/gridtools/rt-spherical-integral
  # regtest/gridtools/rt-uweights-integral
  # regtest/gridtools/rt-weights-integral
  # regtest/isdb/rt-noe-mpi
  # regtest/landmarks/rt-reweight-temp-press
  # regtest/mapping/rt-pca-multi
  # regtest/mapping/rt-rmsd-vector
  # regtest/secondarystructure/rt-coverage-antibetarmsd
  # regtest/secondarystructure/rt-coverage-parabetarmsd-parabeta
  # regtest/secondarystructure/rt33-mpi
  # regtest/secondarystructure/rt33
  # regtest/symfunc/rt-combine
  # regtest/symfunc/rt-simple-center2
  # regtest/symfunc/rt-splane
  # regtest/symfunc/rt-urea
  # regtest/targetGPU/rt-adjmat-eigvecs
  # regtest/targetGPU/rt-adjmat-join-eigvecs
  # regtest/ves/rt-bf-combined
  # regtest/ves/rt-le-1d-chebyshev-welltempered
  # regtest/ves/rt-le-1d-legendre-scaled-welltempered
  # regtest/ves/rt-le-1d-legendre-welltempered
  # regtest/ves/rt-le-1d-outside
  # regtest/ves/rt-le-1d-wavelets-optimization-threshold
  # regtest/ves/rt-le-2d-legendre-scaled-welltempered
  # regtest/ves/rt-le-2d-legendre-welltempered
  # regtest/ves/rt-le-2d-mixedbasis-welltempered
  # regtest/ves/rt-le-2d-outside-mwalkers
  # regtest/ves/rt-le-2d-outside
  # regtest/ves/rt-opt-averaged-sgd-exp-decaying-aver-2
  # regtest/ves/rt-opt-averaged-sgd-exp-decaying-aver
  # regtest/ves/rt-opt-averaged-sgd-monitor-aver-grad-exp-decay
  # regtest/ves/rt-opt-averaged-sgd-monitor-aver-grad
  # regtest/ves/rt-opt-averaged-sgd-monitor-inst-grad
  # regtest/ves/rt-opt-averaged-sgd
  # regtest/volumes/rt-sean-marks

  # regtest/adjmat/rt45
  # regtest/adjmat/rt-adjmat-prod
  # regtest/adjmat/rt-basic-matrix
  # regtest/adjmat/rt-basic-matrix2
  # regtest/adjmat/rt-basic-matrix3
  # regtest/adjmat/rt-basic-matrix4
  # regtest/adjmat/rt-eigvals
  # regtest/adjmat/rt-eigvecs

  # regtest/adjmat/rt-hbond
  # regtest/adjmat/rt-join-eigvals
  # regtest/adjmat/rt-join-eigvecs
  # # regtest/adjmat/rt-link
  # # regtest/adjmat/rt-link-nopbc
  # regtest/adjmat/rt-matmult
  # regtest/adjmat/rt-matrix-func
  # regtest/adjmat/rt-matrixsum
  # regtest/adjmat/rt-perovskite-nucleation
  # regtest/adjmat/rt-sprint
  # regtest/adjmat/rt-topology-derivatives
  # regtest/adjmat/rt-topology-newder
  # regtest/adjmat/rt-topology-newder2
  # regtest/adjmat/rt-topology-symmetric
  # regtest/adjmat/rt-waterbridge
  # regtest/adjmat/rt-waterbridge2
  # regtest/adjmat/rt-waterbridge3

  # regtest/targetGPU/rt-GPU-quaternion-bond-product
  # regtest/targetGPU/rt-GPU-quaternion-bond-product-deriv
  # regtest/targetGPU/rt-GPU-quaternion-bond-product-only

  regtest/targetGPU/rt-multicolvar-filters

  regtest/targetGPU/rt-adjmat-basic-matrix
  regtest/targetGPU/rt-adjmat-basic-matrix2
  regtest/targetGPU/rt-adjmat-basic-matrix3
  regtest/targetGPU/rt-adjmat-basic-matrix4
  regtest/targetGPU/rt-adjmat-eigvals
  regtest/targetGPU/rt-adjmat-eigvecs
  regtest/targetGPU/rt-adjmat-join-eigvals
  regtest/targetGPU/rt-adjmat-join-eigvecs
  regtest/targetGPU/rt-adjmat-matmult
  regtest/targetGPU/rt-adjmat-prod

  regtest/targetGPU/rt-landmarks-calc-dissims

  regtest/targetGPU/rt-coordination-insphere
  regtest/targetGPU/rt-coordination-insphere-nochain
  regtest/targetGPU/rt-cylinder
  regtest/targetGPU/rt-GPU-DISTANCE
  regtest/targetGPU/rt-GPU-DISTANCE2
  regtest/targetGPU/rt-GPU-multicolvar-print
  regtest/targetGPU/rt-GPU-PLANE-bias

  regtest/targetGPU/rt-GPU-secondarystructure
  regtest/targetGPU/rt-GPU-secondarystructure-DRMSD
  regtest/targetGPU/rt-GPU-secondarystructure-DRMSD-play
  regtest/targetGPU/rt-GPU-secondarystructure-DRMSD-play-0
  regtest/targetGPU/rt-GPU-secondarystructure-DRMSD-play-1
  regtest/targetGPU/rt-GPU-secondarystructure-DRMSD-play-2
  regtest/targetGPU/rt-GPU-secondarystructure-DRMSD-play-3
  regtest/targetGPU/rt-GPU-secondarystructure-DRMSD-play-4
  regtest/targetGPU/rt-GPU-secondarystructure-DRMSD-play-5
  regtest/targetGPU/rt-GPU-secondarystructure-DRMSD-play-6
  regtest/targetGPU/rt-GPU-secondarystructure-DRMSD-play-7
  regtest/targetGPU/rt-GPU-TORSION-biasIsSUM
  regtest/targetGPU/rt-volume-around
  # regtest/basic/rt26-high
)

failures=()
if [[ $action == "run" ]]; then
  for dir in "${dirs[@]}"; do
    echo -n "> Running $dir: "
    if (
      cd "$dir" || exit
      make
      grep -vzq FAILURE report.txt
    ) >>/dev/null; then
      echo -e "\e[32mSUCCESS\e[0m"
      # tail "$dir/tmp/out" | grep 'Cycles        Total      Average      Minimum      MaximumCalculating'
      # tail "$dir/tmp/out" | grep Calculating -A1
    else
      echo -e "\e[31mFAILURE\e[0m"
      grep 'Diff for ' "${dir}/report.txt"
      failures[${#failures[@]}]="$dir"
    fi
    # {
    #     grep 'v:' "${dir}"
    #     # grep '4 Calculating' "${dir}/tmp/out"
    #     # grep '5 Applying' "${dir}/tmp/out"
    # }
  done
fi
N=${#dirs[@]}
NF=${#failures[@]}
successes=$((N - NF))
echo -e "\e[32mSUCCESSES\e[0m: ${successes}/${N}"

# for dir in "${failures[@]}"; do
#     echo -e "\e[31mFAILED\e[0m: $dir"
# done
