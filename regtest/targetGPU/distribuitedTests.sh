#!/bin/bash
# checked with shellcheck
# formatted with shfmtv3.36.0

action="run"
initialize="false"
USEGPU="USEGPU"
basedir="rt-GPU-secondarystructure-DRMSD-"
if [[ $1 = "prepare" ]]; then
    action="prepare"
elif [[ $1 = "reset" ]]; then
    action="prepare"
    initialize="true"
    USEGPU=""
elif [[ $1 = "clean" ]]; then
    action="clean"
fi

if [[ $action = "prepare" ]]; then
    inputLine=()
    inputLineNoBias=()
    inputLine[${#inputLine[@]}]="ALPHARMSD RESIDUES=ALL TYPE=DRMSD LESS_THAN={RATIONAL R_0=0.08 NN=8 MM=12 NOSTRETCH} $USEGPU"
    inputLine[${#inputLine[@]}]="ANTIBETARMSD RESIDUES=all TYPE=DRMSD STRANDS_CUTOFF=1.0 LESS_THAN={RATIONAL R_0=0.08 NN=8 MM=12 NOSTRETCH} $USEGPU"
    inputLine[${#inputLine[@]}]="PARABETARMSD RESIDUES=all TYPE=DRMSD STRANDS_CUTOFF=1.0 LESS_THAN={RATIONAL R_0=0.08 NN=8 MM=12 NOSTRETCH} $USEGPU"
    inputLine[${#inputLine[@]}]="ANTIBETARMSD RESIDUES=3-5,8-10 TYPE=DRMSD STYLE=inter LESS_THAN={RATIONAL R_0=0.08 NN=8 MM=12 NOSTRETCH} $USEGPU"
    inputLine[${#inputLine[@]}]="PARABETARMSD RESIDUES=all TYPE=DRMSD LESS_THAN={RATIONAL R_0=0.08 NN=8 MM=12} $USEGPU"

    inputLineNoBias[${#inputLineNoBias[@]}]="ALPHARMSD RESIDUES=2-7 TYPE=DRMSD R_0=0.08 NN=8 MM=12 $USEGPU"
    inputLineNoBias[${#inputLineNoBias[@]}]="ALPHARMSD RESIDUES=ALL TYPE=DRMSD $USEGPU"
    inputLineNoBias[${#inputLineNoBias[@]}]="PARABETARMSD RESIDUES=all TYPE=DRMSD NN=8 MM=12 R_0=0.08 $USEGPU"

    i=0
    for line in "${inputLine[@]}"; do
        dir=${basedir}$i
        mkdir -pv "$dir"

        cat <<EOF >"$dir"/plumed.dat
MOLINFO STRUCTURE=helix.pdb
v: $line
PRINT ARG=v.* STRIDE=1 FILE=colvar FMT=%8.4f
RESTRAINT ARG=v.lessthan AT=2 KAPPA=3 SLOPE=3
EOF
        if [[ $initialize = "true" ]]; then
            (
                cd "$dir" || exit
                touch forces.reference
                touch colvar.reference
                cp ../rt-GPU-secondarystructure-DRMSD/*.pdb .
                cp ../rt-GPU-secondarystructure-DRMSD/*.xyz .
                echo "include ../../scripts/test.make" >Makefile

                cat <<EOF >config
type=driver
arg="--plumed plumed.dat --trajectory-stride 10 --timestep 0.005 --ixyz ala12_trajectory.xyz --dump-forces forces --dump-forces-fmt=%10.6f"

#export NVCOMPILER_ACC_NOTIFY=31
EOF
                make reset
                grep -vzq FAILURE report.txt
            ) >>/dev/null || echo "fail: $i, \"$line\""
        fi
        i=$((i + 1))
    done

    for line in "${inputLineNoBias[@]}"; do
        dir=${basedir}$i
        mkdir -pv "$dir"

        cat <<EOF >"$dir"/plumed.dat
MOLINFO STRUCTURE=helix.pdb
v: $line
PRINT ARG=v.* STRIDE=1 FILE=colvar FMT=%8.4f
EOF
        if [[ $initialize = "true" ]]; then
            (
                cd "$dir" || exit
                touch forces.reference
                touch colvar.reference
                cp ../rt-GPU-secondarystructure-DRMSD/*.pdb .
                cp ../rt-GPU-secondarystructure-DRMSD/*.xyz .
                echo "include ../../scripts/test.make" >Makefile

                cat <<EOF >config
type=driver
arg="--plumed plumed.dat --trajectory-stride 10 --timestep 0.005 --ixyz ala12_trajectory.xyz --dump-forces forces --dump-forces-fmt=%10.6f"

#export NVCOMPILER_ACC_NOTIFY=31
EOF
                make reset
                grep -vzq FAILURE report.txt
            ) >>/dev/null || echo "fail: $i, \"$line\""
        fi
        i=$((i + 1))
    done

fi

if [[ $action = "run" ]]; then
    for dir in rt-GPU-*; do

        (
            cd "$dir" || exit
            make
            grep -vzq FAILURE report.txt
        ) >>/dev/null || {
            echo "####"
            echo "FAILURE in:"
            echo "$dir"
            grep 'v:' "${dir}/plumed.dat"
        } && {
            echo -n ""
            # tail "$dir/tmp/out" | grep 'Cycles        Total      Average      Minimum      MaximumCalculating'
            # tail "$dir/tmp/out" | grep Calculating -A1
        }
    done
fi

if [[ $action = "clean" ]]; then
    for dir in rt-GPU-secondarystructure-DRMSD-play-*; do
        (
            cd "$dir" || exit
            make clean

        )
    done
fi
