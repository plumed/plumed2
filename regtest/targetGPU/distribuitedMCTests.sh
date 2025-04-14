#!/bin/bash
# checked with shellcheck
# formatted with shfmtv3.36.0

action="run"
initialize="false"
USEGPU="USEGPU"
basedir="rt-GPU-MCVT-"
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
    inputLineComponents=()
    inputLine[${#inputLine[@]}]="DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 $USEGPU"
    inputLineComponents[${#inputLineComponents[@]}]="DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 COMPONENTS  $USEGPU"
    inputLineComponents[${#inputLineComponents[@]}]="DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 SCALED_COMPONENTS  $USEGPU"
    inputLine[${#inputLine[@]}]="ANGLE ATOMS1=1,2,3 ATOMS2=4,5,6,7 ATOMS3=8,9,10  $USEGPU"
    inputLine[${#inputLine[@]}]="TORSION ATOMS1=1,2,3,4 VECTORA2=5,6 AXIS2=7,8 VECTORB2=9,10 ATOMS3=11,12,13,14  $USEGPU"
    inputLine[${#inputLine[@]}]="TORSION ATOMS1=1,2,3,4 VECTORA2=5,6 AXIS2=7,8 VECTORB2=9,10 ATOMS3=11,12,13,14 COSINE  $USEGPU"
    #inputLine[${#inputLine[@]}]="TORSION ATOMS1=1,2,3,4 VECTORA2=5,6 AXIS2=7,8 VECTORB2=9,10 ATOMS3=11,12,13,14  $USEGPU"
    #ct5: CUSTOM ARG=tt4 FUNC=cos(x) PERIODIC=NO

    inputLineComponents[${#inputLineComponents[@]}]="POSITION ATOM1=1 ATOM2=3 ATOM3=5  $USEGPU"
    inputLineComponents[${#inputLineComponents[@]}]="POSITION ATOM1=1 ATOM2=3 ATOM3=5 SCALED_COMPONENTS  $USEGPU"
    inputLine[${#inputLine[@]}]="DIPOLE GROUP1=1-6 GROUP2=7-12   $USEGPU"
    inputLineComponents[${#inputLineComponents[@]}]="DIPOLE GROUP1=1-6 GROUP2=7-12 COMPONENTS  $USEGPU"

    i=0
    for line in "${inputLine[@]}"; do
        echo "$i $line"
        dir=${basedir}$i
        mkdir -pv "$dir"

        cat <<EOF >"$dir"/plumed.dat
v: $line
PRINT ARG=v STRIDE=1 FILE=colvar FMT=%8.4f
#RESTRAINT ARG=v AT=2 KAPPA=3 SLOPE=3
EOF
        if [[ $initialize = "true" ]]; then
            (
                cd "$dir" || exit
                #touch forces.reference
                touch colvar.reference

                echo "include ../../scripts/test.make" >Makefile

                cat <<EOF >config
type=driver
arg="--plumed plumed.dat --trajectory-stride 10 --timestep 0.005 --ixyz trajectory.xyz --dump-forces forces --dump-forces-fmt=%10.6f --pdb print_pdb_test.pdb"

extra_files="../../trajectories/trajectory.xyz ../../trajectories/print_pdb_test.pdb"

#export NVCOMPILER_ACC_NOTIFY=31
EOF
                make reset
                grep -vzq FAILURE report.txt
            ) >>/dev/null || echo "fail: $i, \"$line\""
        fi
        i=$((i + 1))
    done

    for line in "${inputLineComponents[@]}"; do
        echo "$i $line"
        dir=${basedir}$i
        mkdir -pv "$dir"

        cat <<EOF >"$dir"/plumed.dat
v: $line
PRINT ARG=v.* STRIDE=1 FILE=colvar FMT=%8.4f
EOF
        if [[ $initialize = "true" ]]; then
            (
                cd "$dir" || exit
                # touch forces.reference
                touch colvar.reference

                echo "include ../../scripts/test.make" >Makefile

                cat <<EOF >config
type=driver
arg="--plumed plumed.dat --trajectory-stride 10 --timestep 0.005 --ixyz trajectory.xyz --dump-forces forces --dump-forces-fmt=%10.6f --pdb print_pdb_test.pdb"

extra_files="../../trajectories/trajectory.xyz ../../trajectories/print_pdb_test.pdb"

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
    for dir in rt-GPU-MCVT*; do

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
