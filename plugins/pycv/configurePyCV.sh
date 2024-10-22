#!/bin/bash

#formatted with shfmtv3.36.0 (https://github.com/mvdan/sh/releases)

#A simple script to configure PyCV if python does not work
PYTHON=""
if which plumed >/dev/null; then
    echo "plumed found"
    PYTHON=$(plumed --no-mpi info --configuration | grep 'python_bin=')
    PYTHON=${PYTHON#python_bin=}
fi

if [ -z "$PYTHON" ]; then
    echo "python not found using plumed"
    echo "serching for  available python with plumed module installed"

    for python_bin in python python3 python3.12 python3.11 python3.10 python3.9 python3.8 python3.7; do
        if $python_bin -c "import plumed" 2>/dev/null; then
            if [ $python_bin != python ]; then
                PYTHON=$python_bin
            fi
            return 1
        fi
    done

fi

if [ -z "$PYTHON" ]; then
    echo "python not found"
    exit 1
fi

cat <<EOF >Makefile.conf
PYTHON=$(which "$PYTHON")
EOF
