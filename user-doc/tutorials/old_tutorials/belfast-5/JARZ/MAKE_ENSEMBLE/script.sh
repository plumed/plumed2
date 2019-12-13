#!/bin/bash
GROMACS_BIN=/usr/local/codes/gromacs/4.6.5/plumed/2/bin/
source $GROMACS_BIN/GMXRC.bash
#which mdrun is this
MDRUN=mdrun_mpi-dp-pl
GROMPP=grompp_mpi-dp-pl 

rm -rf topol.tpr
$GROMPP -c c7eq.gro -p topol.top -f md.mdp
if [ ! -e topol.tpr ] ; then 
	echo " NO TOPOLOGY PRODUCED!!! CHECK YOUR TOP/GRO/MDP"
	exit
fi
#
# do the mdrun
#
$GROMACS_BIN/$MDRUN -plumed plumed.dat
