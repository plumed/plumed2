#!/bin/bash
#
# Note : assumes that plumed is in your path ans it is set correctly
#
#
# wanna the full gromacs suite on your path? Then source the GMXRC
#
GROMACS_BIN=/Users/gareth/MD_code/gromacs-5.1.1/build/bin
GROMACS=$GROMACS_BIN/gmx
source $GROMACS_BIN/GMXRC.bash

ntests=50
for i in `seq 1 $ntests`
do
	#
	# do the topology: this should write a topol.tpr
	#
	rm -rf topol.tpr
	sed s/SEED/$RANDOM/ md.mdp >newmd.mdp
	$GROMACS grompp -c start.gro -p topol.top -f newmd.mdp
	if [ ! -e topol.tpr ] ; then 
		echo " NO TOPOLOGY PRODUCED!!! CHECK YOUR TOP/GRO/MDP"
		exit
	fi
	#
	# do the mdrun
	#
	$GROMACS mdrun -plumed plumed.dat
	mv colvar colvar_$i
	rm -rf \#*
	if [ $i -eq 1 ]; then 
		echo -n "unset key ; set xl \"S(X)\"; set yl \"Z(X)\"; set xr [0:22]; set yr [0:0.01] ; plot  ">script_path.gplt
		echo -n "unset key ; set xl \"S(X)\"; set yl \"Z(X)\"; set xr [-0.15:0.15]; set yr [0:0.2] ; plot  ">script_pca.gplt
	fi
	echo -n "\"colvar_$i\" u 5:4 w l ">>script_path.gplt
	echo -n "\"colvar_$i\" u 2:3 w l ">>script_pca.gplt
	if [ $i -ne $ntests  ]; then 
		echo -n ", ">>script_path.gplt
		echo -n ", ">>script_pca.gplt
	fi
	
done
