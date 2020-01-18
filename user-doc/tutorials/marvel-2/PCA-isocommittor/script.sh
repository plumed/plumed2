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
		echo -n "unset key ; set xl \"t\"; set yl \"Phi\"; set xr [0:0.6]; set yr [-pi:pi] ; plot  ">script_phi.gplt
		echo -n "unset key ; set xl \"t\"; set yl \"Psi\"; set xr [0:0.6]; set yr [-pi:pi] ; plot  ">script_psi.gplt
		echo -n "unset key ; set xl \"t\"; set yl \"PCA\"; set xr [0:0.6]; set yr [-0.15:0.15] ; plot  ">script_pca.gplt
	fi
	echo -n "\"colvar_$i\" u 1:2 w l ">>script_phi.gplt
	echo -n "\"colvar_$i\" u 1:3 w l ">>script_psi.gplt
        echo -n "\"colvar_$i\" u 1:4 w l ">>script_pca.gplt
	if [ $i -ne $ntests  ]; then 
		echo -n ", ">>script_phi.gplt
		echo -n ", ">>script_psi.gplt
		echo -n ", ">>script_pca.gplt
	fi
	
done
echo -n ", -1.424401 w l, 1.197514 w l" >>script_phi.gplt
echo -n ", 1.244763 w l, -1.163951 w l" >>script_psi.gplt 
echo -n ", 0.000000 w l, -0.136399 w l" >>script_pca.gplt
