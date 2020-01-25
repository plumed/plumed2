#!/bin/bash
#
# Note : assumes that plumed is in your path ans it is set correctly
#
#
# wanna the full gromacs suite on your path? Then source the GMXRC
#
GROMACS_BIN=/usr/local/codes/gromacs/4.6.5/plumed/2/bin/
source $GROMACS_BIN/GMXRC.bash
#which mdrun is this
MDRUN=mdrun_mpi-dp-pl
GROMPP=grompp_mpi-dp-pl 

ntests=50
for i in `seq 1 $ntests`
do
	#
	# do the topology: this should write a topol.tpr
	#
	rm -rf topol.tpr
	sed s/SEED/$RANDOM/ md.mdp >newmd.mdp
	$GROMPP -c start.gro -p topol.top -f newmd.mdp
	if [ ! -e topol.tpr ] ; then 
		echo " NO TOPOLOGY PRODUCED!!! CHECK YOUR TOP/GRO/MDP"
		exit
	fi
	#
	# do the mdrun
	#
	$GROMACS_BIN/$MDRUN -plumed plumed.dat
	mv colvar colvar_$i
	rm -rf \#*
	if [ $i -eq 1 ]; then 
		echo -n "unset key ; set xl \"Phi\"; set yl \"Psi\"; set xr [-pi:pi]; set yr [-pi:pi] ; plot  ">script_rama.gplt
		echo -n "set yr [0:0.02] ; unset key ; set xl \"Progress along the path\"; set yl \"Distance from the path\" ; plot  ">script_path_right.gplt
		echo -n "set yr [0:0.02] ; unset key ; set xl \"Progress along the path\"; set yl \"Distance from the path\" ; plot  ">script_path_wrong.gplt
	fi
	echo -n "\"colvar_$i\" u 2:3 w l ">>script_rama.gplt
	echo -n "\"colvar_$i\" u 4:5 w l ">>script_path_right.gplt
	echo -n "\"colvar_$i\" u 6:7 w l ">>script_path_wrong.gplt
	if [ $i -ne $ntests  ]; then 
		echo -n ", ">>script_rama.gplt
		echo -n ", ">>script_path_right.gplt
		echo -n ", ">>script_path_wrong.gplt
	fi
	
done
