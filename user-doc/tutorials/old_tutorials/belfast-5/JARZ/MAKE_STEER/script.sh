#!/bin/bash
GROMACS_BIN=/usr/local/codes/gromacs/4.6.5/plumed/2/bin/
source $GROMACS_BIN/GMXRC.bash
#which mdrun is this
MDRUN=mdrun_mpi-dp-pl
GROMPP=grompp_mpi-dp-pl 

ntests=20
for i in `seq 1 $ntests`
do
	#
	# do the topology: this should write a topol.tpr
	#
	rm -rf topol.tpr
	sed s/SEED/$RANDOM/ md.mdp >newmd.mdp
	$GROMACS_BIN/$GROMPP -c start_$i.gro -p topol.top -f newmd.mdp
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
                echo -n "unset key; set auto ; set xl \"Phi\"; set yl \"Psi\"; set xr [-pi:pi]; set yr [-pi:pi] ; plot  ">script_rama.gplt	
		echo -n "set auto; unset key ; set xl \"Phi\"; set yl \"Work\"; plot  ">script_work.gplt
	fi
	echo -n "\"colvar_$i\" u 2:3 w l ">>script_rama.gplt
	echo -n "\"colvar_$i\" u 6:7 w l ">>script_work.gplt
	if [ $i -ne $ntests  ]; then 
		echo -n ", ">>script_rama.gplt
		echo -n ", ">>script_work.gplt
	fi
	
done
