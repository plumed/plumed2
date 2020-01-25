#!/bin/bash

nrep=4
T1=300.0
T2=400.0


for((i=0;i<nrep;i++))
do 
 ii=$(($i+1))
 temp=`echo " e( ( ${i} / ( ${nrep} - 1 ) ) * l( $T2 / $T1 )) * $T1 " | bc -l | awk '{printf "%6.3f",$1}' `
 #echo $i  $temp
 sed -e "s/TEMP_/${temp}/g" TEMPLATE_md.mdp > md.mdp 
 grompp_mpi -f md.mdp -c conf${i}.gro -p topol.top -o topol${i}.tpr 
 rm md*.mdp
done
