# THIS SCRIPT TAKES IN INPUT:
# the number of 1D free energy profiles to be analysed (fes_#.dat)
# kbt that is the temperature in energy units
# AND GIVES IN OUTPUT:
# the arithmetic averaged difference between successive 1D free energy profiles
# the free energy weighted difference between successive 1D free energy profiles 
# ./2_check_profiles.sh 50 2.49

N_INTERVAL=$(( $1 -1 ))
KT=$2
rm -f average-deviation.dat KL-deviation.dat
 
for i in `seq 0 $(($N_INTERVAL-1))`; do
paste fes_"$(($i+1))".dat fes_"$i".dat | awk '{if($1!="") {if($2<10*'$KT'||$5<10*'$KT') {diff+=($2-$5)^2; c++} }} END {print "'$i'", sqrt(diff/c)}' >> average-deviation.dat 

paste fes_"$(($i+1))".dat fes_"$i".dat | awk '{if($1!="") {norma+=exp(-$2/'$KT'); normb+=exp(-$5/'$KT'); a[NR]=$2; b[NR]=$5 }} 
                                               END {for(i=1;i<=NR;i++) { 
                                                    sumleft+=(log((exp(-a[i]/'$KT')/norma)/(exp(-b[i]/'$KT')/normb)))*(exp(-a[i]/'$KT')/norma); 
                                                    sumright+=(log((exp(-b[i]/'$KT')/normb)/(exp(-a[i]/'$KT')/norma)))*(exp(-b[i]/'$KT')/normb); } 
                                                    print "'$i'", '$KT'*0.5*(sumleft+sumright) }' >> KL-deviation.dat

done

