# THIS SCRIPT TAKES IN INPUT:
# the number of 1D free energy profiles to be analysed (fes_#.dat)
# the number of profiles to be averaged 
# kbt that is the temperature in energy units
# AND GIVES IN OUTPUT:
# the profiles averaged over the defined window and their standard deviation (i.e. if the window is of 5 profiles, fes-mean-sd-0 will have <0-4>, fes-mean-sd-1 <1-5> ...)
# the weigheted error of the profile.
# ./3_check_averaging.sh 50 5 2.49

N_INTERVAL=$(( $1 -1 ))
WINDOW=$(( $2 -1 ))
KT=$3
MUCH=$(( $N_INTERVAL - $WINDOW ))

for i in `seq 0 $MUCH`; do
STRINGA=""
for j in `seq $i $(( $i + $WINDOW))`; do
STRINGA=$STRINGA" fes_$j.dat"
done
paste $STRINGA | awk '{s=0; s2=0; c=0; for(i=2;i<=NF;i+=3) {s+=$i; s2+=$i*$i; c++} print $1, s/c, sqrt(s2/c-(s*s)/(c*c))}'>  fes-mean-sd_$i.dat
done

for i in `seq 0 $MUCH`; do
if [ $(($i+$WINDOW)) -le $MUCH ]; then
paste fes-mean-sd_$i.dat fes-mean-sd_$(($i+$WINDOW)).dat | awk '{diff+=($2-$5)^2; n++; sum+=$3*exp(-$2/'$KT'); c+=exp(-$2/'$KT')} END {print "'$i'", sum/c, diff/n}'
else
paste fes-mean-sd_$i.dat | awk '{diff+=($2-$5)^2; n++; sum+=$3*exp(-$2/'$KT'); c+=exp(-$2/'$KT')} END {print "'$i'", sum/c, 0}'
fi

done > averfes-deviations.dat
