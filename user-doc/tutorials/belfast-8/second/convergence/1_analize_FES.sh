# THIS SCRIPT TAKES IN INPUT:
# the number of 1D free energy profiles to be analysed (fes_#.dat)
# minA and maxA that are the intervals defining the first region
# minB and maxB that are the intervals defining the second region
# kbt that is the temperature in energy units
# AND GIVES IN OUTPUT:
# the difference in free energy between A and B.
# ./1_analize_FES.sh 50 0 1 4 5 2.49

nfes=$(( $1 -1 ))
minA=$2
maxA=$3
minB=$4
maxB=$5
kbt=$6

for i in `seq 0 ${nfes}`
do

 A=`awk 'BEGIN{totA=0.0}{if($1!="#!" && $1>minA && $1<maxA)totA+=exp(-$2/kbt)}END{print -kbt*log(totA)}' minA=${minA} maxA=${maxA} kbt=${kbt} fes_${i}.dat`
 B=`awk 'BEGIN{totB=0.0}{if($1!="#!" && $1>minB && $1<maxB)totB+=exp(-$2/kbt)}END{print -kbt*log(totB)}' minB=${minB} maxB=${maxB} kbt=${kbt} fes_${i}.dat`

 Delta=`echo "${A} - ${B}" | bc -l`
 echo $i $Delta

done > deltaG_AB.dat
