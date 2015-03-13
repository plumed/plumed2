# number of free-energy profiles
nfes=$(( $1 -1 ))
# minimum of basin A
minA=$2
# maximum of basin A
maxA=$3
# minimum of basin B
minB=$4
# maximum of basin B
maxB=$5
# temperature in energy units
kbt=$6

for i in `seq 0 ${nfes}`
do
 # calculate free-energy of basin A
 A=`awk 'BEGIN{tot=0.0}{if($1!="#!" && $1>min && $1<max)tot+=exp(-$2/kbt)}END{print -kbt*log(tot)}' min=${minA} max=${maxA} kbt=${kbt} fes_${i}.dat`
 # and basin B
 B=`awk 'BEGIN{tot=0.0}{if($1!="#!" && $1>min && $1<max)tot+=exp(-$2/kbt)}END{print -kbt*log(tot)}' min=${minB} max=${maxB} kbt=${kbt} fes_${i}.dat`
 # calculate difference
 Delta=`echo "${A} - ${B}" | bc -l`
 # print it
 echo $i $Delta
done
