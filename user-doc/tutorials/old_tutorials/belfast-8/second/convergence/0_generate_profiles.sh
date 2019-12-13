for j in `seq 0 3`; do
HILLS=HILLS.$j
N_INTERVAL=50
KT=2.49

lines=`grep -v \# ../$HILLS | wc -l  | awk '{print $1}'`
stride=`echo $lines/$N_INTERVAL | bc`

mkdir cv$j
plumed sum_hills --hills ../$HILLS --stride $stride --mintozero --outfile cv$j/fes_

done

