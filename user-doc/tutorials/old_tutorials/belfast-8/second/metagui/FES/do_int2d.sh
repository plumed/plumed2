#2d
#PARAMETERS:
colX=3
colY=4
KT=2.49


#START
minX=`sort -k$colX -g FES.4D | head -n 1 | awk '{print $'$colX'}'`
minY=`sort -k$colY -g FES.4D | head -n 1 | awk '{print $'$colY'}'`
maxX=`sort -k$colX -g -r FES.4D | head -n 1 | awk '{print $'$colX'}'`
maxY=`sort -k$colY -g -r FES.4D | head -n 1 | awk '{print $'$colY'}'`
dX=`sort -k$colX -g FES.4D | awk '{print $'$colX'}' | uniq | head -n 2 | awk '{x[NR]=$1} END {print x[2]-x[1]}'`
dY=`sort -k$colY -g FES.4D | awk '{print $'$colY'}' | uniq | head -n 2 | awk '{x[NR]=$1} END {print x[2]-x[1]}'`

echo $colX, $minX, $maxX, $dX
echo $colY, $minY, $maxY, $dY

sort -k$colX -k$colY -g FES.4D | awk 'BEGIN {xold='$minX'; yold='$minY';} 
                                            {if($'$colX'!=xold||$'$colY'!=yold) {print xold, yold, -'$KT'*log(sum); sum=0;} 
                                             sum+=exp(-$7/'$KT'); xold=$'$colX'; yold=$'$colY';}' > tmp."$colX"_"$colY"

awk '{ if($3>800) $3=1000; 
       x[NR]=$1; 
       y[NR]=$2; 
       f[NR]=$3; 
       linee=NR 
     } END { 
       for(i='$minX';i<='$maxX';i+='$dX') { 
         for(j='$minY';j<='$maxY';j+='$dY') {
           fes=1000; 
           for(k=1;k<=linee;k++) { 
             if((x[k]>i-('$dX'/10)&&x[k]<i+('$dX'/10))&&(y[k]>j-('$dY'/10)&&y[k]<j+('$dY'/10))) fes=f[k]; 
           } 
           print i,j,fes;
         } print "";
       }
     }' tmp."$colX"_"$colY" > fes."$colX"_"$colY"

rm -f tmp."$colX"_"$colY" 
#DONE
