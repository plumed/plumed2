# This is the script I used to generate trajectory.xyz
# I leave it here as a reference, but it is not necessary to run it
for((i=-1;i<=10;i++)) ; do
for((j=-1;j<=1;j++)) ; do
for((k=-1;k<=1;k++)) ; do
awk -v n=3 -v a0=2 'BEGIN{
  if(a0=="") a0=1.6796
  if(n=="") n=1
  l=n*a0
  print n*n*n*4
  print l,l,l
  for(i=0;i<n;i++) {
  for(j=0;j<n;j++) {
  for(k=0;k<n;k++) {
    print "Ar",i*a0,j*a0,k*a0
    print "Ar",(i+0.5)*a0,j*a0,(k+0.5)*a0
    print "Ar",(i+0.5)*a0,(j+0.5)*a0,k*a0
    print "Ar",i*a0,(j+0.5)*a0,(k+0.5)*a0
  }
  }
  }
}'  |
  awk -v i=$i -v j=$j -v k=$k 'BEGIN{
      srand(100*i+10*j+k);
    }{
    if(NF!=3) {
       print
    } else {
      if(i<=1){
        for(s=1;s<=3;s++) v[s,s]=$s;
        for(s=1;s<=3;s++) v[1,s]=v[1,s]+i*v[2,s];
        for(s=1;s<=3;s++) v[2,s]=v[2,s]+j*v[3,s];
        for(s=1;s<=3;s++) v[3,s]=v[3,s]+k*v[1,s];
      } else {
        for(s=1;s<=3;s++) for(t=1;t<=3;t++) v[s,t]=rand()-0.5;
      }
      for(s=1;s<=3;s++) for(ss=1;ss<=3;ss++) printf(v[s,ss]" "); print ""
    }
}'
done
done
done > trajectory.xyz
