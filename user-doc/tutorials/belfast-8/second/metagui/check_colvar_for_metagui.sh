ok=0
res=`head -n 2 $1 | grep "ACTIVE" | wc | awk '{print $2}'`
if [ $res != 5 ]; then 
  echo "COLVAR files shoud contain the line #! ACTIVE X Y Z after the FIELDS line"
  echo "X is the total number of biased CV in that replica"
  echo "Y is the index of the collective variables biased in that replica (i.e. 1 for the first, ecc)"
  echo "Z is the label of the replica (i.e. A for the first, B for the second, ...)"
else ok=1
fi
res=`head -n 2 $1 | grep "FIELD" | wc | awk '{print $2}'`
if [ $res -lt 4 ]; then 
  echo "COLVAR or HILLS files should contain a FIELDS field with at least four keywords! (i.e #! FIELDS time cv1)"
else ok=1
fi 
head -n 2 $1 | grep "FIELDS" | awk '{for(i=4;i<=NF;i++) {string="cv"i-3; if($i!=string) print "PROBLEM in '$1':\n FIELDS line should look like:\n #! FIELDS time cv1 ... cvN\n your is:\n",$0}}' 

