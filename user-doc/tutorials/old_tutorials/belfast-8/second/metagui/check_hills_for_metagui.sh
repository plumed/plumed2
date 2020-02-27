ok=0
res=`head -n 1 $1 | grep "ACTIVE" | wc | awk '{print $2}'`
if [ $res != 5 ]; then 
  echo "HILLS files shoud contain the line #! ACTIVE X Y Z as the FIRST line of the file"
  echo "X is the total number of biased CV in that replica"
  echo "Y is the index of the collective variables biased in that replica (i.e. 1 for the first, ecc)"
  echo "Z is the label of the replica (i.e. A for the first, B for the second, ...)"
else ok=1
fi
