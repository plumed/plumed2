BF=bf
out_file=derivs-diff.data
NumFields=`awk '{print NF}' bf.derivs.data | tail -n 1 `
cat $BF.derivs.data | grep -v "^#" | awk '{print $1}' > ${out_file}
for ((i=2; i <= NumFields ; i++))
do
  cat $BF.derivs.data            |  grep -v "^#"  |  awk -v f=$i '{print $f}' > $$.d1  
  cat $BF.derivs-numerical.data  |  grep -v "^#"  |  awk -v f=$i '{print $f}' > $$.d2  
  paste $$.d1 $$.d2 | awk '{printf "%13.5f\n",$1-$2}' > $$.diff
  cp ${out_file} $$.save
  paste $$.save $$.diff > ${out_file}
  rm -f $$.d1 $$.d2 $$.save $$.diff
done
cp ${out_file} $$.data
grep "^#" $BF.derivs.data > $$.header
cat $$.header $$.data > ${out_file}
rm -f $$.data $$.header
