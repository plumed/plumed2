for f in `ls tmp | grep "wavelet_grid.data$"`
do
  cp tmp/${f} ${f}.reference 
done 
