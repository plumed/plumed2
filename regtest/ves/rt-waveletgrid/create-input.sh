for i in `seq 4 20`
do
  echo "BF_WAVELETS ..."
  echo " ORDER=${i}"
  echo " TYPE=DAUBECHIES"
  echo " MINIMUM=0.0 # does not matter"
  echo " MAXIMUM=1.0"
  echo " DUMP_WAVELET_GRID"
  echo " WAVELET_FILE_FMT=%13.3f"
  echo " LABEL=db${i}"
  echo "... BF_WAVELETS" 
  echo " "
done 

for i in `seq 4 20`
do
  echo "BF_WAVELETS ..."
  echo " ORDER=${i}"
  echo " TYPE=SYMLETS"
  echo " MINIMUM=0.0 # does not matter"
  echo " MAXIMUM=1.0"
  echo " DUMP_WAVELET_GRID"
  echo " WAVELET_FILE_FMT=%13.3f"
  echo " LABEL=sym${i}"
  echo "... BF_WAVELETS"
  echo " "
done
