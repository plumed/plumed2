
# paste all biases in a single file
paste "$@" | grep -v \# | awk '{for(i=1;i<=NF/4;i++) printf($(i*4)"  ");printf("\n");}' > ALLBIAS 

# get number of biases
nbias=`head -n 1 ALLBIAS | awk '{print NF}'`

# do wham
python ../SCRIPTS/wham.py ALLBIAS $nbias 2.5

# this is dumping phi on another file
cat "$1" | grep -v \# | awk '{print $2}' > phi 
# this is dumping psi on another file
cat "$1" | grep -v \# | awk '{print $3}' > psi

# the files are pasted 
paste phi weights.dat > phi_weights.dat
paste psi weights.dat > psi_weights.dat

# make free energy
python ../SCRIPTS/do_fes.py phi_weights.dat 1 -3.14 3.14 50 2.5 phi_fes.dat
python ../SCRIPTS/do_fes.py psi_weights.dat 1 -3.14 3.14 50 2.5 psi_fes.dat
