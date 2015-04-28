# count how many simulations have been done
nsym=$(echo "$@" | wc -w)

# feed the wham executable with bias, in kbT units
paste "$@" | grep -v \# | awk '{for(i=1;i<=NF/4;i++) printf($(i*4)/2.5"  ");printf("\n");}' | ./wham.x -n $nsym > weights

# this is dumping time, phi, and psi on another file
cat "$1" | grep -v \# | awk '{print $1,$2,$3}' > first

# the two files are pasted, producing a file with 4 columns
# time, phi, psi, and weight. notice that the latter is written as
# kT*log(w)
paste first weights | awk '{print $1,$2,$3,2.5*log($4)}'
