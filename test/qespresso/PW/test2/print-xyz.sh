grep -A6 "ATOMIC_POSITIONS" $1 | \
awk 'BEGIN{a=0.529177}
{
if($1=="ATOMIC_POSITIONS"){print 6;s++;print "step ",s}
if($1=="Cl"||$1=="C"||$1=="H")printf "%2s %10.4f %10.4f %10.4f\n",$1,$2*a,$3*a,$4*a
}'
