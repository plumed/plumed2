if [ "$1" = --description ] ; then
  echo "create a new collective variable from a template"
  exit 0
fi

if [ "$1" = --help ] || [ "$1" = -h ] ; then
  cat <<EOF
Usage:

  plumed partial_tempering scale < processed.top

where scale is the Hamiltonian scaling factor and
processed.top is a post-processed topology file (i.e. produced with grompp -pp)
where each "hot" atom has a "_" appended to the atom type, e.g.:

     1 amber99_43_     1    RC5    O5'      1    -0.6223         16   ; qtot -0.6223

Notice that the section that should be edited is the [atoms] section for all the
molecules that you wish to affect (typically only for the solute, but you may also
want to change solvent parameters).

Also remember to first produce the processed.top file with grompp -pp. Editing a normal
topol.top file will not work, because it does not contain all the parameters.
The processed.top file should not have any "#include" statement.

# produce a processed topology
grompp -pp
# choose the "hot" atoms
vi processed.top
# generate the actual topology
plumed partial_tempering \$scale < processed.top > topol\$i.top

WARNING: It's not very robust! A suggested check is to compare partial_tempering with scale=1.0
to non-scaled force field. E.g.
grompp -o topol-unscaled.tpr
grompp -pp
vi processed.top # choose the "hot" atoms
plumed partial_tempering 1.0 < processed.top > topol-scaled.top # scale with factor 1
grompp -p topol-scaled.top -o topol-scaled.tpr
# Then do a rerun on a trajectory
mdrun -s topol-unscaled.tpr -rerun rerun.trr
mdrun -s topol-scaled.tpr -rerun rerun.trr
# and compare the resuling energy files. they should be identical

EOF
  exit
fi

awk -v scale=$1 '
function recname()
{
     if($1=="[" && $3=="]") return $2;
     return "";
}
function error(msg)
{
     print "ERROR:",msg > "/dev/stderr" ;
     exit;
}
{
# This is the suffix for "hot" atoms:
  suffix="_";

# format for writing parameters
  CONVFMT="%.12g"

##### PARSING DATABASE #####
# store comments:
  comments="";
  if(a=match($0,";")) comments=substr($0,a);
# remove comments:
  gsub(";.*","");
# echo empty line
  if(NF==0){
    print comments;
    next;
  }
# set name of current block
  if(recname() ) rec=recname();
# if we are in atomtypes section, check which fields are present
# use same heuristics as in src/kernel/toppush.c
  if(rec=="atomtypes" && NF>=4){
    if((length($4)==1 && $4~"[a-zA-Z]")) error("in atomtypes");
    else if((length($6)==1 && $6~"[a-zA-Z]")){
      epsilonfield=8;
    }
    else if((length($5)==1 && $5~"[a-zA-Z]")){
      if(substr($1,0,1) ~"[a-zA-Z]") epsilonfield=7;
      else error("in atomtypes");
    } else error("in atomtypes");
    if(epsilonfield!=NF) error("in atomtypes");
# NOTE: OPLS uses bond types, thus we have to save the bondtype
#   atom[$1]=$2;
# For other force fields (e.g. AMBER) atomtype is used as bondtype
# and column two is ignored (it is just atomic number).
    atom[$1]=$1;
  }

# storing dihedraltypes:
  if(rec=="dihedraltypes" && ($5==9 || $5==4)){
    string=":"$6" "$7" "$8;
    type=$1"-"$2"-"$3"-"$4"-"$5
    if(type in params) params[type]=params[type]string;
    else               params[type]=NR""string;
# parameters are commented since they are used inline
    print "; LINE("NR")",$0,comments
    next;
  }

##### SCANNING #####
# in case new list of atoms for a new molecule, delete the present list
  if(recname()=="atoms"){
    delete list_of_atoms;
    n_of_atoms=0;
  }
# detect amber type of each atom
  if(rec=="atoms" && NF>6){
     name=$2;
     gsub(suffix"$","",name)
     ato[$1]=name;
  }

##### PRINTING #####
# DIHEDRALS
  if(rec=="dihedrals" && ($5==9 || $5==4) ){
    found1=0; found4=0;
    for(j=0;j<n_of_atoms;j++) {
      if($1==list_of_atoms[j]) found1=1;
      if($4==list_of_atoms[j]) found4=1;
    }
    sscale=1.0;
    if(found1)sscale*=sqrt(scale);
    if(found2)sscale*=sqrt(scale);

     if(NF==8){
       printf($1" "$2" "$3" "$4" "$5" "$6" "$7*sscale" "$8 "comments\n");
     } else if(NF==5){
       param="";
       atype[1]=atom[ato[$1]]
       atype[2]=atom[ato[$2]]
       atype[3]=atom[ato[$3]]
       atype[4]=atom[ato[$4]]
       
       progression=NR
       for(iswitch=0;iswitch<32;iswitch++){
         if(iswitch%2==0){
           a1=atype[1]; a2=atype[2]; a3=atype[3]; a4=atype[4];
         } else {
           a1=atype[4]; a2=atype[3]; a3=atype[2]; a4=atype[1];
         }
         if(int(iswitch/2)%2==1) a1="X";
         if(int(iswitch/4)%2==1) a2="X";
         if(int(iswitch/8)%2==1) a3="X";
         if(int(iswitch/16)%2==1) a4="X";
         test=a1"-"a2"-"a3"-"a4"-"$5;
         if(test in params){
           split(params[test],array,":");
           if(array[1]<progression){
             progression=array[1];
             param=params[test];
           }
         }
       }
       
    n=split(param,array,":");
    if(n<=1) error("params not found "$1" "$2" "$3" "$4" "$5" "ato[$1]" "ato[$2]" "ato[$3]" "ato[$4]);
    if($5!=9 && n!=2) error("multiple dihedrals should be type 9");
    
    printf("; parameters for types %s %s %s %s at LINE(%s)\n",ato[$1],ato[$2],ato[$3],ato[$4],array[1]);
    for(i=2;i<=n;i++){
      printf($1" "$2" "$3" "$4" "$5" ");
      split(array[i],array1," ");
      printf(array1[1]" "sscale*array1[2]" "array1[3]);
      printf(comments);
      printf("\n");
    }
   } else error("dihedrals should have 5 or 8 fields");
# ATOMTYPES
  } else if(rec=="atomtypes" && NF>=4){
    for(i=1;i<NF;i++)printf($i" "); print $NF,comments;
    printf($1""suffix" "$1" ");
      for(i=3;i<NF;i++)printf($i" "); print scale*$NF," ; scaled";
# ATOMTYPES (PAIRS)
  } else if(rec=="nonbond_params" && NF>=5){
    print $1,$2,$3,$4,$5,comments
    print $1""suffix,$2,$3,$4,sqrt(scale)*$5," ; scaled";
    print $1,$2""suffix,$3,$4,sqrt(scale)*$5," ; scaled";
    print $1""suffix,$2""suffix,$3,$4,scale*$5," ; scaled";
# ATOMS
  } else if(rec=="atoms" && NF>=7){
     if($2~".*"suffix"$"){
       if(NF>=8) print $1,$2,$3,$4,$5,$6,$7*sqrt(scale),$8,comments;
       if(NF==7) print $1,$2,$3,$4,$5,$6,$7*sqrt(scale),comments;
       list_of_atoms[n_of_atoms]=$1;
       n_of_atoms++;
     }
     else print $0
# EVERYTHING ELSE (just print)
  } else print $0,comments
}
'

