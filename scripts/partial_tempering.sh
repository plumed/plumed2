#! /bin/bash
# vim:ft=awk
if [ "$1" = --description ] ; then
  echo "scale parameters in a gromacs topology to implement solute or partial tempering"
  exit 0
fi

if [ "$1" = --options ] ; then
  echo "--description --gromacs4 --help -h --options"
  exit 0
fi

if [ "$1" = --help ] || [ "$1" = -h ] ; then
  cat <<EOF
Usage:

  plumed partial_tempering [--gromacs4] scale < processed.top

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

WARNING: It's not very robust and there might be force-field dependent issues!
A few tests are strongly suggested.

1. Compare partial_tempering with scale=1.0 to non-scaled force field. E.g.
grompp -o topol-unscaled.tpr
grompp -pp
vi processed.top # choose the "hot" atoms appending "_". You can choose whatever.
plumed partial_tempering 1.0 < processed.top > topol-scaled.top # scale with factor 1
grompp -p topol-scaled.top -o topol-scaled.tpr
# Then do a rerun on a trajectory
mdrun -s topol-unscaled.tpr -rerun rerun.trr
mdrun -s topol-scaled.tpr -rerun rerun.trr
# and compare the resuling energy files. they should be identical

2. Compare partial_tempering with scale=0.5 to non-scaled force field.
Repeat the same procedure but using "plumed partial_tempering 0.5".
Choose all the atoms in all the relevant [atoms] sections (e.g. solute, solvent and ions).
In the two resulting energy files you should see:
long range electrostatics, LJ, and dihedral energy is *half* in the scaled case
all other terms (bonds/bends) are identical.

EOF
  exit
fi

gromacs5=1

if [ "$1" == --gromacs4 ] ; then
  gromacs5=0
  shift
fi

awk -v scale=$1 -v gromacs5=$gromacs5 '
BEGIN{
  combrule=1;
}
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
function warning(msg)
{
     print "WARNING:",msg | "cat 1>&2"
}
function find_matching_torsion(params,atype, a1,a2,a3,a4,iswitch,progression,test,array,param,countX,bestCountX,bestMatch)
{
   progression=NR
   bestCountX=5
   for(iswitch=0;iswitch<32;iswitch++){
     countX=0
     if(iswitch%2==0){
       a1=atype[1]; a2=atype[2]; a3=atype[3]; a4=atype[4];
     } else {
       a1=atype[4]; a2=atype[3]; a3=atype[2]; a4=atype[1];
     }
     if(int(iswitch/2)%2==1){ a1="X"; countX++; }
     if(int(iswitch/4)%2==1){ a2="X"; countX++; }
     if(int(iswitch/8)%2==1){ a3="X"; countX++; }
     if(int(iswitch/16)%2==1){a4="X"; countX++; }
     test=a1"-"a2"-"a3"-"a4"-"$5;
     if(test in params){
       split(params[test],array,":");
       bestMatch=0;
       if(gromacs5) {
         if(countX<bestCountX || (countX==bestCountX && array[1]<progression)) bestMatch=1
       } else {
         if(array[1]<progression) bestMatch=1;
       }
       if(bestMatch){
         progression=array[1];
         bestCountX=countX
         param=params[test];
       }
     }
   }
   return param;
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
# set defaults for nb interactions
  if(rec=="defaults" && NF==5) combrule=$2;
# if we are in atomtypes section, check which fields are present
# use same heuristics as in src/kernel/toppush.c
  if(rec=="atomtypes" && NF>=4){
    if((length($4)==1 && $4~"[a-zA-Z]")){
      bondtypefield=1;
      epsilonfield=6;
    }
    else if((length($6)==1 && $6~"[a-zA-Z]")){
      bondtypefield=2;
      epsilonfield=8;
    }
    else if((length($5)==1 && $5~"[a-zA-Z]")){
      if(substr($2,0,1) ~"[a-zA-Z]"){
        bondtypefield=2;
        epsilonfield=7;
      } else {
        bondtypefield=1;
        epsilonfield=7;
      }
    } else error("in atomtypes");
    if(epsilonfield!=NF) error("in atomtypes");
# NOTE: OPLS uses bond types, thus we have to save the bondtype
# For other force fields (e.g. AMBER) atomtype is used as bondtype
# and column two is ignored (it is just atomic number).
    bondtype[$1]=$bondtypefield;
  }

# storing dihedraltypes:
  if(rec=="dihedraltypes" && ($5==1 || $5==2 || $5==3 || $5==4 || $5==5 || $5==9)){
    if($5==1 || $5==4 || $5==9) string=":"$6" "$7" "$8;
    else if($5==2) string=":"$6" "$7;
    else if($5==3) string=":"$6" "$7" "$8" "$9" "$10" "$11;
    else if($5==5) string=":"$6" "$7" "$8" "$9;
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
  if(rec=="dihedrals" && ($5==1 || $5==2 || $5==3 || $5==4 || $5==5 || $5==9) ){
    found1=0; found4=0;
    for(j=0;j<n_of_atoms;j++) {
      if($1==list_of_atoms[j]) found1=1;
      if($4==list_of_atoms[j]) found4=1;
    }
    sscale=1.0;
    if(found1)sscale*=sqrt(scale);
    if(found4)sscale*=sqrt(scale);

# this is the case in which dihedrals are online:
     if(NF>5){
       printf($1" "$2" "$3" "$4" "$5" ");
       if($5==1 || $5==4 || $5==9){
                                   if(NF!=8) error("dihedrals with type 1,4,9 should have 8 fields");
                                   printf($6" "$7*sscale" "$8);
       } else if($5==2) {
                                   if(NF!=7) error("dihedrals with type 2 should have 7 fields");
                                   printf($6" "$7*sscale);
       } else if($5==3) {
                                   if(NF!=11) error("dihedrals with type 3 should have 11 fields");
                                   printf($6*sscale" "$7*sscale" "$8*sscale" "$9*sscale" "$10*sscale" "$11*sscale);
       } else if($5==5) {
                                   if(NF!=9) error("dihedrals with type 5 should have 9 fields");
                                   printf($6*sscale" "$7*sscale" "$8*sscale" "$9*sscale);
       } else error("dihedrals with more than 5 fields should be 1,2,3,4,5 or 9");
       printf(" "comments"\n");
# this is the case in which we have to search the database
     } else if(NF==5){
       param="";
       atype[1]=bondtype[ato[$1]]
       atype[2]=bondtype[ato[$2]]
       atype[3]=bondtype[ato[$3]]
       atype[4]=bondtype[ato[$4]]

       param=find_matching_torsion(params,atype);
       
       
    n=split(param,array,":");
    if(n<=1) error("params not found "$1" "$2" "$3" "$4" "$5" "atype[1]" "atype[2]" "atype[3]" "atype[4]);
    if($5!=9 && n!=2){
# in case of multiple dihedrals !=9, all parameters should be the same, otherwise I suspect there is some problem
      for(i=3;i<=n;i++){
        if((array[i]-array[2])**2>1e-20) error("multiple dihedrals !=9: parameters "array[i]" and "array[2]" are different\n");
      }
# then, I just take one of the instances
      param=array[1]":"array[2];
      n=split(param,array,":");
    }
    
    printf("; parameters for types %s %s %s %s at LINE(%s)\n",atype[1],atype[2],atype[3],atype[4],array[1]);
    for(i=2;i<=n;i++){
      printf($1" "$2" "$3" "$4" "$5" ");
      split(array[i],array1," ");
      if($5==1 || $5==4 || $5==9){
                                  printf(array1[1]" "array1[2]*sscale" "array1[3]);
      } else if($5==2) {
                                  printf(array1[1]" "array1[2]*sscale);
      } else if($5==3) {
                                  printf(array1[1]*sscale" "array1[2]*sscale" "array1[3]*sscale" "array1[4]*sscale" "array1[5]*sscale" "array1[6]*sscale);
      } else if($5==5) {
                                  printf(array1[1]*sscale" "array1[2]*sscale" "array1[3]*sscale" "array1[4]*sscale);
      } else error("dihedrals with more than 5 fields should be 1,2,3,4,5 or 9");
      printf(comments);
      printf("\n");
    }
   } else error("dihedrals should have at least 5 fields");
# ATOMTYPES
  } else if(rec=="atomtypes" && NF>=4){
    scale2=1.0; # scaling for second to last column
    if(combrule==1) scale2=scale;
    for(i=1;i<NF;i++)printf($i" "); print $NF,comments;
    printf($1""suffix" "bondtype[$1]" ");
    from=3;
    if(NF==6) from=2; # GROMOS does not store bondtype by default, so we should add one column
    for(i=from;i<NF-1;i++)printf($i" "); print scale2*$(NF-1),scale*$NF," ; scaled";
# ATOMTYPES (PAIRS)
  } else if((rec=="pairtypes" || rec=="nonbond_params") && NF>=5){
    scale2=1.0; # scaling for second to last column
    if(combrule==1) scale2=scale;
    print $1,$2,$3,$4,$5,comments
    print $1""suffix,$2,$3,sqrt(scale2)*$4,sqrt(scale)*$5," ; scaled";
    print $1,$2""suffix,$3,sqrt(scale2)*$4,sqrt(scale)*$5," ; scaled";
    print $1""suffix,$2""suffix,$3,scale2*$4,scale*$5," ; scaled";
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

