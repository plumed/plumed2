# This script adds a LGPL header to all the source files
# If the header is already present, it does not touch the file
# Please run it whenever you add a new file so as to keep copyright info there!

test -d ../.git  || {
  echo "ERROR: this script should be run from a git repository!"
  exit 1
}

DIRS=$1

# remove trailing "/"
DIRS=${DIRS%%/*}

test -z "$DIRS" && DIRS=*

for dir in $DIRS
do
test -d "$dir" || continue
cd $dir

FIX_COPYRIGHT=""
test -f COPYRIGHT && FIX_COPYRIGHT="$(<COPYRIGHT)"

for file in *.c *.cpp *.h *.inc.in
do

test -f "$file" || continue

if test -n "$FIX_COPYRIGHT"
then
  COPYRIGHT="$FIX_COPYRIGHT"
  echo -n "Custom header for file: $file"
else
    years=$(git log --follow -M75% --format=%aD $file |
    sort -r -n -k 4 |
    awk -v now="$(date +%Y)" '{if(NR==1)last=$4;}END{
# override last with now
    last=now
# Notice that the script could be simplified, I leave it as is for reference:
    first=$4
    if(first=="") print ""
    else if(first==last) print first;
    else if(first+1==last) print first "," last;
    else print first "-" last}')
  if test -z "$years" ; then
    echo  "Error processing file: $file"
    continue
  fi

  COPYRIGHT="\
   Copyright (c) $years The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>."
  echo -n "LGPL header ($years) for file: $file"
fi


{

plus="+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

cat << EOF 
/* $plus
$COPYRIGHT
$plus */
EOF

awk -v plus=$plus 'BEGIN{
  inheader=0;
}{
  if($1=="/*" && $2==plus) inheader=1;
  if(inheader==0) print $0
  if($1==plus && $2=="*/")  inheader=0;
}' $file

} > $file.tmp

if [[ "$file" == *h && ! "$file" == asmjit_apibegin.h && ! "$file" == asmjit_apiend.h ]] ; then
# if [[ "$file" == *h ]] ; then 
ff="${file//./_}"
guard="__PLUMED_${dir}_${ff}"

awk -v plus=$plus -v guard=$guard '
{
  if(past==1){
    if($1=="#ifndef"){
      line=-10;
    } else if($1=="#define" && line==-9){
      found=1;
    }
    else if(found!=1 && NF==0);
    else past=2;
  }
  if(past==0 || past==2) print $0;
  if(past==0 && $1==plus && $2=="*/") {
    past=1;
    print "#ifndef "guard
    print "#define "guard
  }
  if(NF>0)line++;
}END{
if(!found) print "#endif"
}' $file.tmp > $file.tmp2
mv $file.tmp2 $file.tmp

fi

if cmp -s $file $file.tmp ; then
  echo 
else
  cp $file.tmp $file
  echo " +++ PATCHED"
  git add $file
fi

rm $file.tmp

done

cd -

done


