for file in *.c *.cpp *.h PlumedConfig.h.in
do

if [ $file == "PlumedConfig.h" ] ; then

continue

fi

echo "Applying LGPL header to file $file"

{

plus="+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

cat << EOF 
/* $plus
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
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

cmp -s $file $file.tmp || cp $file.tmp $file

rm $file.tmp

done


