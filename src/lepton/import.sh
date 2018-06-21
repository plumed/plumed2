#! /bin/bash

path=$1

if (($# != 1)) ; then
echo "usage: $0 /path/to/lepton"
echo
echo "We currently use lepton from openmm, so the full procedure is"
echo "cd /some/dir/"
echo "git clone https://github.com/pandegroup/openmm.git"
echo "cd /your/plumed2/src/lepton"
echo "./import /some/dir/openmm/lepton"
exit 0
fi

rm -f *.cpp *.h

cp $path/include/Lepton.h .
cp $path/include/lepton/*.h .
cp $path/src/*.cpp $path/src/*.h .

for file in *.h *.cpp ; do
  sed 's|lepton/||
       s|Lepton|lepton|
       s|LEPTON_USE_JIT|__PLUMED_HAS_ASMJIT|g
       s|asmjit.h|asmjit/asmjit.h|g
      ' $file > tmp
  mv tmp $file
  dos2unix $file
done

for file in *.h ; do
  awk '{
     if($1=="namespace" && $2=="lepton") print "namespace PLMD {"
     print
     if($1=="}" && $2=="//" && $3=="namespace" && $4=="lepton") print "} // namespace PLMD"
  }' $file > tmp
  mv tmp $file
done


for file in *.cpp ; do
  awk '{
     if(!done && $1=="using" && $2=="namespace") {
       print "namespace PLMD {";
       done=1
     }
     print
  }END{
    print "}"
  }' $file > tmp
  mv tmp $file
done

cd ../
./header.sh lepton





