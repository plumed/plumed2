#!/usr/bin/env bash

path=$1

if (($# != 1)) ; then
echo "usage: $0 /path/to/vesin"
echo
echo "All the commands to run:"
echo "    cd /some/dir/"
echo "    git clone https://github.com/luthaf/vesin"
echo "    cd /your/plumed2/src/metatomic"
echo "    ./import.sh /some/dir/vesin"
exit 0
fi

bash -c "$path/create-single-cpp.py > /dev/null"

cp $path/vesin/include/vesin.h vesin.h
mv vesin-single-build.cpp vesin.cpp

# Patch files to follow PLUMED linter
sed 's|#define VESIN_H|#define VESIN_H\n/*INDENT-OFF*/\n|
     s|VESIN_H|__PLUMED_metatomic_vesin_h|g
     s|<stddef.h>|<cstddef>|
     s|<stdint.h>|<cstdint>|
     s|extern "C" {|namespace PLMD {\nnamespace metatomic {\nnamespace vesin {\nextern "C" {|
     s|} // extern "C"|} // extern "C"\n} // namespace vesin\n} // namespace metatomic\n} // namespace PLMD|
    ' vesin.h > tmp
mv tmp vesin.h

sed '1 s|^|/*INDENT-OFF*/\n#include "vesin.h"\n|
     s|<stddef.h>|<cstddef>|
     s|<stdint.h>|<cstdint>|
     s|vesin::|PLMD::metatomic::vesin::|
     s|namespace vesin {|namespace PLMD {\nnamespace metatomic {\nnamespace vesin {|
     s|} // namespace vesin|} // namespace vesin\n} // namespace metatomic\n} // namespace PLMD|
     s|using namespace PLMD::metatomic::vesin::cpu;|using namespace PLMD::metatomic::vesin;\nusing namespace PLMD::metatomic::vesin::cpu;|
    ' vesin.cpp > tmp
mv tmp vesin.cpp

cd ..
./header.sh metatomic
