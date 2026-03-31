#!/usr/bin/env bash

version=$1

if (($# != 1)) ; then
echo "usage: $0 <vesin-version>"
exit 1
fi

rm -f vesin-single-build-v*.tar.gz
wget https://github.com/Luthaf/vesin/releases/download/v$version/vesin-single-build-v$version.tar.gz

tar xf vesin-single-build-v$version.tar.gz
# bash -c "$path/create-single-cpp.py > /dev/null"

# cp $path/vesin/include/vesin.h vesin.h
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

sed '1 s|^|/*INDENT-OFF*/\n|
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
