#! /usr/bin/env bash

path=$1

if (($# != 1)) ; then
echo "usage: $0 /path/to/small_vector.hpp"
echo
echo "Full explanation:"
echo "cd /some/dir/"
echo "git clone https://github.com/gharveymn/small_vector"
echo "cd /your/plumed2/src/small_vector"
echo "./import /some/dir/source/include/gch/small_vector.hpp"
exit 0
fi

  awk '{
     if($1=="namespace" && $2=="gch") {
       print "namespace PLMD {"
     }
     if($1=="noexcept" && $2=="(noexcept" && $3=="(Allocator::operator=" && $4=="(std::move" && $5=="(other))))") fix_gcc=1
     if(fix_gcc) print("#if defined(__clang__) || ! defined(__GNUC__) || __GNUC__>=9")
     print
     if(fix_gcc) print("#endif")
     fix_gcc=0

     if($1=="}" && $2=="//" && $3=="namespace" && $4=="gch") {
       print "} // namespace PLMD"
     }

     if($1=="{" && operator) {
       print("#ifdef _GLIBCXX_DEBUG")
       print(" if (size () <= pos) base::throw_index_error ();")
       print("#endif")
     }
     if($1=="operator[]" && $3=="pos)") operator=1
     else operator=0

  }' $1 |
sed "s/GCH/PLUMED_GCH/g" > small_vector.h

git add small_vector.h
git commit -m update small_vector.h

cd ../
./header.sh small_vector





