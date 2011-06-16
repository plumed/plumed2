#! /bin/bash


rm -f sourceme.sh Makefile.conf

prefix=""
case "$(uname)" in
(Linux)  prefix=linux. ;;
(AIX)    prefix=aix. ;;
(Darwin) prefix=mac. ;;
esac

conflist=$(cd configurations ; echo $prefix*)

PS3="Choose a configuration:"
select conf in $conflist
do
  [[ -n "$conf" ]] && break
done

ln -s configurations/$conf Makefile.conf


case "$conf" in
(aix.*)
  SOEXT=so
  echo "On cineca, you may need to load the following modules:"
  echo "module load make"
  echo "module load sed"
;;
(mac.*)
  SOEXT=dylib  
;;
(*)
  SOEXT=so
esac

echo 'export PATH="$PATH:'"$PWD"'/tools/"' >> sourceme.sh
# this is just for mac:
echo 'export DYLD_LIBRARY_PATH="$DYLD_LIBRARY_PATH:'"$PWD"'/src/"' >> sourceme.sh

cat << EOF >> sourceme.sh
export PLUMED_KERNEL="$PWD/src/libplumedKernel.$SOEXT"
export PLUMED_ROOT="$PWD"
EOF

