#! /bin/bash


rm -f sourceme.sh Makefile.conf

if (($#==1)) ; then

conf=$1

else

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

fi


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

echo 'export PATH="'"$PWD"'/src/lib/:$PATH"' >> sourceme.sh
# this is just for mac:
echo 'export DYLD_LIBRARY_PATH="'"$PWD"'/src/lib:$DYLD_LIBRARY_PATH"' >> sourceme.sh

cat << EOF >> sourceme.sh
export PLUMED_KERNEL="$PWD/src/lib/libplumedKernel.$SOEXT"
EOF

echo "PLUMED will be installed using prefix /usr/local"
echo "If you wish to change this, set PLUMED_PREFIX environment variable before compiling"
echo "Executable will be named 'plumed'"
echo "To add a suffix to this name, set PLUMED_LIBSUFFIX environment variable before compiling"

