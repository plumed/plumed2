#! /bin/bash

# first argument is name of txt file to be updated
file="$1"

# version strings:
version=$(
echo "version short $(
  if test -f ../../VERSION ; then
    grep -v "#" ../../VERSION | sed  's/^\([0-9][0-9]*\.[0-9][0-9]*\).*/\1/'
  else
    echo "Unknown"
  fi
)"

echo "version long $(
  if test -f ../../VERSION ; then
    grep -v "#" ../../VERSION
  else
    echo "Unknown"
  fi
)"

echo "version git $(
  if test -d ../../.git && hash git 2> /dev/null ; then
# in case it does not work, fallback to normal hash (12 char long)
    git describe --long --dirty --always || git rev-parse  --short=12 HEAD
  else
    echo "Unknown"
  fi
)"

)



# first get full list of possible defines from configure.ac and .h
list="$(
  {
  cat ../../configure.ac
  grep -E "^#" ../*/*.h ../*/*.cpp
  } |
  grep -o "__PLUMED_HAS_[A-Z_]*" |
  sed "s/__PLUMED_HAS_//" | sort | uniq | tr A-Z a-z)"

# now get list of -D__PLUMED_HAS flags in Makefile.conf
has="$(
for flag
do
  case "$flag" in
  (-D__PLUMED_HAS_*=0) ;;
  (-D__PLUMED_HAS_*=*) echo ${flag#-D} | sed "s/=.*//" ;;
  (-D__PLUMED_HAS_*)   echo ${flag#-D} ;;
  (*) ;;
  esac
done | sed "s/__PLUMED_HAS_//" | sort | uniq | tr A-Z a-z
)"

# now get a list of other defines in Makefile.conf
def="$(
for flag
do
  case "$flag" in
  (-D__PLUMED_HAS_*) ;;
  (-D*) echo "define ${flag#-D}" ;;
  (*) ;;
  esac
done | sort | uniq
)"

modules=$(
cd ../
for dir in *
do
  if test -f "$dir/module.type"
  then
    mtype="$(cat "$dir/module.type")"
    is=""
    case "$mtype" in
    (always) is=on ;;
    (default-on)
      if test -f $dir.off ; then
         is=off
      else
         is=on
      fi  ;;
    (default-off)
      if test -f $dir.on ; then
         is=on
      else
         is=off
      fi  ;;
    esac
    echo "module $dir $is ($mtype)"
  fi
done
)

{
echo "# version strings"
echo "# syntax: version short/long/git number"
echo "$version"
echo
echo "# python executable"
echo "# syntax: python_bin executable"
echo "# empty string means that python has not been configured"
echo "python_bin $python_bin"
echo
echo "# command to lauch mpi processes"
echo "# syntax: mpiexec command"
echo "# empty string means that mpiexec was not chosen at configure time"
echo "mpiexec $mpiexec"
echo
echo "# list of 'has' options"
echo "# syntax: has name on/of"
echo "# if option xx is on then plumed has beeen compiled with -D__PLUMED_HAS_XX"
{
  for d in $has
  do
    echo "has $d on"
  done
  for u in $list
  do
    found=ko
    for d in $has
    do
      if test "$d" = "$u" ; then
        found=ok
      fi
    done
    if test "$found" = ko ; then
      echo "has $u off"
    fi
  done
} | sort
echo
echo "# other defines"
echo "# syntax: define name=value"
echo "$def"
echo
echo "# list of modules"
echo "# syntax: module name on/off (default-on/default-off/always)"
echo "$modules" | sort
echo
echo "# Makefile.conf file"
echo "# syntax: makefile_conf followed by a single space followed by a line from makefile_conf"
cat "$makefile_conf" | awk '{printf("makefile_conf %s\n",$0)}'

}> $file~

cmp -s $file~ $file || cp $file~ $file
rm $file~



