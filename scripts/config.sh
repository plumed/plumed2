#! /bin/bash

MANUAL="\

Options:
  -h, --help
                    print this help and exit
  -q, --quiet
                    don't write anything, just return true of false
  show
                    dump a full configuration file
  has [word1 [word2]..]
                    check if plumed has features words
  module [word1 [word2]..]
                    check if plumed has enables modules words
  makefile_conf
                    dumps the Makefile.conf file

Examples:

Check if plumed as xdrfile enabled
> plumed config has xdrfile
Check if plumed as xdrfile AND zlib enabled
> plumed config has xdrfile zlib
Check if plumed as module colvar active
> plumed config module colvar
"

# notice that the relative path of config.txt is also
# hardcoded in a comment written in the log from src/core/PlumedMain.cpp
# if you change it here, also change it there!
configfile="$PLUMED_ROOT"/src/config/config.txt

quiet=no
list=no
checklist=""
for opt
do
  case "$opt" in
  (--help|-h) echo "$MANUAL" ; exit ;;
  (--description)
    echo "inquire plumed about how it was configure"
    exit 0
    ;;
  (--options)
    echo "--help -h --description --options --quiet -q --version -v show has module mpiexec makefile_conf python_bin"
    exit 0
    ;;
  (--quiet|-q) quiet=yes ;;
  (--version|-v) action=version ;;
  (show)
    cat "$configfile"
    exit 0
    ;;
  (has) action=has ;;
  (module) action=module ;;
  (python_bin) action=python_bin ;;
  (mpiexec) action=mpiexec ;;
  (makefile_conf)
    cat "$configfile" | awk '{if($1=="makefile_conf") { gsub("^makefile_conf ",""); print} }'
    exit 0
    ;;
  (*)
    checklist="$checklist $opt"
  esac
done

if test -z "$checklist" && ( test "$action" = show  || test  "$action" = has || test "$action" = module)
then
  test "$quiet" = no && echo "nothing to do"
  exit 1
fi

case $action in
(has|module)
  retval=0
  for check in $checklist
  do
  
  ff=$(
    cat "$configfile" |
    awk -v action="$action" -v check="$check" '{ if($1==action && $2==check){ print $3;exit } }'
  )
  
  if test "$ff" = on ; then
    test "$quiet" = no && echo "$check on"
  elif test "$ff" = off ; then
    test "$quiet" = no && echo "$check off"
    retval=1
  else
    test "$quiet" = no && echo "$check not found"
    retval=1
  fi
  
  done
  
  exit $retval
 ;;
(version)
  long=$(cat "$configfile" | grep -v \# | awk '{ if($1=="version" && $2=="long") print $3 }')
  git=$(cat "$configfile" | grep -v \# | awk '{ if($1=="version" && $2=="git") print $3 }')
  echo "Version: $long (git: $git)"
 ;;
(python_bin)
  py=$(cat "$configfile" | grep -v \# | awk '{ if($1=="python_bin") print $2 }')
  if test -n "$py" ; then
    retval=0
    test "$quiet" = no && echo "$py"
  else
    retval=1
    test "$quiet" = no && echo "python not found"
  fi
  exit $retval
;;
(mpiexec)
  mpi=$(cat "$configfile" | grep -v \# | awk '{ if($1=="mpiexec") { sub(" *mpiexec ","",$0);  print} }')
  if test -n "$mpi" ; then
    retval=0
    test "$quiet" = no && echo "$mpi"
  else
    retval=1
    test "$quiet" = no && echo "mpiexec not found"
  fi
  exit $retval
esac

