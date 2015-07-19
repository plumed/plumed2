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

Examples:

Check if plumed as xdrfile enabled
> plumed config has xdrfile
Check if plumed as xdrfile AND zlib enabled
> plumed config has xdrfile zlib
Check if plumed as module colvar active
> plumed config module colvar
"

configfile="$(cat "$PLUMED_ROOT"/src/config/config.txt)"

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
  (--quiet|-q) quiet=yes ;;
  (--version|-v) action=version ;;
  (show)
    echo "$configfile"
    exit 0
    ;;
  (has) action=has ;;
  (module) action=module ;;
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
    echo "$configfile" |
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
  long=$(echo "$configfile" | grep -v \# | awk '{ if($1=="version" && $2=="long") print $3 }')
  git=$(echo "$configfile" | grep -v \# | awk '{ if($1=="version" && $2=="git") print $3 }')
  echo "Version: $long (git: $git)"
 ;;
esac

