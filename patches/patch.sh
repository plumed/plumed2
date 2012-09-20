#! /bin/bash

MANUAL="\

Actions (choose one):
  -h, --help
                    print this help and exit
  -p, --patch
                    patch
  -r, -R, --revert
                    revert
  -l, --list-engines
                    print a list of available MD engines
  -s, --save
                    save, this needs *.preplumed files (*)
  -n NEWENGINE, --new NEWENGINE
                    create a new patch named NEWENGINE (*)
Options:
  -e ENGINE, --engine ENGINE
                    set MD engine to ENGINE (default: choose interactively)
  -m MODE, --mode MODE (default: static)
                    set link mode to MODE, which can be either static, shared or runtime
  --static
                    same as --mode static
  --shared
                    same as --mode shared
  --runtime
                    same as --mode runtime
  -d FILE, --diff FILE
                    set the path to diff file (default: ROOT/patches/ENGINE.diff) (*)
  -f, --force
                    force patching (*)

(*) These options are for developers or for testing only. Be sure to know what
    you are doing before you use them.
"

prefix=""
action=""
engine=""
diff=""
mode=static
force=""
newpatch=

multiple_actions=

for option
do

  prefix_option="$prefix$option"
  prefix=""

  case "$prefix_option" in
    (--help|-h)         echo "$MANUAL" ; exit ;;
    (--patch|-p)        test -n "$action" && multiple_actions=yes ; action=patch ;;
    (--save|-s)         test -n "$action" && multiple_actions=yes ; action=save ;;
    (--revert|-R|-r)    test -n "$action" && multiple_actions=yes ; action=revert ;;
    (--list-engines|-l) test -n "$action" && multiple_actions=yes ; action=list ;;
    (--new=*)           test -n "$action" && multiple_actions=yes ; action=new ; newpatch="${prefix_option#--new=}" ;;
    (--description)     echo "patch an MD engine" ; exit ;;
    (--engine=*) engine="${prefix_option#--engine=}" ;;
    (--mode=*) mode="${prefix_option#--mode=}" ;;
    (--diff=*) diff="${prefix_option#--diff=}" ;;
    (--engine|-e) prefix="--engine=" ;;
    (--root) prefix="--root=" ;;
    (--diff|-d) prefix="--diff=" ;;
    (--mode|-m) prefix="--mode=" ;;
    (--new|-n) prefix="--new=" ;;
    (--static) mode=static ;;
    (--shared) mode=shared ;;
    (--runtime) mode=runtime ;;
    (--force|-f) force=yes ;;
    (*)
      echo "ERROR: Unknown option $prefix_option. Use -h for help."
      exit
  esac
done

if [ -n "$multiple_actions" ] ; then
  echo "Too many actions. -h for help"
  exit
fi

if [ -z "$action" ] ; then
  echo "Nothing to do. -h for help"
  exit
fi

echo "PLUMED patching tool"
echo
if [ -z "$PLUMED_ROOT" ]
then
  echo "ERROR: I cannot find PLUMED"
  echo "Please set PLUMED_ROOT environment variable or use --root"
  exit
fi
if [ ! -d "$PLUMED_ROOT/patches/" ]
then
  echo "ERROR: cannot find $PLUMED_ROOT/patches/ directory"
  echo "Check your PLUMED_ROOT variable or --root option"
  exit
fi

# build MD engines list

mdlist=""
for file in "$PLUMED_ROOT"/patches/*diff
do
  b="${file##*/}"
  b="${b%.diff}"
  mdlist="$mdlist $b"
done

if [ "$action" == list ]
then
  echo "Available MD engines:"
  for file in "$PLUMED_ROOT"/patches/*diff
  do
    b="${file##*/}"
    b="${b%.diff}"
    echo "  $b"
  done
  exit
fi

if [ "$action" == new ]
then
  echo "Creating a new patch"
  if [[ -e "$PLUMED_ROOT"/patches/"$newpatch".diff ]] ; then
      echo "ERROR: patch $newpatch is already defined"
      exit
  fi
  touch "$PLUMED_ROOT"/patches/"$newpatch".diff
  echo "Created file $PLUMED_ROOT/patches/$newpatch.diff"
  echo "Also consider the possibility of adding a $PLUMED_ROOT/patches/$newpatch.config file"
  exit
fi

if [ -z "$engine" ]
then
  PS3="Choose the best matching code/version:"
  select engine in $mdlist; do
    if [[ -n "$engine" ]] ; then
      break
    else
      echo "ERROR: choose in the list above or interrupt (^c)"
    fi
  done
fi

if [ -z "$diff" ]
then
  diff="$PLUMED_ROOT/patches/${engine}.diff"
fi
if [ -z "$config" ]
then
  config="$PLUMED_ROOT/patches/${engine}.config"
fi

echo "MD engine: $engine"
echo "PLUMED location: $PLUMED_ROOT"
echo "diff file: $diff"

if [ -f "$config" ]
then
  echo "sourcing config file: $config"
  source "$config"
fi

case "$mode" in
(static|shared|runtime) ;;
(*)
  echo "I don't understand mode $mode"
  exit
esac


case "$action" in
  (patch)
    if [ ! -f "$diff" ] ; then
      echo "ERROR: MD engine not supported (or mispelled)"
      exit
    fi
    if type -t plumed_preliminary_test 1>/dev/null ; then
      if plumed_preliminary_test || [ "$force" ] ; then
        echo >/dev/null
      else
        echo "ERROR: Preliminary test not passed."
        echo "It seems that this is not $engine, or you are in the wrong directory"
        echo "If you are sure about what you are doing, use -f"
      exit
      fi
    fi
    if [ -L Plumed.h -o -L Plumed.inc ]
    then
      echo "ERROR: you have likely already patched. Revert first (-r)"
      exit
    fi
    if [ ! -f "$PLUMED_ROOT/src/Plumed.inc" ]
    then
      echo "ERROR: cannot find $PLUMED_ROOT/src/Plumed.inc file"
      echo "Compile plumed before patching"
      exit
    fi
    echo "Linking Plumed.h and Plumed.inc ($mode mode)"
    ln -s "$PLUMED_ROOT/src/Plumed.h"
    ln -s "$PLUMED_ROOT/src/Plumed.inc.$mode" Plumed.inc
    bash "$diff"
  ;;
  (save)
    if [ ! -L Plumed.h -o ! -L Plumed.inc ]
    then
      echo "ERROR: I cannot find Plumed.h and Plumed.inc files. You have likely not patched yet."
      exit
    fi
    PREPLUMED=$(find . -name "*.preplumed")
    if ! test "$PREPLUMED" ; then
      echo "ERROR: I cannot find any .preplumed file. You have likely not patched yet."
      exit
    fi
    if type -t plumed_preliminary_test 1>/dev/null ; then
      if plumed_preliminary_test || [ "$force" ] ; then
        echo >/dev/null
      else
        echo "ERROR: Preliminary test not passed."
        echo "It seems that this is not $engine, or you are in the wrong directory"
        echo "If you are sure about what you are doing, use -f"
      exit
      fi
    fi
    echo "Saving your changes to $diff"
    echo "Preplumed files:"
    echo "$PREPLUMED"
    test -e "$diff" && rm "$diff"
    for bckfile in $PREPLUMED ; do
      file="${bckfile%.preplumed}"
      if test -e "$file" ; then
        echo "patch -u -l -b -F 5 --suffix=.preplumed \"${file}\" << \\EOF_EOF" >> "$diff"
        diff -U 5 "${bckfile}" "$file" --label="$bckfile" --label="$file" >> "$diff"
        echo "EOF_EOF"                                                   >> "$diff"
      else
        echo "ERROR: File $file is missing"
      fi
    done
  ;;
  (revert)
    if [ ! -f "$diff" ] ; then
      echo "ERROR: MD engine not supported (or mispelled)"
      exit
    fi
    if [ ! -L Plumed.h -o ! -L Plumed.inc ]
    then
      echo "WARNING: I cannot find Plumed.h and Plumed.inc files. You have likely not patched yet."
    else
    echo "Removing Plumed.h and Plumed.inc"
      rm Plumed.h Plumed.inc
    fi
    PREPLUMED=$(find . -name "*.preplumed")
    if ! test "$PREPLUMED" ; then
      echo "WARNING: I cannot find any .preplumed file. You have likely not patched yet."
    else
      echo "Reverting changes and touching reverted files"
      for bckfile in $PREPLUMED ; do
        file="${bckfile%.preplumed}"
        mv "$bckfile" "$file"
        touch "$file"
      done
    fi
esac



