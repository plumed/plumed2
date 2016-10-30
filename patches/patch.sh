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
  --save-originals
                    same as save, but save also original files (*)
  -n NEWENGINE, --new NEWENGINE
                    create a new patch named NEWENGINE (*)
  -i, --info
                    output information on the patching procedure for a particular code                
Options:
  -e ENGINE, --engine ENGINE
                    set MD engine to ENGINE (default: choose interactively)
  -m MODE, --mode MODE (default: shared)
                    set link mode to MODE, which can be either static, shared or runtime
  --static
                    same as --mode static
  --shared
                    same as --mode shared
  --runtime
                    same as --mode runtime
  -d FILE, --diff FILE
                    set the path to diff file (default: ROOT/patches/ENGINE.diff) (*)
  -q, --quiet
                    do not write loggin information; useful with -i to print just
                    the patching information
  -f, --force
                    force patching (*)

(*) These options are for developers or for testing only. Be sure to know what
    you are doing before you use them.
"

prefix=""
action=""
engine=""
diff=""
mode=shared
force=""
newpatch=

multiple_actions=

otherfiles=
save_originals=
quiet=
mdroot=

for option
do

  prefix_option="$prefix$option"
  prefix=""

  case "$prefix_option" in
    (--help|-h)         echo "$MANUAL" ; exit ;;
    (--patch|-p)        test -n "$action" && multiple_actions=yes ; action=patch ;;
    (--save|-s)         test -n "$action" && multiple_actions=yes ; action=save ;;
    (--save-originals)  test -n "$action" && multiple_actions=yes ; action=save ; save_originals=yes ;;
    (--revert|-R|-r)    test -n "$action" && multiple_actions=yes ; action=revert ;;
    (--list-engines|-l) test -n "$action" && multiple_actions=yes ; action=list ;;
    (--info|-i)         test -n "$action" && multiple_actions=yes ; action=info ;;
    (--new=*)           test -n "$action" && multiple_actions=yes ; action=new ; newpatch="${prefix_option#--new=}" ;;
    (--description)     echo "patch an MD engine" ; exit ;;
    (--engine=*) engine="${prefix_option#--engine=}" ;;
    (--mdroot=*) mdroot="${prefix_option#--mdroot=}" ;;
    (--mode=*) mode="${prefix_option#--mode=}" ;;
    (--diff=*) diff="${prefix_option#--diff=}" ;;
    (--engine|-e) prefix="--engine=" ;;
    (--mdroot) prefix="--mdroot" ;;
    (--root=*) prefix="--root="; PLUMED_ROOT="${prefix_option#--root=}" ;;
    (--diff|-d) prefix="--diff=" ;;
    (--mode|-m) prefix="--mode=" ;;
    (--new|-n) prefix="--new=" ;;
    (--static) mode=static ;;
    (--shared) mode=shared ;;
    (--runtime) mode=runtime ;;
    (--force|-f) force=yes ;;
    (--quiet|-q) quiet=yes ;;
    (*)
      echo "ERROR: Unknown option $prefix_option. Use -h for help."
      exit
  esac
done

if [ -n "$mdroot" ] ; then
  if ! cd "$mdroot" ; then
    echo "Directory $mdroot does not exist"
    exit
  fi
fi

if [ -n "$multiple_actions" ] ; then
  echo "Too many actions. -h for help"
  exit
fi

if [ -z "$action" ] ; then
  echo "Nothing to do. -h for help"
  exit
fi

test -n "$quiet" || echo "PLUMED patching tool"
test -n "$quiet" || echo
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
  test -n "$quiet" || echo "Creating a new patch"
  if [[ -e "$PLUMED_ROOT"/patches/"$newpatch".diff ]] ; then
      echo "ERROR: patch $newpatch is already defined"
      exit
  fi
  touch "$PLUMED_ROOT"/patches/"$newpatch".diff
  test -n "$quiet" || echo "Created file $PLUMED_ROOT/patches/$newpatch.diff"
  test -n "$quiet" || echo "Also consider the possibility of adding a $PLUMED_ROOT/patches/$newpatch.config file"
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
if [ -z "$otherfiles" ]
then
  test -d "$PLUMED_ROOT/patches/${engine}" && otherfiles="$PLUMED_ROOT/patches/${engine}/"
fi

test -n "$quiet" || echo "MD engine: $engine"
test -n "$quiet" || echo "PLUMED location: $PLUMED_ROOT"
test -n "$quiet" || echo "diff file: $diff"

if [ -f "$config" ]
then
  test -n "$quiet" || echo "sourcing config file: $config"
  source "$config"
fi

if [ -d "$otherfiles" ]
then
  test -n "$quiet" || echo "extra files located in: $otherfiles"
fi

case "$mode" in
(static|shared|runtime) ;;
(*)
  echo "I don't understand mode $mode"
  exit
esac


case "$action" in
  (patch)
    if [ ! -e "$diff" ] ; then
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
    if [ ! -f "$PLUMED_ROOT/src/lib/Plumed.inc" ]
    then
      echo "ERROR: cannot find $PLUMED_ROOT/src/lib/Plumed.inc file"
      echo "Compile plumed before patching"
      exit
    fi
    if [ ! -f "$PLUMED_ROOT/src/lib/Plumed.cmake.$mode" ]
    then
      echo "ERROR: cannot find $PLUMED_ROOT/src/lib/Plumed.cmake.$mode file"
      echo "Compile a $mode version of plumed before patching, or change patching mode [static|shared|runtime]"
      exit
    fi
    if type -t plumed_before_patch 1>/dev/null ; then
      test -n "$quiet" || echo "Executing plumed_before_patch function"
      plumed_before_patch
    fi
    test -n "$quiet" || echo "Linking Plumed.h and Plumed.inc ($mode mode)"
    ln -s "$PLUMED_INCLUDEDIR/$PLUMED_PROGRAM_NAME/wrapper/Plumed.h" Plumed.h
    ln -s "$PLUMED_ROOT/src/lib/Plumed.inc.$mode" Plumed.inc
    ln -s "$PLUMED_ROOT/src/lib/Plumed.cmake.$mode" Plumed.cmake

    if [ -d "$diff" ]; then
      test -n "$quiet" || echo "Patching with on-the-fly diff from stored originals"
      PREPLUMED=$(cd "$diff" ; find . -name "*.preplumed" | sort)
      for bckfile in $PREPLUMED ; do
        file="${bckfile%.preplumed}"
        if test -e "$file" ; then
          diff -U 5 "$diff/$bckfile" "$diff/$file" --label="$bckfile" --label="$file" |
          patch -u -l -b -F 5 --suffix=.preplumed "$file"
        else
          echo "ERROR: File $file is missing"
        fi
      done
    else
      test -n "$quiet" || echo "Patching with stored diff"
      bash "$diff"
    fi

    if [ "$PLUMED_IS_INSTALLED" = no  ] && [ "$mode" = shared ] ; then
      echo ""
      echo "You are patching in shared mode from a non installed PLUMED"
      echo "Be warned that if you 'make clean' PLUMED the patched code won't work anymore"
    fi

    if [ "$mode" = runtime ] ; then
      echo ""
      echo "You are patching in runtime mode"
      echo "Be warned that when you will run MD you will use the PLUMED version pointed at"
      echo "by the PLUMED_KERNEL environment variable"
    fi

    echo ""
    if grep -q "D__PLUMED_HAS_MPI=1" "$PLUMED_ROOT"/src/config/compile_options.sh ; then
      echo "PLUMED is compiled with MPI support so you can configure $engine with MPI" 
    else
      echo "PLUMED is compiled WITHOUT MPI support so you CANNOT configure $engine with MPI"
    fi

    
    if type -t plumed_after_patch 1>/dev/null ; then
      test -n "$quiet" || echo "Executing plumed_after_patch function"
      plumed_after_patch
    fi

  ;;
  (info)
    if type -t plumed_patch_info 1>/dev/null ; then
      test -n "$quiet" || echo "Executing plumed_patch_info function"
      plumed_patch_info
    else
      echo "No special info for this MD engine"
    fi
  ;;
  (save)
    if [ ! -L Plumed.h -o ! -L Plumed.inc ]
    then
      echo "ERROR: I cannot find Plumed.h and Plumed.inc files. You have likely not patched yet."
      exit
    fi
    PREPLUMED=$(find . -name "*.preplumed" | sort)
    if ! test "$PREPLUMED" ; then
      echo "ERROR: I cannot find any .preplumed file. There is nothing to save."
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
    test -n "$quiet" || echo "Saving your changes to $diff"
    test -n "$quiet" || echo "Preplumed files:"
    test -n "$quiet" || echo "$PREPLUMED"
    if [ -d "$diff" ] && [ -z "$save_originals" ]; then
      echo "This patch uses the dir format (originals are saved)"
      echo "Are you sure you want to save the single diff?"
      echo "Possible reasons to do it are:"
      echo "* because of licence you do not want to store originals in plumed"
      echo "* you do not want to do a big change on the diff files now"
      answer=
      PS3="Choose:"
      select answer in single-diff originals ; do
        if [[ -n "$answer" ]] ; then
          break
        else
          echo "ERROR: choose in the list above or interrupt (^c)"
        fi
      done
      if [ "$answer" = originals ] ; then
        test -n "$quiet" || echo "saving originals"
        save_originals=yes
      else
        test -n "$quiet" || echo "saving single diff file"
      fi
    fi
    test -e "$diff" && rm -r "$diff"
    for bckfile in $PREPLUMED ; do
      file="${bckfile%.preplumed}"
      if test -e "$file" ; then
        if [ -n "$save_originals" ] ; then
          xx="$diff/$file"
          mkdir -p "${xx%/*}"
          cp "$file" "$diff/$file"
          cp "$bckfile" "$diff/$bckfile"
        else
          echo "patch -u -l -b -F 5 --suffix=.preplumed \"${file}\" << \\EOF_EOF" >> "$diff"
          diff -U 5 "${bckfile}" "$file" --label="$bckfile" --label="$file" >> "$diff"
          echo "EOF_EOF"                                                   >> "$diff"
        fi
      else
        echo "ERROR: File $file is missing"
      fi
    done
cat <<EOF
* If you want your patch to perform some arbitrary action before/after
* patching the diff files, just add a function named
* plumed_before_patch/plumed_after_patch to the ${diff%diff}config file.
* Do not forget to also add an equivalent plumed_before_revert/plumed_after_revert
* function
EOF
  ;;
  (revert)
    if [ ! -e "$diff" ] ; then
      echo "ERROR: MD engine not supported (or mispelled)"
      exit
    fi
    if type -t plumed_before_revert 1>/dev/null ; then
      test -n "$quiet" || echo "Executing plumed_before_revert function"
      plumed_before_revert
    fi
    if [ ! -L Plumed.h -o ! -L Plumed.inc ]
    then
      echo "WARNING: I cannot find Plumed.h and Plumed.inc files. You have likely not patched yet."
    else
    test -n "$quiet" || echo "Removing Plumed.h and Plumed.inc"
      rm Plumed.h Plumed.inc Plumed.cmake
    fi
    PREPLUMED=$(find . -name "*.preplumed")
    if ! test "$PREPLUMED" ; then
      echo "No .preplumed file found, nothing to restore."
    else
      test -n "$quiet" || echo "Reverting changes and touching reverted files"
      for bckfile in $PREPLUMED ; do
        file="${bckfile%.preplumed}"
        mv "$bckfile" "$file"
        touch "$file"
      done
    fi
    if type -t plumed_after_revert 1>/dev/null ; then
      test -n "$quiet" || echo "Executing plumed_after_revert function"
      plumed_after_revert
    fi
esac



