#! /bin/bash

VIMFILE="$PLUMED_ROOT"/src/lib/plumed.vim

if [ "$PLUMED_IS_INSTALLED" = yes ] ; then
  VIMFILE="$PLUMED_ROOT/vim/syntax/$PLUMED_PROGRAM_NAME".vim
fi

if [ "$1" = --description ] ; then
  echo "convert plumed input file to colored html using vim syntax"
  exit 0
fi

if [ "$1" = --help ] ; then
cat << EOF
Usage:
  plumed vim2html [options] < input > output
  plumed vim2html [options] input > output
  plumed vim2html [options] input output
(one of the three)

Options can be:
--annotate-syntax Reports annotated syntax on output.
                  Also reports non-annotated regions on stderr

EOF
  exit 0
fi

if [ ! -f "$VIMFILE" ] ; then
  echo "Cannot find vimfile at $VIMFILE"
  exit 1
fi

annotate_syntax=no
inputfile=""
outputfile=""

for opt
do
  case "$opt" in
    (--annotate-syntax)
      annotate_syntax=yes ;;
    (-*)
      echo "ERROR: Unknown option $opt. Use -h for help."
      exit 1 ;;
    (*)
      if [ -z "$inputfile" ] ; then
        inputfile="$opt"
      elif [ -z "$outputfile" ] ; then
        outputfile="$opt" ;
      else
        echo "ERROR: too many arguments"
        exit 1
      fi
  esac
done

temp="$( mktemp -dt plumed.XXXXX)"

{
if [ -n "$inputfile" ] ; then
  cat "$inputfile"
else
  cat
# this is to remove \r characters
fi | tr '\r' ' '
} | {
cd $temp
{
if [ "$annotate_syntax" = "yes" ]; then

vim - -c ":syntax off" \
      -c ":source $VIMFILE" \
      -c ":call PlumedAnnotateSyntax()" \
      -c ":wq! syntax.log" -c ":q!" 1> /dev/null 2> /dev/null

else
# this is a workaround for a bug in TOhtml
# that does not color properly when there are empty lines
sed "s/^$/ /" |
vim - -c ":syntax off" \
      -c ":source $VIMFILE" \
      -c ":let g:html_use_css = 0" \
      -c ":let g:html_no_progress" \
      -c ":let g:html_number_lines = 1" \
      -c ":TOhtml" \
      -c ":wq! out.html" -c ":q!" 1> /dev/null 2> /dev/null
fi
}
}

if [ "$annotate_syntax" = "yes" ] ; then

output="$(awk '{if($1=="ANNOTATION") print}' $temp/syntax.log 2>/dev/null)"
error="$(awk '{if($1=="ERROR") print}' $temp/syntax.log >&2 2>/dev/null)"

else

output="$(cat $temp/out.html | sed 's|^<title>.*</title>$||' | sed 's|^<meta.*$||')"

fi

if [ -n "$outputfile" ] ; then
  echo "$output" > "$outputfile"
else
  echo "$output"
fi

echo "$error" >&2

if [ -n "$error" ] ; then
  exit 1
fi





