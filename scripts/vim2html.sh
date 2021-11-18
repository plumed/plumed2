#! /bin/bash

VIMFILE="$PLUMED_ROOT/vim/syntax/plumed.vim"

if [ "$1" = --description ] ; then
  echo "convert plumed input file to colored html using vim syntax"
  exit 0
fi

if [ "$1" = --options ] ; then
  echo "--description --options --annotate-syntax --pdf --crop --fs --colors"
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
--pdf             To produce a pdf file with syntax highlighting.
--crop            Crop the pdf file to the only written part. Useful to insert the pdf in a LaTex file as image.
--fs              Specify the fontsize of the pdf output.
--colors          Specify the color palette. Allowed values are: default/ac
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
pdf=no
crop=no
fontsize=7.5
prefix=""
colors=""
for opt
do
  prefixopt="$prefix$opt"
  prefix=""
  case "$prefixopt" in
    (--annotate-syntax)
      annotate_syntax=yes ;;
    (--pdf)
      pdf=yes;;
    (--crop)
      crop=yes;;
    (--fs)
      prefix="--fs=";;
    (--fs=*)
      fontsize="${prefixopt#--fs=}";;
    (--colors)
      prefix="--colors=";;
    (--colors=*)
      colors="${prefixopt#--colors=}";;
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

case "$colors" in
(ac)
  cat > $temp/colors.vim << EOF
: set bg=dark
: hi Statement  ctermfg=5 cterm=bold guifg=#ffff60 gui=none
: highlight Type ctermfg=2 cterm=bold
: hi Comment  guifg=#80a0ff ctermfg=darkblue
: hi Constant ctermfg=red guifg=#ffa0a0 cterm=bold
EOF
;;
(*)
  echo > $temp/colors.vim
;;
esac

if [ "$pdf" = "yes" ]; then

# this is a workaround for a bug in TOhtml
# that does not color properly when there are empty lines
  sed "s/^$/ /" |
  vim - -c ":syntax off" \
        -c ":source $VIMFILE" \
        -c ":source $temp/colors.vim" \
        -c ":set printoptions=header:0" \
        -c ":set printfont=courier:h$fontsize" \
        -c ":hardcopy > out.ps" \
        -c ":q!" 1> /dev/null 2> /dev/null
else
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
      -c ":source $temp/colors.vim" \
      -c ":let g:html_use_css = 0" \
      -c ":let g:html_no_progress" \
      -c ":let g:html_number_lines = 1" \
      -c ":TOhtml" \
      -c ":wq! out.html" -c ":q!" 1> /dev/null 2> /dev/null
  fi 
fi


}
}

if [ "$pdf" = "yes" ]; then
  ps2pdf $temp/out.ps $temp/out.pdf
  output="$temp/out.pdf"
  if [ "$crop" = "yes" ]; then
    (cd $temp;     pdfcrop $temp/out.pdf $temp/outcrop.pdf > /dev/null   )
    output="$temp/outcrop.pdf"
  fi
else
  if [ "$annotate_syntax" = "yes" ] ; then
    output="$(awk '{if($1=="ANNOTATION") print}' $temp/syntax.log 2>/dev/null)"
    error="$(awk '{if($1=="ERROR") print}' $temp/syntax.log >&2 2>/dev/null)"
  else
    output="$(cat $temp/out.html | sed 's|^<title>.*</title>$||' | sed 's|^<meta.*$||')"
  fi
fi


if [ "$pdf" = "yes" ]; then
  if [ -n "$outputfile" ] ; then
    cp $output $outputfile
  else
    echo "You should specify an output file if using pdf option!"
    exit 1;
  fi
else

  if [ -n "$outputfile" ] ; then
    echo "$output" > "$outputfile"
  else
    echo "$output"
  fi

fi

echo "$error" >&2

if [ -n "$error" ] ; then
  exit 1
fi





