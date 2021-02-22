#! /bin/bash

# bck.meup.sh creates a plumed-like backup of the given file(s)
# Options:
#  -i informs if the file is backed up
#  -v verbose (-i + tells if files do not exist)
#  -h shows usage

#check for options
info=false
verbose=false
if [ "$1" = -i ]
then
  info=true
elif [ "$1" = -v ]
then
  info=true
  verbose=true
elif [ "$1" = -h ]
then
  echo "  USAGE: $0 [-i/-v] [files to be backed up]"
fi

#do the backup
for file in "$@"
do
  if [ -f "$file" ] || [ -d "$file" ]
  then
    file=${file%/}
    name=${file##*/}
    path=${file%${name}}
    i=0
    bck_file=${path}bck.${i}.${name}
    while [ -f "$bck_file" ] || [ -d "$bck_file" ]
    do
      i=$[$i+1]
      bck_file=${path}bck.${i}.${name}
    done
    mv $file $bck_file
    if [ "$info" = true ]
    then 
      echo "--> ${0##*/}:  mv $file $bck_file"
    fi
  else
    if [ "$verbose" = true ] && [ "$file" != -v ]
    then 
      echo "--> ${0##*/}:  file \"$file\" not found"
    fi
  fi
done
