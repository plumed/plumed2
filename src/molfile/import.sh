#!/bin/bash

# A script to import a subset of molfile plugins into plumed2's tree

if (( $# != 1 ))
then
  cat <<EOF
Usage: ./import.sh <vmd_plugins_dir>  

 <vmd_plugins_dir>  should contain e.g. molfile_plugin/src/Gromacs.h
EOF
  exit 1
fi



function mycp {
	f=$1 
	t=$2
	echo "Copying $f into $t  ..."
	cp $f $t
}

# Note that a directory name is not allowed as second argument
function mycp_wrap {
	i=$1
	f=$2 
	t=$3
	echo "Copying $f into $t (wrapping into #ifndef $i) ..."
	(echo "#ifndef $i";
         cat $f;
	 echo "#endif";) > $t
}



#      <plugins_dir> should be the one containing dcdplugin.cpp
#                    usually .../plugins/molfile_plugin/src
#      <header_dir>  should contain molfile_plugin.h, vmdplugin.h
#		    usually .../plugins/include

PD="$1"
IFNDEF=__PLUMED_HAS_EXTERNAL_MOLFILE_PLUGINS

# List of files to import from molfile_plugin/src
for i in endianswap.h fastio.h Gromacs.h largefiles.h periodic_table.h readpdb.h; do
	mycp $PD/molfile_plugin/src/$i .
done

# List of files to import from include/
for i in molfile_plugin.h vmdplugin.h; do
	mycp $PD/include/$i .
done

mycp $PD/molfile_plugin/LICENSE COPYRIGHT

# List of "known-good" plugins. Some renaming is necessary
mycp_wrap $IFNDEF $PD/molfile_plugin/src/dcdplugin.c dcdplugin.cpp
mycp_wrap $IFNDEF $PD/molfile_plugin/src/gromacsplugin.C gromacsplugin.cpp
mycp_wrap $IFNDEF $PD/molfile_plugin/src/pdbplugin.c pdbplugin.cpp


# Generate static header
plugins="dcdplugin gromacsplugin pdbplugin"
echo "Generating libmolfile_plugin.h: $plugins"
rm libmolfile_plugin.h
./create_static_header.sh MOLFILE molfile libmolfile_plugin.h $plugins

# Add copyright notice at top
echo "Adding COPYRIGHT headers"
(cd ..; ./header.sh molfile > /dev/null)


