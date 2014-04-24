#!/bin/sh

# Script for creating header file for static libraries

libname=$1
shift
prefix=$1
shift
target=$1
shift
plugins=$*

#
# boilerplate header
#
echo "#ifndef LIB${libname}_PLUGIN_H" >> $target
echo "#define LIB${libname}_PLUGIN_H" >> $target

#
# include the generic vmd plugin header
#
echo "#include \"vmdplugin.h\"" >> $target
echo "" >> $target

#
# all plugin API's are extern "C".
#
echo "#ifdef __cplusplus" >> $target
echo "extern \"C\" {" >> $target
echo "#endif" >> $target
echo "" >> $target

#
# Function declarations
#
for p in $plugins
do
##  XXX not legal in standard Bourne shell, so we do it in the Makefile now
##  name=${p%%.so}
  echo "extern int ${prefix}_${p}_init(void);" >> $target
  echo "extern int ${prefix}_${p}_register(void *, vmdplugin_register_cb);" >> $target
  echo "extern int ${prefix}_${p}_fini(void);" >> $target
done
echo "" >> $target

# macros for init, register, fini
echo "#define ${libname}_INIT_ALL \\" >> $target
for p in $plugins
do
##  XXX not legal in standard Bourne shell, so we do it in the Makefile now
##  name=${p%%.so}
  echo "    ${prefix}_${p}_init(); \\" >> $target
done
echo "" >> $target

echo "#define ${libname}_REGISTER_ALL(v, cb) \\" >> $target
for p in $plugins
do
##  XXX not legal in standard Bourne shell, so we do it in the Makefile now
##  name=${p%%.so}
  echo "    ${prefix}_${p}_register(v, cb); \\" >> $target
done
echo "" >> $target

echo "#define ${libname}_FINI_ALL \\" >> $target
for p in $plugins
do
##  XXX not legal in standard Bourne shell, so we do it in the Makefile now
##  name=${p%%.so}
  echo "    ${prefix}_${p}_fini(); \\" >> $target
done
echo "" >> $target

#
# wrap it up
#
echo "#ifdef __cplusplus" >> $target
echo "}" >> $target
echo "#endif" >> $target
echo "#endif" >> $target
