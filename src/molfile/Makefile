# Adding a minimal subset of molfile plugins from VMD plugins
# These should hopefully build on all supported platforms.

# Instead of disabling code via #ifdefs, I disable the make rules
# ifdef $(findstring __PLUMED_HAS_INTERNAL_MOLFILE_PLUGINS,$(CPPFLAGS))) 

#
# ./configure --host=i686-w64-mingw32 --disable-shared LDFLAGS="-static"

USE=

include ../maketools/make.module

# Local rules for molfile plugins - replacement is delayed 
CPPFLAGS+=-DSTATIC_PLUGIN  -DVMDPLUGIN=molfile_$(patsubst %.o,%,$@)

# -fpermissive works around a missing cast in pdbplugin 
#          char *newstr=realloc(pdb->meta->remarks, newlen + 1);


