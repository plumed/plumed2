# Adding a minimal subset of molfile plugins from VMD plugins
# These should hopefully build on all supported platforms.

#
# ./configure --host=i686-w64-mingw32 --disable-shared LDFLAGS="-static"

USE=



# Instead of disabling code via #ifdefs, I disable the make rules
ifdef $(findstring __PLUMED_INTERNAL_MOLFILE_PLUGINS,$(CPPFLAGS))) 

include ../maketools/make.module

#OBJS=gromacsplugin.o dcdplugin.o

CPPFLAGS+=-DSTATIC_PLUGIN

%.o: %.c
	$(CC) -o $@ -c $(CFLAGS) $(CPPFLAGS)  $< -DVMDPLUGIN=molfile_$(patsubst %.o,%,$@)

%.o: %.cpp
	$(CXX) -o $@ -c $(CFLAGS) $(CPPFLAGS)  $< -DVMDPLUGIN=molfile_$(patsubst %.o,%,$@)

else

# Make links and obj targets do nothing

$(info  Not building internal molfile plugins)

links:
	true

obj:
	true

endif



