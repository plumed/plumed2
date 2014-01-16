USE=core wrapper config tools

# generic makefile
include ../maketools/make.module

# Use the correct (internal/external) libmolfile_plugin.h header. 
# Done to avoid modifying VMD's headers with #ifdefs. TG
ifneq (,$(findstring __PLUMED_INTERNAL_MOLFILE_PLUGINS,$(CPPFLAGS)))
$(info Adding internal molfile include path)
CPPFLAGS+=-I../molfile
endif

