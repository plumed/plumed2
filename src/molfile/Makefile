# Adding a minimal subset of molfile plugins from catdcd-5.1
# These should hopefully build on all supported platforms.

#
# ./configure --host=i686-w64-mingw32 --disable-shared LDFLAGS="-static"

USE=

include ../maketools/make.module

OBJS=gromacsplugin.o dcdplugin.o

CFLAGS+=-DSTATIC_PLUGIN

default: $(OBJS) libmolfile_plugin.h

libmolfile_plugin.h: $(OBJS)
	./create_static_header.sh MOLFILE molfile $@ $(patsubst %.o,%,$^)


%.o: %.c
	$(CC) -o $@ -c $(CFLAGS) $< -DVMDPLUGIN=molfile_$(patsubst %.o,%,$@)

%.o: %.cpp
	$(CXX) -o $@ -c $(CFLAGS) $< -DVMDPLUGIN=molfile_$(patsubst %.o,%,$@)



# generic makefile

