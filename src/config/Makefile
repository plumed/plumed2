
# include the machine dependent configuration
-include ../../Makefile.conf

# sed script to edit program name defaults to do nothing
ifndef program_transform_name
  program_transform_name=s,x,x,
endif

.PHONY: all obj
ifdef PLUMED_FORCE
.PHONY: ConfigInstall.cpp Config.cpp version.h
endif

# default target
all:
	$(MAKE) PLUMED_FORCE=yes ConfigInstall.cpp Config.cpp version.h
	$(MAKE) Config.o ConfigInstall.o compile_options.sh

obj: all

version.h:
	@echo "Updating version.h"
	@ ../maketools/update-version.sh version.h

ConfigInstall.cpp: Config.cpp.in
	@echo "Updating ConfigInstall.cpp"
	@sed "s/@SOEXT@/$(SOEXT)/g" Config.cpp.in | \
         sed "s/@ISINSTALLED@/true/g" | \
         sed "s|@PLUMED_ROOT@|$(libdir)/$(program_name)|g" | \
         sed "s|@htmldir@|$(htmldir)|g" | \
         sed "s|@includedir@|$(includedir)|g" | \
         sed "s|@program_name@|$(program_name)|g" | \
         sed "s|@libdir@|$(libdir)|g" > $@~
	cmp -s $@~ $@ || cp $@~ $@
	rm $@~

# file to import compilation options inside c++
Config.cpp: Config.cpp.in
	@echo "Updating Config.cpp"
	@sed "s/@SOEXT@/$(SOEXT)/g" $< | \
         sed "s/@ISINSTALLED@/false/g" | \
         sed "s|@htmldir@|xxxxNAxxxx|g" | \
         sed "s|@includedir@|xxxxNAxxxx|g" | \
         sed "s|@program_name@|xxxxNAxxxx|g" | \
         sed "s|@PLUMED_ROOT@|$(CURDIR:src/config=)|g" > $@~
	cmp -s $@~ $@ || cp $@~ $@
	rm $@~

links: 
	@echo Nothing to do in $(CURDIR)


Config.o: Config.cpp Makefile.conf.xxd version.h
ConfigInstall.o: ConfigInstall.cpp Makefile.conf.xxd version.h

Makefile.conf.xxd: ../../Makefile.conf
	cat ../../Makefile.conf | xxd -i > Makefile.conf.xxd

# file to import compilation options inside a bash script
compile_options.sh:
	@echo Building compile_options.sh, a sourceable bash script with some compilation options
	@echo "compile='"$(CXX) -c -I'"$$PLUMED_INCLUDEDIR/$$PLUMED_PROGRAM_NAME"'/ $(CPPFLAGS) $(CXXFLAGS) -o"'" > $@
	@echo "link='"$(LDSO) -o"'" >> $@
	@echo "soext='"$(SOEXT)"'" >> $@

.PHONY: clean
clean:
	rm -f Config.cpp compile_options.sh Makefile_conf version.h
	rm -fr deps 
	rm -f Config.cpp ConfigInstall.cpp compile_options.sh install.conf
	rm -f $(CLEANLIST)

# generic makefile rules
include ../maketools/make.rules
