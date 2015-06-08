
# include the machine dependent configuration
-include ../../Makefile.conf

# sed script to edit program name defaults to do nothing
ifndef program_transform_name
  program_transform_name=s,x,x,
endif

# default target
.PHONY: all obj
all:
	@echo "Updating ConfigInstall.cpp"
	@sed "s/@SOEXT@/$(SOEXT)/g" Config.cpp.in | \
         sed "s/@ISINSTALLED@/true/g" | \
	env prefix="${prefix}" env program_transform_name="${program_transform_name}" ../maketools/update-plumedroot.sh ConfigInstall.cpp
	@echo "Updating version file"
	@ ../maketools/update-version.sh version.h
	$(MAKE) do-all

obj: all

do-all: Config.o ConfigInstall.o compile_options.sh 
	make ConfigInstall.o

ConfigInstall.cpp:
	@echo "Updating ConfigInstall.cpp"
	@sed "s/@SOEXT@/$(SOEXT)/g" Config.cpp.in | \
         sed "s/@ISINSTALLED@/true/g" | \
         ../maketools/update-plumedroot.sh ConfigInstall.cpp
	@echo "Updating version file"
	@ ../maketools/update-version.sh version.h

# file to import compilation options inside c++
Config.cpp: Config.cpp.in Config.h
	@echo Building Config.cpp, a object containing configuration infos
	@sed "s/@SOEXT@/$(SOEXT)/g" $< | \
         sed "s/@ISINSTALLED@/false/g" | \
         sed "s|@PLUMED_ROOT@|$(CURDIR:src/config=)|g" > $@

links: 
	@echo Nothing to do in $(CURDIR)


Config.o: Config.cpp Makefile.conf.xxd version.h
ConfigInstall.o: ConfigInstall.cpp Makefile.conf.xxd version.h

Makefile.conf.xxd: ../../Makefile.conf
	cat ../../Makefile.conf | xxd -i > Makefile.conf.xxd

# file to import compilation options inside a bash script
compile_options.sh:
	@echo Building compile_options.sh, a sourceable bash script with some compilation options
	@echo "compile='"$(CXX) -c -I'$$PLUMED_ROOT'/src $(CPPFLAGS) $(CXXFLAGS) -o"'" > $@
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
