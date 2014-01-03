
# include the machine dependent configuration
-include ../../Makefile.conf

# default target
.PHONY: all obj
all:
	@echo "Updating ConfigInstall.cpp"
	@sed "s/@SOEXT@/$(SOEXT)/g" Config.cpp.in | \
         sed "s/@ISINSTALLED@/true/g" | \
	env prefix="${prefix}" ../maketools/update-plumedroot.sh ConfigInstall.cpp
	$(MAKE) do-all

obj: all

do-all: Config.o ConfigInstall.o compile_options.sh 
	make ConfigInstall.o

ConfigInstall.cpp:
	@echo "Updating ConfigInstall.cpp"
	@sed "s/@SOEXT@/$(SOEXT)/g" Config.cpp.in | \
         sed "s/@ISINSTALLED@/true/g" | \
         ../maketools/update-plumedroot.sh ConfigInstall.cpp

# file to import compilation options inside c++
Config.cpp: Config.cpp.in Config.h
	@echo Building Config.cpp, a object containing configuration infos
	@sed "s/@SOEXT@/$(SOEXT)/g" $< | \
         sed "s/@ISINSTALLED@/false/g" | \
         sed "s|@PLUMED_ROOT@|$(CURDIR:src/config=)|g" > $@

links: 
	@echo Nothing to do in $(CURDIR)
	

Config.o: Config.cpp Makefile.conf.xxd
ConfigInstall.o: ConfigInstall.cpp Makefile.conf.xxd

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
	rm -f Config.cpp compile_options.sh Makefile_conf
	rm -fr deps 
	rm -f Config.cpp ConfigInstall.cpp compile_options.sh install.conf
	rm -f $(CLEANLIST)

# generic makefile rules
include ../maketools/make.rules
