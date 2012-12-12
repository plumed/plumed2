
# include the machine dependent configuration
-include ../../Makefile.conf

# default target
.PHONY: all
all: Config.o compile_options.sh

# file to import compilation options inside c++
Config.cpp: Config.cpp.in Config.h
	@echo Building Config.cpp, a object containing configuration infos
	@sed "s/@SOEXT@/$(SOEXT)/g" $< | \
         sed "s|@PLUMED_ROOT@|$(CURDIR:src/config=)|g" > $@

Config.o: Config.cpp Makefile.conf.xxd

Makefile.conf.xxd: ../../Makefile.conf
	cat ../../Makefile.conf | xxd -i > Makefile.conf.xxd

# file to import compilation options inside a bash script
compile_options.sh:
	@echo Building compile_options.sh, a sourceable bash script with some compilation options
	@echo "compile='"$(CXX) -c -I$(realpath .) $(CPPFLAGS) $(CXXFLAGS) -o"'" > $@
	@echo "link='"$(LDSO) -o"'" >> $@
	@echo "soext='"$(SOEXT)"'" >> $@

.PHONY: clean
clean:
	rm -f Config.cpp compile_options.sh Makefile_conf
	rm -fr deps
	rm -f $(CLEANLIST)

# generic makefile rules
include ../maketools/make.rules
