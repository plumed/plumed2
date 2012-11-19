
# include the machine dependent configuration
-include ../../Makefile.conf

# default target
.PHONY: all
all: PlumedConfig.h Makefile.conf.xxd compile_options.sh

# file to import compilation options inside c++
PlumedConfig.h: PlumedConfig.h.in ../../Makefile.conf
	@echo Building PlumedConfig.h, an includable c++ header with some compilation options
	@sed "s/@SOEXT@/$(SOEXT)/g" $< | \
         sed "s|@PLUMED_ROOT@|$(CURDIR:src/config=)|g" > $@

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
	rm -f PlumedConfig.h compile_options.sh
	rm -f $(CLEANLIST)

# generic makefile rules
include ../maketools/make.rules
