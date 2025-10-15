USE=bias colvar core tools function

#generic makefile
include ../maketools/make.module
CXXFLAGS:= $(CXXFLAGS) -Wshadow -Werror

ifneq (,$(findstring -D__PLUMED_HAS_LIBTORCH,$(CPPFLAGS)))
#torch gives some problems with some unnamed warning in variadic macros
#suppress all warnings
EMMIVox.o: CXXFLAGS+=-w 
endif
