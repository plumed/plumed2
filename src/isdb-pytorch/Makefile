USE=colvar core tools

include ../maketools/make.module
CXXFLAGS:= $(CXXFLAGS) -Wshadow -Werror

ifneq (,$(findstring -D__PLUMED_HAS_LIBTORCH,$(CPPFLAGS)))
# LibTorch headers require suppression of variadic macro and other warnings.
EMMIVox.o: CXXFLAGS+=-w
endif
