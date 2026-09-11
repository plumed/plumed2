# I want to include the configuration only if the target is not only astyle or clean or a combination of astyle and clean
# And I want an error if the configuration does not exists
#no targets
ifeq ($(MAKECMDGOALS),)
  include ./Makefile.conf
# otherwise include if the user asked (also) for something other than clean and astyle
else ifneq ($(MAKECMDGOALS),$(filter clean astyle,$(MAKECMDGOALS)))
  include ./Makefile.conf
endif

#Dependency tracking based on https://make.mad-scientist.net/papers/advanced-auto-dependency-generation/#tldr
#this assumes gcc
DEPDIR := .deps
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.d
NVCC = nvcc
#tested with nvcc with :"Build cuda_11.7.r11.7/compiler.31442593_0"
#-dc adds relocatable device code 
#-dlto Perform link-time optimization of device code. 
NVCCCFLAGS = -dc -dlto
NVCCLDFLAGS = -shared -dlto
OBJS = Coordination.o cudaHelpers.o

ADDCPPFLAGS=$(PLUMED_INCLUDE)
ADDCLDFLAGS=

ifeq ($(SOEXT),dylib)
  SONAME_OPTION:=-Xlinker -install_name
else
  SONAME_OPTION:=-Xlinker -soname
endif

#this makes make search for files in src
VPATH = src

.PHONY: clean check all astyle

all: CudaCoordination.$(SOEXT)


%.o: %.cu $(DEPDIR)/%.d | $(DEPDIR)
	@echo Compiling object $@
	@$(NVCC) $(NVCCCFLAGS) -ccbin $(CXX) -c $(DEPFLAGS) $(CPPFLAGS) $(ADDCPPFLAGS) $(CXXFLAGS) $< -o $@

$(DEPDIR): ; @mkdir -p $@

DEPFILES := $(OBJS:%.o=$(DEPDIR)/%.d)
$(DEPFILES):
include $(wildcard $(DEPFILES))

CudaCoordination.$(SOEXT): $(OBJS)
	@echo Linking $@
	@$(NVCC) $(NVCCLDFLAGS) -ccbin $(LDSHARED)  $(ADDCLDFLAGS) $(SONAME_OPTION),"$(notdir $@)" $(DYNAMIC_LIBS) $(PLUMED_KERNEL) -o $@ $^

#clean should work even if Makefile.conf is not here
clean:
	@rm -fv $(OBJS) CudaCoordination.so CudaCoordination.dylib
	@rm -rf .deps

check: all
	$(MAKE) -C regtest
	$(MAKE) -C regtest checkfail

astyle:
	@../astyle.sh

