include ./Makefile.conf

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

clean:
	@rm -fv $(OBJS) CudaCoordination.$(SOEXT)
	@rm -rf .deps

check: all
	$(MAKE) -C regtest
	$(MAKE) -C regtest checkfail
