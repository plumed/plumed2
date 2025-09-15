include ./Makefile.conf

#Dependency tracking based on https://make.mad-scientist.net/papers/advanced-auto-dependency-generation/#tldr
#this assumes gcc
DEPDIR := .deps
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.d
#tested with nvcc with :"Build cuda_11.7.r11.7/compiler.31442593_0"
#-dc adds relocatable device code 
#-dlto Perform link-time optimization of device code. 
CVS_OBJS = $(wildcard *.cpp)
TOOLS_OBJS = extra/colvar/ColvarInput.o extra/tools/Pbc.o extra/tools/SwitchingFunction.o extra/tools/HistogramBead.o

OBJS = $(CVS_OBJS:.cpp=.o) $(TOOLS_OBJS)
ADDCPPFLAGS=$(PLUMED_INCLUDE) -D__PLUMED_HAS_OPENACC=1 -acc -acc="gpu"
ADDCLDFLAGS=-acc -static-nvidia

ifeq ($(SOEXT),dylib)
  SONAME_OPTION:=-Xlinker -install_name
else
  SONAME_OPTION:=-Xlinker -soname
endif

PLUGIN=plumedOpenACC.$(SOEXT)

NVCXX ?= nvc++
NVLDSHARED = $(LDSHARED:$(CXX)=$(NVCXX))
all: $(PLUGIN)

%.o: %.cpp $(DEPDIR)/%.d | $(DEPDIR)
	@echo Compiling object $@
	$(NVCXX) -c $(DEPFLAGS) $(CPPFLAGS) $(ADDCPPFLAGS) $(CXXFLAGS) $< -o $@

extra/colvar/%.o: plumed/colvar/%.cpp $(DEPDIR)/extra/colvar/%.d | $(DEPDIR)/extra/colvar extra/colvar
	@echo Compiling object $@
	$(NVCXX) -c $(DEPFLAGS) $(CPPFLAGS) $(ADDCPPFLAGS) $(CXXFLAGS) $< -o $@

extra/tools/%.o: plumed/tools/%.cpp $(DEPDIR)/extra/tools/%.d | $(DEPDIR)/extra/tools extra/tools
	@echo Compiling object $@
	$(NVCXX) -c $(DEPFLAGS) $(CPPFLAGS) $(ADDCPPFLAGS) $(CXXFLAGS) $< -o $@


$(DEPDIR): ; @mkdir -p $@
$(DEPDIR)/extra/colvar: ; @mkdir -p $@
$(DEPDIR)/extra/tools: ; @mkdir -p $@
extra/colvar: ; mkdir -p $@
extra/tools: ; mkdir -p $@

DEPFILES := $(OBJS:%.o=$(DEPDIR)/%.d)
$(DEPFILES):
include $(wildcard $(DEPFILES))

$(PLUGIN): $(OBJS)
	@echo Linking $@
	$(NVLDSHARED)  $(ADDCLDFLAGS) $(SONAME_OPTION),"$(notdir $@)" $(DYNAMIC_LIBS) $(PLUMED_KERNEL) -o $@ $^

clean:
	@rm -fv $(OBJS) $(PLUGIN)
	@rm -rf .deps

check: all
	$(MAKE) -C regtest
	$(MAKE) -C regtest checkfail
