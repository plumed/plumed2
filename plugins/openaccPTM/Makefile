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
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$(@:.o=.d)
#tested with nvcc with :"Build cuda_11.7.r11.7/compiler.31442593_0"
#-dc adds relocatable device code 
#-dlto Perform link-time optimization of device code. 
CVS_OBJS = $(patsubst src/%,%, $(wildcard src/*.cpp))
TOOLS_OBJS = extra/colvar/ColvarInput.o extra/tools/Pbc.o extra/tools/SwitchingFunction.o extra/tools/HistogramBead.o extra/gridtools/SumOfKernels.o

#this makes make search for files in src
VPATH = src
OBJS = $(CVS_OBJS:.cpp=.o) $(TOOLS_OBJS) $(EXTRA_OBJS)
# ACCTARGET selects what the OpenACC code is compiled for: "gpu" (default), "multicore" or "host".
# "host" runs the kernels serially on the CPU and needs no GPU nor OpenMP runtime (used by the CI to run the regtests).
ACCTARGET ?= gpu
# __PLUMED_HAS_OPENACC is not used by the code (which uses the standard _OPENACC): it declares that
# the plugin provides the openacc feature, for the module tools that report the __PLUMED_HAS_ macros
# of a compiled module (e.g. "plumed modules --has-feature openacc")
ADDCPPFLAGS=$(PLUMED_INCLUDE) -I$(CURDIR)/plumed -D__PLUMED_HAS_OPENACC=1 -acc=$(ACCTARGET)
# -static-nvidia is needed for every target: it keeps the OpenACC runtime private to the plugin. Otherwise
# a PLUMED built with GCC and OpenMP resolves the acc_* calls of the plugin to libgomp ("libgomp: TODO").
ADDCLDFLAGS=-acc=$(ACCTARGET) -static-nvidia
# code that is only needed with GPUs (e.g. keeping the CUDA context of the MD code, see OpenACCDebye.cpp)
ifneq (,$(findstring gpu,$(ACCTARGET)))
  ADDCPPFLAGS+=-DPLUMED_OPENACC_GPU
endif

# The soname must be passed to the linker as a single token: nvc++ splits
# "-Xlinker -soname,name" at the comma and then treats "name" as an input
# file, which makes ld fail with "input file 'name' is the same as output file".
ifeq ($(SOEXT),dylib)
  SONAME_FLAGS=-Xlinker -install_name -Xlinker $(notdir $@)
else
  SONAME_FLAGS=-Xlinker -soname=$(notdir $@)
# The NVIDIA runtime is linked statically in the plugin (-static-nvidia), so it
# is finalized when the plugin is unloaded. PLUMED unloads the libraries
# of LOAD with dlclose when it is destroyed, and an MD code destroys PLUMED
# before freeing its own GPU memory (GROMACS >= 2025 does it just before
# destroying the forcerec): the runtime shutdown then invalidates the MD code's
# allocations, which fails with "Freeing of the device buffer failed ...
# cudaErrorInvalidValue". With nodelete dlclose does not unload the plugin, and
# the runtime is finalized at process exit, as for a statically linked code.
  ADDCLDFLAGS+=-Xlinker -z -Xlinker nodelete
endif

PLUGIN=plumedOpenACC.$(SOEXT)

NVCXX ?= nvc++
NVLDSHARED = $(LDSHARED:$(CXX)=$(NVCXX))
# configure sets NVCXX=mpic++ when PLUMED has MPI: the wrapper must then call
# nvc++ and not the compiler that was used to build MPI. These variables are
# read by the OpenMPI, MPICH and Intel MPI wrappers respectively (and ignored
# by nvc++ itself); they are exported to the compile and link commands, and a
# value already set in the environment takes precedence.
export OMPI_CXX ?= nvc++
export MPICH_CXX ?= nvc++
export I_MPI_CXX ?= nvc++
all: $(PLUGIN)

%.o: %.cpp $(DEPDIR)/%.d | $(DEPDIR)
	@echo Compiling object $@
	@$(NVCXX) -c $(DEPFLAGS) $(CPPFLAGS) $(ADDCPPFLAGS) $(CXXFLAGS) $< -o $@

extra/colvar/%.o: plumed/colvar/%.cpp $(DEPDIR)/extra/colvar/%.d | $(DEPDIR)/extra/colvar extra/colvar
	@echo Compiling object $@
	@$(NVCXX) -c $(DEPFLAGS) $(CPPFLAGS) $(ADDCPPFLAGS) $(CXXFLAGS) $< -o $@

extra/tools/%.o: plumed/tools/%.cpp $(DEPDIR)/extra/tools/%.d | $(DEPDIR)/extra/tools extra/tools
	@echo Compiling object $@
	@$(NVCXX) -c $(DEPFLAGS) $(CPPFLAGS) $(ADDCPPFLAGS) $(CXXFLAGS) $< -o $@

extra/gridtools/%.o: plumed/gridtools/%.cpp $(DEPDIR)/extra/gridtools/%.d | $(DEPDIR)/extra/gridtools extra/gridtools
	@echo Compiling object $@
	@$(NVCXX) -c $(DEPFLAGS) $(CPPFLAGS) $(ADDCPPFLAGS) $(CXXFLAGS) $< -o $@

$(DEPDIR): ; @mkdir -p $@
$(DEPDIR)/extra/colvar: ; @mkdir -p $@
$(DEPDIR)/extra/tools: ; @mkdir -p $@
$(DEPDIR)/extra/gridtools: ; @mkdir -p $@
extra/colvar: ; mkdir -p $@
extra/tools: ; mkdir -p $@
extra/gridtools: ; mkdir -p $@

DEPFILES := $(OBJS:%.o=$(DEPDIR)/%.d)
$(DEPFILES):
include $(wildcard $(DEPFILES))

$(PLUGIN): $(OBJS)
	@echo Linking $@
	@$(NVLDSHARED)  $(ADDCLDFLAGS) $(SONAME_FLAGS) -o $@ $^ $(DYNAMIC_LIBS) $(PLUMED_KERNEL)

clean:
	@rm -fv $(OBJS) $(PLUGIN)
	@rm -rf .deps
	@$(MAKE) -C gpucheck clean

check: all
	$(MAKE) -C regtest
	$(MAKE) -C regtest checkfail

# on a node with GPUs: checks that the plugin leaves the GPU of the MD code alone (see gpucheck/)
check-gpu: all
	$(MAKE) -C gpucheck

astyle:
	@../astyle.sh

