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
DEPDIR := .deps
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$(@:.o=.d)

#this makes make search for files in src
VPATH = src
OBJS = $(patsubst src/%.cpp,%.o,$(wildcard src/*.cpp))
# ACCTARGET selects what the OpenACC code is compiled for: "gpu" (default), "multicore" or "host".
# "host" runs the kernels serially on the CPU and needs no GPU (used by the CI to run the regtests).
ACCTARGET ?= gpu
ADDCPPFLAGS=$(PLUMED_INCLUDE) -D__PLUMED_HAS_OPENACC=1 -acc=$(ACCTARGET)
# -static-nvidia keeps the OpenACC runtime private to the plugin (see plugins/openaccPTM/Makefile)
ADDCLDFLAGS=-acc=$(ACCTARGET) -static-nvidia

# The soname must be passed to the linker as a single token: nvc++ splits
# "-Xlinker -soname,name" at the comma (see plugins/openaccPTM/Makefile).
ifeq ($(SOEXT),dylib)
  SONAME_FLAGS=-Xlinker -install_name -Xlinker $(notdir $@)
else
  SONAME_FLAGS=-Xlinker -soname=$(notdir $@)
# The NVIDIA runtime is linked statically: with nodelete the dlclose done by PLUMED
# does not finalize it before the MD code frees its GPU memory (see plugins/openaccPTM/Makefile).
  ADDCLDFLAGS+=-Xlinker -z -Xlinker nodelete
endif

PLUGIN=plumedOpenACCSAXS.$(SOEXT)

NVCXX ?= nvc++
NVLDSHARED = $(LDSHARED:$(CXX)=$(NVCXX))
# configure sets NVCXX=mpic++ when PLUMED has MPI: the wrapper must then call nvc++
export OMPI_CXX ?= nvc++
export MPICH_CXX ?= nvc++
export I_MPI_CXX ?= nvc++
all: $(PLUGIN)

%.o: %.cpp $(DEPDIR)/%.d | $(DEPDIR)
	@echo Compiling object $@
	@$(NVCXX) -c $(DEPFLAGS) $(CPPFLAGS) $(ADDCPPFLAGS) $(CXXFLAGS) $< -o $@

$(DEPDIR): ; @mkdir -p $@

DEPFILES := $(OBJS:%.o=$(DEPDIR)/%.d)
$(DEPFILES):
include $(wildcard $(DEPFILES))

$(PLUGIN): $(OBJS)
	@echo Linking $@
	@$(NVLDSHARED) $(ADDCLDFLAGS) $(SONAME_FLAGS) $(DYNAMIC_LIBS) $(PLUMED_KERNEL) -o $@ $^

clean:
	@rm -fv $(OBJS) $(PLUGIN)
	@rm -rf .deps

check: all
	$(MAKE) -C regtest
	$(MAKE) -C regtest checkfail

astyle:
	@../astyle.sh
