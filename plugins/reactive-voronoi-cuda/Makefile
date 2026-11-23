PLUMED ?= plumed
NVCC ?= nvcc
CUDA_ARCH ?= native
NVCCFLAGS ?= -O3 -std=c++17

.PHONY: all clean check check-nogpu astyle
all:

ifneq ($(filter-out clean astyle,$(or $(MAKECMDGOALS),all)),)
  SOEXT := $(shell $(PLUMED) --no-mpi info --soext)
  ifeq ($(SOEXT),)
    $(error Set PLUMED to an installed plumed executable)
  endif
  TARGET = libReactiveVoronoiCuda.$(SOEXT)
  PLUMED_INCLUDE ?= $(shell $(PLUMED) --no-mpi info --include-dir)
  PLUMED_KERNEL ?= $(shell $(PLUMED) --no-mpi info --lib-dir)/libplumedKernel.$(SOEXT)

  # Keep the runtime alive across PLUMED destruction inside an MD process.
  ifeq ($(shell uname -s),Linux)
    PLATFORM_LDFLAGS = -Xlinker -z -Xlinker nodelete -Xlinker -rpath -Xlinker "$(dir $(PLUMED_KERNEL))"
  endif

all: $(TARGET)

$(TARGET): src/ReactiveVoronoiCuda.cu $(PLUMED_KERNEL)
	$(NVCC) $(NVCCFLAGS) -arch=$(CUDA_ARCH) --cudart=shared -Xcompiler=-fPIC -shared -I"$(PLUMED_INCLUDE)" $(PLATFORM_LDFLAGS) $(LDFLAGS) $< "$(PLUMED_KERNEL)" -o $@

check: all
	$(MAKE) -C regtest check PLUMED="$(PLUMED)" LIBRARY="$(abspath $(TARGET))"
check-nogpu: all
	CUDA_VISIBLE_DEVICES="" $(MAKE) -C regtest check-nogpu PLUMED="$(PLUMED)" LIBRARY="$(abspath $(TARGET))"
endif

astyle:
	@../astyle.sh

clean:
	rm -f libReactiveVoronoiCuda.so libReactiveVoronoiCuda.dylib
