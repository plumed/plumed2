NVCC ?= nvcc
CUDA_ARCH ?= native
NVCCFLAGS ?= -O3 -std=c++17
TARGET = libReactiveVoronoiCuda.so

# Keep the CUDA runtime alive across PLUMED destruction inside an MD process.
ifeq ($(shell uname -s),Linux)
  PLATFORM_LDFLAGS = -Xlinker -z -Xlinker nodelete
endif

.PHONY: all clean check
all: $(TARGET)

$(TARGET): ReactiveVoronoiCuda.cu
	$(NVCC) $(NVCCFLAGS) -arch=$(CUDA_ARCH) --cudart=shared -Xcompiler=-fPIC -shared $(PLATFORM_LDFLAGS) $(LDFLAGS) $< -o $@

check: all
	$(MAKE) -C regtest check LIBRARY=$(abspath $(TARGET))

clean:
	rm -f $(TARGET)
