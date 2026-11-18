# Checks, on a node with GPUs, that the GPU plugins leave the GPU of the MD code alone
# (see device_check.cpp). Run with "make check-gpu" from the openaccPTM directory (after
# "make", ACCTARGET=gpu) or from the cudaCoord one. Not part of the regtests: the CI has no GPU.
# PLUGINS lists the plugins to check (by default the openaccPTM one).

CUDACXX ?= nvcc
PLUGINS ?= $(CURDIR)/../plumedOpenACC.so
PLUMED_INCLUDE_DIR = $(shell plumed --no-mpi info --include-dir)

all: check

device_check: device_check.cpp
	$(CUDACXX) -std=c++17 -I$(PLUMED_INCLUDE_DIR) -D__PLUMED_WRAPPER_IMPLEMENTATION=1 -D__PLUMED_HAS_DLOPEN -D__PLUMED_WRAPPER_LINK_RUNTIME=1 $< -o $@ -ldl

# The kernel is loaded at runtime (PLUMED_KERNEL), as the MD codes do; the plugin
# then needs the symbols of the kernel to be global.
# Called by "make check-gpu" of a plugin, PLUMED_KERNEL is not the one of the environment
# but the link flags of the Makefile.conf of the plugin ("-L/path/libplumedKernel.so" or
# "/path/libplumedKernel.so -Xlinker ..."): the library is the kernel the plugin was built
# against, that is the one to use.
KERNEL_LIB = $(patsubst -L%,%,$(firstword $(PLUMED_KERNEL)))
check: device_check
	PLUMED_KERNEL="$(KERNEL_LIB)" PLUMED_LOAD_NAMESPACE=GLOBAL ./device_check follow $(PLUGINS)
	PLUMED_KERNEL="$(KERNEL_LIB)" PLUMED_LOAD_NAMESPACE=GLOBAL ./device_check md $(PLUGINS)

clean:
	rm -f device_check device_check_*.log

.PHONY: all check clean
