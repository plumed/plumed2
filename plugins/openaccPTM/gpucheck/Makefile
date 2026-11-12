# Checks, on a node with GPUs, that the plugin leaves the GPU of the MD code alone
# (see device_check.cpp). Run with "make check-gpu" from the plugin directory, after
# "make" (ACCTARGET=gpu). Not part of the regtests: the CI has no GPU.

NVCXX ?= nvc++
PLUMED_INCLUDE_DIR = $(shell plumed --no-mpi info --include-dir)

all: check

device_check: device_check.cpp
	$(NVCXX) -cuda -I$(PLUMED_INCLUDE_DIR) -D__PLUMED_WRAPPER_IMPLEMENTATION=1 -D__PLUMED_HAS_DLOPEN -D__PLUMED_WRAPPER_LINK_RUNTIME=1 $< -o $@ -ldl

# the kernel is loaded at runtime (PLUMED_KERNEL), as the MD codes do; the plugin
# then needs the symbols of the kernel to be global
check: device_check
	PLUMED_LOAD_NAMESPACE=GLOBAL ./device_check ../plumedOpenACC.so follow
	PLUMED_LOAD_NAMESPACE=GLOBAL ./device_check ../plumedOpenACC.so md

clean:
	rm -f device_check device_check_*.log

.PHONY: all check clean
