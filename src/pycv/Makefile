USE=core tools bias colvar
ADDCPPFLAGS += $(shell python3-config --cflags)
ADDCPPFLAGS += -Ipybind11-2.3.0/include

# generic makefile
include ../maketools/make.module


# Workaround conda bugs on OSX
#  - stack size not allowed when making dylib
#  - need to pass export_dynamic to linker  https://bugs.python.org/issue21122
remove_me := -Wl,-stack_size,1000000
PYTHON_LDFLAGS_SO := $(subst $(remove_me),,$(shell python3-config --ldflags)) -fno-lto -Wl,-export_dynamic -undefined dynamic_lookup

# Uninstalled plumed dylib location (use $(libdir) for the installed one)
#LIBPLUMEDKERNEL_SO := $(abspath ../lib/libplumedKernel.$(SOEXT))
LIBPLUMEDKERNEL_SO := $(HOME)/compile/plumed2/src/lib/libplumedKernel.$(SOEXT)

PythonCV.$(SOEXT): PythonCV.o
	@$(warning Warning: linking $@ to uninstalled kernel $(LIBPLUMEDKERNEL_SO))
	@$(info Flags being added to usual mklib: $(PYTHON_LDFLAGS_SO))
	$(LDSHARED) $(PYTHON_LDFLAGS_SO) $(DYNAMIC_LIBS) $(LIBPLUMEDKERNEL_SO) $^ -o $@


# https://www.cmcrossroads.com/article/dumping-every-makefile-variable
.PHONY: printvars
printvars:
	@$(foreach V,$(sort $(.VARIABLES)), $(if $(filter-out environment% default automatic,  $(origin $V)),$(warning $V=$($V) ($(value $V)))))


.PHONY: testpy
testpy:
	@$(info You may need to set e.g. PYTHONHOME=$$CONDA_PREFIX)
	../lib/plumed-shared driver --ixyz trajectory.xyz --plumed plumed.dat
