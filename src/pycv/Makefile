USE=core tools bias colvar function
ADDCPPFLAGS += $(shell python3-config --cflags)
ADDCPPFLAGS += -Ipybind11-2.3.0/include

# generic makefile
include ../maketools/make.module

# Note to myself: to run regtests,
# export PLUMED_PROGRAM_NAME=/Users/toni/compile/plumed2-pycv/src/lib/plumed


# Workaround conda bugs on OSX
#  - stack size not allowed when making dylib
#  - need to pass export_dynamic to linker  https://bugs.python.org/issue21122
remove_flags := -Wl,-stack_size,1000000
PYTHON_LDFLAGS_SO := $(subst $(remove_flags),,$(shell python3-config --ldflags)) -fno-lto # -Wl,-export_dynamic -undefined dynamic_lookup


# Uninstalled plumed dylib location (use $(libdir) for the installed one)
LIBPLUMEDKERNEL_SO := $(abspath ../lib/libplumedKernel.$(SOEXT))

PYCV.$(SOEXT): PythonCV.o PythonFunction.o PythonPlumedBase.o
	@$(warning Warning: linking $@ to uninstalled kernel $(LIBPLUMEDKERNEL_SO))
	@$(info Flags being added to usual mklib: $(PYTHON_LDFLAGS_SO))
	$(LDSHARED) $(PYTHON_LDFLAGS_SO) $(DYNAMIC_LIBS) $(LIBPLUMEDKERNEL_SO) $^ -o $@


# https://www.cmcrossroads.com/article/dumping-every-makefile-variable
.PHONY: printvars
printvars:
	@$(foreach V,$(sort $(.VARIABLES)), $(if $(filter-out environment% default automatic,  $(origin $V)),$(warning $V=$($V) ($(value $V)))))


