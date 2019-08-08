USE=core tools bias colvar
ADDCPPFLAGS += $(shell python3-config --cflags)
ADDCPPFLAGS += -Ipybind11-2.3.0/include

# generic makefile
include ../maketools/make.module


# Workaround conda bug on OSX - stack size not allowed when making dylib
remove_me := -Wl,-stack_size,1000000
PYTHON_LDFLAGS_SO := $(subst $(remove_me),,$(shell python3-config --ldflags))

# Uninstalled plumed dylib location (use $(libdir) for the installed one)
LIBPLUMEDKERNEL_SO := $(realpath ../lib/libplumedKernel.$(SOEXT))

PythonCV.$(SOEXT): PythonCV.o
	echo $(warning Warning: linking to uninstalled kernel $(LIBPLUMEDKERNEL_SO))
	$(LDSHARED) $(PYTHON_LDFLAGS_SO) $(DYNAMIC_LIBS) $(LIBPLUMEDKERNEL_SO) $^ -o $@


# https://www.cmcrossroads.com/article/dumping-every-makefile-variable
.PHONY: printvars
printvars:
	@$(foreach V,$(sort $(.VARIABLES)), $(if $(filter-out environment% default automatic,  $(origin $V)),$(warning $V=$($V) ($(value $V)))))
