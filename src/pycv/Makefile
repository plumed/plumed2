USE=core tools bias colvar function
ADDCPPFLAGS += -Ipybind11-2.5.0/include

# generic makefile
include ../maketools/make.module


# ========================================
# Stuff for building the dynamic object

# Workaround conda bugs on OSX
#  - stack size not allowed when making dylib
#  - need to pass export_dynamic to linker  https://bugs.python.org/issue21122
remove_flags := -Wl,-stack_size,1000000
add_flags := -fno-lto  -Wl,-export_dynamic -undefined dynamic_lookup 
PYTHON_LDFLAGS_SO := $(subst $(remove_flags),,$(shell python3-config --ldflags)) $(add_flags)

# You may or may not want to pass this to the linking command. If you
# do, you are giving an hint as to where to find the kernel at
# runtime. Having it hard-coded may not be desirable, especially if
# the kernel has PYCV already built-in (in which case you get an
# "interpreter is already running" error)

# LIBPLUMEDKERNEL_SO := $(abspath ../lib/libplumedKernel.$(SOEXT))

PYCV.$(SOEXT): PythonCV.o PythonFunction.o PythonPlumedBase.o
	@$(info  )
	@$(info Adding the following flags to the usual PLUMED mklib: $(PYTHON_LDFLAGS_SO))
	@$(info  )
	$(LDSHARED) $(PYTHON_LDFLAGS_SO) $(DYNAMIC_LIBS) $(LIBPLUMEDKERNEL_SO) $^ -o $@



# ========================================
# For debugging

# https://www.cmcrossroads.com/article/dumping-every-makefile-variable
.PHONY: printvars
printvars:
	@$(foreach V,$(sort $(.VARIABLES)), $(if $(filter-out environment% default automatic,  $(origin $V)),$(warning $V=$($V) ($(value $V)))))



# Note: to run regtests, use e.g.
# export PLUMED_PROGRAM_NAME=/Users/toni/compile/plumed2-pycv/src/lib/plumed
