
# include the machine dependent configuration
ifneq ($(MAKECMDGOALS),clean)
  -include ../../Makefile.conf
endif

# openmp flags should not be used here
ifdef CXXFLAGS_NOOPENMP
  CXXFLAGS := $(CXXFLAGS_NOOPENMP)
endif

# if machine dependent configuration has been found:
ifdef GCCDEP

# source for plumed wrapper
SRC_WRAPPER=PlumedStatic.cpp

# source for plumed dynamic wrapper
SRC_DYNAMIC_WRAPPER=Plumed.c

# all sources
ALL_SRC=$(SRC_WRAPPER) $(SRC_DYNAMIC_WRAPPER)

# corresponding objects:
OBJ_WRAPPER=PlumedStatic.o
OBJ_DYNAMIC_WRAPPER=Plumed.o

# dependency files for all c++ sources
CPP_SRC=$(wildcard *.cpp)
ALL_DEP=$(addprefix deps/,$(CPP_SRC:.cpp=.d))

# default target
.PHONY: all obj
all: $(OBJ_WRAPPER) $(OBJ_DYNAMIC_WRAPPER)

obj: all


# there is a single C file here, and this is its dependency:
Plumed.o: Plumed.c Plumed.h
PlumedStatic.o: PlumedStatic.cpp Plumed.h

# Individual dependencies

# if machine dependent configuration has not been found:
else

.PHONY: error
error:
	@echo No configuration available
	@echo First run ./configure in the root directory
endif

# this target is available anyway

.PHONY: clean links
clean:
	rm -fr deps 
	rm -fr $(CLEANLIST)
links:
	@echo Nothing to do in $(CURDIR)

# generic makefile rules
include ../maketools/make.rules

