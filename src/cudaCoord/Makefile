USE=colvar core tools

# generic makefile
include ../maketools/make.module
OBJ=$(SRC:.cu=.o)

Coordination.o:Coordination.cu
	nvcc $(CXXFLAGS) $< -o $@