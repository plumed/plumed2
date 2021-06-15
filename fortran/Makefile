-include ../Makefile.conf
-include ../src/lib/Plumed.inc

SRC=plumed.f03 main.f03
OBJ=plumed.o main.o

EXE=test-runtime test

all: $(EXE)

gen: 
	python3 gen.py > plumed.f03

test-runtime: $(OBJ) $(PLUMED_RUNTIME_DEPENDENCIES)
	$(FC) $(OBJ) -o $@ $(PLUMED_RUNTIME_LOAD)

test: $(OBJ) $(PLUMED_SHARED_DEPENDENCIES)
	$(FC) $(OBJ) -o $@ $(PLUMED_SHARED_LOAD)

clean:
	rm -f *.o
	rm -f $(EXE)
	rm -f *.mod
	rm -f plumed.f03

# this is due to module inclusion:
main.o: plumed.o

%.o : %.f03
	$(FC) -c $(FFLAGS) $(DOUBLEPRECISION) $< -o $@
