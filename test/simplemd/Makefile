
-include ../../Makefile.conf
-include ../../src/Plumed.inc

SRC=routines.f90 simplemd.f90 ran1.f gasdev.f
OBJ=routines.o simplemd.o ran1.o gasdev.o

EXE=simplemd-dummy simplemd-static simplemd

all: $(EXE)

simplemd-dummy: $(OBJ) $(PLUMED_RUNTIME_DEPENDENCIES)
	$(LDF90) $(OBJ) -o $@ $(PLUMED_RUNTIME_LOAD)

simplemd-static: $(OBJ) $(PLUMED_STATIC_DEPENDENCIES)
	$(LDF90) $(OBJ) -o $@ $(PLUMED_STATIC_LOAD)

simplemd: $(OBJ) $(PLUMED_SHARED_DEPENDENCIES)
	$(LDF90) $(OBJ) -o $@ $(PLUMED_SHARED_LOAD)

clean:
	rm -f *.o
	rm -f $(EXE)
	rm -f *.mod
	rm -f energies.dat output.xyz trajectory.xyz

# this is due to module inclusion:
simplemd.o: routines.o

%.o : %.f90
	$(FC) -c $(FFLAGS) $(DOUBLEPRECISION) $< -o $@

%.o : %.f
	$(FC) -c $(FFLAGS) $(DOUBLEPRECISION) $< -o $@



