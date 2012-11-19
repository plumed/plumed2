
-include ../../Makefile.conf
-include ../../src/lib/Plumed.inc

EXE=plumed plumedf plumedc \
    plumed-static plumedf-static plumedc-static \
    plumed-dummy plumedf-dummy plumedc-dummy

all: $(EXE)

plumed-dummy: main.o $(PLUMED_RUNTIME_DEPENDENCIES)
	$(LD) main.o -o $@ $(PLUMED_RUNTIME_LOAD)

plumed-static: main.o $(PLUMED_STATIC_DEPENDENCIES)
	$(LD) main.o -o $@ $(PLUMED_STATIC_LOAD)

plumed: main.o $(PLUMED_SHARED_DEPENDENCIES)
	$(LD) main.o -o $@ $(PLUMED_SHARED_LOAD)

plumedf-dummy: mainf.o $(PLUMED_RUNTIME_DEPENDENCIES)
	$(LDF90) mainf.o -o $@ $(PLUMED_RUNTIME_LOAD)

plumedf-static: mainf.o $(PLUMED_STATIC_DEPENDENCIES)
	$(LDF90) mainf.o -o $@ $(PLUMED_STATIC_LOAD)

plumedf: mainf.o $(PLUMED_SHARED_DEPENDENCIES)
	$(LDF90) mainf.o -o $@ $(PLUMED_SHARED_LOAD)

plumedc-dummy: mainc.o $(PLUMED_RUNTIME_DEPENDENCIES)
	$(LD) mainc.o -o $@ $(PLUMED_RUNTIME_LOAD)

plumedc-static: mainc.o $(PLUMED_STATIC_DEPENDENCIES)
	$(LD) mainc.o -o $@ $(PLUMED_STATIC_LOAD)

plumedc: mainc.o $(PLUMED_SHARED_DEPENDENCIES)
	$(LD) mainc.o -o $@ $(PLUMED_SHARED_LOAD)

clean:
	rm -f *.o
	rm -f $(EXE)
	rm -f makefile.dep

makefile.dep : *.c *.cpp *.h
	for i in *.c *.cpp; do g++  $(CPPFLAGS) -MM "$${i}"; done > $@

-include makefile.dep

%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@



