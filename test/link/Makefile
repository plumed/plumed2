
-include ../../Makefile.conf
-include ../../src/Plumed.inc

EXE=plumed plumedf plumedc \
    plumed-static plumedf-static plumedc-static \
    plumed-dummy plumedf-dummy plumedc-dummy

all: $(EXE)

plumed-dummy: main.o $(PLUMED_WRAPPER)
	$(LD) $(LDFLAGS) $^ -o $@ $(PLUMED_LIBS)

plumed-static: main.o $(PLUMED_OBJ)
	$(LD) $(LDFLAGS) $^ -o $@ $(PLUMED_LIBS) $(PLUMED_DYNAMIC_LIBS)

plumed: main.o $(PLUMED_SHARED_OBJ)
	$(LD) $(LDFLAGS) $^ -o $@ $(PLUMED_LIBS)

plumedf-dummy: mainf.o $(PLUMED_WRAPPER)
	$(LDF90) $(LDFLAGS) $^ -o $@ $(PLUMED_LIBS)

plumedf-static: mainf.o $(PLUMED_OBJ)
	$(LDF90) $(LDFLAGS) $^ -o $@ $(PLUMED_LIBS) $(PLUMED_DYNAMIC_LIBS)

plumedf: mainf.o $(PLUMED_SHARED_OBJ)
	$(LDF90) $(LDFLAGS) $^ -o $@ $(PLUMED_LIBS)

plumedc-dummy: mainc.o $(PLUMED_WRAPPER)
	$(LD) $(LDFLAGS) $^ -o $@ $(PLUMED_LIBS)

plumedc-static: mainc.o $(PLUMED_OBJ)
	$(LD) $(LDFLAGS) $^ -o $@ $(PLUMED_LIBS) $(PLUMED_DYNAMIC_LIBS)

plumedc: mainc.o $(PLUMED_SHARED_OBJ)
	$(LD) $(LDFLAGS) $^ -o $@ $(PLUMED_LIBS)

clean:
	rm -f *.o
	rm -f $(EXE)
	rm -f makefile.dep

makefile.dep : *.c *.cpp *.h
	for i in *.c *.cpp; do g++  $(CPPFLAGS) -MM "$${i}"; done > $@

-include makefile.dep

%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@



