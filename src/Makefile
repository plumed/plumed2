
-include ../Makefile.conf
-include ../imd/Makefile.conf

ifndef srcdir
srcdir=.
endif

SRC_KERNEL=$(subst $(srcdir)/PlumedStatic.cpp,,$(wildcard $(srcdir)/*.cpp))
OBJ_KERNEL=$(subst $(srcdir)/,,$(SRC_KERNEL:.cpp=.o))

SRC_WRAPPER=$(srcdir)/PlumedStatic.cpp
OBJ_WRAPPER=PlumedStatic.o

SRC_DYNAMIC_WRAPPER=$(srcdir)/Plumed.c
OBJ_DYNAMIC_WRAPPER=Plumed.o

ALL_SRC=$(SRC_KERNEL) $(SRC_WRAPPER) $(SRC_DYNAMIC_WRAPPER)

.PHONY: all clean

ifdef GCCDEP
all: libplumedKernel.$(SOEXT) libplumed.$(SOEXT) Plumed.o Plumed.inc Plumed.inc.runtime Plumed.inc.shared Plumed.inc.static
else
all:
	@echo No configuration available
	@echo First run ./configure.sh in the root directory
endif

PlumedConfig.h: PlumedConfig.h.in
	sed "s/@SOEXT@/$(SOEXT)/g" $< > $@

check:
	echo $(OBJ_KERNEL) $(SRC_KERNEL)

libplumedKernel.$(SOEXT): $(OBJ_KERNEL)
	$(LDSO) -o $@ $^ $(DYNAMIC_LIBS)

libplumed.$(SOEXT): $(OBJ_KERNEL) $(OBJ_WRAPPER)
	$(LDSO) -o $@ $^ $(DYNAMIC_LIBS)

clean:
	rm -f *.o
	rm -f makefile.dep
	rm -f PlumedConfig.h
	rm -f Plumed.inc*
	rm -f *~
ifdef SOEXT
	rm -f *.$(SOEXT)
endif

Plumed.inc: $(OBJ_KERNEL) $(OBJ_WRAPPER) libplumed.$(SOEXT) Plumed.o
	echo "PLUMED_OBJ=" $(realpath $(OBJ_KERNEL)) $(realpath $(OBJ_WRAPPER))  > $@
	echo "PLUMED_SHARED_OBJ=" $(realpath libplumed.$(SOEXT)) >> $@
	echo "PLUMED_WRAPPER=" $(realpath Plumed.o) >> $@
	echo "PLUMED_INCLUDE="-I$(realpath ../include) >> $@
	echo "PLUMED_LIBS="$(LIBS) >> $@
	echo "PLUMED_DYNAMIC_LIBS="$(DYNAMIC_LIBS) >> $@
	echo "PLUMED_LDFLAGS="$(LDFLAGS) >> $@

Plumed.inc.runtime: Plumed.inc
	echo "# PLUMED: runtime installation" > $@
	cat $< >> $@
	echo PLUMED_LOAD= '$$(PLUMED_WRAPPER) $$(PLUMED_LIBS) $$(PLUMED_LDFLAGS)' >> $@
	echo PLUMED_DEPENDENCIES= >> $@

Plumed.inc.static: Plumed.inc
	echo "# PLUMED: static installation" > $@
	cat $< >> $@
	echo PLUMED_LOAD= '$$(PLUMED_OBJ) $$(PLUMED_LIBS) $$(PLUMED_DYNAMIC_LIBS) $$(PLUMED_LDFLAGS)' >> $@
	echo PLUMED_DEPENDENCIES='$$(PLUMED_OBJ)' >> $@

Plumed.inc.shared: Plumed.inc
	echo "# PLUMED: shared installation" > $@
	cat $< >> $@
	echo PLUMED_LOAD= '$$(PLUMED_SHARED_OBJ) $$(PLUMED_LIBS) $$(PLUMED_LDFLAGS)' >> $@
	echo PLUMED_DEPENDENCIES='$$(PLUMED_SHARED_OBJ)' >> $@

makefile.dep : $(ALL_SRC) $(srcdir)/*.h $(srcdir)/PlumedConfig.h
ifdef GCCDEP
	for i in $(ALL_SRC) ; do $(GCCDEP) $(CPPFLAGS) -MM "$${i}"; done > $@
endif

-include makefile.dep

