include Make.inc

OBJS = PythonCVInterface.o PythonPlumedBase.o PlumedPythonEmbeddedModule.o

ADDCPPFLAGS=$(shell python3-config --cflags --embed) $(shell python -m pybind11 --includes) $(PLUMED_INCLUDE)
ADDCLDFLAGS=$(shell python3-config --ldflags --embed) 
ifeq ($(SOEXT),dylib)
  SONAME_OPTION:=-Wl,-install_name
else
  SONAME_OPTION:=-Wl,-soname
endif

all: PythonCVInterface.$(SOEXT)

%.o: %.cpp
	$(CXX) -c $(CPPFLAGS) $(ADDCPPFLAGS) $(CXXFLAGS) $< -o $@

PythonCVInterface.$(SOEXT): $(OBJS)
	$(LDSHARED) $(ADDCLDFLAGS) $(SONAME_OPTION),"$(notdir $@)" $(DYNAMIC_LIBS) $(PLUMED_KERNEL) -o $@ $^

clean:
	rm -f $(OBJS) PythonCVInterface.$(SOEXT)