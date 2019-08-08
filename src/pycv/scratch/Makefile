PYBINDINC := -I../pybind11-2.3.0/include
CXXFLAGS := $(shell python3-config --cflags) $(PYBINDINC)  -std=c++14 -fPIC -fno-lto
LDFLAGS := $(shell python3-config --ldflags) -fno-lto # -Wl,-export_dynamic

# The export dynamic thing was required in OSX to avoid obscure
# "symbol not found" runtime errors, which may be related to
# LTO. https://bugs.python.org/issue21122

# LTO also seems to slow down linking to the extreme, and sometimes
# "internal errors" the compiler, but can't find a way to disable it
# in OSX.


default: mybind py


mybind.o: mybind.cpp
	$(CXX) -std=c++14 $(CXXFLAGS) $(PYBINDINC) -c $< -o $@

mybind: mybind.o
	$(CXX)  $^ -o $@ $(LDFLAGS)

clean:
	rm mybind.o mybind


test:
	PYTHONHOME=$(CONDA_PREFIX) ./py simple mul 3 4

