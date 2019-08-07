CXXFLAGS := $(shell python3-config --cflags)
LDFLAGS := $(shell python3-config --ldflags)
PYBINDINC := -I../pybind11/include

default: mybind py


mybind.o: mybind.cpp
	$(CXX) -std=c++14 $(CXXFLAGS) $(PYBINDINC) -c $< -o $@

mybind: mybind.o
	$(CXX) $(LDFLAGS) $^ -o $@


py: py.cpp

test:
	PYTHONHOME=$(CONDA_PREFIX) ./py simple mul 3 4

