CXXFLAGS := $(shell python3-config --cflags) -I../pybind11/include  -std=c++14 -fPIC
LDFLAGS := $(shell python3-config --ldflags) 
PYBINDINC := -I../pybind11/include
CXX=x86_64-conda_cos6-linux-gnu-g++

default: mybind py


mybind.o: mybind.cpp
	$(CXX) -std=c++14 $(CXXFLAGS) $(PYBINDINC) -c $< -o $@

mybind: mybind.o
	$(CXX)  $^ -o $@ $(LDFLAGS)




test:
	PYTHONHOME=$(CONDA_PREFIX) ./py simple mul 3 4

