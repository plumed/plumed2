PYBINDINC := -I../pybind11-2.3.0/include
#CXXFLAGS := $(shell python3.7m-config --cflags) $(PYBINDINC)  -std=c++14 -fPIC -fno-lto
#LDFLAGS := $(shell python3.7m-config --ldflags) -fno-lto -Wl,-export_dynamic

CXXFLAGS := -std=c++14 -I/Users/toni/anaconda3/envs/jax/include/python3.7m -I/Users/toni/anaconda3/envs/jax/include/python3.7m  -g  -Wall -Wstrict-prototypes -I/Users/toni/anaconda3/envs/jax/include -arch x86_64 -I/Users/toni/anaconda3/envs/jax/include -arch x86_64 -I../pybind11-2.3.0/include  -std=c++14 -fPIC -fno-lto -I../pybind11-2.3.0/include
LDFLAGS := -L/Users/toni/anaconda3/envs/jax/lib/python3.7/config-3.7m-darwin -lpython3.7m -ldl -framework CoreFoundation -Wl,-stack_size,1000000 -framework CoreFoundation -fno-lto -Wl,-export_dynamic

# The export dynamic thing was required in OSX to avoid obscure
# "symbol not found" runtime errors, which may be related to
# LTO. https://bugs.python.org/issue21122

# LTO also seems to slow down linking to the extreme, and sometimes
# "internal errors" the compiler, but can't find a way to disable it
# in OSX.

#CXX=x86_64-conda_cos6-linux-gnu-g++

default: mybind py


mybind.o: mybind.cpp
	$(CXX) -std=c++14 $(CXXFLAGS) $(PYBINDINC) -c $< -o $@

mybind: mybind.o
	$(CXX)  $^ -o $@ $(LDFLAGS)

clean:
	rm mybind.o mybind


test:
	PYTHONHOME=$(CONDA_PREFIX) ./py simple mul 3 4

