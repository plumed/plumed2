PYTHON ?= python3
PLUMED ?= plumed
LIBRARY ?= ../libReactiveVoronoiCuda.so

.PHONY: check clean
check:
	rm -rf tmp
	$(PYTHON) backend_checks.py --plumed $(PLUMED) --library $(abspath $(LIBRARY)) --work $(CURDIR)/tmp

clean:
	rm -rf tmp
