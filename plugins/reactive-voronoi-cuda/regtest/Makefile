PYTHON ?= python3
PLUMED ?= plumed
LIBRARY ?= ../libReactiveVoronoiCuda.so

.PHONY: check check-nogpu clean
check:
	rm -rf tmp
	$(PYTHON) backend_checks.py --plumed $(PLUMED) --library $(abspath $(LIBRARY)) --work $(CURDIR)/tmp

check-nogpu:
	rm -rf tmp-nogpu
	$(PYTHON) backend_checks.py --plumed "$(PLUMED)" --cpu-only --library "$(abspath $(LIBRARY))" --work "$(CURDIR)/tmp-nogpu"

clean:
	rm -rf tmp tmp-nogpu
