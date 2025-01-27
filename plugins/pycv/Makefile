
PYTHON=python
#optional, just in case you want to override python with ./configurePyCV.sh
-include Makefile.conf
#this makefiles assume that pip and pytest are installed
.PHONY: clean check check_standalone check_python all

all: pycv_here

pycv_here: src/*.cpp src/*.h src/pycv/*.py
	@echo installing pycv
	$(PYTHON) -m pip install .
	@touch $@

clean:
	@$(PYTHON) -m  pip uninstall pycv -y
	@rm -fv pycv_here

check_standalone: pycv_here
	$(MAKE) -C regtest testclean
	$(MAKE) -C regtest checkfail

#just in case  pytest is still not installed we install it before the tests
check_python: pycv_here
	@$(PYTHON) -m pip install pytest
	@$(PYTHON) -m pytest -v

check: check_standalone check_python
