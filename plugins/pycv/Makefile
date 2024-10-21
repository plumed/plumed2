#this makefiles assume that pip and pytest are installed
.PHONY: clean check check_standalone check_python all

all: pycv_here

pycv_here: src/*.cpp src/*.h src/pycv/*.py
	@echo installing pycv
	pip3 install .
	touch $@

clean:
	pip3 uninstall pycv
	rm -fv pycv_here

check_standalone: pycv_here
	$(MAKE) -C regtest testclean
	$(MAKE) -C regtest checkfail

#just in case  pytest is still not installed we install it before the tests
check_python: pycv_here
	pip3 install pytest
	pytest -v

check: check_standalone check_python
