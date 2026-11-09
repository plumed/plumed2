.PHONY: all install clean
all:
	$(EXTENSION_PYTHON) extensions/build.py build
install: all
	$(EXTENSION_PYTHON) extensions/build.py install
clean:
	rm -rf _extensions
