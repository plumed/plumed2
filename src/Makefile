
.PHONY: all clean lib install uninstall

# by default we compile the lib directory, which in turn requests all the needeed modules
all: lib

lib: 
	$(MAKE) -C lib

install:
	$(MAKE) -C lib install

uninstall:
	$(MAKE) -C lib uninstall

# list of all the contained directories which can be cleaned
SUBDIRS :=$(subst /Makefile,,$(wildcard */Makefile))
SUBDIRSCLEAN :=$(addsuffix .clean,$(SUBDIRS))
.PHONY: $(SUBDIRS)

# this is cleaning all subdirectories
clean: $(SUBDIRSCLEAN)
	rm -f cppcheck.log cppcheck.fatal

$(SUBDIRSCLEAN): %.clean:
	$(MAKE) -C $* clean

.PHONY: cppcheck
cppcheck:
	maketools/cppcheck

.PHONY: mod-light mod-heavy mod-reset
mod-light:
	maketools/mod.sh light
mod-heavy:
	maketools/mod.sh heavy
mod-reset:
	maketools/mod.sh reset
	

.PHONY: help
help:
	@echo "Available targets are"
	@echo "  all: compile everything (default)"
	@echo "  clean: clean everything"
	@echo "  cppcheck: check source (requires cppcheck installed)"
	@echo "  mod-light: disable all optional modules"
	@echo "  mod-heavy: enable all optional modules"
	@echo "  mod-reset: reset optional modules to their default"


