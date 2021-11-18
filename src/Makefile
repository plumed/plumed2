
.PHONY: all clean lib install uninstall install-html astyle nmcheck

# by default we compile the lib directory, which in turn requests all the needeed modules
all: lib

lib: 
	$(MAKE) -C lib

install:
	$(MAKE) -C lib install

install-html:
	$(MAKE) -C lib install-html

uninstall:
	$(MAKE) -C lib uninstall

# list of all the contained directories which can be cleaned
SUBDIRS :=$(subst /Makefile,,$(wildcard */Makefile))
SUBDIRSCLEAN :=$(addsuffix .clean,$(SUBDIRS))
.PHONY: $(SUBDIRS)

# this is cleaning all subdirectories
clean: $(SUBDIRSCLEAN)
	rm -f codecheck.log codecheck.fatal

$(SUBDIRSCLEAN): %.clean:
	$(MAKE) -C $* clean

.PHONY: cppcheck
cppcheck:
	maketools/codecheck --cppcheck

.PHONY: plumedcheck
plumedcheck:
	maketools/codecheck --plumedcheck

.PHONY: codecheck
codecheck:
	maketools/codecheck

.PHONY: help
help:
	@echo "Available targets are"
	@echo "  all: compile everything (default)"
	@echo "  clean: clean everything"
	@echo "  codecheck: check source (requires cppcheck  and gawk installed)"
	@echo "  cppcheck: check source (requires cppcheck and gawk installed)"
	@echo "  plumedcheck: check source (requires gawk installed)"

astyle:
	$(MAKE) -C ../astyle
	./astyle.sh

nmcheck:
	maketools/nmcheck */*.o


