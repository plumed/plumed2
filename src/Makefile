
.PHONY: all clean lib

# by default we compile the lib directory, which in turn requests all the needeed modules
all: lib lib-static

lib: 
	$(MAKE) -C lib

lib-static: 
	$(MAKE) -C lib all-static
     
# list of all the contained directories which can be cleaned
SUBDIRS =$(subst /Makefile,,$(wildcard */Makefile))
SUBDIRSCLEAN=$(addsuffix .clean,$(SUBDIRS))
.PHONY: $(SUBDIRS)

# this is cleaning all subdirectories
clean: $(SUBDIRSCLEAN)

$(SUBDIRSCLEAN): %.clean:
	$(MAKE) -C $* clean

.PHONY: cppcheck
cppcheck:
	cppcheck --std=c++03 --std=posix -j 4 --platform=unix64 -q --enable=all */*.{h,cpp}

.PHONY: help
help:
	@echo "Available targets are"
	@echo "  all: compile everything (default)"
	@echo "  clean: clean everything"
	@echo "  cppcheck: check source (requires cppcheck installed)"


