
.PHONY: all clean lib

# by default we compile the lib directory, which in turn requests all the needeed modules
all: lib

lib: 
	$(MAKE) -C lib
     
# list of all the contained directories which can be cleaned
SUBDIRS =$(subst /Makefile,,$(wildcard */Makefile))
SUBDIRSCLEAN=$(addsuffix .clean,$(SUBDIRS))
.PHONY: $(SUBDIRS)

# this is cleaning all subdirectories
clean: $(SUBDIRSCLEAN)

$(SUBDIRSCLEAN): %.clean:
	$(MAKE) -C $* clean

.PHONY: help
help:
	@echo "Available targets are"
	@echo "  all: compile everything (default)"
	@echo "  clean: clean everything"

