ACTIVE_MODULES := $(shell ./scripts/find_modules.sh) basic
#SUBDIRS = $(subst /Makefile,,$(wildcard rt*/Makefile))
SUBDIRS := $(shell for dir in $(ACTIVE_MODULES) ; do ls $$dir/rt*/Makefile | sed 's/\Makefile//' 2>/dev/null; done)
SUBDIRSCLEAN := $(addsuffix .clean,$(SUBDIRS))
SUBDIRSVALGRIND := $(addsuffix .valgrind,$(SUBDIRS))

.PHONY: all clean valgrind $(SUBDIRS) $(SUBDIRSCLEAN) $(SUBDIRSVALGRIND)

all: $(SUBDIRS)
	scripts/check
	scripts/what

clean: $(SUBDIRSCLEAN)

valgrind: $(SUBDIRSVALGRIND)
	scripts/check

$(SUBDIRS):
	$(MAKE) -C $@

$(SUBDIRSCLEAN): %.clean:
	$(MAKE) -C $* clean

$(SUBDIRSVALGRIND): %.valgrind:
	$(MAKE) -C $* valgrind
