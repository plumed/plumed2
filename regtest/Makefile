SUBDIRS := $(subst /,,$(dir $(shell ls */Makefile)))

SUBDIRSCLEAN := $(addsuffix .clean,$(SUBDIRS))
SUBDIRSVALGRIND := $(addsuffix .valgrind,$(SUBDIRS))
SUBDIRSTESTCLEAN := $(addsuffix .testclean,$(SUBDIRS))

.PHONY: all checkfail copytodoc  clean valgrind $(SUBDIRS) $(SUBDIRSCLEAN) $(SUBDIRSVALGRIND) 

all: $(SUBDIRS)
	scripts/check
	scripts/what

checkfail:
	scripts/check --fail

copytodoc:
	scripts/copytodoc

clean: $(SUBDIRSCLEAN)

valgrind: $(SUBDIRSVALGRIND)
	scripts/check

testclean: $(SUBDIRSTESTCLEAN)

$(SUBDIRS):
	$(MAKE) -C $@

$(SUBDIRSCLEAN): %.clean:
	$(MAKE) -C $* clean

$(SUBDIRSVALGRIND): %.valgrind:
	$(MAKE) -C $* valgrind

$(SUBDIRSTESTCLEAN): %.testclean:
	$(MAKE) -C $* testclean

