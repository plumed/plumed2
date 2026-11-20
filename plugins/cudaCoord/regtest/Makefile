SUBDIRS := $(subst /,,$(dir $(shell ls */Makefile)))

SUBDIRSCLEAN := $(addsuffix .clean,$(SUBDIRS))
SUBDIRSVALGRIND := $(addsuffix .valgrind,$(SUBDIRS))
SUBDIRSTESTCLEAN := $(addsuffix .testclean,$(SUBDIRS))

.PHONY: all checkfail checknogpu copytodoc clean valgrind $(SUBDIRS) $(SUBDIRSCLEAN) $(SUBDIRSVALGRIND) 

all: $(SUBDIRS)
	scripts/check


checkfail:
	scripts/check --fail

# runs only the tests in nogpu/, and fails if one of them fails
# (scripts/check would count the tests of the other directories, which are not run, as passed)
checknogpu:
	$(MAKE) -C nogpu
	@fail= ; for f in nogpu/rt*/report.txt ; do \
	  if ! test -f $$f || grep -q FAILURE $$f ; then echo "+ ERROR in $$f" ; cat $$f ; fail=1 ; fi ; \
	done ; test -z "$$fail" && echo "+ the tests that do not need a GPU passed"

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
