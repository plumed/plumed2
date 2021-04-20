SUBDIRS := $(subst /Makefile,,$(wildcard rt*/Makefile))
SUBDIRSCLEAN := $(addsuffix .clean,$(SUBDIRS))
SUBDIRSTESTCLEAN := $(addsuffix .testclean,$(SUBDIRS))
SUBDIRSVALGRIND := $(addsuffix .valgrind,$(SUBDIRS))

.PHONY: all clean valgrind $(SUBDIRS) $(SUBDIRSCLEAN) $(SUBDIRSVALGRIND)

all: $(SUBDIRS)

clean: $(SUBDIRSCLEAN)

valgrind: $(SUBDIRSVALGRIND)

testclean: $(SUBDIRSTESTCLEAN)

$(SUBDIRS):
	$(MAKE) -C $@

$(SUBDIRSCLEAN): %.clean:
	$(MAKE) -C $* clean

$(SUBDIRSVALGRIND): %.valgrind:
	$(MAKE) -C $* valgrind

$(SUBDIRSTESTCLEAN): %.testclean:
	$(MAKE) -C $* testclean

