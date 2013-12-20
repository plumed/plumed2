SUBDIRS := $(subst /Makefile,,$(wildcard rt*/Makefile))
SUBDIRSCLEAN := $(addsuffix .clean,$(SUBDIRS))
SUBDIRSVALGRIND := $(addsuffix .valgrind,$(SUBDIRS))

.PHONY: all clean valgrind $(SUBDIRS) $(SUBDIRSCLEAN) $(SUBDIRSVALGRIND)

all: $(SUBDIRS)

clean: $(SUBDIRSCLEAN)

valgrind: $(SUBDIRSVALGRIND)

$(SUBDIRS):
	$(MAKE) -C $@

$(SUBDIRSCLEAN): %.clean:
	$(MAKE) -C $* clean

$(SUBDIRSVALGRIND): %.valgrind:
	$(MAKE) -C $* valgrind
