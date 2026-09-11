SUBDIRS := $(subst /,,$(dir $(shell ls */Makefile)))
SUBDIRSCLEAN := $(addsuffix .clean,$(SUBDIRS))
SUBDIRSASTYLE := $(addsuffix .astyle,$(SUBDIRS))


.PHONY: all clean $(SUBDIRSCLEAN) astyle $(SUBDIRSASTYLE)

all:
	@echo The plugins should be compiled and installed with their own rules, please cd in the directroies you are interested in and follow the instructions

clean:$(SUBDIRSCLEAN)

$(SUBDIRSCLEAN): %.clean:
	$(MAKE) -C $* clean

astyle:$(SUBDIRSASTYLE)

$(SUBDIRSASTYLE): %.astyle:
	$(MAKE) -C $* astyle

