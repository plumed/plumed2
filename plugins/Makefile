SUBDIRS := $(subst /,,$(dir $(shell ls */Makefile)))
SUBDIRSCLEAN := $(addsuffix .clean,$(SUBDIRS))


.PHONY: all clean $(SUBDIRSCLEAN)

all:

clean:$(SUBDIRSCLEAN)

$(SUBDIRSCLEAN): %.clean:
	$(MAKE) -C $* clean

