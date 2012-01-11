SUBDIRS = $(subst /Makefile,,$(wildcard rt*/Makefile))
SUBDIRSCLEAN = $(addsuffix .clean,$(SUBDIRS))

.PHONY: all clean $(SUBDIRS) $(SUBDIRSCLEAN)

all: $(SUBDIRS)
	scripts/check

clean: $(SUBDIRSCLEAN)

$(SUBDIRS):
	$(MAKE) -C $@

$(SUBDIRSCLEAN): %.clean:
	$(MAKE) -C $* clean
