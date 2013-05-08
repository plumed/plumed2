-include ../Makefile.conf

SUBDIRS = link

.PHONY: all clean $(SUBDIRS)

all: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

clean:
	cd link && make clean

