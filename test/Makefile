-include ../Makefile.conf

SUBDIRS = link simplemd

.PHONY: all clean $(SUBDIRS)

all: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

clean:
	cd link && make clean
	cd simplemd && make clean

