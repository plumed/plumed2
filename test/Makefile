# include the machine dependent configuration
ifneq ($(MAKECMDGOALS),clean)
  -include ../Makefile.conf
endif

SUBDIRS = link

.PHONY: all clean $(SUBDIRS)

all: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

clean:
	cd link && make clean

