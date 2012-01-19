-include Makefile.conf


SRCDIRS = src extensions test
SUBDIRS = $(SRCDIRS) user-doc developer-doc imd regtest

SUBDIRSCLEAN=$(addsuffix .clean,$(SUBDIRS))

     
.PHONY: all clean $(SRCDIRS) doc docclean check

# if machine dependent configuration has been found:
ifdef GCCDEP
all: $(SRCDIRS)
     
$(SRCDIRS):
	$(MAKE) -C $@
     
# compile plumed before tests:
test: src

# doxygen
doc:
	$(MAKE) -C user-doc
	$(MAKE) -C developer-doc

# regtests
check: src test
	$(MAKE) -C regtest

else

all:
	@echo No configuration available
	@echo First run ./configure.sh
endif

# these targets are available also without configuration

clean: $(SUBDIRSCLEAN)
	rm -f *~ */*~ */*/*~

$(SUBDIRSCLEAN): %.clean:
	$(MAKE) -C $* clean

fullclean:
	make clean
	rm -f Makefile.conf
	rm -f sourceme.sh


docclean:
	cd user-doc && make clean
	cd developer-doc && make clean



