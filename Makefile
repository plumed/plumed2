-include Makefile.conf


SRCDIRS := src test
SUBDIRS := $(SRCDIRS) user-doc developer-doc regtest

SUBDIRSCLEAN:=$(addsuffix .clean,$(SUBDIRS))

     
.PHONY: all lib clean $(SRCDIRS) doc docclean check

# if machine dependent configuration has been found:
ifdef GCCDEP
all:
	$(MAKE) lib

lib:
	$(MAKE)	-C src

install:
	$(MAKE) -C src install

uninstall:
	$(MAKE) -C src uninstall

$(SRCDIRS):
	$(MAKE) -C $@
     
# compile plumed before tests:
test: src

# doxygen
doc:
	$(MAKE) -C user-doc
	$(MAKE) -C developer-doc

docs:
	$(MAKE) doc

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
	rm -fr autoconf/auto*
	rm -f autoconf/Makefile.conf
	rm -f autoconf/sourceme.sh
	rm -f autoconf/config.*


docclean:
	cd user-doc && make clean
	cd developer-doc && make clean



