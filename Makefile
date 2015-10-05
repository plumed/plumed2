-include Makefile.conf


SRCDIRS := src test
SUBDIRS := $(SRCDIRS) user-doc developer-doc regtest

SUBDIRSCLEAN:=$(addsuffix .clean,$(SUBDIRS))

     
.PHONY: all lib clean $(SRCDIRS) doc docclean check cppcheck distclean

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

# standard target (according to GNU doc)
html:
	$(MAKE) doc

# standard target (according to GNU doc)
install-html: doc
	$(MAKE) -C src install-html

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
	@echo First run ./configure
endif

# these targets are available also without configuration

clean: $(SUBDIRSCLEAN)
	rm -f *~ */*~ */*/*~

$(SUBDIRSCLEAN): %.clean:
	$(MAKE) -C $* clean

distclean: fullclean

fullclean:
	make clean
	rm -f Makefile.conf
	rm -f sourceme.sh
	rm -f config.log 
	rm -f */*.on */*.off


docclean:
	cd user-doc && make clean
	cd developer-doc && make clean

cppcheck:
	$(MAKE) -C src cppcheck



