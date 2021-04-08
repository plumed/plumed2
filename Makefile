ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),fullclean)
 -include Makefile.conf
endif
endif
endif


SRCDIRS := src test
SUBDIRS := $(SRCDIRS) user-doc developer-doc regtest macports vim astyle python

SUBDIRSCLEAN:=$(addsuffix .clean,$(SUBDIRS))

     
.PHONY: all lib clean $(SRCDIRS) doc docclean check installcheck cppcheck distclean all_plus_docs macports codecheck plumedcheck astyle nmcheck

# if machine dependent configuration has been found:
ifdef GCCDEP
all:
	$(MAKE) lib
	$(MAKE) -C vim
	$(MAKE) -C python

# target useful for macports
# it builds the code then the documentation
all_plus_docs:
	$(MAKE) all
	$(MAKE) docs

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
# perform tests using non-installed plumed
check:
	LD_LIBRARY_PATH="$(realpath .)/src/lib":${LD_LIBRARY_PATH} PLUMED_PREPEND_PATH="$(realpath .)/src/lib" $(MAKE) -C regtest
	$(MAKE) -C regtest checkfail

# perform tests using the installed version of plumed
installcheck:
	PLUMED_PREPEND_PATH="$(bindir)" PLUMED_PROGRAM_NAME="$(program_name)" $(MAKE) -C regtest
	$(MAKE) -C regtest checkfail

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
	rm -f stamp-h
	rm -f */*.on */*.off


docclean:
	cd user-doc && make clean
	cd developer-doc && make clean

cppcheck:
	$(MAKE) -C src cppcheck

codecheck:
	$(MAKE) -C src codecheck

plumedcheck:
	$(MAKE) -C src plumedcheck

macports:
	$(MAKE) -C macports

# stamp-h file keeps track of when ./configure was last applied
# the procedure below is taken from:
# https://www.gnu.org/software/autoconf/manual/autoconf-2.69/html_node/Automatic-Remaking.html#Automatic-Remaking
Makefile.conf: stamp-h

stamp-h: sourceme.sh.in Makefile.conf.in config.status
	./config.status

config.status: configure
	./config.status --recheck

astyle:
	$(MAKE) -C astyle
	$(MAKE) -C src astyle

ifeq ($(use_debug_glibcxx),yes)
nmcheck:
	@echo SKIPPING nmcheck due to debug flags
else
nmcheck:
	$(MAKE) -C src nmcheck
endif


