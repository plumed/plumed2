-include Makefile.conf

SUBDIRS = src extensions test
     
.PHONY: all clean $(SUBDIRS) doc docclean

ifdef GCCDEP
all: $(SUBDIRS)
else
all:
	@echo No configuration available
	@echo First run ./configure.sh
endif
     
$(SUBDIRS):
	$(MAKE) -C $@
     
# compile plumed before tests:
test: src

clean:
	cd src && make clean
	cd extensions && make clean
	cd test && make clean
	cd imd && make clean
	rm -f *~ */*~ */*/*~

fullclean:
	make clean
	rm -f Makefile.conf
	rm -f sourceme.sh

doc:
	cd user-doc && make
	cd developer-doc && make

docclean:
	cd user-doc && make clean
	cd developer-doc && make clean


