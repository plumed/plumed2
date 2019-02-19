# include the machine dependent configuration
ifneq ($(MAKECMDGOALS),clean)
  -include ../Makefile.conf
endif

.PHONY: all clean

plumed_compiled := $(wildcard ../src/lib/plumed)

ifeq ($(strip $(plumed_compiled)),)

all:
	@echo You must compile plumed before building the user manual

else

all:
ifeq ($(make_doc),yes)
	make_pdfdoc=$(make_pdfdoc) ./extract
	# temporarily removed for travis to build:
	#./go-spelling
	make_pdfdoc=$(make_pdfdoc) ./go-doxygen
	cp *png html/
	./report_errors
endif
endif

clean:
	rm -fr html latex automatic *~ *PP.txt *PP.md errors manual.pdf regtests
