# include the machine dependent configuration
ifneq ($(MAKECMDGOALS),clean)
  -include ../Makefile.conf
endif

.PHONY: all clean spelling

plumed_compiled := $(wildcard ../src/lib/plumed)

ifeq ($(strip $(plumed_compiled)),)

all:
	@echo You must compile plumed before building the user manual

else

all:
ifeq ($(make_doc),yes)
	make_pdfdoc=$(make_pdfdoc) ./extract
	# check examples in input file:
	#./go-example-check
	make_pdfdoc=$(make_pdfdoc) ./go-doxygen
	cp *png html/
	./report_errors
endif
endif

spelling:
	./go-spelling

clean:
	rm -fr html latex automatic *~ *PP.txt *PP.md errors manual.pdf regtests
