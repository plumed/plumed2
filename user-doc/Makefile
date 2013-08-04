# include the machine dependent configuration
-include ../Makefile.conf

.PHONY: all clean

plumed_compiled := $(wildcard ../src/lib/plumed)

ifeq ($(strip $(plumed_compiled)),)

all:
	@echo You must compile plumed before building the user manual

else

all:
	./extract
	doxygen
	cp Octocat.png html/
	./report_errors

endif

clean:
	rm -fr html latex automatic *~ *PP.txt errors 
