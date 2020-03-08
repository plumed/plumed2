# include the machine dependent configuration
ifneq ($(MAKECMDGOALS),clean)
  -include ../Makefile.conf
endif

.PHONY: all clean

plumed_compiled := $(wildcard ../src/lib/plumed)

ifeq ($(strip $(plumed_compiled)),)

all:
	@echo You must compile plumed before building the VIM syntax file

else

ifeq ($(program_can_run),no)

all:
	@echo PLUMED cannot run on this machine, so you cannot generate VIM syntax file

else

all:
	./vimsyntax.sh

endif

endif

clean:
	rm -fr syntax help
