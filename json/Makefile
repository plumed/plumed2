# include the machine dependent configuration
ifneq ($(MAKECMDGOALS),clean)
  -include ../Makefile.conf
endif

.PHONY: all clean

plumed_compiled := $(wildcard ../src/lib/plumed)

ifeq ($(strip $(plumed_compiled)),)

all:
	@echo You must compile plumed before building the JSON syntax file

else

ifeq ($(program_can_run),no)

all:
	@echo PLUMED cannot run on this machine, so you cannot generate JSON syntax file

else

all:
	./get_action_list.sh > action_list
	../src/lib/plumed --no-mpi gen_json --actions action_list > syntax.json
ifeq (, $(shell which jq))
	@echo Could not find jq executable so unable to test if syntax.json is valid json.  Consider running sudo apt-get install jq
else
	cat syntax.json | jq empty
endif

endif

endif

clean:
	rm -fr syntax.json action_list
