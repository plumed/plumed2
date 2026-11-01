# include the machine dependent configuration
ifneq ($(MAKECMDGOALS),clean)
  -include ../Makefile.conf
endif

.PHONY: all clean

plumed_compiled := $(wildcard ../src/lib/plumed)

# syntax.json is generated with "plumed gen_json", which is in the cltools module.
# cltools is a default-on module, so it is disabled when src/cltools.off exists
# (see src/maketools/find_modules.sh), e.g. with ./configure --enable-modules=none
cltools_disabled := $(wildcard ../src/cltools.off)

ifeq ($(strip $(plumed_compiled)),)

all:
	@echo You must compile plumed before building the JSON syntax file

else

ifneq ($(strip $(cltools_disabled)),)

# a file left by a previous build with cltools enabled would be installed otherwise
all:
	@echo The cltools module is disabled, so you cannot generate JSON syntax file
	rm -f syntax.json action_list

else

ifeq ($(program_can_run),no)

all:
	@echo PLUMED cannot run on this machine, so you cannot generate JSON syntax file

else

# syntax.json is only regenerated if the executable (which is relinked whenever the
# kernel changes), the scripts or the sources scanned to build the list of actions have changed.
# A partially written syntax.json is removed if something goes wrong.
.DELETE_ON_ERROR:

all: syntax.json

syntax.json: ../src/lib/plumed get_action_list.sh $(wildcard ../*/*/*.cpp)
	./get_action_list.sh > action_list
	../src/lib/plumed --no-mpi gen_json --actions action_list > syntax.json
ifeq (, $(shell which jq))
	@echo Could not find jq executable so unable to test if syntax.json is valid json.  Consider running sudo apt-get install jq
else
	cat syntax.json | jq empty
endif

endif

endif

endif

clean:
	rm -fr syntax.json action_list
