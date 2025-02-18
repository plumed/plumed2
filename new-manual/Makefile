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
	rm -rf docs site
	python build_manual.py
	mkdocs build
endif
endif

clean:
	rm -rf docs site

