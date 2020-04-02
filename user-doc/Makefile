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
	make_pdfdoc=$(make_pdfdoc) ./go-doxygen
	cp *png html/
	./report_errors
endif
endif

spelling:
	./go-spelling

showerrors:
	./report_errors

clean:
	rm -fr example-check html latex automatic *~ *PP.txt *PP.md errors manual.pdf regtests
