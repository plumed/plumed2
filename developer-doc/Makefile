# include the machine dependent configuration
ifneq ($(MAKECMDGOALS),clean)
  -include ../Makefile.conf
endif

.PHONY: all clean coverage

all:
ifeq ($(make_doc),yes)
	./extract
	./go-doxygen
	cp *png html/
endif

clean:
	rm -fr automatic html latex *~ links-to-user-manual.txt coverage coverage.info

coverage:
	lcov -c -d ../src/ --output-file coverage.info
	genhtml -o coverage -t "plumed test coverage" coverage.info
