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
# remove system headers:
	lcov -r coverage.info "/opt/*" "/usr/*" -o coverage-tmp.info && mv coverage-tmp.info coverage.info
# remove internal blas:
	lcov -r coverage.info "*/blas/*.cpp" "*/blas/*.h" -o coverage-tmp.info && mv coverage-tmp.info coverage.info
# remove internal lapack:
	lcov -r coverage.info "*/lapack/*.cpp" "*/lapack/*.h" -o coverage-tmp.info && mv coverage-tmp.info coverage.info
# remove internal molfile:
	lcov -r coverage.info "*/molfile/*.cpp" "*/molfile/*.h" -o coverage-tmp.info && mv coverage-tmp.info coverage.info
	genhtml -o coverage -t "plumed test coverage" coverage.info
