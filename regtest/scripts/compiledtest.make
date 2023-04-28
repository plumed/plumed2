tests := $(patsubst %.cpp,%,$(wildcard *.cpp))

.PHONY: dotest test reset clean valgrind testclean

test:
	@env PLUMED_MAKE=$(MAKE) ../../scripts/run

-include Makefile.conf


dotest: $(tests)
	@for test in $^ ;do \
	printf  "%-20s: " "$${test}"; \
	if ! ./$$test > /dev/null; then echo "fail, run ./$${test} to see the output"; \
	else echo "success"; fi ;\
	done

%:%.cpp
	$(CXX) $(CPPFLAGS) $(ADDCPPFLAGS) $(CXXFLAGS) ${PLUMED_LIB} ${PLUMED_KERNEL} $< -o $@

reset:
	../../scripts/reset

clean:
	rm -fr tmp/ report.txt

valgrind:
	../../scripts/run --valgrind

testclean:
	$(MAKE) test
	rm -fr tmp/

