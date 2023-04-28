tests := $(patsubst %.cpp,%,$(wildcard *.cpp))

.PHONY: test reset clean valgrind testclean

-include Makefile.conf

test: $(tests)
	@{ for test in $^ ;do \
	printf  "%-20s: " "$${test}"; \
	if ! ./$$test > /dev/null; then echo "fail, run ./$${test} to see the output"; \
	else echo "success"; fi ;\
	done } > output

#-w deactivate the warnings in gcc, these tests are not aimed at catching warnings
%:%.cpp
	@$(CXX) $(CPPFLAGS) $(ADDCPPFLAGS) $(CXXFLAGS) ${PLUMED_LIB} ${PLUMED_KERNEL} -w $< -o $@

reset:
	../../scripts/reset

clean:
	rm -fr tmp/ report.txt

valgrind:
	../../scripts/run --valgrind

testclean:
	$(MAKE) test
	rm -fr tmp/

