
tests := $(patsubst %.cpp,%,$(wildcard *.cpp))
#needs plumed installed with the environmnet variables set up to work (eg with a module)
.PHONY: all $(run)

all: $(tests)
	@for test in $^ ;do \
	printf  "%-20s: " "$${test}"; \
	if ! ./$$test > /dev/null; then echo "fail"; \
	else echo "success"; fi ;\
	done

%:%.cpp
	$(CXX) -lplumed -lplumedKernel $< -o $@

clean:
	rm $(tests)
	rm *.log