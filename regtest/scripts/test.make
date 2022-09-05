
test:
	@env PLUMED_MAKE=$(MAKE) ../../scripts/run

reset:
	../../scripts/reset

clean:
	rm -fr tmp/ report.txt

valgrind:
	../../scripts/run --valgrind

testclean:
	$(MAKE) test
	rm -fr tmp/

