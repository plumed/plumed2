.PHONY: astyle

astyle:
	$(MAKE) -C build/gcc
	ln -fs build/gcc/bin/astyle .

clean:
	$(MAKE) -C build/gcc clean
	rm -f astyle
