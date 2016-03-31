# include the machine dependent configuration
-include ../Makefile.conf

.PHONY: all clean

all:
ifeq ($(make_doc),yes)
	./extract
	./go-doxygen
	cp *png html/
endif

clean:
	rm -fr automatic html latex *~ links-to-user-manual.txt
