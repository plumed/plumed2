# include the machine dependent configuration
ifneq ($(MAKECMDGOALS),clean)
  -include ../Makefile.conf
endif

.PHONY: all clean

all:
ifeq ($(make_doc),yes)
	./extract
	./go-doxygen
	cp *png html/
endif

clean:
	rm -fr automatic html latex *~ links-to-user-manual.txt
