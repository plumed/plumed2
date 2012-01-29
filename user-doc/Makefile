# include the machine dependent configuration
-include ../Makefile.conf

.PHONY: all clean

all: 
	./extract
	doxygen
	./report_errors  

clean:
	rm -fr html latex automatic *~ *PP.txt errors 
