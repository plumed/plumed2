# include the machine dependent configuration
-include ../Makefile.conf

.PHONY: all clean

all: GENMANUAL
	./extract
	doxygen
	./report_errors  

GENMANUAL: GenManual.cpp ../src/libplumed.$(SOEXT)
	$(CC) GenManual.cpp -ldl -lstdc++ ../src/libplumed.$(SOEXT) -o GENMANUAL 

clean:
	rm -fr html latex automatic *~ *PP.txt GENMANUAL
