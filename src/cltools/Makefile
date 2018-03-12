USE=core wrapper config tools molfile

# generic makefile
include ../maketools/make.module

Completion.o: completion.xxd

completion.xxd: completion.sh
	cat completion.sh | xxd -i > completion.xxd

