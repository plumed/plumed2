USE=core config tools molfile xdrfile

# generic makefile
include ../maketools/make.module

Completion.o: completion.xxd

completion.xxd: completion.sh ../../scripts/*.sh ../../patches/*.sh
	{ ../maketools/make-scripts-options.sh ; cat completion.sh ; }| ../maketools/xxd  > completion.xxd

