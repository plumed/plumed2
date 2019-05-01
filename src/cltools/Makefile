USE=core config tools molfile

# generic makefile
include ../maketools/make.module

Completion.o: completion.xxd

completion.xxd: completion.sh ../../scripts/*.sh ../../patches/*.sh
	{ ../maketools/make-scripts-options.sh ; cat completion.sh ; }| xxd -i > completion.xxd

