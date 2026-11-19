USE=core config tools xdrfile

# generic makefile
include ../maketools/make.module

Completion.o: completion.xxd completion-scripts.xxd

completion.xxd: completion.sh
	../maketools/xxd < completion.sh > completion.xxd

completion-scripts.xxd: ../maketools/make-scripts-options.sh ../../scripts/*.sh ../../patches/*.sh
	../maketools/make-scripts-options.sh | ../maketools/xxd > completion-scripts.xxd

CXXFLAGS:= $(CXXFLAGS) -Wshadow -Werror
