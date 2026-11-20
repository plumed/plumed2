USE=core config tools xdrfile

# generic makefile
include ../maketools/make.module

CLEANLIST += completion-scripts.h

Completion.o: completion.xxd completion-scripts.h

completion.xxd: completion.sh
	../maketools/xxd < completion.sh > completion.xxd

completion-scripts.h: ../maketools/make-scripts-options.sh ../../scripts/*.sh ../../patches/*.sh
	../maketools/make-scripts-options.sh | ../maketools/xxd > completion-scripts.h

CXXFLAGS:= $(CXXFLAGS) -Wshadow -Werror
