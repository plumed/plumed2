# generic makefile
include ../maketools/make.module

CLEANLIST += completion.h completion-scripts.h

Completion.o: completion.h completion-scripts.h

completion.h: completion.sh
	../maketools/xxd < completion.sh > completion.h

completion-scripts.h: ../maketools/make-scripts-options.sh ../../scripts/*.sh ../../patches/*.sh
	../maketools/make-scripts-options.sh | ../maketools/xxd > completion-scripts.h

CXXFLAGS:= $(CXXFLAGS) -Wshadow -Werror
