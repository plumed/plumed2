#! /usr/bin/env bash

source "$PLUMED_ROOT"/src/config/compile_options.sh

#nvcc "$2" -Xcompiler -fPIC -c -o "$kernel"
compile="nvcc -g -ccbin ${compile}"
#compile=${compile/-O3/-g}

if [[ ${SILENT_CUDA_COMPILATION} ]]; then
  #echo "disabled warning"
  compile=${compile//-Wall/}
  compile=${compile//-pedantic/}
  #-w suppress the warnings
  compile=${compile/-c /-w -c }
fi

for opt in -W -pedantic -f; do
  compile=${compile//${opt}/-Xcompiler ${opt}}
done

link_command=$link_uninstalled

if [[ -z ${link_command:+x} ]]; then
  link_command=$link_installed
fi

link_command="nvcc -shared${link_command#*-shared}"
link_command=${link_command/-rdynamic/-Xcompiler -rdynamic}
link_command=${link_command/-Wl,/-Xlinker }
#link_command=${link_command/-fopenmp/-Xcompiler -fopenmp}
for opt in -f; do
  link_command=${link_command//${opt}/-Xcompiler ${opt}}
done

compile=${compile// -o/}
link_command=${link_command// -o/}

cat <<EOF
compile := $compile
link    := $link_command
all: Coordination.so

%.o: %.cu ndReduction.h cudaHelpers.cuh
	\$(compile) -o \$@  $<

Coordination.so: ndReduction.o Coordination.o CoordinationFloat.o
	\$(link) -lcusparse -o \$@ $^

clean:
	rm Coordination.so ndReduction.o Coordination.o CoordinationFloat.o
EOF