#! /bin/bash

set -e
set -x

echo "installing latest gawk"
wget http://git.savannah.gnu.org/cgit/gawk.git/snapshot/gawk-4.1.4.tar.gz
tar xzf gawk-4.1.4.tar.gz
cd gawk-4.1.4
./configure --prefix="$HOME/opt"
make -j 4
make install
cd ../

