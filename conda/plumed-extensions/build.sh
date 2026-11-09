#!/usr/bin/env bash
set -euo pipefail
./configure --with-extension-base="$PREFIX" --enable-extension-modules=all --prefix="$PREFIX"
make
make install
