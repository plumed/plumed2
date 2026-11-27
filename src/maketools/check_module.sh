#!/bin/sh
set -e
# Build support directories are not selectable modules.
case "$1" in config|wrapper) exit 0 ;; esac
src=$(CDPATH= cd "$(dirname "$0")/.." && pwd)
active=$("$src/maketools/module-selection" "$src" --config "$src/../Makefile.conf" single "$1")
if test -n "$active"; then exit 0; fi
echo "ERROR: this module depends on $1, which is presently disabled; rerun ./configure"
exit 1
