#!/bin/sh
# Defaults remain live: modules added after configuration follow module.yml.
set -e
src=$(CDPATH= cd "$(dirname "$0")/.." && pwd)
exec "$src/maketools/module-selection" "$src" --config "$src/../Makefile.conf"
