# Library publication check

Run `bash regtest/mklib/check-mklib-publication.sh` to exercise atomic library
publication, preservation of the last usable library after a failed link, and
the up-to-date no-clobber path in `scripts/mklib.sh`. The checker builds a
temporary fixture with stub compiler and linker commands and does not require a
configured PLUMED build.
