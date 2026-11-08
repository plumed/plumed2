# Library publication check

Run `bash regtest/mklib/check-mklib-publication.sh` to verify that the older
single-source `scripts/mklib.sh` interface keeps its previous output when a
replacement link fails. The checker uses stub compiler and linker commands and
does not require a configured PLUMED build.
