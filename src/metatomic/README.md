# Metatomic module for PLUMED


## Building the code

See [the main documentation](../../user-doc/METATOMICMOD.md) for more
information on how to compile and use this module.


This module uses [vesin](https://github.com/Luthaf/vesin) to compute neighbor
lists. If you need to update it, you should use the `import.sh` script:

```bash
git clone https://github.com/Luthaf/vesin
cd /path/to/plumed2/src/metaomic
./import.sh /path/to/vesin
```
