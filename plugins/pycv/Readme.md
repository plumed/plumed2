## Preparation
```bash
python3 -m venv pycvenv
source ./pycvenv/bin/activate
pip install -U pip
pip install -r requirements
```
jax is not required for installation, but for some test
### Install jax:

Go to the original guide in the [jax documenation](https://jax.readthedocs.io/en/latest/installation.html)

The command shoudl be similar to:
 - example if you have a cuda12 compatible device (a wheel for cuda will be installed alognside jax):
`pip install "jax[cuda12_pip]" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html`
 - example if you have a cuda12 compatible device, and cuda already installed on your system:
`pip install "jax[cuda12_local]" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html`

## Standard compilation

If you have a plumed that supports plumed mklib with multiple files
```bash
./standaloneCompile.sh
```

## Develop

This version is slightly more complex, but if are often modifying the cpp files compiles only the modified files. (You'll need to call `./prepareMakeForDevelop.sh` only the first time)
```bash
./prepareMakeForDevelop.sh
make
```

## Set up tests

Create a symbolic link to the `scripts` directory in a plumed source and you can execute the tests. Plumed must be runnable to execute tests
```bash
ln -s path/to/plumed/source/regtest/scripts .
make
```