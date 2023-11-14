## Preparation
```bash
python3 -m venv pycvenv
source ./pycvenv/bin/activate
pip install -U pip
pip install -r requirements
```
jax is not required for installation, but for some test

## Install

```bash
./tandaloneCompile.sh
```

## Develop

```bash
./prepareMakeForDevelop.sh
make
```