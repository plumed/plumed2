#!/usr/bin/env bash
set -euo pipefail
export PLUMED_AUTOLOAD=0
python extensions/register.py --plumed "$PREFIX/bin/plumed" --module pycv \
  --library "../../python${PY_VER}/site-packages/pycv/PythonCVInterface${SHLIB_EXT}" \
  --variant "python${PY_VER}" --features PYTHON_EMBED --output "$SRC_DIR/pycv-registration"
python extensions/conda_config.py --module pycv -- python -m pip install ./plugins/pycv --no-deps --no-build-isolation -v \
  --config-settings=cmake.define.Python_EXECUTABLE="$PYTHON" \
  --config-settings=cmake.define.PLUMED_MODULE_IDENTITY="$SRC_DIR/pycv-registration/module_identity.cpp"
mkdir -p "$PREFIX/lib/plumed/plugins"
cp "$SRC_DIR/pycv-registration/pycv.plumed" "$PREFIX/lib/plumed/plugins/"

python extensions/conda_config.py --module pycv --check-features --prefix "$PREFIX"
