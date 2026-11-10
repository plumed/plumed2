#!/usr/bin/env bash
set -euo pipefail
export PLUMED_AUTOLOAD=0
requires=()
runtime_setting=()
if [[ $(uname) == Darwin ]]; then
  python extensions/register.py --plumed "$PREFIX/bin/plumed" --module pycv_runtime \
    --library "../../python${PY_VER}/site-packages/pycv/PythonCVRuntime${SHLIB_EXT}" \
    --variant "python${PY_VER}" --output "$SRC_DIR/pycv-runtime-registration"
  requires=(--requires pycv_runtime)
  runtime_setting=(--config-settings=cmake.define.PLUMED_RUNTIME_IDENTITY="$SRC_DIR/pycv-runtime-registration/module_identity.cpp")
fi
python extensions/register.py --plumed "$PREFIX/bin/plumed" --module pycv \
  --library "../../python${PY_VER}/site-packages/pycv/PythonCVInterface${SHLIB_EXT}" \
  "${requires[@]}" --variant "python${PY_VER}" --features PYTHON_EMBED --output "$SRC_DIR/pycv-registration"
python -m pip install ./plugins/pycv --no-deps --no-build-isolation -v \
  --config-settings=cmake.define.Python_EXECUTABLE="$PYTHON" "${runtime_setting[@]}" \
  --config-settings=cmake.define.PLUMED_MODULE_IDENTITY="$SRC_DIR/pycv-registration/module_identity.cpp"
mkdir -p "$PREFIX/lib/plumed/plugins"
cp "$SRC_DIR/pycv-registration/pycv.plumed" "$PREFIX/lib/plumed/plugins/"
if [[ $(uname) == Darwin ]]; then
  cp "$SRC_DIR/pycv-runtime-registration/pycv_runtime.plumed" "$PREFIX/lib/plumed/plugins/"
fi
