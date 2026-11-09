#!/usr/bin/env bash
set -euo pipefail
modules=$(python - <<'PY'
import json
with open('extensions/catalog.json') as stream:
    catalog = json.load(stream)
print(':'.join(name for name, module in catalog['modules'].items()
               if module['type'] == 'default-off' and not module['required_features']))
PY
)
./configure --with-extension-base="$PREFIX" --enable-extension-modules="$modules" --prefix="$PREFIX"
make
make install
