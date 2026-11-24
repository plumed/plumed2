#!/usr/bin/env bash
set -euo pipefail
install -m 755 scripts/mklib.sh "$PREFIX/lib/plumed/scripts/mklib.sh"
# The alternate shell-only entry point uses the selected kernel's environment.
cat > "$PREFIX/lib/plumed/plumed-mklib" <<'LAUNCHER'
#!/usr/bin/env bash
exec "$(dirname "$0")/../../bin/plumed" mklib "$@"
LAUNCHER
chmod +x "$PREFIX/lib/plumed/plumed-mklib"
