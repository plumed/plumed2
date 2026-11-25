#!/usr/bin/env bash
set -euo pipefail
repo=$(cd "$(dirname "$0")/../.." && pwd)
fixture=$(mktemp -d "${TMPDIR:-/tmp}/plumed-mklib-publication.XXXXXX")
trap 'rm -rf "$fixture"' EXIT
mkdir -p "$fixture/src/config" "$fixture/output"
export PLUMED_ROOT="$fixture" FIXTURE="$fixture"
cat > "$fixture/compiler" <<'COMPILER'
#!/usr/bin/env bash
set -eu
if [[ " $* " == *" -o "* ]]; then
  output=
  input=
  while [[ $# -gt 0 ]]; do
    case "$1" in
      -o) output=$2; shift 2 ;;
      *) [[ -f $1 ]] && input=$1; shift ;;
    esac
  done
else
  output=$1
  input=$2
fi
cp "$input" "$output"
COMPILER
cat > "$fixture/linker" <<'LINKER'
#!/usr/bin/env bash
set -eu
[[ ${FAIL_LINK:-no} != yes ]] || exit 1
if [[ " $* " == *" -o "* ]]; then
  output=
  input=
  while [[ $# -gt 0 ]]; do
    case "$1" in
      -o) output=$2; shift 2 ;;
      *) [[ -f $1 ]] && input=$1; shift ;;
    esac
  done
else
  output=$1
  input=$2
fi
cp "$input" "$output"
LINKER
chmod +x "$fixture/compiler" "$fixture/linker"
cat > "$fixture/src/config/compile_options.sh" <<'OPTIONS'
compile='"$FIXTURE/compiler"'
link_uninstalled='"$FIXTURE/linker"'
link_installed='"$FIXTURE/linker"'
link_shared='"$FIXTURE/linker"'
link_libs_uninstalled=''
link_libs_installed=''
PLUMED_IS_INSTALLED=no
soext=so
OPTIONS
cd "$fixture"
printf 'new library\n' > input.cpp
printf 'old library\n' > input.so
bash "$repo/scripts/mklib.sh" input.cpp
[[ $(cat input.so) == 'new library' ]]
printf 'preserved library\n' > input.so
if FAIL_LINK=yes bash "$repo/scripts/mklib.sh" input.cpp; then
  echo 'mklib accepted a failed link' >&2
  exit 1
fi
[[ $(cat input.so) == 'preserved library' ]]
echo 'mklib failed-link preservation passed'
