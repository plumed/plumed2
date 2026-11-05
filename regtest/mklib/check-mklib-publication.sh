#!/usr/bin/env bash
set -euo pipefail

repo=$(cd "$(dirname "$0")/../.." && pwd)
mklib=${MKLIB_SCRIPT:-$repo/scripts/mklib.sh}
fixture=$(mktemp -d "${TMPDIR:-/tmp}/plumed-mklib-publication.XXXXXX")
child=
cleanup() {
  touch "$fixture/release" 2>/dev/null || true
  if [[ -n "$child" ]]; then wait "$child" 2>/dev/null || true; fi
  rm -rf "$fixture"
}
trap cleanup EXIT
assert_contents() {
  if [[ ! -f $1 ]] || [[ $(cat "$1") != "$2" ]]; then
    echo "unexpected contents or missing library: $1 (expected $2)" >&2
    exit 1
  fi
}
mkdir -p "$fixture/src/config" "$fixture/output"
export PLUMED_ROOT="$fixture" FIXTURE="$fixture"

cat > "$fixture/compiler" <<'EOF'
#!/usr/bin/env bash
set -eu
while [[ $# -gt 0 ]]; do
  case "$1" in
    -o) output=$2; shift 2 ;;
    *) input=$1; shift ;;
  esac
done
cp "$input" "$output"
EOF
cat > "$fixture/linker" <<'EOF'
#!/usr/bin/env bash
set -eu
while [[ $# -gt 0 ]]; do
  case "$1" in
    -o) output=$2; shift 2 ;;
    *) input=$1; shift ;;
  esac
done
touch "$FIXTURE/link-started"
for ((attempt=0; attempt<500; ++attempt)); do
  [[ -f "$FIXTURE/release" ]] && break
  sleep 0.01
done
[[ -f "$FIXTURE/release" ]]
[[ ${FAIL_LINK:-no} != yes ]] || exit 1
cp "$input" "$output"
EOF
chmod +x "$fixture/compiler" "$fixture/linker"
cat > "$fixture/src/config/compile_options.sh" <<'EOF'
compile='"$FIXTURE/compiler"'
link_shared='"$FIXTURE/linker"'
link_libs_uninstalled=''
soext=so
EOF

cd "$fixture"
printf 'new library\n' > input.cpp
printf 'old library\n' > output/library.so
touch -t 202001010000 output/library.so

bash "$mklib" -n -o output/library.so input.cpp > build.log 2>&1 &
child=$!
for ((attempt=0; attempt<500; ++attempt)); do
  [[ -f link-started ]] && break
  sleep 0.01
done
if [[ ! -f link-started ]]; then
  echo 'linker did not start' >&2
  exit 1
fi
# A reader must still see the previous complete library during the link.
assert_contents output/library.so 'old library'
touch release
wait "$child"
child=
assert_contents output/library.so 'new library'

# A failed replacement must preserve the last usable library as well.
printf 'preserved library\n' > output/library.so
if FAIL_LINK=yes bash "$mklib" -o output/library.so input.cpp > failure.log 2>&1; then
  echo 'mklib accepted a failed link' >&2
  exit 1
fi
assert_contents output/library.so 'preserved library'

# The no-clobber fast path must not invoke the linker for up-to-date inputs.
touch -t 202001010000 input.cpp
rm link-started
bash "$mklib" -n -o output/library.so input.cpp > current.log 2>&1
if [[ -f link-started ]]; then
  echo 'mklib relinked an up-to-date library' >&2
  exit 1
fi
echo 'mklib atomic publication and failed-link preservation passed'
