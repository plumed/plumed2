#!/bin/bash
#formatted with shfmt v3.6.0 https://github.com/mvdan/sh/releases
function usage() {
    echo "${0#*/} Runs cppcheck ONLY on the files in src that have changed from the selected base branch"
    echo 'By default it uses "master" as the base branch, but you can change it with --rev=REV:'
    echo ""
    echo "usage: ${0#*/} [--base=REV] change the base branch to use as base for the diff"
    echo "       ${0#*/} -h|--help   print this help"
    echo ""
    echo "for example \`${0#*/} --base=upstream/master\`"
}

baserev=${PLMD_MYCCPCHECK_REV:-master}

for arg; do
    case $arg in
    --base=*) baserev=${arg#*=} ;;
    -h | --help)
        usage
        exit 0
        ;;
    *)
        usage
        exit 1
        ;;
    esac
done
srcdir=$(git rev-parse --show-toplevel)

cd "$srcdir" >/dev/null || {
    echo "Could not cd to $srcdir"
    exit 1
}

files=$(
    git diff-index --diff-filter=ACMR "$baserev" | awk '/src.*(\.h|\.cpp|\.inc\.in)/{print $6}'
)
# shellcheck disable=SC2086
cppcheck --std=c++17 -j 4 --platform=unix64 --language=c++ \
    -U__GNUG__ -U__PLUMED_HAS_EXTERNAL_LAPACK -U__PLUMED_HAS_EXTERNAL_BLAS \
    -UGMX_CYGWIN -UF77_NO_UNDERSCORE -U_GLIBCXX_DEBUG -DNDEBUG -U__PLUMED_PBC_WHILE \
    -U__PLUMED_HAS_ASMJIT \
    -D__PLUMED_WRAPPER_CXX_EXPLICIT=explicit \
    --config-exclude=small_vector/ \
    --template='[{file}:{line}] ({severity}) :{id}: {message}' --enable=all --suppress=missingIncludeSystem --inline-suppr --force \
    $files
