#! /bin/bash

# This script generates a Portfile in science/plumed
# Currently the portfile is aimed at testing the currect git hash
# TODO:
# - Analyze configure.ac to generate the list of variants
# - Allow for a proper release port.
#   This would require a portfile that is based on a tag (e.g. v2.3.0),
#   optionally including patches, that can be then uploaded to macports

prefix=
if git describe --exact-match --tags HEAD 2>/dev/null 1>/dev/null
then
  version=$(git describe --exact-match --tags HEAD)
  case "$version" in
  (v*)
    version=${version#v}
    prefix=v
  esac
else
  version=$( git log -1 --format="%h")
fi

if test -n "$plumed_repository" ; then
# get this from environment
  repository="$plumed_repository"
else
# parent directory:
  repository="${PWD%/*}"
fi

mkdir -p science/plumed

cat Portfile.in |
sed "
  s/@_VERSION_@/$version/
  s/@_REVISION_@/0/
" | awk '{
  if($1=="@_FETCH_@"){
    print "fetch.type          git"
    print "git.url             '$repository'"
# notice that if instead of hashtag we want to put a version, then it should be
# git.branch          v${version}
    print "git.branch          '$prefix'${version}"
  } else print
}'  > science/plumed/Portfile




