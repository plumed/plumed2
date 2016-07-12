#! /bin/bash

# This script generates a Portfile in science/plumed
# Currently the portfile is aimed at testing the currect git hash
# TODO:
# - Analyze configure.ac to generate the list of variants
# - Allow for a proper release port.
#   This would require a portfile that is based on a tag (e.g. v2.3.0),
#   optionally including patches, that can be then uploaded to macports

hash=$( git log -1 --format="%h")

mkdir -p science/plumed

# parent directory:
repository="${PWD%/*}"

cat Portfile.in |
sed "
  s/@_VERSION_@/$hash/
  s/@_REVISION_@/0/
" | awk '{
  if($1=="@_FETCH_@"){
    print "fetch.type          git"
    print "git.url             '$repository'"
# notice that if instead of hashtag we want to put a version, then it should be
# git.branch          v${version}
    print "git.branch          ${version}"
  } else print
}'  > science/plumed/Portfile




