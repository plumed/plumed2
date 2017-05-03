#! /bin/bash

# This script generates a Portfile in science/plumed
# Currently the portfile is aimed at testing the currect git hash

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

modules=$(
  grep default-off ../src/*/module.type | sed "s|.*/src/||" | sed "s|/.*||" | awk '{printf("%s ",$1)}'
)


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
}' | awk -v modules="$modules" '{
  if($1=="@_MODULES_@"){
    var=$2
    split(modules,list);
    for(mod in list){
      print "variant mod_"list[mod]" description {Enable "list[mod]" module} {"
      print "  set "var" ${"var"}+"list[mod]
      print "}"
      print ""
    }
    print "variant allmodules description {Enable all optional modules} {"
    for(mod in list){
      print "  set "var" ${"var"}+"list[mod]
    }
    print "}"
    print ""
  } else print
}'  > science/plumed/Portfile




