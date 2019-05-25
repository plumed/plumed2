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
mkdir -p python/py-plumed

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
  } else if($1=="@_PYTHON_@"){
    ver[0]="26"
    ver[1]="27"
    ver[2]="33"
    ver[3]="34"
    ver[4]="35"
    ver[5]="36"
    nver=6;
    for(i=0;i<nver;i++){
# number with dot (e.g. 2.7):
      verdot= substr(ver[i],1,1)"."substr(ver[i],2,1)
      printf("%s", "variant python" ver[i] " description {Bindings for python" ver[i] "} ")
      printf("%s", "conflicts");
      for(j=0;j<nver;j++) if(i!=j) printf("%s"," python" ver[j]);
      print " {"
      print "  depends_lib-append port:python" ver[i];
      print "  depends_lib-append port:py" ver[i] "-numpy";
      print "  configure.args-replace --disable-python --enable-python";
      print "  configure.args-append PYTHON_BIN=${prefix}/bin/python" verdot
      print "}"
      print ""
    }

  } else print
}'  > science/plumed/Portfile


cat PortfilePython.in |
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
}' > python/py-plumed/Portfile





