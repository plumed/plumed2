#! /bin/bash

# use this script to make a new release
# just type ./release.sh
# and follow instructions

confirm () {
    # call with a prompt string or use a default
    read -r -p "${1:-Are you sure? [y/N]} " response
    case $response in
        [yY][eE][sS]|[yY]) 
            true
            ;;
        *)
            false
            ;;
    esac
}

VALIDATED=
if [ "$1" = --validated ] ; then
  VALIDATED=yes
elif test -n "$1" ; then
  echo ERROR
  echo "cannot understand $1 option"
  exit 1
fi


echo "*** Only use this script if you want to release a new PLUMED version ***"
echo "*** follow instructions below ***"

ls src README 1>/dev/null 2>/dev/null ||  {
  echo "Launch from root directory"
  exit 1
}

read -r -p "Type the version number (e.g. 2.1.3 or 2.2b): " version

case "$version" in
(2.?*)
  version="${version#2.}"
  ;;
(*)
  echo "ERROR"
  echo "Use a version in form 2.*"
  exit 1
  ;;
esac

shortversion=$(echo "$version" | sed  's/^\([0-9][0-9]*\).*/\1/' )

version=2.$version
shortversion=2.$shortversion

echo "major version v$shortversion"
echo "minor version v$version"

if test $version == $shortversion ; then
  echo "ERROR"
  echo "please add a patch identified so that tag and branch have different names"
  exit 1
fi

if test "$(git tag -l | awk '{if($1=="'v$version'")print "yes"}')" == yes ; then
  echo "tag v$version already exists"
  exit 1
fi

set -e
if ! test "$VALIDATED" ; then
  echo "checking out branch v$shortversion"
  git checkout v$shortversion
  echo 
  msg="Travis tests for v$version

[makedoc]"
  echo "Now I will add an empty commit and push the result to origin"
  echo "I will use the following commands:"
  echo "***"
  echo "git commit --allow-empty -m \"$msg\""
  echo "git push origin v$shortversion"
  echo "***"
  confirm || exit
  git commit --allow-empty -m "$msg"
  git push origin v$shortversion
  echo
  echo "Now you should go at this link:"
  echo "  http://travis-ci.org/plumed/plumed2/builds"
  echo "and wait for travis to finish the tests"
  echo "In case of failures, fix and repeat the procedure"
  echo "Also check the online manual, that should be here:"
  echo "  http://plumed.github.io/doc-v$shortversion"
  echo "In case of success, relaunch this script as \"./release.sh --validated\""
else
  {
    grep  \# VERSION 
    echo $version
  } > VERSION.tmp
  mv VERSION.tmp VERSION
  echo "Here's the new VERSION file:"
  echo "***"
  cat VERSION
  echo "***"
  msg="Release v$version

[makedoc]"
  echo "Now I will add it, prepare a release commit, add a tag named v$version"
  echo "push it to origin and create a tgz file"
  echo "I will use the following commands:"
  echo "***"
  echo "git add VERSION"
  echo "git commit --allow-empty -m \"$msg\""
  echo "git tag v$version"
  echo "git push origin v$version"
  echo "git archive -o plumed-$version.tgz --prefix plumed-$version/ v$version"
  echo "***"
  confirm || exit
  git add VERSION
  git commit --allow-empty -m "$msg"
  git tag v$version
  git push origin v$shortversion v$version
  git archive -o plumed-$version.tgz --prefix plumed-$version/ v$version
  echo
  echo "Done!"
  echo
  echo "Please upload the file plumed-$version.tgz on the download directory"
  echo "Remember to notify the mailing list"
fi







