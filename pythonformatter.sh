#!/bin/bash
#formatted with shfmt -ci -s -bn -i 2
#checked with shellcheck

# this script imitates the behaviour of the src/astile.sh:
# if there is need to reformat the file reformat it and git add it
# it respect the setting in the ruff.toml file stored in the base directory of the repo

# I decided to use ruff since it is quite fast and should not impact much the perfomance of astyle

filelist=$(
  ruff format --check \
    | grep 'Would reformat' \
    | awk '{ printf "%s ", $3 }'
)
IFS=' ' read -r -a files <<<"$filelist"

# making ruff work on all file at the same time is MUCH faster
ruff format -q "${files[@]}"
# then creating a report in style of src/astyle.sh
for file in "${files[@]}"; do
  echo "astyle-ruff $file +++ PATCHED"
  git add "$file"
done
