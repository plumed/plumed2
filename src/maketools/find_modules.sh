#!/bin/bash
# sort ../basic ../tools $(shell for d in $(MODULES) ; do test -f ../$$d.off || echo ../$$d ; done ; for dir in ../* ; do test -d $$dir && test -e $$dir.on && echo $$dir ; done

cd ../
for dir in *
do
  
  if test -f $dir/module.type
  then
    case "$(cat "$dir/module.type")" in
    (always) echo $dir ;;
    (default-on) test -f $dir.off || echo $dir ;;
    (default-off) test -f $dir.on && echo $dir ;;
    esac
  fi
done
