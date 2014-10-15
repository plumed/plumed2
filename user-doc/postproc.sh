for file in html/*.html
do
awk -v version=$(plumed info --version) '{
  if(match($0,"<span>Main&#160;Page</span>")){
    sub("Main","Home",$0);
    sub("Page","(v"version")",$0);
    print
  } else print
}' $file > $file.tmp
mv $file.tmp $file
done

