# Breadth-first closure, including cycles. Input: name enabled dependencies...
function fail(message) {
  print "module dependencies: " message > "/dev/stderr"
  bad=1
}
{
  known[$1]=1
  if ($2) {
    selected[$1]=1
    queue[++tail]=$1
    path[$1]=$1
  }
  for (i=3; i<=NF; ++i) use[$1]=use[$1] " " $i
}
END {
  split(excluded, words, " ")
  for (i in words) forbidden[words[i]]=1
  for (head=1; head<=tail; ++head) {
    name=queue[head]
    if (forbidden[name]) fail(path[name] " is explicitly disabled")
    n=split(use[name], deps, " ")
    for (i=1; i<=n; ++i) {
      dep=deps[i]
      if (dep == "" || dep == "config") continue
      reason=path[name] " -> " dep
      if (!(dep in known)) {
        fail(reason ": unknown module")
      } else if (forbidden[dep]) {
        fail(reason ": required module " dep " is explicitly disabled")
      } else if (!(dep in selected)) {
        if (mode == "check") {
          fail(reason ": required module " dep " is unavailable after feature checks")
        } else {
          selected[dep]=1
          queue[++tail]=dep
          path[dep]=reason
          additions[++count]=dep " " reason
        }
      }
    }
  }
  if (bad) exit 1
  for (i=1; i<=count; ++i) print additions[i]
}
