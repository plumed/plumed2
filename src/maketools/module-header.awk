# Read only the deliberately restricted, top-level build fields in module.yml.
# Keep the grammar in sync with module_header.py and PlumedModules.cmake.
function fail(message) {
  print FILENAME ": " message > "/dev/stderr"
  bad=1
}
{ sub(/\r$/, "") }
/^["\047]?(type|use)["\047]?[ \t]*:/ {
  if ($0 ~ /^type: (always|default-on|default-off)$/) {
    if (++nt != 1) fail("duplicate type")
    type=$2
  } else if ($0 ~ /^use: \[([a-z][a-z0-9_-]*(, [a-z][a-z0-9_-]*)*)?\]$/) {
    if (++nu != 1) fail("duplicate use")
    use=$0
    sub(/^use: \[/,"",use)
    sub(/\]$/,"",use)
    gsub(/, /," ",use)
  } else fail("expected type: always/default-on/default-off or use: [core, tools]")
}
END {
  if (nt != 1 || nu != 1) fail("exactly one type and use field required")
  if (bad) exit 1
  if (field == "type") print type
  else if (field == "use") print use
  else if (field == "table") print FILENAME, type, use
  else { fail("unknown header field " field); exit 1 }
}
