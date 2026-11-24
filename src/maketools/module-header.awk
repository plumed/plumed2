# Read only the deliberately restricted, top-level build fields in module.yml.
# Keep the grammar in sync with module_header.py and PlumedModules.cmake.
# type and use are required; toolchain (default cxx) and placement (default any)
# are optional.
function fail(message) {
  print FILENAME ": " message > "/dev/stderr"
  bad=1
}
BEGIN { toolchain="cxx"; placement="any" }
{ sub(/\r$/, "") }
/^["\047]?(type|use|toolchain|placement)["\047]?[ \t]*:/ {
  if ($0 ~ /^type: (always|default-on|default-off)$/) {
    if (++nt != 1) fail("duplicate type")
    type=$2
  } else if ($0 ~ /^use: \[([a-z][a-z0-9_-]*(, [a-z][a-z0-9_-]*)*)?\]$/) {
    if (++nu != 1) fail("duplicate use")
    use=$0
    sub(/^use: \[/,"",use)
    sub(/\]$/,"",use)
    gsub(/, /," ",use)
  } else if ($0 ~ /^toolchain: (cxx|cuda)$/) {
    if (++nc != 1) fail("duplicate toolchain")
    toolchain=$2
  } else if ($0 ~ /^placement: (any|extension)$/) {
    if (++np != 1) fail("duplicate placement")
    placement=$2
  } else fail("expected type: always/default-on/default-off, use: [core, tools], toolchain: cxx/cuda or placement: any/extension")
}
END {
  if (nt != 1 || nu != 1) fail("exactly one type and use field required")
  if (bad) exit 1
  if (field == "type") print type
  else if (field == "use") print use
  else if (field == "toolchain") print toolchain
  else if (field == "placement") print placement
  else if (field == "table") print FILENAME, type, use, toolchain, placement
  else { fail("unknown header field " field); exit 1 }
}
