# Runtime prerequisite queries shared by the regression runner and its tests.
# Probe command support separately: a false feature result must never trigger
# fallback to the kernel-wide build flags.
plumed_module_help=$(PLUMED_AUTOLOAD=0 $plumed_nompi modules --help 2>/dev/null) || plumed_module_help=
plumed_has_modules=no
plumed_has_scoped_modules=no
if [[ $plumed_module_help == *'plumed modules'* ]]; then
  plumed_has_modules=yes
  if [[ $plumed_module_help == *'--module '* ]]; then
    plumed_has_scoped_modules=yes
  fi
fi

plumed_check_module() {
  if [[ $plumed_has_modules == yes ]]; then
    $plumed_nompi modules -q --has "$1"
  else
    $plumed_nompi config module "$1" -q
  fi
}

plumed_check_feature() {
  local requirement=$1 module feature
  if [[ $requirement == *:* ]]; then
    module=${requirement%%:*}
    feature=${requirement#*:}
    [[ -n $module && -n $feature && $feature != *:* ]] || return 1
    if [[ $plumed_has_scoped_modules == yes ]]; then
      $plumed_nompi modules --module "$module" --has-feature "$feature" -q
    else
      # Legacy monolithic builds need both the module and the global feature.
      $plumed_nompi config module "$module" -q &&
        $plumed_nompi config has "$feature" -q
    fi
  elif [[ $plumed_has_modules == yes ]]; then
    $plumed_nompi modules -q --has-feature "$requirement"
  else
    $plumed_nompi config has "$requirement" -q
  fi
}

# A device available at runtime, unlike the features above, that are compiled in:
# "plumed <device>-devices -q" (e.g. cuda-devices, provided by the gpu_cuda_hip
# module) returns 0 when it sees a device. The query runs in its own process, so
# it leaves the devices of the test alone; a missing tool means no device.
plumed_check_device() {
  [[ $1 =~ ^[a-z][a-z0-9_]*$ ]] || return 1
  $plumed_nompi "$1-devices" -q >/dev/null 2>&1
}
