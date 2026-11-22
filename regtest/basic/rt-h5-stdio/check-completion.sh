# Source this after check.py has created the filename fixtures.
function h5_stdio_complete() {
  eval "$($plumed_nompi completion)"
}

function h5_complete_case() {
  local label=$1 line=$2 expected=$3 candidate
  shift 3
  COMP_WORDS=("${PLUMED_PROGRAM_NAME:-plumed}" "$@")
  COMP_CWORD=$((${#COMP_WORDS[@]}-1))
  COMP_LINE=$line
  COMP_POINT=${#COMP_LINE}
  h5_stdio_complete
  for candidate in "${COMPREPLY[@]}"; do
    if [[ "$candidate" == "$expected" ]]; then
      echo "$label: OK"
      return 0
    fi
  done
  echo "$label: missing $expected"
  return 1
}

h5_complete_case 'h5cat command' 'plumed h5c' h5cat h5c
h5_complete_case 'h5cat option' 'plumed h5cat --h' --help h5cat --h
h5_complete_case 'h5cat first file' 'plumed h5cat fir' first.h5 h5cat fir
h5_complete_case 'h5cat second file' 'plumed h5cat first.h5 sec' second.h5 h5cat first.h5 sec
h5_complete_case 'h5cat quoted filename' "plumed h5cat 'file with sp" 'file with spaces.h5' h5cat 'file with sp'
h5_complete_case 'h5cat after delimiter' 'plumed h5cat -- -l' -leading.h5 h5cat -- -l
h5_complete_case 'h5convert input stream' 'plumed h5convert --input -' - h5convert --input -
h5_complete_case 'h5convert output stream' 'plumed h5convert --input first.h5 --output -' - h5convert --input first.h5 --output -
h5_complete_case 'h5convert input file' 'plumed h5convert --input fir' first.h5 h5convert --input fir
h5_complete_case 'h5convert equals split' 'plumed h5convert --input=fir' first.h5 h5convert --input = fir
h5_complete_case 'h5convert equals unsplit' 'plumed h5convert --input=fir' --input=first.h5 h5convert --input=fir
h5_complete_case 'h5convert equals stream' 'plumed h5convert --output=-' - h5convert --output = -
h5_complete_case 'h5convert equals stream unsplit' 'plumed h5convert --output=-' --output=- h5convert --output=-

h5_complete_case 'h5cat TAR file' 'plumed h5cat first.h5.t' first.h5.tar h5cat first.h5.t
h5_complete_case 'h5convert batch option' 'plumed h5convert --b' --batch-size h5convert --b
h5_complete_case 'h5convert batch count' 'plumed h5convert --batch-size 4' 4096 h5convert --batch-size 4
h5_complete_case 'h5convert equals batch count' 'plumed h5convert --batch-size=4' --batch-size=4096 h5convert --batch-size=4
