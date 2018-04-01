    local cur prev opts cmd_found i cmd_test comp1 comp2 l
    shopt -s extglob
    COMPREPLY=()

    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"
# these options can be combined with other commands
    opts="--load --no-mpi --mpi"

# check if one of the previous keywords is a command
    for((i=1;i<COMP_CWORD;i++)); do
      cmd_found=""
      for cmd_test in $cmds ; do
        if [[ "$cmd_test" == "${COMP_WORDS[i]}" ]] ; then
          eval "local comp=\"\$cmd_keys_${cmd_test//-/_}\""
          comp1=""
          comp2=""
          for l in $comp ; do
            case $l in
            (-*) comp1="$comp1 $l" ;;
            (*) comp2="$comp2 $l" ;;
            esac
          done
          case "$cur" in 
            (-*) COMPREPLY=( $(compgen -W "$comp1" -- $cur ) ) ;;
            (*)  COMPREPLY=( $(compgen -W "$comp2" -o bashdefault -- ${cur}) ) ;;
          esac
          return 0
        fi
      done
      if [[ "$cmd_found" == 1 ]]; then
        break
      fi
    done

# if previous is --load, autocomplete with dynamic library
    if [[ "${prev}" == --load ]] ; then
      COMPREPLY=( $(compgen -X '!*.@(dylib|so)' -- $cur ) )
      return 0
    fi

    comp1=""
    comp2=""
    for l in $opts $cmds ; do
      case $l in
      (-*) comp1="$comp1 $l" ;;
      (*) comp2="$comp2 $l" ;;
      esac
    done

# complete with options or commands
    case "${cur}" in
# show options only if completing a "-" 
    (-*) COMPREPLY=( $(compgen -W "$comp1" -- ${cur}) ) ;;
    (*)  COMPREPLY=( $(compgen -W "$comp2" -- ${cur}) ) ;;
    esac
    return 0
