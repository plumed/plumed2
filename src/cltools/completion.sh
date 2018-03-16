    local cur prev opts
    shopt -s extglob
    COMPREPLY=()

    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"
# these options can be combined with other commands
    opts="--load"
# this option should be the first
    if((COMP_CWORD==1)) ; then
      opts="$opts --no-mpi --mpi"
    fi

    local cmd_found i cmd_test

# check if one of the previous keywords is a command
    for((i=1;i<COMP_CWORD;i++)); do
      cmd_found=""
      for cmd_test in $cmds ; do
        if [[ "$cmd_test" == "${COMP_WORDS[i]}" ]] ; then
          eval "local comp=\"\$cmd_keys_${cmd_test//-/_}\""
          case "$cur" in 
            (-*) COMPREPLY=( $(compgen -W "$comp" -- $cur ) ) ;;
            (*)  COMPREPLY=( $(compgen -o bashdefault -- ${cur}) ) ;;
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

# complete with options or commands
    case "${cur}" in
# show options only if completing a "-" 
    (-*) COMPREPLY=( $(compgen -W "${opts} ${cmds}" -- ${cur}) ) ;;
    (*)  COMPREPLY=( $(compgen -W "${cmds}" -- ${cur}) ) ;;
    esac
    return 0
