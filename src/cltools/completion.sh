    local -a plumed_words=() plumed_load_options=() plumed_candidates=()
    local plumed_executable=${COMP_WORDS[0]} plumed_cur plumed_prev=""
    local plumed_word="" plumed_quote="" plumed_char plumed_next plumed_line
    local plumed_index plumed_command="" plumed_value="" plumed_kind=commands
    # Bash passes the actual Readline replacement fragment as $2. In Bash 3
    # COMP_WORDS can still contain a larger word such as key=value.
    local plumed_trim="" plumed_original=${2-${COMP_WORDS[COMP_CWORD]}}
    local plumed_command_index=0 plumed_file_prefix="" plumed_context_prefix=""
    local plumed_output plumed_header=yes plumed_start plumed_duplicate
    local LC_ALL=C
    COMPREPLY=()

    # Readline splits unquoted '=' and ':' inside words. Tokenize only the text
    # before the cursor, without eval or expansion, preserving quoted spaces.
    if [[ -n ${COMP_LINE+x} ]]; then
      plumed_line=${COMP_LINE:0:COMP_POINT}
      for ((plumed_index=0; plumed_index<${#plumed_line}; plumed_index++)); do
        plumed_char=${plumed_line:plumed_index:1}
        plumed_next=${plumed_line:plumed_index+1:1}
        if [[ "$plumed_quote" == "'" ]]; then
          if [[ "$plumed_char" == "'" ]]; then plumed_quote=""; else plumed_word+=$plumed_char; fi
        elif [[ "$plumed_char" == \\ && "$plumed_quote" != "'" ]]; then
          if [[ -z "$plumed_quote" || "$plumed_next" == [\$\`\"\\] ]]; then
            plumed_word+=$plumed_next
            ((plumed_index+=1))
          else
            plumed_word+=$plumed_char
          fi
        elif [[ "$plumed_quote" == '"' ]]; then
          if [[ "$plumed_char" == '"' ]]; then plumed_quote=""; else plumed_word+=$plumed_char; fi
        else
          case "$plumed_char" in
            "'"|'"') plumed_quote=$plumed_char ;;
            ' '|$'\t'|$'\n')
              if [[ -n "$plumed_word" ]]; then plumed_words+=("$plumed_word"); plumed_word=""; fi ;;
            ';'|'|'|'&') plumed_words=(); plumed_word="" ;;
            *) plumed_word+=$plumed_char ;;
          esac
        fi
      done
      plumed_words+=("$plumed_word")
    else
      # Also support completion callers which provide COMP_WORDS alone.
      for ((plumed_index=0; plumed_index<=COMP_CWORD; plumed_index++)); do
        plumed_word=${COMP_WORDS[plumed_index]}
        if ((${#plumed_words[@]})) && [[ "$plumed_word" == = || "$plumed_word" == : || ${plumed_words[${#plumed_words[@]}-1]} == *= || ${plumed_words[${#plumed_words[@]}-1]} == *: ]]; then
          plumed_words[${#plumed_words[@]}-1]+=$plumed_word
        else
          plumed_words+=("$plumed_word")
        fi
      done
    fi
    plumed_cur=${plumed_words[${#plumed_words[@]}-1]}
    for ((plumed_index=1; plumed_index<${#plumed_words[@]}-1; plumed_index++)); do
      plumed_word=${plumed_words[plumed_index]}
      if [[ -z "$plumed_command" ]]; then
        case "$plumed_word" in
          --load|--load-module)
            if ((plumed_index+1==${#plumed_words[@]}-1)); then
              plumed_prev=$plumed_word; break
            fi
            ((plumed_index+=1))
            plumed_load_options+=("$plumed_word" "${plumed_words[plumed_index]}") ;;
          --load=*|--load-module=*) plumed_load_options+=("$plumed_word") ;;
          --tmpi-threads) ((plumed_index+=1)) ;;
          --*) ;;
          *) plumed_command=$plumed_word; plumed_command_index=$plumed_index ;;
        esac
      else
        plumed_prev=$plumed_word
      fi
    done
    if [[ -z "$plumed_command" ]]; then
      case "$plumed_cur" in
        --load=*|--load-module=*) plumed_prev=${plumed_cur%%=*}; plumed_cur=${plumed_cur#*=} ;;
      esac
      case "$plumed_prev" in
        --load) plumed_kind=libraries ;;
        --load-module) plumed_kind=installed-modules ;;
        *) if [[ "$plumed_cur" == -* ]]; then plumed_kind=global-options; fi ;;
      esac
    else
      # Tools own argument semantics. Pass the complete tool-local context,
      # including the unfinished word, without interpreting it in the shell.
      plumed_start=$((plumed_command_index+1))
      plumed_output=$("$plumed_executable" --no-mpi "${plumed_load_options[@]}" completion --list context --command "$plumed_command" -- "${plumed_words[@]:plumed_start}" 2>/dev/null) || return 0
      while IFS= read -r plumed_value; do
        if [[ "$plumed_header" == yes ]]; then
          plumed_kind=${plumed_value%%$'\t'*}
          plumed_context_prefix=${plumed_value#*$'\t'}
          case "$plumed_kind" in files|words) ;; *) return 0 ;; esac
          plumed_header=no
        else
          plumed_candidates+=("$plumed_value")
        fi
      done <<< "$plumed_output"
      if [[ -n "$plumed_context_prefix" && "$plumed_cur" == "$plumed_context_prefix"* ]]; then
        if [[ "$plumed_original" == "$plumed_context_prefix"* ]]; then
          plumed_file_prefix=$plumed_context_prefix
        fi
        plumed_cur=${plumed_cur#"$plumed_context_prefix"}
      fi
    fi

    # Return only the Readline fragment being replaced, for literal tool values
    # as well as filenames (e.g. the value after '=' or ':').
    if [[ "$plumed_cur" == *"$plumed_original" ]]; then
      plumed_trim=${plumed_cur:0:${#plumed_cur}-${#plumed_original}}
      if [[ -z ${2+x} && ( -z "$plumed_trim" || "${COMP_WORDBREAKS-=:}" != *"${plumed_trim: -1}"* ) ]]; then
        plumed_trim=""
      fi
    fi
    for plumed_value in "${plumed_candidates[@]}"; do
      if [[ -n "$plumed_value" && "$plumed_value" == "$plumed_cur"* ]]; then
        COMPREPLY+=("${plumed_file_prefix}${plumed_value#"$plumed_trim"}")
      fi
    done
    if [[ "$plumed_kind" == words ]]; then return 0; fi

    if [[ "$plumed_kind" == libraries || "$plumed_kind" == files ]]; then
      # Mark filenames for Readline quoting without changing shell options.
      if type compopt &>/dev/null; then
        compopt -o filenames 2>/dev/null || :
      else
        # Bash 3 has no per-invocation compopt. Upgrade an old registration
        # while preserving its flags; its default filename fallback handles
        # this first invocation, and later invocations use the new setting.
        local plumed_registration
        plumed_registration=$(complete -p -- "$plumed_executable" 2>/dev/null) || plumed_registration=""
        if [[ -n "$plumed_registration" && "$plumed_registration" != *'-o filenames'* ]]; then
          eval "${plumed_registration/complete /complete -o filenames }"
          return 0
        fi
      fi
      while IFS= read -r plumed_value; do
        if [[ -d "$plumed_value" ]]; then
          plumed_value+=/
        elif [[ "$plumed_kind" == libraries ]]; then
          case "$plumed_value" in *.so|*.so.*|*.dylib|*.dll) ;; *) continue ;; esac
        fi
        plumed_duplicate=no
        for plumed_word in "${plumed_candidates[@]}"; do
          if [[ "$plumed_word" == "$plumed_value" ]]; then plumed_duplicate=yes; break; fi
        done
        if [[ "$plumed_duplicate" == yes ]]; then continue; fi
        COMPREPLY+=("${plumed_file_prefix}${plumed_value#"$plumed_trim"}")
      done < <(compgen -f -- "$plumed_cur")
      return 0
    fi
    plumed_output=$("$plumed_executable" --no-mpi "${plumed_load_options[@]}" completion --list "$plumed_kind" --command "$plumed_command" 2>/dev/null) || return 0
    while IFS= read -r plumed_value; do
      if [[ -n "$plumed_value" && "$plumed_value" == "$plumed_cur"* ]]; then COMPREPLY+=("${plumed_value#"$plumed_trim"}"); fi
    done <<< "$plumed_output"
    return 0
