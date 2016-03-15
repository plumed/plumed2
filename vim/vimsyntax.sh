#! /bin/bash

# Making sure that plumed executable is available
echo -n "Searching for plumed ..."
if plumed --no-mpi 2>/dev/null 1>/dev/null ; then
  echo " found"
else
  echo " not found"
  echo -n "Sourcing sourceme.sh and searching again ..."
  if source ../sourceme.sh && plumed --no-mpi 2>/dev/null 1>/dev/null ; then
    echo " found"
  else
    echo "ERROR: you should compile plumed first!"
    exit 1
  fi
fi

mkdir -p syntax help

cat > syntax/plumedf.vim << \EOF
syn match plumedFCol2S   /\v\s+/    nextgroup=plumedFCol1 contained
syn match plumedFCol2    /\v\S+/    nextgroup=plumedFCol2S contained
syn match plumedFCol1S   /\v\s+/    nextgroup=plumedFCol2 contained
syn match plumedFCol1    /\v\S+/    nextgroup=plumedFCol1S

syn match plumedNothing  /\v.*/    contained
syn match plumedSCol2    /\v\S+/    nextgroup=plumedNothing contained
syn match plumedSCol1S   /\v\s+/    nextgroup=plumedSCol2 contained
syn match plumedSCol1    /\v\S+/    nextgroup=plumedSCol1S contained

highlight link plumedFCol1 Constant
highlight link plumedFCol2 Keyword
highlight link plumedSCol1 Type
highlight link plumedSCol2 Constant

syntax match   plumedFComment excludenl /\v#.*$/
highlight link plumedFComment Comment

syntax region plumedFieldLine  matchgroup=plumedFieldLine start=/\v[ \t]*\#\![ \t]+FIELDS[ \t]+/ excludenl end=/$/ contains=plumedFCol1
syntax region plumedSetLine    matchgroup=plumedSetLine start=/\v[ \t]*\#\![ \t]+SET[ \t]+/ excludenl end=/$/ contains=plumedSCol1
highlight link plumedFieldLine Type
highlight link plumedSetLine   Type
EOF



actions=$(
plumed --no-mpi manual --action 2>&1 | awk '{
  if(NR==1) next;
  if(NF!=1) exit;
  print $1
}'
)


actions="$(
for a in $actions
do

plumed --no-mpi manual --action $a --vim 2>/dev/null | head -n 1

done
)"

{

cat << \EOF
" Vim syntax file
" Language: PLUMED

if exists("b:current_syntax")
  finish
endif

let b:current_syntax="plumed"

let s:path=expand('<sfile>:p:h') . "/../"

" All except space and hash are in word
set iskeyword=33,34,36-126

" Matching dots, possibly followed by a comment
" Only dots are part of the match
syntax match plumedDots /\v\.\.\.(\s*(#.*)*$)@=/ contained
highlight link plumedDots Type

let b:plumedActions=[]
let b:plumedDictionary={}

EOF
for a in $actions ; do
action_name="${a%%,*}" 
action_name_=$(echo $action_name | sed s/-/_/g)

echo 'call add(b:plumedActions,{"word":"'"$action_name"'"})'

dictionary='{"word":"LABEL=","menu":"add a label"}'



#echo "let b:plumedFullMan[\"$action_name\"]=\"$(
#eval plumed" --no-mpi manual --action $action_name --vim 2>/dev/null | awk '{if(NR>1)printf("%s\\n",$0)}' )"\"

{
echo "****************************************"
echo "Short helpfile for action $action_name"
echo "****************************************"
plumed --no-mpi manual --action $action_name --vim 2>/dev/null | awk '{if(NR>1) print}'
} > help/$action_name.txt


for l in $(echo "$a" | sed 's/,/ /g')
do
  string=
  case "$l" in
  (*:LABEL)
# this is treated differently
    string=
  ;;
  (flag:*)
    dictionary="$dictionary"'
{"word":"'${l#flag:}'","menu":"(flag)"}'
# syntax match   plumedKeywordsDISTANCE "\v<COMPONENTS>" contained
    string='"\v<'${l#flag:}'>"' ;;
  (numbered:*)
    dictionary="$dictionary"'
{"word":"'${l#*:}'","menu":"(numbered)"}'
# syntax match   plumedKeywordsMOVINGRESTRAINT "\v<KAPPA[0-9]*\=[^{ ]*" contained
    string='"\v<'${l#*:}'[0-9]*\=[^{ #]*"' ;;
  (*:*)
    dictionary="$dictionary"'
{"word":"'${l#*:}'="}'
# syntax match   plumedKeywordsMOVINGRESTRAINT "\v<KAPPA[0-9]*\=[^{ ]*" contained
    string='"\v<'${l#*:}'[0-9]*\=[^{ #]*"' ;;
  esac
  test -n "$string" && echo "syntax match   plumedKeywords$action_name_ $string contained contains=plumedStringInKeyword"
done

dictionary="$(
  echo "$dictionary" | sort | tr '\n' ',' | sed 's/,$//'
)"
echo "let b:plumedDictionary[\"plumedLine$action_name\"]=[$dictionary]"

cat << \EOF | sed s/ACTION/$action_name/g | sed s/ACTNAME/$action_name_/g
" single line, with explicit LABEL
" matching action at beginning of line, till the end of the line
" can contain all the keywords associated with this action, plus strings, label, and comments
syntax region plumedLineACTNAME matchgroup=plumedActionACTNAME start=/\v^\s*ACTION>/ excludenl end=/$/ contains=plumedComment,plumedKeywordsACTNAME,plumedLabel,plumedStringOneline fold
" multiple line, with explicit LABEL
" first row might contain extra words before arriving at the dots
" thus continuation dots are matched by plumedDots
" matching action at beginning of line, followed by dots and possibly by a comment
" ends on dots, possibly followed by the same action name and possibly a comment
" comments and initial dots are not part of the match
" can contain all the keywords associated with this action, plus strings, label, and comments
syntax region plumedCLineACTNAME matchgroup=plumedActionACTNAME start=/\v^\s*ACTION>(.+\.\.\.\s*(#.*)*$)@=/ end=/\v^\s*\.\.\.(\s+ACTION)?\s*((#.*)*$)@=/ contains=plumedComment,plumedKeywordsACTNAME,plumedLabel,plumedString,plumedDots fold
" single line, with label: syntax
" matching label followed by action
" can contain all the keywords associated with this action, plus strings and comments
" labels are not allwed
syntax region plumedLLineACTNAME matchgroup=plumedActionACTNAME start=/\v^\s*[^ #@][^ #]*:\s+ACTION/ excludenl end=/$/ contains=plumedComment,plumedKeywordsACTNAME,plumedStringOneline fold
" multiple line, with label: syntax
" first row might contain extra words before arriving at the dots
" thus continuation dots are matched by plumedDots
" matching label, action, dots, and possibly comment
" comments and dots are not part of the match
" ends on dots, possibly followed by the same label and possibly a comment
" comments and initial dots are not part of the match
syntax region plumedLCLineACTNAME matchgroup=plumedActionACTNAME start=/\v^\s*\z([^ #@][^ #]*\:)\s+ACTION>(.+\.\.\.\s*(#.*)*$)@=/ end=/\v^\s*\.\.\.(\s+\z1)?\s*((#.*)*$)@=/ contains=plumedComment,plumedKeywordsACTNAME,plumedString,plumedDots fold
" this is a hack required to match the ACTION when it is in the second line
syntax match plumedSpecialACTNAME /\v(\.\.\.\s*(#.*)*\_s*)@<=ACTION>/ contained
highlight link plumedSpecialACTNAME Type
" multiple line, with label: syntax
" here ACTION is on the second line
" matching label, dots, possibly comments, newline, then action name
" comments, dots, and action are not part of the match
" ends on dots possibly followed by the same label and possibly a comment
syntax region plumedLCLineACTNAME matchgroup=plumedActionACTNAME start=/\v^\s*\z([^ #@][^ #]*\:)\s+(\.\.\.\s*(#.*)*\_s*ACTION)@=/ end=/\v^\s*\.\.\.(\s+\z1)?\s*((#.*)*$)@=/ contains=plumedComment,plumedKeywordsACTNAME,plumedString,plumedSpecialACTNAME,plumedDots fold
highlight link plumedActionACTNAME Type
highlight link plumedKeywordsACTNAME Statement
EOF

done
cat << \EOF
" comments and strings last, with highest priority
syntax region  plumedString start=/\v\{/  end=/\v\}/ contained contains=plumedString fold
syntax region  plumedStringOneline start=/\v\{/  end=/\v\}/ oneline contained contains=plumedStringOneline fold
highlight link plumedString String
highlight link plumedStringOneline String
syntax match   plumedStringInKeyword /\v(<[^ #]+\=)@<=[^ #]+/ contained
highlight link plumedStringInKeyword String

" Matching label
syntax match   plumedLabel "\v<LABEL\=[^ #]*" contained contains=plumedLabelWrong
highlight link plumedLabel Type

" Errors
syntax match   plumedLabelWrong "\v<LABEL\=\@[^ #]*" contained
highlight link plumedLabelWrong Error

syntax region  plumedComment start="\v^\s*ENDPLUMED>" end="\%$" fold
syntax match   plumedComment excludenl "\v#.*$"
highlight link plumedComment Comment


fun! PlumedGuessRegion()
" this is to find the match
" first, sync syntax
            syn sync fromstart
" find the syntactic attribute of the present region
            let col=col(".")
            let line=line(".")
            let key=""
            let stack=synstack(line,col)
            if(len(stack)>0)
              let key = synIDattr(stack[0], "name")
            endif
            if(key=~"^plumed[LC]*Line.*")
              return substitute(key,"^plumed[LC]*Line","","")
            endif
            return ""
endfunction

fun! PlumedContextManual()
  if(exists("b:plumed_helpfile"))
    quit
    return
  endif
  let m=PlumedGuessRegion()
  if(m=="")
    return
  else
    let name=s:path . "/help/" . m . ".txt"
    if(exists("b:plumed_helpfile_vertical"))
      execute 'rightbelow vsplit | view ' name
    else
      execute 'rightbelow split | view ' name
    endif
    let b:plumed_helpfile=1
  endif
endfunction

fun! PlumedManualV()
  let b:plumed_helpfile_vertical=1
endfunction

fun! PlumedManualH()
  unlet b:plumed_helpfile_vertical
endfunction

command! -nargs=0 PHelp call PlumedContextManual()

" autocomplete function
fun! PlumedComplete(findstart, base)
" this is to find the start of the word to be completed
          if a:findstart
            " locate the start of the word
            let line = getline('.')
            let start = col('.') - 1
            while start > 0 && line[start - 1] =~ '[a-zA-Z\_]'
              let start -= 1
            endwhile
            return start
          else
" this is to find the match
" first, sync syntax
            syn sync fromstart
" find the syntactic attribute of the present region
            let col=col(".")
            let line=line(".")
            let key=""
            if col!=1
              let key = synIDattr(synID(line, col-1, 1), "name")
            else
              let stack=synstack(line,col)
              if(len(stack)>0)
                let key = synIDattr(stack[0], "name")
              endif
            endif
            let comp=[]
" normalize key removing L/C/LC indication
            let key1=substitute(key,"^plumed[LC]*Line","plumedLine","")
            if key ==""
" if outside of any region, complete with list of actions
              let comp=b:plumedActions
            elseif has_key(b:plumedDictionary,key1)
" if inside a region in the form "plumedLineXXX"
" complete with keywords associated to action XXX
              let comp=b:plumedDictionary[key1]
            endif
            " find months matching with "a:base"
            let res = []
            for m in comp
" this is to allow m to be a dictionary
" with a word and a one-liner
              if(type(m)==type({}))
                let n=m["word"]
              else
                let n=m
              endif
              if n =~ '^' . a:base
                if(n!="LABEL=" || key =~ "^plumedLine.*" || key =~ "^plumedCLine.*")
                call add(res, m)
                endif
              endif
" in principle comp could be a heterogeneous list
" so it should be unlet to iterate the loop
              unlet m
            endfor
            return res
          endif
        endfun
setlocal omnifunc=PlumedComplete

" inspect the entire file to find lines containing
" non highlighted characters
fun! PlumedAnnotateSyntax()
" buffer where errors are written
  let buffer=[]
  let l=1
" loop over lines
  while l <= line("$")
    let line=getline(l)
    let p=0
" line is assumed right a priori
    let wrong=0
" position the cursor and redraw the screen
    call cursor(l,1)
    redraw! "! is required for some reason
    while p <len(line)
      let stack=synstack(l,p+1)
      if line[p] !~ "[ \t]"
        if(len(stack)==0)
          let wrong=1
        elseif(synIDattr(stack[len(stack)-1],"name")=~"^plumed[LC]*Line.*")
          let wrong=1
        endif
      endif
      let annotation=""
      for s in stack
        let annotation=annotation."+".synIDattr(s,"name")
      endfor
      call add(buffer,printf("ANNOTATION %5d %3d %s %s",l,p,line[p],annotation))
      let p=p+1
    endwhile
    
    if(wrong)
      call add(buffer,"ERROR AT LINE ".l." : ".line)
    endif
    let l=l+1
  endwhile
" dump the buffer on a new window
  new
  for l in buffer
    put=l
  endfor
endfun

EOF

} > syntax/plumed.vim


# colors:
# Constant
# Identifier
# Statement
# PreProc
# Type
# Special
# Underlined
# Ignore
# Error
# Todo



