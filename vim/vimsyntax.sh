#! /bin/bash

source ../sourceme.sh

mkdir -p syntax help

cat > syntax/plumedf.vim << \EOF

if exists("b:current_syntax")
  finish
endif

let b:current_syntax="plumedf"

if exists("g:plumed_shortcuts")
  noremap <buffer> <F3> :PMinus<CR>
  noremap <buffer> <F4> :PPlus<CR>
endif

command! -nargs=0 -count=1 PPlus call PlumedMove(<count>)
command! -nargs=0 -count=1 PMinus call PlumedMove(-<count>)
command! -nargs=0 -count=0 PCol call PlumedColumn(<count>)

" move highlight by count columns
" count can be negative
function! PlumedMove(count)
  if !exists("b:plumedf_column")
    let b:plumedf_column=0
  endif
  let b:plumedf_column=b:plumedf_column+a:count
  if(b:plumedf_column<0)
    let b:plumedf_column=0
  endif
  call PlumedColumn(b:plumedf_column)
endfunction

" highlight column col
function! PlumedColumn(col)
 let b:plumedf_column=a:col
 let s=""
 let c=0
 while c<=a:col+2
   let c=c+1
   if(c==a:col)
     let s=s.'X'
   else
     let s=s.'.'
   endif
 endwhile
 call PlumedColumnPattern(s)
endfunction

" show a given highlight pattern
" I keep this separate to allow future implementation
" of multiple column highlight
function! PlumedColumnPattern(col)
 syntax clear
 let coll=split(a:col,'\zs') " splits in characters
 let cmax=len(coll)
 let c=cmax
 while c >0
  let cc=c+1
  if(cc>cmax)
    let cc=cmax-1
  endif
  execute 'syn match plumedFCol' . c . 'S /\v\s+/             nextgroup=plumedFCol' . cc . ' contained'
  let contained='contained'
  if(c==1)
     let contained=''
  endif
  execute 'syn match plumedFCol'  . c . ' /\v\S+/             nextgroup=plumedFCol' . c  . 'S ' . contained
  if(coll[c-1]=='B' || (coll[c-1]=='.' && (c-1)%2!=0) )
    let type="Keyword"
  elseif(coll[c-1]=='A' || (coll[c-1]=='.' && (c-1)%2==0) )
    let type="Constant"
  elseif(coll[c-1]=='X')
    let type="Todo"
  endif
  execute 'highlight link plumedFCol' .c. ' ' . type
  let c -= 1
 endwhile

 syn match plumedNothing  /\v.*/    contained
 syn match plumedSCol2    /\v\S+/    nextgroup=plumedNothing contained
 syn match plumedSCol1S   /\v\s+/    nextgroup=plumedSCol2 contained
 syn match plumedSCol1    /\v\S+/    nextgroup=plumedSCol1S contained
 highlight link plumedSCol1 Type
 highlight link plumedSCol2 Constant
 syntax match   plumedFComment excludenl /\v#.*$/
 highlight link plumedFComment Comment
 syntax region plumedFieldLine  matchgroup=plumedFieldLine start=/\v[ \t]*\#\![ \t]+FIELDS[ \t]+/ excludenl end=/$/ contains=plumedFCol1
 syntax region plumedSetLine    matchgroup=plumedSetLine start=/\v[ \t]*\#\![ \t]+SET[ \t]+/ excludenl end=/$/ contains=plumedSCol1
 highlight link plumedFieldLine Type
 highlight link plumedSetLine   Type
endfunction

" initialize to column zero
call PlumedColumn(0)

EOF



actions=$(
../src/lib/plumed --no-mpi manual --action 2>&1 | awk '{
  if(NR==1) next;
  if(NF!=1) exit;
  print $1
}'
)


actions="$(
for a in $actions
do

../src/lib/plumed --no-mpi manual --action $a --vim 2>/dev/null | awk -v a=$a 'BEGIN{
  help="help/" a ".txt"
  print "****************************************" > help
  print "Short helpfile for action " a > help
  print "****************************************" > help
}{
  if(NR==1){ print}
  else print > help
}'

done
)"

{

cat << \EOF
" Vim syntax file
" Language: PLUMED

if exists("b:current_syntax")
  finish
endif

syntax case ignore

let b:current_syntax="plumed"

if exists("g:plumed_shortcuts")
  noremap  <buffer> <F2> :PHelp<CR>
  inoremap <buffer> <F2> <Esc>:PHelp<CR>
endif

" path for plumed plugin
let s:path=expand('<sfile>:p:h:h')

" All except space and hash are in word
setlocal iskeyword=33,34,36-126

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

dictionary='{"word":"LABEL=","menu":"(label)"}'

for l in $(echo "$a" | sed 's/,/ /g')
do
  string=
  case "$l" in
  (*:LABEL)
# this is treated differently
  ;;
  (flag:*)
    dictionary="$dictionary"'
{"word":"'${l#flag:}'","menu":"(flag)"}'
  ;;
  (numbered:*)
    dictionary="$dictionary"'
{"word":"'${l#*:}'","menu":"(numbered)"}'
  ;;
  (*:*)
    dictionary="$dictionary"'
{"word":"'${l#*:}'=","menu":"(option)"}'
  ;;
  esac
done

dictionary="$(
  echo "$dictionary" | sort | tr '\n' ',' | sed 's/,$//'
)"
echo "let b:plumedDictionary[\"$action_name\"]=[$dictionary]"

done

cat << \EOF
function! PlumedDefineSyntax()

  for key in sort(keys(b:plumedDictionary))
    call add(b:plumedActions,{"word":key})
  endfor

for a in b:plumedActions
  let action=a["word"]
" vim variables cannot contain -
" we convert it to triple ___
  let action_=substitute(action,"-","___","g")

  for b in b:plumedDictionary[action]
    if(b["menu"]=="(numbered)")
      let string='"\v<' . b["word"] . '[0-9]*\=[^{ #]*"'
    elseif(b["menu"]=="(option)")
" this is necessary since word for option is e.g ."RESTART="
      let string='"\v<' . substitute(b["word"],"=","","") . '[0-9]*\=[^{ #]*"'
    elseif(b["menu"]=="(flag)")
      let string='"\v<' . b["word"] . '>"'
    endif
    execute 'syntax match   plumedKeywords' . action_ . ' ' . string . ' contained contains=plumedStringInKeyword,plumedFillTodo'
  endfor

" single line, with explicit LABEL
" matching action at beginning of line, till the end of the line
" can contain all the keywords associated with this action, plus strings, label, and comments
execute 'syntax region plumedLine' . action_ . ' matchgroup=plumedAction' . action_ . ' start=/\v^\s*' . action . '>/ excludenl end=/$/ contains=plumedComment,plumedKeywords' . action_ . ',plumedLabel,plumedStringOneline,plumedFillTodo fold'
" multiple line, with explicit LABEL
" first row might contain extra words before arriving at the dots
" thus continuation dots are matched by plumedDots
" matching action at beginning of line, followed by dots and possibly by a comment
" ends on dots, possibly followed by the same action name and possibly a comment
" comments and initial dots are not part of the match
" can contain all the keywords associated with this action, plus strings, label, and comments
execute 'syntax region plumedCLine' . action_ . ' matchgroup=plumedAction' . action_ . ' start=/\v^\s*' . action . '>(.+\.\.\.\s*(#.*)*$)@=/ end=/\v^\s*\.\.\.(\s+' . action . ')?\s*((#.*)*$)@=/ contains=plumedComment,plumedKeywords' . action_ . ',plumedLabel,plumedString,plumedDots,plumedFillTodo fold'
" single line, with label: syntax
" matching label followed by action
" can contain all the keywords associated with this action, plus strings and comments
" labels are not allwed
execute 'syntax region plumedLLine' . action_ . ' matchgroup=plumedAction' . action_ . ' start=/\v^\s*[^ #@][^ #]*:\s+' . action . '/ excludenl end=/$/ contains=plumedComment,plumedKeywords' . action_ . ',plumedStringOneline,plumedFillTodo fold'
" multiple line, with label: syntax
" first row might contain extra words before arriving at the dots
" thus continuation dots are matched by plumedDots
" matching label, action, dots, and possibly comment
" comments and dots are not part of the match
" ends on dots, possibly followed by the same label and possibly a comment
" comments and initial dots are not part of the match
execute 'syntax region plumedLCLine' . action_ . ' matchgroup=plumedAction' . action_ . ' start=/\v^\s*\z([^ #@][^ #]*\:)\s+' . action . '>(.+\.\.\.\s*(#.*)*$)@=/ end=/\v^\s*\.\.\.(\s+\z1)?\s*((#.*)*$)@=/ contains=plumedComment,plumedKeywords' . action_ . ',plumedString,plumedDots,plumedFillTodo fold'
" this is a hack required to match the ACTION when it is in the second line
execute 'syntax match plumedSpecial' . action_ . ' /\v(\.\.\.\s*(#.*)*\_s*)@<=' . action . '>/ contained'
execute 'highlight link plumedSpecial' . action_ . ' Type'
" multiple line, with label: syntax
" here ACTION is on the second line
" matching label, dots, possibly comments, newline, then action name
" comments, dots, and action are not part of the match
" ends on dots possibly followed by the same label and possibly a comment
execute 'syntax region plumedLCLine' . action_ . ' matchgroup=plumedAction' . action_ . ' start=/\v^\s*\z([^ #@][^ #]*\:)\s+(\.\.\.\s*(#.*)*\_s*' . action . ')@=/ end=/\v^\s*\.\.\.(\s+\z1)?\s*((#.*)*$)@=/ contains=plumedComment,plumedKeywords' . action_ . ',plumedString,plumedSpecial' . action_ . ',plumedDots,plumedFillTodo fold'
execute 'highlight link plumedAction' . action_ . ' Type'
execute 'highlight link plumedKeywords' . action_ . ' Statement'
endfor

" comments and strings last, with highest priority
syntax region  plumedString start=/\v\{/  end=/\v\}/ contained contains=ALL fold
syntax region  plumedStringOneline start=/\v\{/  end=/\v\}/ oneline contained contains=ALL fold
highlight link plumedString String
highlight link plumedStringOneline String
syntax match   plumedStringInKeyword /\v(<[^ #]+\=)@<=[^ #]+/ contained
highlight link plumedStringInKeyword String

" Matching label
syntax match   plumedLabel "\v<LABEL\=[^ #]*" contained contains=plumedLabelWrong,plumedFillTodo
highlight link plumedLabel Type

" Errors
syntax match   plumedLabelWrong "\v<LABEL\=\@[^ #]*" contained
highlight link plumedLabelWrong Error

" Todo
" This is used since the Trieste tutorials to indicate
" parts that should be filled by the user
syntax match   plumedFillTodo "__FILL__"
highlight link plumedFillTodo Todo

syntax region  plumedComment start="\v^\s*ENDPLUMED>" end="\%$" fold
syntax match   plumedComment excludenl "\v#.*$"
highlight link plumedComment Comment
endfunction

call PlumedDefineSyntax()


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
" this is to allow closing the window with F2
    if exists("g:plumed_shortcuts")
      noremap  <buffer> <F2> :PHelp<CR>
    endif
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
            while start > 0 && line[start - 1] =~ '[a-zA-Z\_\=\-\.]'
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
" retrieve action name
" normalize ___ to -
            let key1=substitute(substitute(key,"^plumed[LC]*Line","",""),"___","-","g")
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
"           if("..." =~ '^' . a:base && (key=~"^plumedLLine.*" || key=~"^plumedLine.*"))
"             call add(res,{"word":"...","menu":"(start multiline statement)"})
"           endif
"           if("..." =~ '^' . a:base && (key=~"^plumedLCLine.*" || key=~"^plumedCLine.*") && getline('.')=~'^\s*$')
"              call add(res,{"word":"...","menu":"(end multiline statement)"})
"           endif
            if("#" =~ '^' . a:base && key!="plumedComment") 
               call add(res,{"word":"#","menu":"(add comment)"})
            endif
            if("ENDPLUMED" =~ '^' . a:base && key =="")
               call add(res,{"word":"ENDPLUMED","menu":"(end input)"})
            endif
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



