
if exists("b:current_syntax")
  finish
endif

let b:current_syntax="plumedf"

if get(g:, 'plumed_shortcuts', 0)
  noremap <buffer> <F3> :PMinus<CR>
  noremap <buffer> <F4> :PPlus<CR>
endif

command! -buffer -nargs=0 -count=1 PPlus call PlumedMove(<count>)
command! -buffer -nargs=0 -count=1 PMinus call PlumedMove(-<count>)
command! -buffer -nargs=? -count=0 PCol call PlumedSelectColumn(<count>, <q-args>)
command! -buffer -nargs=0 PField call PlumedField()

if exists('g:plumedf_synmaxcol')
  let &l:synmaxcol = g:plumedf_synmaxcol
endif

" Resolve names against the nearest preceding header, without moving the cursor.
function! PlumedFields() abort
  let row = line('.')
  while row >= 1
    if getline(row) =~# '^\s*#!\s\+FIELDS\s'
      return split(substitute(getline(row), '^\s*#!\s\+FIELDS\s*', '', ''))
    endif
    let row -= 1
  endwhile
  return []
endfunction

function! PlumedSelectColumn(count, name) abort
  if empty(a:name)
    call PlumedColumn(a:count)
  else
    let column = index(PlumedFields(), a:name) + 1
    if !column
      echoerr 'No such field in the preceding FIELDS header: ' . a:name
      return
    endif
    call PlumedColumn(column)
  endif
endfunction

function! PlumedField() abort
  let text = getline('.')
  let prefix = strpart(text, 0, col('.'))
  let column = len(split(prefix))
  if text =~# '^\s*#!\s\+FIELDS\s'
    let column -= 2
  elseif text =~# '^\s*#'
    echo 'Not a data field'
    return
  endif
  let fields = PlumedFields()
  if column > 0 && column <= len(fields)
    echo column . ': ' . fields[column - 1]
  else
    echo 'No field name at cursor'
  endif
endfunction

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
 let b:current_syntax="plumedf"
 let coll=split(a:col,'\zs') " splits in characters
 let cmax=len(coll)
 let groups=[]
 for index in range(cmax)
   let mark=coll[index]
   let role = mark ==# 'X' ? 'Selected' : (mark ==# 'B' || (mark ==# '.' && index%2) ? 'Even' : 'Odd')
   call add(groups, 'plumedFCol' . (index+1) . role)
   execute 'highlight default link ' . groups[-1] . ' plumedF' . role
 endfor
 highlight default link plumedFOdd Constant
 highlight default link plumedFEven Keyword
 highlight default link plumedFSelected Todo
 let c=cmax
 while c > 0
   let cc = c == cmax ? cmax-1 : c+1
   execute 'syn match plumedFCol' . c . 'S /\v\s+/ nextgroup=' . groups[cc-1] . ' contained'
   let contained = c == 1 ? '' : 'contained'
   execute 'syn match ' . groups[c-1] . ' /\v\S+/ nextgroup=plumedFCol' . c . 'S ' . contained
   let c -= 1
 endwhile

 syn match plumedNothing  /\v.*/    contained
 syn match plumedSCol2    /\v\S+/    nextgroup=plumedNothing contained
 syn match plumedSCol1S   /\v\s+/    nextgroup=plumedSCol2 contained
 syn match plumedSCol1    /\v\S+/    nextgroup=plumedSCol1S contained
 highlight default link plumedSCol1 Type
 highlight default link plumedSCol2 Constant
 syntax match   plumedFComment excludenl /\v#.*$/
 highlight default link plumedFComment Comment
 execute 'syntax region plumedFieldLine matchgroup=plumedFieldLine start=/\v^[ \t]*\#\![ \t]+FIELDS[ \t]+/ excludenl end=/$/ contains=' . groups[0]
 syntax region plumedSetLine    matchgroup=plumedSetLine start=/\v[ \t]*\#\![ \t]+SET[ \t]+/ excludenl end=/$/ contains=plumedSCol1
 highlight default link plumedFieldLine Type
 highlight default link plumedSetLine   Type
endfunction

" initialize to column zero
call PlumedColumn(get(b:, 'plumedf_column', 0))

