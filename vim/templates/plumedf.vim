
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

