" Runtime syntax for the selected PLUMED installation (Vim 8 or later).
if exists('g:loaded_plumed_cache')
  finish
endif
let g:loaded_plumed_cache = 1
let s:pending = {}
let s:failed = {}
let s:temporary = ''
let s:serial = 0

function! s:Absolute(path) abort
  return a:path =~# '^/' || a:path =~# '^\a:[/\\]'
endfunction

function! s:CacheRoot() abort
  if exists('g:plumed_cache_dir')
    return fnamemodify(expand(g:plumed_cache_dir), ':p')
  elseif s:Absolute($XDG_CACHE_HOME)
    return $XDG_CACHE_HOME . '/plumed/vim'
  elseif has('macunix')
    return expand('~/Library/Caches/plumed/vim')
  elseif has('win32') && !empty($LOCALAPPDATA)
    return $LOCALAPPDATA . '/plumed/vim'
  endif
  return expand('~/.cache/plumed/vim')
endfunction

function! s:TemporaryRoot() abort
  if empty(s:temporary)
    let s:temporary = tempname()
    call mkdir(s:temporary, 'p', 0700)
  endif
  return s:temporary
endfunction

function! s:Stamp(path) abort
  let path = resolve(fnamemodify(a:path, ':p'))
  return [path, getftime(path), getfsize(path)]
endfunction

" Describe files without loading libraries or hashing large binaries. Read the
" small registration files too, so new/removed packages invalidate the cache.
function! s:Describe(root) abort
  let executable = get(g:, 'plumed_executable', '')
  if empty(executable)
    for candidate in [a:root . '/../src/lib/plumed',
          \ a:root . '/../../../bin/plumed', a:root . '/../../../bin/plumed.exe']
      if executable(candidate)
        let executable = candidate
        break
      endif
    endfor
  endif
  let executable = exepath(empty(executable) ? 'plumed' : executable)
  let executable = empty(executable) ? '' : resolve(fnamemodify(executable, ':p'))
  let metadata = ['plumed-vim-cache-1', a:root, s:Stamp(executable),
        \ $PLUMED_ROOT, $PLUMED_KERNEL, $PLUMED_PLUGIN_PATH, $PLUMED_AUTOLOAD,
        \ $PLUMED_MODULE_OVERRIDE, $LD_LIBRARY_PATH, $DYLD_LIBRARY_PATH,
        \ $DYLD_FALLBACK_LIBRARY_PATH, get(g:, 'plumed_modules', [])]
  for template in ['plumed.vim.in', 'plumedf.vim']
    call add(metadata, s:Stamp(a:root . '/templates/' . template))
  endfor
  call add(metadata, s:Stamp(a:root . '/autoload/plumed.vim'))
  let kernels = glob(escape(fnamemodify(executable, ':h:h'), '*?[]{}') . '/lib/libplumedKernel.*', 0, 1)
        \ + glob(escape(a:root, '*?[]{}') . '/../src/lib/libplumedKernel.*', 0, 1)
        \ + glob(escape(fnamemodify(executable, ':h'), '*?[]{}') . '/plumedKernel.dll', 0, 1)
  if !empty($PLUMED_KERNEL)
    call add(kernels, $PLUMED_KERNEL)
  endif
  for selection in split($PLUMED_MODULE_OVERRIDE, ';')
    let separator = stridx(selection, '=')
    if separator >= 0
      call add(metadata, s:Stamp(strpart(selection, separator + 1)))
    endif
  endfor
  let directories = [a:root . '/../plugins']
        \ + split($PLUMED_PLUGIN_PATH, has('win32') ? ';' : ':')
  for kernel in uniq(sort(kernels))
    call add(metadata, s:Stamp(kernel))
    call add(directories, fnamemodify(resolve(kernel), ':h') . '/plumed/plugins')
  endfor
  for directory in uniq(sort(map(directories, 'resolve(fnamemodify(v:val, ":p"))')))
    call add(metadata, directory)
    for registration in sort(globpath(escape(directory, ',*?[]{}'), '*.plumed', 0, 1))
      let lines = readfile(registration)
      call add(metadata, [s:Stamp(registration), lines])
      for line in lines
        if line =~# '^library '
          let library = strpart(line, 8)
          if !s:Absolute(library)
            let library = directory . '/' . library
          endif
          call add(metadata, s:Stamp(library))
        endif
      endfor
    endfor
  endfor
  return {'root': a:root, 'executable': executable, 'key': sha256(string(metadata))}
endfunction

function! s:Current(entry) abort
  let pointer = a:entry . '/current'
  if filereadable(pointer)
    try
      let lines = readfile(pointer, '', 1)
    catch /^Vim\%((\a\+)\)\=:E484:/
      " Another editor can replace the hint between the check and the read.
      let lines = []
    endtry
    if !empty(lines) && lines[0] =~# '^generation-[0-9-]\+$'
      let path = a:entry . '/' . lines[0]
      if filereadable(path . '/syntax/plumed.vim') && isdirectory(path . '/help')
        return path
      endif
    endif
  endif
  " Vim's rename() unlinks an existing destination before renaming, even on
  " Unix. During that gap use an immutable completed generation. Never pick a
  " directory merely because its generator has already written syntax/help.
  for marker in reverse(sort(globpath(escape(a:entry, ',*?[]{}'), 'generation-*/complete', 0, 1)))
    let path = fnamemodify(marker, ':h')
    if filereadable(path . '/syntax/plumed.vim') && isdirectory(path . '/help')
      return path
    endif
  endfor
  return ''
endfunction

function! s:Writable(entry) abort
  try
    if !isdirectory(a:entry)
      call mkdir(a:entry, 'p', 0700)
    endif
    return filewritable(a:entry) == 2
  catch
    return 0
  endtry
endfunction

function! s:Source(path) abort
  if get(b:, 'current_syntax', '') ==# 'plumed'
    syntax clear
    unlet! b:current_syntax
  endif
  execute 'source ' . fnameescape(a:path . '/syntax/plumed.vim')
  let b:plumed_cache_path = a:path
endfunction

" Refresh visible buffers without moving the cursor, changing windows or
" touching their contents. Hidden buffers pick up the cache on BufEnter.
function! s:RefreshBuffers(state) abort
  if !exists('*win_findbuf') || !exists('*win_execute')
    return
  endif
  for buffer in getbufinfo()
    if getbufvar(buffer.bufnr, 'plumed_cache_key', '') ==# a:state.key
      for window in win_findbuf(buffer.bufnr)
        " This editor just completed this generation. Do not reread a shared
        " hint that another publisher can temporarily remove or replace.
        let command = 'call ' . expand('<SID>') . 'Source(' . string(a:state.stage) . ')'
        call win_execute(window, command)
      endfor
    endif
  endfor
endfunction

function! s:Finished(state, job, code) abort
  let state = a:state
  let published = 0
  call remove(s:pending, state.key)
  if has_key(state, 'job')
    call remove(state, 'job')
  endif
  try
    if a:code != 0 || !filereadable(state.stage . '/syntax/plumed.vim')
      throw 'gen_vim failed (exit ' . a:code . ')'
    endif
    if s:Describe(state.root).key !=# state.key
      throw 'installation changed while gen_vim was running; use :PlumedRefresh'
    endif
    " Publish the immutable completion marker first. The shared pointer is a
    " fast lookup hint: readers can use completed generations while rename()
    " temporarily removes it. Once marked, never delete this generation.
    call writefile([], state.stage . '/complete')
    let published = 1
    let pointer = state.stage . '/current.new'
    call writefile([fnamemodify(state.stage, ':t')], pointer)
    if rename(pointer, state.entry . '/current') != 0
      throw 'cannot publish generated Vim files'
    endif
    call s:RefreshBuffers(state)
  catch
    let s:failed[state.key] = 1
    let g:plumed_cache_error = v:exception
    if !published
      call delete(state.stage, 'rf')
    endif
  endtry
endfunction

function! s:Generate(state) abort
  let state = a:state
  let s:serial += 1
  let state.stage = state.entry . '/generation-' . getpid() . '-' . localtime() . '-' . s:serial
  call mkdir(state.stage, '', 0700)
  let command = [state.executable, '--no-mpi']
  for module in get(g:, 'plumed_modules', [])
    call extend(command, ['--load-module', module])
  endfor
  call extend(command, ['gen_vim', '--outdir', state.stage])
  if isdirectory(state.root . '/templates')
    call extend(command, ['--template-dir', state.root . '/templates'])
  endif
  let s:pending[state.key] = state
  if has('job') && exists('*win_execute') && get(g:, 'plumed_vim_async', 1)
    let state.job = job_start(command, {'in_io': 'null', 'out_io': 'null',
          \ 'err_io': 'file', 'err_name': state.stage . '/generation.log',
          \ 'exit_cb': function('s:Finished', [state])})
    if job_status(state.job) ==# 'fail'
      call s:Finished(state, state.job, -1)
    endif
  else
    " Minimal Vim builds have a synchronous first-use fallback. A warm cache
    " never starts a process, on either path.
    call system(join(map(copy(command), 'shellescape(v:val)'), ' ') . ' 2>&1')
    call s:Finished(state, 0, v:shell_error)
  endif
endfunction

function! s:Update(root, force) abort
  let state = s:Describe(a:root)
  let b:plumed_cache_root = a:root
  let b:plumed_cache_key = state.key
  let state.entry = s:CacheRoot() . '/' . state.key
  let current = s:Current(state.entry)
  if empty(current) && !empty(s:temporary)
    let state.entry = s:temporary . '/' . state.key
    let current = s:Current(state.entry)
  endif
  if !empty(current) && !a:force
    if get(b:, 'plumed_cache_path', '') !=# current || get(b:, 'current_syntax', '') !=# 'plumed'
      call s:Source(current)
    endif
    return 1
  endif
  if empty(state.executable) || has_key(s:pending, state.key)
        \ || (!a:force && has_key(s:failed, state.key))
    return get(b:, 'current_syntax', '') ==# 'plumed' && !empty(get(b:, 'plumed_cache_path', ''))
  endif
  if !s:Writable(state.entry)
    let state.entry = s:TemporaryRoot() . '/' . state.key
    if !s:Writable(state.entry)
      return 0
    endif
  endif
  call s:Generate(state)
  if !has_key(s:pending, state.key)
    let current = s:Current(state.entry)
    if !empty(current)
      call s:Source(current)
      return 1
    endif
  endif
  return get(b:, 'current_syntax', '') ==# 'plumed' && !empty(get(b:, 'plumed_cache_path', ''))
endfunction

" Called by the installed/generated syntax file before its static dictionary.
function! plumed#Load(root) abort
  if !get(g:, 'plumed_vim_dynamic', 1) || !exists('*sha256')
    return 0
  endif
  try
    return s:Update(resolve(fnamemodify(a:root, ':p')), 0)
  catch
    let g:plumed_cache_error = v:exception
    return 0
  endtry
endfunction

function! plumed#Update(force) abort
  if &l:filetype !=# 'plumed' || !exists('b:plumed_cache_root')
        \ || !get(g:, 'plumed_vim_dynamic', 1)
    return
  endif
  try
    call s:Update(b:plumed_cache_root, a:force)
  catch
    let g:plumed_cache_error = v:exception
  endtry
endfunction

command! -bar PlumedRefresh call plumed#Update(1)
augroup plumed_cache
  autocmd!
  autocmd BufEnter * call plumed#Update(0)
augroup END
