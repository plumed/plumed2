# Using the VIM syntax file

For the impatient use:
- Add the following to your .vimrc file:

````
" Enable syntax
:syntax on
" This allows including the proper PLUMED syntax file:
:let &runtimepath.=','.$PLUMED_VIMPATH
" The former command requires PLUMED_VIMPATH to be set. Alternatively, use this:
" let &runtimepath.=',/usr/local/lib/plumed/vim'
" properly adjusted to the path where PLUMED is installed.
" This makes autocompletion work in the expected way:
:set completeopt=longest,menuone
" This enables bindings of F2/F3/F4 to plumed specific commands:
:let plumed_shortcuts=1
````

- When you open a PLUMED input file, you can enable syntax highlighting with:

````
:set ft=plumed
````

This will also enable autocompletion. Use `<CTRL-X><CTRL-O>` to autocomplete a word.

- If you want to fold multi-line statements, type

````
:setlocal foldmethod=syntax
````

- While editing a plumed input file, you can use command `:PHelp` (or shortcut `<F2>`)
  to show in a split window a short help about the action defined in the line where the cursor is.
  Typing `:PHelp` again (or pushing `<F2>`) you will
  close that window. With `<CTRL-W><CTRL-W>` you go back and forth between the two windows.
- When you open a file starting with `#! FIELDS`, VIM will automatically understand it
  is a PLUMED output file (VIM filetype = plumedf) and will color fields and data columns with
  alternating colors. Typing `:PPlus` and `:PMinus` (or pushing `<F3>` and `<F4>`)
  you can move a highlighted column.

See below for more detailed instructions.

## Configuration

When PLUMED is compiled, directories `help` and `syntax` will appear in `builddir/vim`.
They contain a VIM plugin that can be used to highlight proper PLUMED instructions in a PLUMED
input file and to quickly retrieve help.
There is also a file `builddir/vim/scripts.vim` that helps VIM in recognizing PLUMED output files.

> [!WARNING]
> Notice that these file do not appear if you are cross compiling.
> In this case, you must copy the plugin files from another machine.

To make VIM aware of these files, you should copy them to your `$HOME/.vim` directory.
Later you can
enable plumed syntax with the command

````
:set ft=plumed
````

If you work in an environment where several PLUMED versions are installed (e.g. using env modules),
we recommend the following procedure:

- Install PLUMED
- Add to your `.vimrc` file the following line:

````
:let &runtimepath.=','.$PLUMED_VIMPATH
````

The modulefile provided with PLUMED should set the PLUMED_VIMPATH environment variable
to the proper path.
Thus, when working with a given PLUMED module loaded, you should be able to
enable to proper syntax by just typing

````
:set ft=plumed
````

in VIM.
Notice that the variable `PLUMED_VIMPATH` is also set in the `sourceme.sh` script in the build directory.
Thus, if you modify your `.vimrc` file as suggested, you will be able to use the correct syntax both
when using an installed PLUMED and when running from a just compiled copy.
Finally, in case you have both a pre-installed PLUMED **and** you have your development version
the following command would give you the optimal flexibility:

````
:let &runtimepath.=','.$PLUMED_VIMPATH.',/opt/local/lib/plumed/vim/'
````

The environment variable `PLUMED_VIMPATH`, if set, will take the precedence.
Otherwise, vim will resort to the hard coded path.
In this case we assumed that there is a PLUMED installed in `/opt/local/` (e.g. using MacPorts),
but you can override it sourcing a `sourceme.sh` file in the compilation directory
or loading a PLUMED module with `module load plumed`.

If you are tired of typing `:set ft=plumed`, you can use a modeline.
Add to your `.vimrc` file the following commands

````
:set modeline
:set modelines=5
````

Then, at the beginning of your PLUMED input file, put the following comment:

```plumed
# vim:ft=plumed
d: DISTANCE ATOMS=1,2
RESTRAINT ARG=d AT=0.0 KAPPA=1.0
```

Now, every time you open this file, you will see it highlighted.

## Syntax highlighting

The syntax file contains a definition of all possible PLUMED actions and keywords.
It is designed to allow for a quick validation of the PLUMED input file before running it.
As such, all the meaningful words in the input should be highlighted:
- Valid action names (such as `METAD`) and labels (such as `m:` or `LABEL=m`) will be
  highlighted in the brightest way (`Type` in VIM). Those are the most important words.
- Keyword and flag names (such as `ATOMS=` or `COMPONENTS` when part of the action [DISTANCE](DISTANCE.md)) will be highlighted with a different color
  (`Statement` in VIM).
- Values provided by users (such as the number of the atoms following `ATOMS=`) will be highlighted with a different color
  (`String` in VIM).
- Comments will be highlighted as comments (`Comment` in VIM).
- String `__FILL__` (extensively used in tutorials to indicate parts to be completed) is highlighted (`Todo` in VIM).

If you see something that is not highlighted and appears in black, this is likely going to result in an error at runtime.
Think of this as a sort of preliminary spell-check.
For this checks to be effective, we recommend to use a syntax file generated with
exactly the same version of PLUMED that you are using.
In case you find that parts of an input file that is valid are not highlighted, then please
report it as a bug.
On the contrary, you cannot expect the VIM syntax file to recognize all possible errors
in a PLUMED input. Thus, a file for  which the highlighting looks correct might still contain errors.

## Multi-line folding

Notice that syntax highlighting also allow VIM to properly fold multi-line actions.
Try to do the following:
- Open a PLUMED input file
- Enable PLUMED syntax

````
:set ft=plumed
````

- Enable syntax-based folding

````
:setlocal foldmethod=syntax
````

Now look at what happened to all the multi-line statements in PLUMED (i.e. those using
continuation lines).  As you can see, they will be folded into single lines.
Folded lines can be expanded with `zo` and folded with `zc`. Look at VIM documentation to
learn more.
In case you want to use this feature, we suggest you to put both label
and action type on the first line of multi-line statements. E.g.

```plumed
d: DISTANCE ATOMS=1,2

m: METAD ...
  ARG=d
  HEIGHT=1.0
  SIGMA=0.5
  PACE=100
...
```

will be folded to

````
d: DISTANCE ATOMS=1,2

+--  6 lines: m: METAD ...------------------------------------------------------
````

and

```plumed
d: DISTANCE ATOMS=1,2

METAD LABEL=m ...
  ARG=d
  HEIGHT=1.0
  SIGMA=0.5
  PACE=100
...
```

will be folded to

````
d: DISTANCE ATOMS=1,2

+--  6 lines: METAD LABEL=m ...-------------------------------------------------
````

This will allow you to easily identify the folded lines by seeing the most important information,
that is the action type (`METAD`) and its label (`m`). This feature is convenient if
you want to browse files that contain a lot of actions defined on multiple lines.

## Autocompletion

Another VIM feature that comes when you load PLUMED syntax is autocompletion of PLUMED
actions and keywords. Open your favorite PLUMED input file and set it to PLUMED syntax highlighting with

````
:set ft=plumed
````

Now go into insert mode pressing `i` and type `DU` followed by `<CTRL+X><CTRL+O>`.
Here `<CTRL+X>` stands for autocompletion and `<CTRL+O>` for omnifunc autocompletion. You will see a short
menu listing the following actions

````
DUMPATOMS       
DUMPDERIVATIVES 
DUMPFORCES      
DUMPMASSCHARGE  
DUMPMULTICOLVAR 
DUMPPROJECTIONS 
````

That is, all the actions starting with `DU`.
You can navigate it with up and down arrows so as to choose the
best match.

Notice that the default behavior of VIM is to use the first match by default.
In the first example (`DU<CTRL+X><CTRL+O`), it would be `DUMPATOMS`.
The following settings make it work as most of the people expect:

````
:set completeopt=longest,menuone
````

With these settings, in the first example (`DU<CTRL+X><CTRL+O`) VIM will only complete up to the longest common part (`DUMP`).

As you can imagine,
if you use autocompletion after you have typed the word `DISTANCE` followed by a space you will see
a menu listing `LABEL=`, `COMPONENTS`, etc. Basically, all the keywords that are possibly used within a `DISTANCE` line
will be shown. This is very useful if you do not remember the exact name of the keywords associated with
a given action.

## Quick help

You can also retrieve quick explanation of the input options for a specific action.
Try to do the following. Enable plumed syntax:

````
:set ft=plumed
````

Then add the following line

````
DISTANCE
````

Now, in normal mode, go with the cursor on the `DISTANCE` line and type

````
:PHelp
````

A new split window should appear containing some documentation about the [DISTANCE](DISTANCE.md) collective variable.
You can go back and forth between the two windows with `<CTRL+W><CTRL+W>`, as usually in vim.
Notice that if you are in the help window and type `:PHelp` this window will be closed.

To make the navigation easier, you can add a shortcut in your .vimrc file. For example, adding:

````
: nmap <F2> : PHelp<CR>
````

you should be able to open and close the manual hitting the F2 key.
This is done automatically in the PLUMED syntax file if you add `let plumed_shortcuts=1` to your
vimrc file.

## Displaying output files

Most of the PLUMED output files look like this

````
#! FIELDS A B C
1 2 3
````

This is useful since in the header you can see the name of the quantities that are printed in the
data lines. However, when you have an output file with many columns it might be a bit error prone to count them.
To simplify this, when PLUMED syntax for VIM is configured properly VIM should be able to:
- Detect that this file is a PLUMED output file with fields, automatically setting its type to `plumedf`. If not, just type
  `:set ft=plumedf`.
- Show this file with syntax highlighting to increase its readability.

Notice that the syntax file for the output files (`plumedf.vim`) is not the same one that is used for the PLUMED
input file (`plumed.vim`).

To make output files more readable, vim will show `FIELDS` and `SET` words in a different color,
and data columns with alternating colors (e.g. dark/light/dark/light).
The colors in the columns are consistent with those shown in the FIELD line.
In the example above, 1, 2, and 3 will be of the same color as A, B, and C respectively.
This should make it much easier to find which columns correspond to a given quantity.

It is also possible to highlight a specific field of the file. Typing

````
:5PCol
````

you will highlight the fifth field. Notice that in the `FIELDS` line (the first line of the file)
the seventh word of the line will be highlighted, which is the one containing the name of the field.
This allows for easy matching of values shown
in the file and tags provided in the `FIELDS` line.
The highlighted column can be moved back and forth using `:PPlus` and `:PMinus`.
Adding a count to the command will move the highlighted column more. E.g. `:2PPlus` will move
the column to the right twice.

If you have a long output file, it might be convenient to split it with
`:split` so that one of the two windows will only show the header. The other window
can be used to navigate the file.


To make the navigation easier, you can add a shortcut in your .vimrc file. For example, adding:

````
: map <F3> :PMinus<CR>
: map <F4> :PPlus<CR>
````

you should be able to move the highlight column using F3 and F4 buttons.
This is done automatically in the PLUMED syntax file if you add `let plumed_shortcuts=1` to your
vimrc file.
