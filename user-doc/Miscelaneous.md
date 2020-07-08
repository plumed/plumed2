\page Miscellaneous Miscellaneous

- \subpage comments
- \subpage ContinuationLines
- \subpage VimSyntax
- \subpage BashAutocompletion
- \subpage includes
- \subpage load
- \subpage embed
- \subpage degub
- \subpage exchange-patterns
- \subpage mymodules
- \subpage special-replica-syntax
- \subpage parsing-constants
- \subpage misc

\page comments Comments

If you are an organized sort of person who likes to remember what the hell you were trying to do when you ran a 
particular simulation you might find it useful to put comments in your input file.  In PLUMED you can do this as 
comments can be added using a # sign.  On any given line everything after the # sign is ignored so 
erm... yes add lines of comments or trailing comments to your hearts content as shown below (using Shakespeare is optional):

\plumedfile
# This is the distance between two atoms:
d1: DISTANCE ATOMS=1,2 
UPPER_WALLS ARG=d1 AT=3.0 KAPPA=3.0 LABEL=Snout # In this same interlude it doth befall.
# That I, one Snout by name, present a wall.
\endplumedfile
(see \ref DISTANCE and \ref UPPER_WALLS)

An alternative to including comments in this way is to use the command \subpage ENDPLUMED.  Everything in the PLUMED input after this
keyword will be ignored.

\page ContinuationLines Continuation lines

If your input lines get very long then editing them using vi and other such text editors becomes a massive pain in the arse.  
We at PLUMED are aware of this fact and thus have provided a way of doing line continuations so as to make your life that much 
easier - aren't we kind?  Well no not really, we have to use this code too.  Anyway, you can do continuations by using the "..." syntax
as this makes this: 

\plumedfile
DISTANCES ATOMS1=1,300 ATOMS2=1,400 ATOMS3=1,500 LABEL=dist
\endplumedfile
(see \ref DISTANCES)

equivalent to this:

\plumedfile
DISTANCES ...
  LABEL=dist
# we can also insert comments here
  ATOMS1=1,300
# multiple kewords per line are allowed
  ATOMS2=1,400 ATOMS3=1,500
#empty lines are also allowed

... DISTANCES
\endplumedfile

Notice that the closing `...` is followed by the word `DISTANCES`. This is optional, but might be
useful to find more easily which is the matching start of the statement. The following is equally correct
\plumedfile
dist: DISTANCES ...
# we can also insert comments here
  ATOMS1=1,300
# multiple kewords per line are allowed
  ATOMS2=1,400 ATOMS3=1,500
#empty lines are also allowed

...
\endplumedfile

Notice that PLUMED makes a check that the word following the closing `...` is actually identical to
the first word in the line with the first `...`. If not, it will throw an error.
Also notice that you might put more than one word in the first line. E.g.
\plumedfile
DISTANCES LABEL=dist ...
# we can also insert comments here
  ATOMS1=1,300
# multiple kewords per line are allowed
  ATOMS2=1,400 ATOMS3=1,500
#empty lines are also allowed
...
\endplumedfile
or, equivalently,
\plumedfile
dist: DISTANCES ...
# we can also insert comments here
  ATOMS1=1,300
# multiple kewords per line are allowed
  ATOMS2=1,400 ATOMS3=1,500
#empty lines are also allowed
...  
\endplumedfile

\page BashAutocompletion Using bash autocompletion

When possible, PLUMED tries to install bash autocompletion so that
you do not have to do anything. Just use the `<TAB>` key to complete
plumed commands (e.g. `plumed dr<TAB>`) or even options (e.g. `plumed driver --i<TAB>`).
In case this does not work, you might have to add the following lines to your .bashrc file:
\verbatim
_plumed() { eval "$(plumed --no-mpi completion 2>/dev/null)";}
complete -F _plumed -o default plumed
\endverbatim

\par Effect

When typing on the the shell you should observe the following behavior.
\verbatim
> plumed <TAB>
\endverbatim
will autocomplete with the names of the available PLUMED commands (e.g. `driver`, `help`, etc).
\verbatim
> plumed -<TAB>
\endverbatim
will autocomplete with the available PLUMED options (e.g. `--no-mpi`, etc).

PLUMED also knows which are the options available for each command
(e.g. `plumed driver --natoms`). So, the following
\verbatim
> plumed driver -<TAB>
\endverbatim
(notice the `-`) will autocomplete to the options of `plumed driver`. On the contrary
\verbatim
> plumed driver --ixtc <TAB>
\endverbatim
(notice the there is no `-` before `<TAB>`) will autocomplete to the files in the current directory.

Also notice that every time you use the `<TAB>` key to autocomplete the command `plumed` will be invoked.
This should allow the correct commands and options to be reported depending on the exact `plumed` command
in the current execution path. For instance, if you have multiple PLUMED versions installed with
env modules, you should be able to see the commands available in the currently loaded version.
Clearly, this feature will only be available if `plumed` can run on this machine (that is: will not work
if you are cross compiling). This is not a problem since you are not expecting to run the `plumed` command
in this specific case.

\par Technicalities

At configure time if the variable `BASH_COMPLETION_DIR` is defined it will be used to decide where PLUMED
autocompletion should be installed. Otherwise, configure will look for the presence of the `bash-completion` package
and, in case it is installed on the same prefix as PLUMED, also PLUMED autocompletion will be installed.
Finally, if none of these two conditions are satisfied, autocompletion will not be enabled. You will
have to change your bashrc file once adding the following lines:
\verbatim
_plumed() { eval "$(plumed --no-mpi completion 2>/dev/null)";}
complete -F _plumed -o default plumed
\endverbatim
The command `plumed completion` just writes on its standard output the body of a bash function that is
then used by bash to construct the autocompletion.
The `--no-mpi` flag makes it more likely that the command can be executed correctly e.g. when you are on the login node of a cluster and
PLUMED was compiled with MPI but the login node does not support MPI. In other cases, it is harmless.
The `-o default` options will make sure that if `plumed --no-mpi completion` returns an error the default bash completion
will be used. This is what will happen if you load an older PLUMED version for which the `completion` command is not available yet.
In future PLUMED versions the `plumed completion` command might return more sophisticated functions. You should
be able to benefit of these features without ever changing your bash configuration file again.

\par Multiple versions and suffixes

In case you have multiple versions of PLUMED installed in separate env modules there is nothing more to do.
However, if you have have multiple versions of PLUMED installed with different suffixes you should
consistently add more lines to your profile file. For instance, if you installed two executables named
`plumed` and `plumed_mpi` your configuration file should look like:
\verbatim
_plumed() { eval "$(plumed --no-mpi completion 2>/dev/null)";}
complete -F _plumed -o default plumed
_plumed_mpi() { eval "$(plumed_mpi --no-mpi completion 2>/dev/null)";}
complete -F _plumed_mpi -o default plumed_mpi
\endverbatim

\page VimSyntax Using VIM syntax file

For the impatient use:
- Add the following to your .vimrc file:
\verbatim
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
\endverbatim
- When you open a PLUMED input file, you can enable syntax highlighting with:
\verbatim
:set ft=plumed
\endverbatim
This will also enable autocompletion. Use `<CTRL-X><CTRL-O>` to autocomplete a word.
- If you want to fold multi-line statements, type
\verbatim
:setlocal foldmethod=syntax
\endverbatim
- While editing a plumed input file, you can use command `:PHelp` (or shortcut `<F2>`)
  to show in a split window a short help about the action defined in the line where the cursor is.
  Typing `:PHelp` again (or pushing `<F2>`) you will
  close that window. With `<CTRL-W><CTRL-W>` you go back and forth between the two windows.
- When you open a file starting with `#! FIELDS`, VIM will automatically understand it
  is a PLUMED output file (VIM filetype = plumedf) and will color fields and data columns with
  alternating colors. Typing `:PPlus` and `:PMinus` (or pushing `<F3>` and `<F4>`)
  you can move a highlighted column.

See below for more detailed instructions.

\par Configuration

When PLUMED is compiled, directories `help` and `syntax` will appear in `builddir/vim`.
They contain a VIM plugin that can be used to highlight proper PLUMED instructions in a PLUMED
input file and to quickly retrieve help.
There is also a file `builddir/vim/scripts.vim` that helps VIM in recognizing PLUMED output files.

\warning
Notice that these file do not appear if you are cross compiling.
In this case, you must copy the plugin files from another machine.

To make VIM aware of these files, you should copy them to your `$HOME/.vim` directory.
Later you can
enable plumed syntax with the command
\verbatim
:set ft=plumed
\endverbatim

If you work in an environment where several PLUMED versions are installed (e.g. using env modules),
we recommend the following procedure:
- Install PLUMED
- Add to your `.vimrc` file the following line:
\verbatim
:let &runtimepath.=','.$PLUMED_VIMPATH
\endverbatim

The modulefile provided with PLUMED should set the PLUMED_VIMPATH environment variable
to the proper path.
Thus, when working with a given PLUMED module loaded, you should be able to
enable to proper syntax by just typing
\verbatim
:set ft=plumed
\endverbatim
in VIM.
Notice that the variable `PLUMED_VIMPATH` is also set in the `sourceme.sh` script in the build directory.
Thus, if you modify your `.vimrc` file as suggested, you will be able to use the correct syntax both
when using an installed PLUMED and when running from a just compiled copy.
Finally, in case you have both a pre-installed PLUMED **and** you have your development version
the following command would give you the optimal flexibility:
\verbatim
:let &runtimepath.=','.$PLUMED_VIMPATH.',/opt/local/lib/plumed/vim/'
\endverbatim
The environment variable `PLUMED_VIMPATH`, if set, will take the precedence.
Otherwise, vim will resort to the hard coded path.
In this case we assumed that there is a PLUMED installed in `/opt/local/` (e.g. using MacPorts),
but you can override it sourcing a `sourceme.sh` file in the compilation directory
or loading a PLUMED module with `module load plumed`.

If you are tired of typing `:set ft=plumed`, you can use a modeline.
Add to your `.vimrc` file the following commands
\verbatim
:set modeline
:set modelines=5
\endverbatim
Then, at the beginning of your PLUMED input file, put the following comment:
\plumedfile
# vim:ft=plumed
d: DISTANCE ATOMS=1,2
RESTRAINT ARG=d AT=0.0 KAPPA=1.0
\endplumedfile
Now, every time you open this file, you will see it highlighted.

\par Syntax highlighting

The syntax file contains a definition of all possible PLUMED actions and keywords.
It is designed to allow for a quick validation of the PLUMED input file before running it.
As such, all the meaningful words in the input should be highlighted:
- Valid action names (such as `METAD`) and labels (such as `m:` or `LABEL=m`) will be
  highlighted in the brightest way (`Type` in VIM). Those are the most important words.
- Keyword and flag names (such as `ATOMS=` or `COMPONENTS` when part of the action \ref DISTANCE) will be highlighted with a different color
  (`Statement` in VIM).
- Values provided by users (such as the number of the atoms following `ATOMS=`) will be highlighted with a different color
  (`String` in VIM).
- Comments (see \ref comments) will be highlighted as comments (`Comment` in VIM).
- String `__FILL__` (extensively used in tutorials to indicate parts to be completed) is highlighted (`Todo` in VIM).

If you see something that is not highlighted and appears in black, this is likely going to result in an error at runtime.
Think of this as a sort of preliminary spell-check.
For this checks to be effective, we recommend to use a syntax file generated with
exactly the same version of PLUMED that you are using.
In case you find that parts of an input file that is valid are not highlighted, then please
report it as a bug.
On the contrary, you cannot expect the VIM syntax file to recognize all possible errors
in a PLUMED input. Thus, a file for  which the highlighting looks correct might still contain errors.

\par Multi-line folding

Notice that syntax highlighting also allow VIM to properly fold multi-line actions.
Try to do the following:
- Open a PLUMED input file
- Enable PLUMED syntax
\verbatim
:set ft=plumed
\endverbatim
- Enable syntax-based folding
\verbatim
:setlocal foldmethod=syntax
\endverbatim

Now look at what happened to all the multi-line statements in PLUMED (i.e. those using
\ref ContinuationLines).  As you can see, they will be folded into single lines.
Folded lines can be expanded with `zo` and folded with `zc`. Look at VIM documentation to
learn more.
In case you want to use this feature, we suggest you to put both label
and action type on the first line of multi-line statements. E.g.
\plumedfile
d: DISTANCE ATOMS=1,2

m: METAD ...
  ARG=d
  HEIGHT=1.0
  SIGMA=0.5
  PACE=100
...
\endplumedfile
will be folded to
\verbatim
d: DISTANCE ATOMS=1,2

+--  6 lines: m: METAD ...------------------------------------------------------
\endverbatim
and
\plumedfile
d: DISTANCE ATOMS=1,2

METAD LABEL=m ...
  ARG=d
  HEIGHT=1.0
  SIGMA=0.5
  PACE=100
...
\endplumedfile
will be folded to
\verbatim
d: DISTANCE ATOMS=1,2

+--  6 lines: METAD LABEL=m ...-------------------------------------------------
\endverbatim
This will allow you to easily identify the folded lines by seeing the most important information,
that is the action type (`METAD`) and its label (`m`). This feature is convenient if
you want to browse files that contain a lot of actions defined on multiple lines.

\par Autocompletion

Another VIM feature that comes when you load PLUMED syntax is autocompletion of PLUMED
actions and keywords. Open your favorite PLUMED input file and set it to PLUMED syntax highlighting with
\verbatim
:set ft=plumed
\endverbatim
Now go into insert mode pressing `i` and type `DU` followed by `<CTRL+X><CTRL+O>`.
Here `<CTRL+X>` stands for autocompletion and `<CTRL+O>` for omnifunc autocompletion. You will see a short
menu listing the following actions
\verbatim
DUMPATOMS       
DUMPDERIVATIVES 
DUMPFORCES      
DUMPMASSCHARGE  
DUMPMULTICOLVAR 
DUMPPROJECTIONS 
\endverbatim
That is, all the actions starting with `DU`.
You can navigate it with up and down arrows so as to choose the
best match.

Notice that the default behavior of VIM is to use the first match by default.
In the first example (`DU<CTRL+X><CTRL+O`), it would be `DUMPATOMS`.
The following settings make it work as most of the people expect:
\verbatim
:set completeopt=longest,menuone
\endverbatim
With these settings, in the first example (`DU<CTRL+X><CTRL+O`) VIM will only complete up to the longest common part (`DUMP`).

As you can imagine,
if you use autocompletion after you have typed the word `DISTANCE` followed by a space you will see
a menu listing `LABEL=`, `COMPONENTS`, etc. Basically, all the keywords that are possibly used within a `DISTANCE` line
will be shown. This is very useful if you do not remember the exact name of the keywords associated with
a given action.

\par Quick help

You can also retrieve quick explanation of the input options for a specific action.
Try to do the following. Enable plumed syntax:
\verbatim
:set ft=plumed
\endverbatim
Then add the following line
\verbatim
DISTANCE
\endverbatim
Now, in normal mode, go with the cursor on the `DISTANCE` line and type
\verbatim
:PHelp
\endverbatim
A new split window should appear containing some documentation about the \ref DISTANCE collective variable.
You can go back and forth between the two windows with `<CTRL+W><CTRL+W>`, as usually in vim.
Notice that if you are in the help window and type `:PHelp` this window will be closed.

To make the navigation easier, you can add a shortcut in your .vimrc file. For example, adding:
\verbatim
: nmap <F2> : PHelp<CR>
\endverbatim
you should be able to open and close the manual hitting the F2 key.
This is done automatically in the PLUMED syntax file if you add `let plumed_shortcuts=1` to your
vimrc file.

\par Displaying output files

Most of the PLUMED output files look like this
\verbatim
#! FIELDS A B C
1 2 3
\endverbatim
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
\verbatim
:5PCol
\endverbatim
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
\verbatim
: map <F3> :PMinus<CR>
: map <F4> :PPlus<CR>
\endverbatim
you should be able to move the highlight column using F3 and F4 buttons.
This is done automatically in the PLUMED syntax file if you add `let plumed_shortcuts=1` to your
vimrc file.

\page includes Including other files 

If, for some reason, you want to spread your PLUMED input over a number of files you can use \subpage INCLUDE as shown below:

\plumedfile
INCLUDE FILE=filename
\endplumedfile

So, for example, a single "plumed.dat" file:

\plumedfile
DISTANCE ATOMS=1,2 LABEL=dist
RESTRAINT ARG=dist AT=2.0 KAPPA=1.0
\endplumedfile

could be split up into two files as shown below:
 
\plumedfile
DISTANCE ATOMS=1,2 LABEL=dist
INCLUDE FILE=toBeIncluded.inc
\endplumedfile
plus a "toBeIncluded.inc" file
\plumedfile
#SETTINGS FILENAME=toBeIncluded.inc
# this is toBeIncluded.inc
RESTRAINT ARG=dist AT=2.0 KAPPA=1.0
\endplumedfile

However, when you do this it is important to recognize that \ref INCLUDE is a real directive that is only resolved
after all the \ref comments have been stripped and the \ref ContinuationLines have been unrolled.  This means it
is not possible to do things like:

\plumedfile
# this is wrong:
DISTANCE INCLUDE FILE=options.dat
RESTRAINT ARG=dist AT=2.0 KAPPA=1.0
\endplumedfile

\page load Loading shared libraries

You can introduce new functionality into PLUMED by placing it directly into the src directory and recompiling the 
PLUMED libraries.  Alternatively, if you want to keep your code independent from the rest of PLUMED (perhaps
so you can release it independently - we won't be offended), then you can create your own dynamic library.  To use this 
in conjunction with PLUMED you can then load it at runtime by using the \subpage LOAD keyword as shown below:

\plumedfile
LOAD FILE=library.so
\endplumedfile
 
N.B.  If your system uses a different suffix for dynamic libraries (e.g. macs use .dylib) then PLUMED will try to 
automatically adjust the suffix accordingly.

\page embed Embed a separate PLUMED instance

\subpage PLUMED

\page degub Debugging the code

The \subpage DEBUG action provides some functionality for debugging the code that may be useful if you are doing 
very intensive development of the code of if you are running on a computer with a strange architecture.

\page exchange-patterns Changing exchange patterns in replica exchange

Using the \subpage RANDOM_EXCHANGES keyword it is possible to make exchanges between randomly
chosen replicas. This is useful e.g. for bias exchange metadynamics \cite piana.

\page special-replica-syntax Special replica syntax

(this part of the manual is based on \ref trieste-5-replica-special-syntax).

In many cases, we need to run multiple replicas with almost identical PLUMED files.
These files might be prepared with cut-and-paste, which is very error prone,
or could be set up with some smart bash or python script. Additionally,
one can take advantage of the \ref INCLUDE keyword so as to have a shared input
file with common definitions and specific input files with replica-dependent keywords.
However, as of PLUMED 2.4, we introduced a simpler manner to manipulate multiple replica
inputs with tiny differences. Look at the following example:

\plumedfile 
#SETTINGS NREPLICAS=3
# Compute a distance
d: DISTANCE ATOMS=1,2

# Apply a restraint.
RESTRAINT ARG=d AT=@replicas:1.0,1.1,1.2 KAPPA=1.0
# On replica 0, this means:
#   RESTRAINT ARG=d AT=1.0 KAPPA=1.0
# On replica 1, this means:
#   RESTRAINT ARG=d AT=1.1 KAPPA=1.0
# On replica 2, this means:
#   RESTRAINT ARG=d AT=1.2 KAPPA=1.0
\endplumedfile

If you prepare a single `plumed.dat` file like this one and feeds it to PLUMED while using 3 replicas,
the 3 replicas will see the very same input except for the `AT` keyword, that sets the position of the restraint.
Replica 0 will see a restraint centered at 1.0, replica 1 centered at 1.1, and replica 2 centered at 1.2.

The `@replicas:` keyword is not special for \ref RESTRAINT or for the `AT` keyword. Any keyword in PLUMED can accept that syntax.
For instance, the following single input file can be used to setup a bias exchange metadynamics \cite piana simulations:
\plumedfile
#SETTINGS NREPLICAS=2
# Compute distance between atoms 1 and 2
d: DISTANCE ATOMS=1,2

# Compute a torsional angle
t: TORSION ATOMS=30,31,32,33

# Metadynamics.
METAD ...
  ARG=@replicas:d,t
  HEIGHT=1.0
  PACE=100
  SIGMA=@replicas:0.1,0.3
  GRID_MIN=@replicas:0.0,-pi
  GRID_MAX=@replicas:2.0,pi
...
# On replica 0, this means:
#  METAD ARG=d HEIGHT=1.0 PACE=100 SIGMA=0.1 GRID_MIN=0.0 GRID_MAX=2.0
# On replica 1, this means:
#  METAD ARG=t HEIGHT=1.0 PACE=100 SIGMA=0.3 GRID_MIN=-pi GRID_MAX=pi
\endplumedfile

This would be a typical setup for a bias exchange simulation.
Notice that even though variables `d` and `t` are both read in both replicas,
`d` is only computed on replica 0 (and `t` is only computed on replica 1).
This is because variables that are defined but not used are never actually calculated by PLUMED.

If the value that should be provided for each replica is a vector, you should use curly braces as delimiters.
For instance, if the restraint acts on two variables, you can use the following input:

\plumedfile
#SETTINGS NREPLICAS=3
# Compute distance between atoms 1 and 2
d: DISTANCE ATOMS=10,20

# Compute a torsional angle
t: TORSION ATOMS=30,31,32,33

# Apply a restraint:
RESTRAINT ...
  ARG=d,t
  AT=@replicas:{{1.0,2.0} {3.0,4.0} {5.0,6.0}}
  KAPPA=1.0,3.0
...
# On replica 0 this means:
#  RESTRAINT ARG=d AT=1.0,2.0 KAPPA=1.0,3.0
# On replica 1 this means:
#  RESTRAINT ARG=d AT=3.0,4.0 KAPPA=1.0,3.0
# On replica 2 this means:
#  RESTRAINT ARG=d AT=5.0,6.0 KAPPA=1.0,3.0
\endplumedfile

Notice the double curly braces. The outer ones are used by PLUMED to know there the argument of the `AT` keyword ends,
whereas the inner ones are used to group the values corresponding to each replica.
Also notice that the last example can be split in multiple lines exploiting the fact that
within multi-line statements (enclosed by pairs of `...`) newlines are replaced with simple spaces:
\plumedfile
#SETTINGS NREPLICAS=3
d: DISTANCE ATOMS=10,20
t: TORSION ATOMS=30,31,32,33
RESTRAINT ...
  ARG=d,t
# indentation is not required (this is not python!)
# but makes the input easier to read
  AT=@replicas:{
    {1.0,2.0}
    {3.0,4.0}
    {5.0,6.0}
  }
  KAPPA=1.0,3.0
...
\endplumedfile

In short, whenever there are keywords that should vary across replicas, you should set them using the `@replicas:` keyword.
As mentioned above, you can always use the old syntax with separate input file, and this is recommended when the
number of keywords that are different is large.

\page parsing-constants Parsing constants

You might have noticed that from time to time constants are specified using strings rather than numbers.
An example is the following

\plumedfile
#SETTINGS MOLFILE=regtest/basic/rt65/AA.pdb
MOLINFO STRUCTURE=AA.pdb  MOLTYPE=rna
e1: TORSION ATOMS=@epsilon-1
t: METAD ARG=e1 SIGMA=0.15 PACE=10 HEIGHT=2 GRID_MIN=-pi GRID_MAX=pi GRID_BIN=200
\endplumedfile

Notice that the boundaries for `GRID_MIN` and `GRID_MAX` are `-pi` and `pi`. Until PLUMED 2.3,
we used a very dummy parses that could recognize only `pi` as a special string, plus strings such
as `0.5pi` and `-pi`. However, as of version 2.4, we use the Lepton library in order to parse every constant
that we read. This means that you can also employ more complicated expressions such as `1+2` or `exp(10)`:

\plumedfile
#SETTINGS MOLFILE=regtest/basic/rt65/AA.pdb
MOLINFO STRUCTURE=AA.pdb  MOLTYPE=rna
e1: TORSION ATOMS=@epsilon-1
RESTRAINT ARG=e1 AT=1+0.5
\endplumedfile

Notice that this applies to any quantity read by plumed as a real number, but does not apply
yet to integer numbers (e.g.: the PACE argument of \ref METAD).

\page misc Frequently used tools

@DICTIONARY@
<TABLE ALIGN="center" FRAME="void" WIDTH="95%%" CELLPADDING="5%%">
<TR>
<TD WIDTH="5%"> 
\subpage Regex </TD><TD> </TD><TD> POSIX regular expressions can be used to select multiple actions when using ARG (i.e. \ref PRINT).
</TD>
</TR>
<TR>
<TD WIDTH="5%"> 
\subpage Files </TD><TD> </TD><TD> Dealing with Input/Output
</TD>
</TR>
</TABLE>

