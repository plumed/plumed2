\page Syntax Getting Started 

To run PLUMED you need to provide one input file.  In this file you specify what it
is that PLUMED should do during the course of the run.  Typically this will involve calculating 
one or more collective variables, perhaps calculating a function of these CVs
 and then doing some analysis of values of your collective variables/functions or running
some free energy method. A very brief introduction to the syntax used in the PLUMED input file
is provided in this <a href="http://www.youtube.com/watch?v=PxJP16qNCYs"> 10-minute video </a>.

Within this input file every line is an instruction for PLUMED to perform some particular action.  This could be
 the calculation of a colvar, an occasional analysis of the trajectory or a biasing of the dynamics.  The first
word in these lines specify what particular action is to be performed.  This is then followed by a number of keywords
which provide PLUMED with more details as to how the action is to be performed.  These keywords are either single words
(in which they tell PLUMED to do the calculation in a particular way - for example NOPBC tells PLUMED to not use the periodic
boundary conditions when calculating a particular colvar) or they can be words followed by an equals sign and a comma separated 
list _with no spaces_ of numbers or characters (so for example ATOMS=1,2,3,4 tells PLUMED to use atom numbers 1,2,3 and 4 in 
the calculation of a particular colvar).
The reason why spaces are not admitted is that PLUMED should be able to understand when the list of atoms
ended and a new keyword should be expected. 
Space separated lists can be used instead of comma separated list if the entire list
is enclosed in curly braces (e.g. ATOMS={1 2 3 4}).  Please note that you can split commands over multiple lines by using
\ref ContinuationLines. 

The most important of these keywords is the label keyword as it is only by using these labels that we can pass data 
from one action to another.  As an example if you do:

\plumedfile
DISTANCE ATOMS=1,2
\endplumedfile

Then PLUMED will do nothing other than read in your input file.  In contrast if you do:

\plumedfile
DISTANCE ATOMS=1,2 LABEL=d1
PRINT ARG=d1 FILE=colvar STRIDE=10
\endplumedfile

then PLUMED will print out the value of the distance between atoms 1 and 2 every 10 steps to the file colvar as you have told
PLUMED to take the value calculated by the action d1 and to print it. You can use any character string to label your actions
as long as it does not begin with the symbol \@.  Strings beginning with \@ are used by within PLUMED to reference special, 
code-generated groups of atoms and to give labels to any Actions for which the user does not provide a label in the input. 

Notice that if a word followed by a column is added at the beginning of the line (e.g. pippo:), PLUMED automatically
removes it and adds an equivalent label (LABEL=pippo).
Thus, a completely equivalent result can be obtained with the following shortcut:
\plumedfile
d1: DISTANCE ATOMS=1,2
PRINT ARG=d1 FILE=colvar STRIDE=10
\endplumedfile

Also notice that all the actions can be labeled, and that many actions besides normal collective variables can define
one or more value, which can be then referred using the corresponding label.

Actions can be referred also with POSIX regular expressions (see \ref Regex) if regex library is available on your system
and detected at configure time.
You can also add \ref comments to the input or set up your input over multiple files and then create a composite input by
\ref includes.

More information on the input syntax as well as details on the the various trajectory
analysis tools that come with PLUMED are given in: 

- \ref colvarintro tells you about the ways that you can calculate functions of the positions of the atoms.
- \ref Analysis tells you about the various forms of analysis you can run on trajectories using PLUMED.
- \ref Bias tells you about the methods that you can use to bias molecular dynamics simulations with PLUMED.

\section units Plumed units
By default the PLUMED inputs and outputs quantities in the following units:

- Energy - kJ/mol
- Length - nanometers
- Time - picoseconds

Unlike PLUMED 1 the units used are independent of the MD engine you are using.  If you want to change these units you can do this using the 
\subpage UNITS keyword. 

