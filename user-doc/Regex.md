\page Regex Regular Expressions

When you use need to pass many arguments to a PLUMED action, being them
components of a few collective variables or also multiple collective variables,
you might find it convenient to use [regular expressions](https://en.wikipedia.org/wiki/Regular_expression).

Since version 2.1, plumed takes advantage of a configuration scripts that
detects libraries installed on your system. If regex library is found,
then you will be able to use regular expressions to refer to collective variables
or function names.

Regular expressions are enclosed in round braces and must not contain spaces (the components 
names have no spaces indeed, so why use them?).

As an example the command:
\plumedfile
d1: DISTANCE ATOMS=1,2 COMPONENTS
PRINT ARG=(d1\.[xy])   STRIDE=100 FILE=colvar FMT=%8.4f
\endplumedfile
will cause both the d1.x and d1.y components of the DISTANCE action to be printed.

Notice that selection does not happen in alphabetic order, nor in the order in which `[xy]` are listed, but rather in the order in which
the two variables have been created by PLUMED.
Also notice that the
`.` character must be escaped as `\.` in order to interpret it as a literal `.`. An un-escaped dot is a wildcard which is matched by any character,
So as an example
\plumedfile
d1: DISTANCE ATOMS=1,2 COMPONENTS
dxy: DISTANCE ATOMS=1,3

# this will match d1.x,d1.y,dxy
PRINT ARG=(d1.[xy])   STRIDE=100 FILE=colvar FMT=%8.4f

# while this will match d1.x,d1.y only
PRINT ARG=(d1\.[xy])   STRIDE=100 FILE=colvar FMT=%8.4f
\endplumedfile

You can concatenate more than one regular expression by using comma separated regular expressions.
The resulting matches will be concatenated:
\plumedfile
t1: TORSION ATOMS=5,7,9,15
t2: TORSION ATOMS=7,9,15,17
d1: DISTANCE ATOMS=7,17 COMPONENTS

# The first expression matches d1.x and d1.y
# The second expression matches t1 and t2
PRINT ARG=(d1\.[xy]),(t[0-9]) STRIDE=100 FILE=colvar FMT=%8.4f
# Thus this is the same as ARG=d1.x,d1.y,t1,t2
\endplumedfile

Be aware that if you have overlapping selections they will be duplicated.
As an alternative you could use the "or" operator `|`:
\plumedfile
t1: TORSION ATOMS=5,7,9,15
t2: TORSION ATOMS=7,9,15,17
d1: DISTANCE ATOMS=7,17 COMPONENTS

# Here is a single regular expression
PRINT ARG=(d1\.[xy]|t[0-9]) STRIDE=100 FILE=colvar FMT=%8.4f
# Thus this is the same as ARG=t1,t2,d1.x,d1.y
\endplumedfile

this selects the same set of arguments as the previous example.

\note
Be careful you do not confuse regular expressions, which are triggered by the parenthesis `()` and only available when
PLUMED has been compiled with the regex library, with the capability of PLUMED to use `*` as a wildcard in arguments:
\plumedfile
d1: DISTANCE ATOMS=1,2 COMPONENTS
# this is a regular expression that selects all components of d1
# i.e. d1.x d1.y and d1.z
PRINT ARG=(d1\..*)   STRIDE=100 FILE=colvar_reg FMT=%8.4f

# this is a wildcard that selects all the components of d1 as well
PRINT ARG=d1.*   STRIDE=100 FILE=colvar_wild FMT=%8.4f
\endplumedfile
Regular expressions are way more flexible than wildcards!

You can check the log to see whether or not your regular expression is picking the set of components you desire.

For more information on regular expressions visit http://www.regular-expressions.info/reference.html.

