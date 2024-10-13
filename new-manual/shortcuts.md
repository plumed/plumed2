Shortcuts
---------

Many of the commands in a PLUMED input file are shortcuts. When PLUMED reads one of these shortcuts replaces it with multiple lines of PLUMED input and then reads 
in all the [actions](actions.md) that are defiend in this longer and more complicated input.  The reason for implementing these commands in this way is twofold:

1. It ensures that there is less code within PLUMED to maintain.
2. It allows users of PLUMED to understand how various methods are implemented within PLUMED.   

We believe that the second item in this list is a particularly important objective.  We cannot guarantee that PLUMED contains the fastest implementations of these methods 
on all the various architectures that are available to the modern simulator.  We hope, however, that users can read the manual, understand our shortcut implementations and use that understanding when implementing
faster versions of these methods that may be required to run a particular system with a particular computer architecture. In this way we hope that PLUMED serves as a textbook about 
the various methods that are implemented within it as well as a simulation code. 

## Visualising how shortcuts work

You can see an example of a shortcut in the following input:

```plumed
c1: COORDINATIONNUMBER SPECIES=1-100 SWITCH={RATIONAL R_0=0.2} MEAN
```

Whenever shortcuts are used in example inputs in this manual, the [nest](www.plumed-nest.org) or the [tutorials](www.plumed-tutorials.org) site you have the option to expand the shortcut and see all the various constituent actions
that are used in the expanded input that evlauates the final quantity of interest. To expand the input you simply hover over the name of the action and use the approriate link
from the tooltip that appears. 

## Values from shortcuts

In the manual pages for shortcuts you will see that these commands define [values](specifying_arguments.md) just like [actions](actions.md).  The values calculated by shortcuts that are 
specified in the manual can be referenced in other parts of the PLUMED input in the [usual ways](specifying_arguments.md).  Notice, however, that the name of the value that appears in the
output file may differ from the name of the component.  For example, if you use the following input:

```plumed
# This is a shortcut
d1: DISTANCES ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8 ATOMS5=9,10 LESS_THAN={RATIONAL R_0=0.1}
# This is not a shortcut
PRINT ARG=d1.lessthan FILE=colvar1
``` 

The output values in the output colvar1 file will be in a column with the heading `d1_lessthan` rather than `d1.lessthan` as `d1_lessthan` is the label of the action that actually calculates the 
quantity of interest.

PLUMED only allows you access the subset of the values that defined in the manual from the shortcut input using the `label.component` syntax.  So if we take the input above:

```plumed 
d1: DISTANCES ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8 ATOMS5=9,10 LESS_THAN={RATIONAL R_0=0.1}

# You can output d1.lessthan like this
PRINT ARG=d1.lessthan FILE=colvar1

# or you can output it like this
PRINT ARG=d1_lessthan FILE=colvar2

# You cannot output d1_lt like this as the lt component is not defined in the manual
# PRINT ARG=d1.lt FILE=colvar3
# But you can output it like this
PRINT ARG=d1_lt FILE=colvar4
```

## Shortcuts and wildcards

If you have the following input:

```plumed 
d1: DISTANCES ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8 ATOMS5=9,10 LESS_THAN={RATIONAL R_0=0.1} MEAN
PRINT ARG=d1.* FILE=colvar2
```

The values `d1_lessthan` and `d1_mean` will be output in the colvar2 file as these are the components of the DISTANCES shortcut that are defined in the manual. `d1_lt` is not output as this is an intermediate value that is not defined in the manual.

The same result can be obtained using the following input:

```plumed 
d1: DISTANCES ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8 ATOMS5=9,10 LESS_THAN={RATIONAL R_0=0.1} MEAN
PRINT ARG=* FILE=colvar
```

as the * wildcard instructs PLUMED to output values from all the shortcuts in the input file that are described in the manual.  However, if you use a *.* wildcard everything (including the value of `d1_lt` and all the other intermediate values) will 
be output. 
