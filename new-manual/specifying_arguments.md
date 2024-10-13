Specifying arguments
--------------------

Many of the actions in PLUMED output one or more values that can be used in the input to other actions.  These values can be
scalars, vectors, matrices or functions evaluated on a grid.

## Single component actions

If an action outputs only one value it is referenced later in the input using the label of the action that calculated it.
For example, in the input below, the DISTANCE action calculates one scalar every MD time step.  This scalar is then reused 
by the PRINT action, which the time-series of scalars to a file called colvar1.

```plumed
d1: DISTANCE ATOMS=1,2
PRINT ARG=d1 FILE=colvar1
```

This next input does something similar. However, the DISTANCE action now calculates the distances between 4 separate pairs of atoms.
The PRINT action will thus output a vector with four elements to the colvar2 file every MD time step.

```plumed
d2: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8
PRINT ARG=d2 FILE=colvar2
```

Lastly, consider the following input, the first here calculates a 10x10 matrix of switching functions that act on the distances between
each pair of atoms.  The PRINT action then outputs the 100 elements of this vector to the colvar3 file on every MD time step.

```plumed
c1: CONTACT_MATRIX GROUP=1-10 SWITCH={RATIONAL R_0=0.2}
PRINT ARG=c1 FILE=colvar3
```

If the label refers to a function evaluated on a grid you can output the data using PRINT. However, we would recommend using DUMPGRID instead in these cases.

## Multi-component actions

If an action outputs multiple values the component values are referenced later in the input using the string that contains the label of the action 
followed by a dot and the required component name. The action pages in the manual will tell you the names of the component values that are output
by each action.  

The following example input illustrates this idea.  

```plumed
d3: DISTANCE ATOMS=1,2 COMPONENTS
PRINT ARG=d3.x,d3.y.d3.z FILE=colvar
```

For this input the DISTANCE action now outputs 3 values because the COMPONENTS keyword is present. This flag tells PLUMED that you would like to 
evaluate the vector connecting the two atoms rather than the modulus of the vector.  The PRINT command thus outputs

- d3.x - a scalar that tells us the x component of the vector connecting atoms 1 and 2 
- d3.y - a scalar that tells us the y component of the vector connecting atoms 1 and 2
- d3.z - a scalar that tells us the z component of the vector connecting atoms 1 and 2

When actions output multiple values these values can be vectors rather than scalars as illustrated by the input below:

```plumed
d4: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8 COMPONENTS
PRINT ARG=d4.x,d4.y,d4.z FILE=colvar
```

The COMPONENTS flag again tells PLUMED that you would like to evaluate the vectors connecting each pair of atoms.  The PRINT command thus outputs:

- d4.x - a vector with four components that tells us the x components of the vectors connecting each pair of atoms
- d4.y - a vector with four components that tells us the y components of the vectors connecting each pair of atoms
- d4.z - a vector with four components that tells us the z components of the vectors connecting each pair of atoms

Actions that output multiple components that are all matrices also exist and work in a similar way.

## Wildcards

When you have many values in the input to a PLUMED action it can be useful to use wildcards.  As a first example, you can create 
an input that is equivalent to this one:

```plumed
d3: DISTANCE ATOMS=1,2 COMPONENTS
PRINT ARG=d3.x,d3.y.d3.z FILE=colvar
```

by using the * widlcard to print all the values output by the action with label d3 in the PRINT action as follows:

```plumed
d3: DISTANCE ATOMS=1,2 COMPONENTS
PRINT ARG=d3.* FILE=colvar
```

Notice that the * will also work in the same way if the values that are output by the action with label d3 are vectors, matrices or functions on grids.

This second example uses a wildcard to output all the x components for the DISTANCE commands:

```plumed
dd1: DISTANCE ATOMS=1,2 COMPONENTS
dd2: DISTANCE ATOMS=3,4 COMPONENTS
dd3: DISTANCE ATOMS=5,6 COMPONENTS
dd4: DISTANCE ATOMS=7,8 COMPONENTS
PRINT ARG=*.x FILE=colvar
```

We would note, however, that a similar result can be produced by using the following input:

```plumed
dd1: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8 COMPONENTS
PRINT ARG=dd1.x FILE=colvar
```

and that in this second case the calculation of the distances will be parallelized.  The `*.component` wildcard is nevertheless still useful in cases where
you want to use all the simulation bias that are acting upon the system in the input to an action.

If you want all the values in the PLUMED input to be used in the input for an action you can do this:

```plumed
dd1: DISTANCE ATOMS=1,2 COMPONENTS
PRINT ARG=* FILE=colvar
```

or this:

```plumed
dd1: DISTANCE ATOMS=1,2 COMPONENTS
PRINT ARG=*.* FILE=colvar
```

__We recommend not using the * and *.* wildcards as there are often many values defined in a PLUMED input file and you will rarely want to use all them in the input to an action.  
Furthermore, if you are using [shortcuts](shortcuts.md) you may not even be aware of all the values that are defined in your input as your shortcut actions work by creating intermediate actions 
that calculate intermediate values.__

## Regular expressions

Another way to use many values in the input to a PLUMED action is to use [regular expressions](https://en.wikipedia.org/wiki/Regular_expression).
which are enabled when PLUMED is compiled with the regex library.  To learn more about regular expression you can read the following [site](http://www.regular-expressions.info/reference.html).

In PLUMED Regular expressions are enclosed in round braces and must not contain spaces.  You can duplicate the `label.*` and `*.component` wildcards that were described in the previous section
by using regular expressions as shown below:

```plumed
dd1: DISTANCE ATOMS=1,2 COMPONENTS
dd2: DISTANCE ATOMS=3,4 COMPONENTS
dd3: DISTANCE ATOMS=5,6 COMPONENTS
dd4: DISTANCE ATOMS=7,8 COMPONENTS

# This is equivalent to dd1.* in regex
PRINT ARG=(dd1\..*) FILE=colvar1

# And this equivalent to *.x in regex
PRINT ARG=(.*\.x) FILE=colvar2
```

Perhaps, more usefully another input that used regular expressions is as follows:

```plumed
d1: DISTANCE ATOMS=1,2 COMPONENTS
PRINT ARG=(d1\.[xy]) FILE=colvar5 
```

The PRINT command here will output the d1.x and d1.y components of the vector connecting atoms 1 and 2 to the colvar5 file. Notice that in your regular expressions the
`.` character must be escaped as `\.` in order to interpret it as a literal `.`. An un-escaped dot is a wildcard which is matched by any character,
so for example

```plumed
d1: DISTANCE ATOMS=1,2 COMPONENTS
dxy: DISTANCE ATOMS=1,3

# this will match d1.x,d1.y,dxy
PRINT ARG=(d1.[xy]) FILE=colvar

# while this will match d1.x,d1.y only
PRINT ARG=(d1\.[xy]) FILE=colvar
```

You can concatenate more than one regular expression by using comma separated regular expressions as illustrated below

```plumed
t1: TORSION ATOMS=5,7,9,15
t2: TORSION ATOMS=7,9,15,17
d1: DISTANCE ATOMS=7,17 COMPONENTS

# The first expression matches d1.x and d1.y
# The second expression matches t1 and t2
PRINT ARG=(d1\.[xy]),(t[0-9]) FILE=colvar 
# So the command here is the same as ARG=d1.x,d1.y,t1,t2
```

It is perhaps better, however, to the use "or" operator `|` in a single regular expression as shown below:

```plumed
t1: TORSION ATOMS=5,7,9,15
t2: TORSION ATOMS=7,9,15,17
d1: DISTANCE ATOMS=7,17 COMPONENTS

# Here is a single regular expression
PRINT ARG=(d1\.[xy]|t[0-9]) FILE=colvar
# Thus this is the same as ARG=t1,t2,d1.x,d1.y
```

Constructing the regular expression this way is preferable as each value will only be generated once in the input. In the first example with the comma
separated expressions, duplicate values can can appear if there are values that match more than one of the constituent regular expressions.

**Please note that [shortcuts](shortcuts.md) work by creating additional actions to calculate intermediate values. You need to take care when you are designing regular expressions 
for inputs that use shortcuts to not have regular expressions that match the labels of the intermediate values.**

## Reshaping values

There may be cases where it is useful to access a particular element of a vector or matrix or to combine multiple scalars or vectors into a new value.  The following
actions allow you to complete operations of this type:

- SELECT_COMPONENTS - take $n$ of the elements from the input vector/matrix and output an $n$-dimensional vector that contains these components only.
- CONCATENATE - combine all the input scalar/vectors/matrices into a single output vector/matrix.
- FLATTEN - vectorise a matrix
- VSTACK - construct a matrix by stacking multiple input vectors together 

You can then use these outputs from these actions in the input to later actions in your input file.


