__THIS MODULE HAS BEEN INCLUDED FOR BACK COMPATIBILITY.  ALL THE ACTIONS WITHIN IT ARE NOW SHORTCUTS.  IF YOU USE THESE 
FUNCTIONALITIES WE WOULD STRONGLY ENCOURAGE YOU TO LOOK AT THE EXPANDED VERSIONS OF THESE INPUTS AND TO LEARN HOW THE NEWER
AND MORE POWERFUL SYNTAX OPERATES.  THE DOCUMENTATION THAT FOLLOWS IS COPIED FROM THE OLD MANUAL WHERE THIS NEW SYNTAX WAS NOT 
AVAILABLE__

Oftentimes, when you do not need one of the collective variables described elsewhere in the manual, what you want instead is a
function of a distribution of collective variables of a particular type.  In other words, you would like to calculate a
function something like this:

$$
s = \sum_i g[f(\{X\}_i)]
$$

In this expression $g$ is a function that takes in one argument and $f$ is a function that takes a set of atomic positions
as argument. The symbol $\{X\}_i$ is used to indicate the fact that the function $f$ is evaluated for a number of different
sets of atoms.  If you would just like to output the values of all the various $f$ functions you can use the command [DUMPMULTICOLVAR](DUMPMULTICOLVAR.md)

This functionality is useful if you need to calculate a minimum distance or the number of coordination numbers greater than a 3.0.
To avoid duplicating the code to calculate an angle or distance many times and to make it easier to implement very complex collective
variables PLUMED provides these sort of collective variables using so-called MultiColvars.  MultiColvars are named in this way because a single
PLUMED action can be used to calculate a number of different collective variables.  For instance the [DISTANCES](DISTANCES.md)
action can be used to calculate the minimum distance, the number of distances less than a certain value, the number of
distances within a certain range... A more detailed introduction to multicolvars is provided in this
<a href="http://www.youtube.com/watch?v=iDvZmbWE5ps">10-minute video</a>. 

To instruct PLUMED to calculate a multicolvar you give an instruction that looks something like this:

````
NAME <atoms involved> <parameters> <what am I calculating> TOL=0.001 LABEL=label
````

Oftentimes the simplest way to specify the atoms involved is to use multiple instances of the ATOMS keyword
i.e. ATOMS1, ATOMS2, ATOMS3,...  Separate instances of the quantity specified by NAME are then calculated for
each of the sets of atoms.  For example if the command issued contains the following:

```plumed
d: DISTANCES ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6
```

The distances between atoms 1 and 2, atoms 3 and 4, and atoms 5 and 6 are calculated. Obviously, generating
this sort of input is rather tedious so short cuts are also available many of the collective variables.
These are described on the manual pages for the actions.

After specifying the atoms involved you sometimes need to specify some parameters that required in the
calculation.  For instance, for [COORDINATIONNUMBER](COORDINATIONNUMBER.md) - the number of atoms in the first coordination
sphere of each of the atoms in the system - you need to specify the parameters for a switching function (see [LESS_THAN](LESS_THAN.md))
that will tell us whether or not an atom is in the first coordination sphere.  Details as to how to do this
are provided on the manual pages.

Once you have specified the base quantities that are to be calculated from the atoms involved and any parameters
you need to specify what function of these base quantities is to be calculated.  For most multicolvars you can calculate
the minimum, the number less than a target value, the number within a certain range, the number more than a target
value and the average value directly.

## Multicolvar functions and biases

In older versions of PLUMED there were actions that took the label of a multicolvar as input and computed things from the vector
of CVs directly from them. In newer versions this functionality is easier to use as the multicolvars are shortcuts.  You thus use 
the outputs from the commmands that output the vectors of CV values in the input to other actions directly. 
