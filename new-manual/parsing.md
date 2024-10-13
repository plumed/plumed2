Parsing action input 
--------------------

The input to any PLUMED [action](actions.md) consists of a series of keyword value pairs and flags that turn on or off various options.  
In the example input below:

```plumed
d1: DISTANCE ATOMS=1,2 NOPBC
```

`d1` is the action's [label](specifying_arguments.md) and `DISTANCE` is the action's name.  The input to the action is then the keyword value pair ATOMS and the flag NOPBC.
Every example input you encounter in this manual, in the [nest](www.plumed-nest.org) and the [tutorials](www.plumed-tutorials.org) has tooltips
that explain what options are turned on by flags and what the values provided in each keyword pair represents.  You will also find details of
all the keywords and flags that are available for an action on the manual page for that action.

## Reading constants

In the input to keywords that read real numbers you can use constants that are specified using strings rather than numbers.
An example is the following

```plumed
#SETTINGS MOLFILE=regtest/basic/rt65/AA.pdb
MOLINFO STRUCTURE=AA.pdb  MOLTYPE=rna
e1: TORSION ATOMS=@epsilon-1
t: METAD ARG=e1 SIGMA=0.15 PACE=10 HEIGHT=2 GRID_MIN=-pi GRID_MAX=pi GRID_BIN=200
```

Notice here that the boundaries for `GRID_MIN` and `GRID_MAX` are `-pi` and `pi`. 

Any real numbers that are read in from input are interpretted using the [Lepton library](Custom.md) so you can thus employ
complicated expressions such as `1+2` or `exp(10)` as shown in the in the input as shown below:

```plumed
#SETTINGS MOLFILE=regtest/basic/rt65/AA.pdb
MOLINFO STRUCTURE=AA.pdb  MOLTYPE=rna
e1: TORSION ATOMS=@epsilon-1
RESTRAINT ARG=e1 AT=1+0.5
```

The lepton library is able to interpret any of thse following constants:

| Symbol | Description | 
| :----- |:------------|
| `e` | Euler's number - the base for the natural logarithm |
| `log2e` | $1 / \log(2)$ |
| `log10e` | $1 / \log(10)$ |
| `ln2` | $\log(2)$ |
| `ln10` | $\log(10)$ |
| `pi`   | the circle constant $\pi$ |
| `pi_2` | $\pi / 2$ |
| `pi_4` | $\pi / 4$ |
| `sqrt2 | $\sqrt(2)$ |
| `sqrt1_2 ` | $\sqrt(0.5)$ |

Notice that this functionality cannot be used when the keyword takes integer numbers in input 
(e.g.: the PACE argument for METAD).

## Special replica syntax

PLUMED provides a number of ways to prepare the multiple replicas with almost identical PLUMED files that are required to run a multiple replica simulation:

* You can repare input files for such calculation using cut-and-paste but that is is very error prone.
* You can write a smart bash or python script to generate all the inputs.
* You can use different inputs for the various replicas that all contain an INCLUDE directive to include a file that contains the parts of the input that are common to all replicas 

We think, however, the best option is to use features that have been available from PLUMED 2.4 onwards that allow you 
manipulate multiple replica inputs that have only tiny differences in the the input.  The following example illustrates how the syntax 
for this option operates:  

```plumed
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
````

If you prepare a single `plumed.dat` file like this one and feeds it to PLUMED while using 3 replicas,
the 3 replicas will see PLUMED input files that are the same except for the `AT` keyword, that sets the position of the restraint.
Replica 0 will see a restraint centered at 1.0, replica 1 centered at 1.1, and replica 2 centered at 1.2.

The `@replicas:` keyword is not specific to RESTRAINT and the `AT` keyword. Any keyword in PLUMED can accept that syntax.
For instance, the following single input file can be used to setup a bias exchange metadynamics simulations:

```plumed
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
```

This input contains a typical setup for a bias exchange simulation.
Notice that even though the actions with labels `d` and `t` are read for both replicas,
`d` is only computed on replica 0 (and `t` is only computed on replica 1).
This is because variables that are defined but not used are never actually calculated by PLUMED.

If the value that should be provided for each replica is a vector, you should use curly braces as delimiters.
For instance, if you have a restraint that acts on two variables, you can use the following input:

```plumed
#SETTINGS NREPLICAS=3
# Compute distance between atoms 1 and 2
d: DISTANCE ATOMS=10,20

# Compute a torsional angle
t: TORSION ATOMS=30,31,32,33

# Apply a restraint:
# RESTRAINT ...
#  ARG=d,t
#  AT=@replicas:{{1.0,2.0} {3.0,4.0} {5.0,6.0}}
#  KAPPA=1.0,3.0
# ...
# On replica 0 this means:
#  RESTRAINT ARG=d AT=1.0,2.0 KAPPA=1.0,3.0
# On replica 1 this means:
#  RESTRAINT ARG=d AT=3.0,4.0 KAPPA=1.0,3.0
# On replica 2 this means:
#  RESTRAINT ARG=d AT=5.0,6.0 KAPPA=1.0,3.0
```

Notice the double curly braces. The outer ones are used by PLUMED to know where the argument of the `AT` keyword ends,
whereas the inner ones are used to group the values corresponding to each replica.
Also notice that the last example can be split in multiple lines exploiting the fact that
within multi-line statements (enclosed by pairs of `...`) newlines are replaced with simple spaces:

```plumed
#SETTINGS NREPLICAS=3
d: DISTANCE ATOMS=10,20
t: TORSION ATOMS=30,31,32,33
# RESTRAINT ...
#  ARG=d,t
# indentation is not required (this is not python!)
# but makes the input easier to read
# NEED TO WORK ON THIS INPUT WITH PLUMED2HTML
#  AT=@replicas:{
#    {1.0,2.0}
#    {3.0,4.0}
#    {5.0,6.0}
#  }
#  KAPPA=1.0,3.0
# ...
```

In short, whenever there are keywords that should vary across replicas, you should set them using the `@replicas:` keyword.
As mentioned above, you can always use the old syntax with separate input files and this is fact recommended when the
differences in the inputs for the various replicas are substantial.

