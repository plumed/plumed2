Actions
-------

Every line in a PLUMED input file gives the input for an action or a [shortcut](shortcuts.md). The input for an action 
can all be on a single line as shown below:

```plumed
d1: DISTANCE ATOMS=1,2 COMPONENTS
```

Alternatively, you can split the input for an action over multiple lines by using a continuation as shown below:

```plumed
d1: DISTANCE ...
   ATOMS=1,2 
   COMPONENTS
...
```

The input for this action can also be written as follows:

```plumed
DISTANCE ...
   ATOMS=1,2
   COMPONENTS
   LABEL=d1
... DISTANCE
```

Notice that the closing `...` is followed by the word `DISTANCE` here. This might be
useful when it comes to matching the start and ends of multi-line statements.  However, 
it is important to note that PLUMED checks that the word following the closing `...` is identical to
the first word in the line with the first `...`. If these are not, it will throw an error. This is why we used 
`LABEL=d1` here in place of `d1: DISTANCE`. If you want to have matching action names at the start and end of your 
multiline statements you are perhaps best off doing:

```plumed
DISTANCE LABEL=d1 ...
   ATOMS=1,2
   COMPONENTS
... DISTANCE
```

so as to ensure that the label appears at the start of the multiline input.

## Adding comments

If you are an organized sort of person who likes to remember what the hell you were trying to do when you ran a
particular simulation you might find it useful to put comments in your input file.  In PLUMED you can do this as
comments can be added using a # sign.  On any given line everything after the # sign is ignored so
erm... yes add lines of comments or trailing comments to your hearts content as shown below (using Shakespeare is optional):

```plumed
# This is the distance between two atoms:
d1: DISTANCE ATOMS=1,2
Snout: UPPER_WALLS ARG=d1 AT=3.0 KAPPA=3.0   # In this same interlude it doth befall.
# That I, one Snout by name, present a wall.
```

PLUMED ignores any text in comments.  Consequently, if you provide the input for an action in a comment like this:

```plumed
# d1: DISTANCE ATOMS=1,2 COMPONENTS
```

a [DISTANCE](DISTANCE.md) command is not created. Similar behaviour is observed when you use the [ENDPLUMED](ENDPLUMED.md) comand.  For example, if your 
input is as follows:

```plumed
d1: DISTANCE ATOMS=1,2 COMPONENTS
ENDPLUMED
d2: DISTANCE ATOMS=3,4
```

PLUMED will evaluate the distance between atom 1 and 2 but will not evaluate the distance between atoms 3 and 4.

Lastly, note that you can include comments in the input for an action or [shortcut](shortcuts.md) that is split over multiple lines as shown below:

```plumed
dist: DISTANCE ...
# we can also insert comments here
  ATOMS1=1,300
# multiple kewords per line are allowed
  ATOMS2=1,400 ATOMS3=1,500
#empty lines are also allowed

... 
```

## Using INCLUDE files

If, for some reason, you want to spread your PLUMED input over a number of files you can use [INCLUDE](INCLUDE.md) as shown below:

````
INCLUDE FILE=filename
````

So, for example, instead of having a single "plumed.dat" file:

```plumed
DISTANCE ATOMS=1,2 LABEL=dist
RESTRAINT ARG=dist AT=2.0 KAPPA=1.0
```

you can split the input over a file like the one below and a file called `extras/toBeIncluded.inc` 

```plumed
DISTANCE ATOMS=1,2 LABEL=dist
INCLUDE FILE=extras/toBeIncluded.inc
```

However, when you do this it is important to recognize that INCLUDE is a real directive that is only resolved
after all the comments have been stripped and the continuation lines have been unrolled.  This means it
is not possible to do things like:

```plumed
# this is wrong:
DISTANCE INCLUDE FILE=options.dat
RESTRAINT ARG=dist AT=2.0 KAPPA=1.0
```

## Parallelism

PLUMED uses parallelism to improve the performance of many actions.  Within PLUMED parallelisation is done using a combination of
[MPI](https://en.wikipedia.org/wiki/Message_Passing_Interface) and [OpenMP](https://en.wikipedia.org/wiki/OpenMP).  In addition, 
we are starting to use [openACC](https://www.openacc.org) in some actions in order to write code than can run on GPUs.

Importantly, all parallism within PLUMED is implemented within actions.  __We do not parallelise over actions.__  Consequently, if you run 
this input with MPI or openMP:

```plumed
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=3,4
d3: DISTANCE ATOMS=5,6
PRINT ARG=d1,d2,d3 FILE=colvar
``` 

all the three distances will be calculated on "the same processor/thread."  If, however, you run the following input with MPI or openMP:

```plumed
d: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6
PRINT ARG=d FILE=colvar
```

the three distances are calculated by different processors/threads.

### MPI

Unless you disable it by using the option `--disable-mpi` during configuration, PLUMED will search for an MPI installation during configuration.
If an MPI installation is found then PLUMED will build with MPI when you run make.  To run a PLUMED [driver](driver.md) calculation with MPI you 
will run a command something like this:

```plumed
mpirun -np 4 plumed driver --ixyz traj.xyz
``` 

When this command is run, any actions that can use MPI will run on four MPI processes.

### OpenMP

If your compiler supports openMP then the features of PLUMED that use it have it enabled by default.  However, if you want to disable 
these features you can do so by configuring PLUMED with the option  `--disable-openmp`.

When you use PLUMED with some MD codes the number of OpenMD threads is set automatically by the MD code. If, however, you are using [driver](driver.md),
another [command line tool](module_cltools.md) or an MD code that does not set the number of OpenMD threads to be used you will need to set the environment 
variable PLUMED\_NUM\_THREADS equal to the number of threads you wish PLUMED to use at runtime.  To be clear, you can use the PLUMED\_NUM\_THREADS to set 
the number of threads to use with PLUMED even if your MD code has an option for setting the number of threads.  With gromacs for instance you can ensure
that PLUMED runs on 8 threads by using the following commands: 

```bash
export PLUMED_NUM_THREADS=8
mdrun -plumed
```

or by using this single command:

```bash
mdrun -plumed -ntomp 8
```

In the first case PLLUMED uses 8 OpenMP threads while gromacs only uses 1 (this is usually sub optimal).
In the second case GROMACS and plumed will the same number of OpenMP threads.

Notice that:

- Using OpenMP is likely to improve the performance, but could also slow down the code in some case.
- Using OpenMP can give results that are slightly different because of numerical round off and different order in summations. This should be harmless.
- The optimum number of threads is not necessary "all of them", nor should be equal to the number of threads used to parallelize the MD calculation. 
- Not all CVs are parallelized with openMP. 
- You might also want to tune the environmental variable PLUMED\_CACHELINE\_SIZE.
  The default of 512 is the size of cache lines on your machine. The PLUMED\_CACHELINE\_SIZE variable is used
  by PLUMED to decrease the number of threads to be used in each loop so as to avoid clashes in memory access. This variable is expected to affect performance only, not results.




