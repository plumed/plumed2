Actions
-------

Every command in a PLUMED input file gives the input for an action or a [shortcut](shortcuts.md). The input for an action 
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
