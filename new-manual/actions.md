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

a DISTANCE command is not created. Similar behaviour is observed when you use the ENDPLUMED comand.  For example, if your 
input is as follows:

```plumed
d1: DISTANCE ATOMS=1,2 COMPONENTS
ENDPLUMED
d2: DISTANCE ATOMS=3,4
```

PLUMED will evaluate the distance between atom 1 and 2 but will not evaluate the distance between atoms 3 and 4.
