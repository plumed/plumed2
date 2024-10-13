Specifying atoms
-----------------

Many of the actions in PLUMED take a list of atom positions in input.  Within PLUMED
atoms are specified using their numerical indices in the molecular dynamics input file.

In PLUMED lists of atoms can be either provided directly inside the definition of each action, or
predefined as a GROUP that can be reused multiple times. Lists of atoms can be written as:

- comma separated lists of numbers (`g1: GROUP ATOMS=10,11,15,20`)
- numerical ranges.  So `g2: GROUP ATOMS=10-20` is equivalent to `g2: GROUP ATOMS=10,11,12,13,14,15,16,17,18,19,20`
- numerical ranges with a stride. So `g3: GROUP ATOMS=10-100:10` is equivalent to `g3: GROUP ATOMS=10,20,30,40,50,60,70,80,90,100`
- atom ranges with a negative stride. So `g4: GROUP ATOMS=100-10:-10` is equivalent to `g4: GROUP ATOMS=100,90,80,70,60,50,40,30,20,10`

If you want to use the atoms in a group as the input for a action you use the label of the group as shown in the following input:

```plumed
g5: GROUP ATOMS=1,2
d1: DISTANCE ATOMS=g5
```

A few other wasys of using groups in the input to actions are also available:

- `@mdatoms` indicate all the physical atoms present in the MD engine (e.g. `DUMPATOMS ATOMS=@mdatoms`).
- `@allatoms` indicates all atoms, including the virtual atoms that are only defined in PLUMED (e.g. `DUMPATOMS ATOMS=@allatoms`).
- `@ndx` uses a [GROMACS index file](https://manual.gromacs.org/archive/5.0.2/online/ndx.html). `@ndx:index.ndx` picks the first group in the file. `{@ndx:{index.ndx protein}}` picks
   the group named `protein`.

## The MOLINFO command

You can access many useful shortcuts for specifying atomic positions to PLUMED by using the MOLINFO, which takes 
a PDB structure file in input. This command allows you to access the following list of shortcuts.

{:#browse-table .display}
| Name | Description |
|:---------|:---------|
{% for item in site.data.grouplist %}| {{ item.name }} | {{ item.description }} |
{% endfor %}

## Virtual atoms 

Sometimes, when calculating a collective variable, you may not want to use the positions of a number of atoms directly. Instead
 you may wish to use the position of a virtual atom whose position is generated based on the positions of a collection
of other atoms.  For example you might want to use the center of mass of a group of atoms. PLUMED has a number of routines
for calculating the positions of these virtual atoms from lists of atoms that are in the vatom module.

To specify to an action that you want to use the position of a virtual atom to calculate an action rather than one of the atoms
in your system you simply use the label for your virtual atom in place of the usual numerical index. Virtual
atoms and normal atoms can be mixed together in the input to actions as shown below:

```plumed
com1: COM ATOMS=1,10 
d1: DISTANCE ATOMS=11,com1
```

If you don't want to calculate CVs from the virtual atom.  That is to say you just want to monitor the position of a virtual atom
(or any set of atoms) over the course of your trajectory you can do this using DUMPATOMS.

The list of the virtual atoms defined in PLUMED can be obtained by using the command `GROUP ATOMS=@allatoms REMOVE=@mdatoms`.

## Broken molecules and Periodic Boundary Conditions

PLUMED is designed so that for the majority of the CVs implemented the periodic boundary conditions are treated
in the same manner as they would be treated in the host code.  In some codes this can be problematic when the actions 
you are using involve some property of a molecule.  These codes allow the atoms in the molecules to become separated by
periodic boundaries, a fact which PLUMED could only deal with if the topology is passed from the MD code to PLUMED.  Doing this
work would involve a lot laborious programming and goes against our original aim of having a general patch that can be implemented
in a wide variety of MD codes.  Consequentially, we have implemented a more pragmatic solution to this problem - the user specifies
in input any molecules (or parts of molecules) that must be kept in tact throughout the simulation run using the WHOLEMOLECULES command.

The following input computes the end-to-end distance for a polymer of 100 atoms and keeps it at a value around 5.

```plumed
WHOLEMOLECULES ENTITY0=1-100
e2e: DISTANCE ATOMS=1,100 NOPBC
RESTRAINT ARG=e2e KAPPA=1 AT=5
```

Notice that NOPBC is used to in the DISTANCE action so as to ensure that if the end-to-end distance is larger than half the simulation box the distance
is compute properly. Also notice that, since many MD codes break molecules across cell boundary, it might be necessary to use the
WHOLEMOLECULES keyword (also notice that it should be before distance).

Notice that most expressions are invariant with respect to a change in the order of the atoms,
but some of them depend on that order. E.g., with WHOLEMOLECULES it could be useful to
specify atom lists in a reversed order.

```plumed
# to see the effect, one could dump the atoms as they were before molecule reconstruction:
# DUMPATOMS FILE=dump-broken.xyz ATOMS=1-20
WHOLEMOLECULES STRIDE=1 ENTITY0=1-20
DUMPATOMS FILE=dump.xyz ATOMS=1-20
```

Notice that there are other ways to manipulate the coordinates stored within PLUMED:

- FIT_TO_TEMPLATE aligns atoms to a template structure.
- WRAPAROUND brings a set of atom as close as possible to another set of atoms.
- RESET_CELL rotates the periodic cell.

<script>
$(document).ready(function() {
var table = $('#browse-table').DataTable({
  "dom": '<"search"f><"top"il>rt<"bottom"Bp><"clear">',
  language: { search: '', searchPlaceholder: "Search project..." },
  buttons: [
        'copy', 'excel', 'pdf'
  ],
  "order": [[ 0, "desc" ]]
  });
$('#browse-table-searchbar').keyup(function () {
  table.search( this.value ).draw();
  });
  hu = window.location.search.substring(1);
  searchfor = hu.split("=");
  if( searchfor[0]=="search" ) {
      table.search( searchfor[1] ).draw();
  }
});
</script>
