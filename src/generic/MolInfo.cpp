/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "core/ActionRegister.h"
#include "core/GenericMolInfo.h"

namespace PLMD {
namespace generic {

//+PLUMEDOC TOPOLOGY MOLINFO
/*
This command is used to provide information on the molecules that are present in your system.
A simple examples where this command is used to provide the information on which atoms
are in the backbone of a protein to the [ALPHARMSD](ALPHARMSD.md) action is shown below.

```plumed
#SETTINGS MOLFILE=regtest/basic/rt32/helix.pdb
MOLINFO STRUCTURE=regtest/basic/rt32/helix.pdb
ALPHARMSD RESIDUES=all TYPE=DRMSD LESS_THAN={RATIONAL R_0=0.08 NN=8 MM=12} LABEL=a
```

The following example prints the distance corresponding to the hydrogen bonds
in a GC Watson-Crick pair.

```plumed
#SETTINGS MOLFILE=regtest/basic/rt-ermsd/ref.pdb
MOLINFO STRUCTURE=regtest/basic/rt-ermsd/ref.pdb MOLTYPE=dna
hb1: DISTANCE ATOMS=@N2-2,@O2-15
hb2: DISTANCE ATOMS=@N1-2,@N3-15
hb3: DISTANCE ATOMS=@O6-2,@N4-15
PRINT ARG=hb1,hb2,hb3
```

The last example shows you how you can use MOLINFO to calculate torsion angles

```plumed
#SETTINGS MOLFILE=regtest/basic/rt32/helix.pdb
MOLINFO MOLTYPE=protein STRUCTURE=regtest/basic/rt32/helix.pdb
t1: TORSION ATOMS=@phi-3
t2: TORSION ATOMS=@psi-4
PRINT ARG=t1,t2 FILE=colvar STRIDE=10
```

The information on the molecules in your system can either be provided in the form of a pdb file
or as a set of lists of atoms that describe the various chains in your system. If a pdb file
is used plumed the MOLINFO command will endeavor to recognize the various chains and residues that
make up the molecules in your system using the chainIDs and resnumbers from the pdb file. You can
then use this information in later commands to specify atom lists in terms residues.  For example
using this command you can find the backbone atoms in your structure automatically.
Since the introduction of PLUMED 2.7 you have been able to use multiple MOLINFO actions. Every time you perform an atom
selection, the last available MOLINFO action will be used. This allows you to provide multiple PDB
files, for instance using different naming conventions (see issue 134).

!!! caution "check PDB parser"

    Please be aware that the PDB parser in plumed is far from perfect. You should thus check the log file
    and examine what plumed is actually doing whenever you use the MOLINFO action.
    Also make sure that the atoms are listed in the pdb with the correct order.
    If you are using gromacs, the safest way is to use reference pdb file
    generated with `gmx editconf -f topol.tpr -o reference.pdb`.

More information of the PDB parser implemented in PLUMED can be found in the documentation for the [PDB2CONSTANT](PDB2CONSTANT.md) command.

If the flag `WHOLE` is used, the reference structure is assumed to be whole, i.e. not broken by PBC.
**This will impact PBC calculations in all the actions following the MOLINFO action!**.
In particular, actions that reconstruct molecules locally will do the reconstruction using a
minimum spanning tree rather than the order in which atoms are provided.
Notice that in PLUMED 2.8, 2.9, and 2.10 this would only have affected the behavior of \ref WHOLEMOLECULES
actions, whereas starting with PLUMED 2.11 all actions that reconstruct PBCs are affected.
If you want a variable to ignore this, you can just repeate the MOLINFO action without
the `WHOLE` flag: all the following actions will use the traditional reconstruction algorithm,
based on the order in which atoms are listed.

Providing `MOLTYPE=protein`, `MOLTYPE=rna`, or `MOLTYPE=dna` will instruct plumed to look
for known residues from these three types of molecule. In other words, this keyword is available for
historical reasons and to allow future extensions where alternative lists will be provided.
As of now, you can just ignore this keyword.

Using MOLINFO extends the possibility of atoms selection using the @ special
symbol. All the various ways that the @ symbol can be used are listed in the following table:

@MOLINFOGROUPS@

In the table above you will find the following shortcuts that allow you to select atoms that are not part of one specific residue:

````
@nucleic : all atoms that are part of a DNA or RNA molecule
@protein : all atoms that are part of a protein
@water : all water molecules
@ions : all the ions
@hydrogens : all hydrogen atoms (those for which the first non-number in the name is a H)
@nonhydrogens : all non hydrogen atoms (those for which the first non-number in the name is not a H)
````

!!! caution ""

    Be careful since these choices are based on common names used in PDB files. Always check if
    the selected atoms are correct.

When you would like atoms from a specific residue you will use one of the atom selection commands from the above table with the following form:

````
@"definition"-chain_residuenum
@"definition"-chainresiduenum
@"definition"-residuenum
````

So for example

````
@psi-1 will select the atoms defining the psi torsion of residue 1
@psi-C1  or @psi-C_1 will define the same torsion for residue 1 of chain C.
@psi-3_1 will define the same torsion for residue 1 of chain 3.
````

Using the underscore to separate chain and residue is available as of PLUMED 2.5 and allows you to select chains
with a numeric id.

Breaking the table above down we see that for protein residues, the following groups are available:

````
# quadruplets for dihedral angles
@phi-#
@psi-#
@omega-#
@chi1-#
@chi2-#
@chi3-#
@chi4-#
@chi5-#

# all sidechain atoms (excluding glycine, including all hydrogens)
@sidechain-#
# all backbone atoms (including hydrogens)
@back-#
````

while for DNA or RNA residues, the following groups are available:

````
# quadruplets for backbone dihedral angles
@alpha-#
@beta-#
@gamma-#
@delta-#
@epsilon-#
@zeta-#

# quadruplets for sugar dihedral angles
@v0-#
@v1-#
@v2-#
@v3-#
@v4-#

# quadruplet corresponding to the chi torsional angle
@chi-#

# backbone, sugar, and base heavy atoms
@back-#
@sugar-#
@base-#

# ordered triplets of atoms on the 6-membered ring of nucleobases
# namely:
#  C2/C4/C6 for pyrimidines
#  C2/C6/C4 for purines
@lcs-#
````

Notice that `zeta` and `epsilon` groups should not be used on 3' end residue and `alpha` and `beta`
should not be used on 5' end residue.

Furthermore it is also possible to pick single atoms using the syntax
`atom-chain_residuenum`, `@atom-chainresiduenum` or `@atom-residuenum`.
As of PLUMED 2.5, this also works when the residue is not a protein/rna/dna residue.
For instance, `@OW-100` will select oxygen of water molecule with residue number 100.

Finally, notice that some shortcuts are available even when you not using the MOLINFO command.

!!! caution ""

    If a residue-chain is repeated twice in the reference pdb only the first entry will be selected.

!!! caution ""

    At the moment the HA1 atoms in a GLY residues are treated as if they are the CB atoms. This may or
    may not be true - GLY is problematic for secondary structure residues as it is achiral.

!!! caution ""

    If you use WHOLEMOLECULES RESIDUES=1-10 for a 18 amino acid protein
    ( 18 amino acids + 2 terminal groups = 20 residues ) the code will fail as it will not be able to
    interpret terminal residue 1.

## Advanced atom selection with mdtraj or MDAnalysis

Since PLUMED 2.6 it is possible to use the expressive selection syntax of [mdtraj](http://mdtraj.org/latest/atom_selection.html) and/or [MDAnalysis](https://www.mdanalysis.org/docs/documentation_pages/selections.html):

```plumed
#SETTINGS MOLFILE=regtest/basic/rt32/helix.pdb
MOLINFO STRUCTURE=regtest/basic/rt32/helix.pdb PYTHON_BIN=python
g1: GROUP ATOMS=@mda:backbone
g2: GROUP ATOMS={@mda:{resnum 1 or resid 3:5}}
g3: GROUP ATOMS={@mda:{resid 3:5} @mda:{resnum 1}}
g4: GROUP ATOMS={@mdt:{protein and (backbone or resname ALA)}}
g5: GROUP ATOMS={@mdt:{mass 5.5 to 20}} # masses guessed by mdtraj based on atom type!
g6: GROUP ATOMS={@mda:{resid 3:5} 1-10 @mda:{resnum 1}}
```

Here `@mda:` indicates that the `MDAnalysis` language is used, whereas `@mdt:` indicates that the `mdtraj` language is used. Notice that these languages typically select atoms in order. If you want to specify a different order, you can chain definitions as in `g3` above (compare with `g2`). Selections can be also chained with standard PLUMED selections (see `g6`).

The double braces are required due to the way PLUMED parses atom lists. In particular:

- The outer braces are needed to show PLUMED where the `ATOMS=...` option ends.
- The inner braces are needed to show PLUMED where each selector ends.

MDAnalysis also supports geometric selectors based on atomic coordinates. These selectors **are static** and return lists computed using the coordinates stored in the `MOLINFO` pdb file.

In order to use this syntax you should check the following points at runtime:

1. `plumed --no-mpi config has subprocess` prints `subprocess on` (should be ok on most UNIX systems).
2. You have a python interpreter with mdtraj and/or MDAnalysis installed. You can check using:
   - `python -c "import mdtraj"`
   - `python -c "import MDAnalysis"`

   In order to install these packages refer to their documentation. Pip or conda install should be ok, provided you make sure the correct python interpreter is in the execution PATH at runtime. Notice that you will only need the package(s) related to the syntax that you want to use.
3. In case you installed these modules on a python with a different name (e.g. `python3.6`), the correct check is:
   - `python3.6 -c "import mdtraj"`
   - `python3.6 -c "import MDAnalysis"`

   If this is the case, you should set the environment variable `export PYTHON_BIN=python3.6` or `export PLUMED_PYTHON_BIN=python3.6` (higher priority). Alternatively, directly provide the interpreter in the PLUMED input file
   using `MOLINFO PYTHON_BIN=python3.6` (even higher priority).
4. The PDB file that you provide to `MOLINFO` should have consecutive atom numbers starting from 1. This is currently enforced since reading atom numbers out of order (as PLUMED does) is not supported by other packages.

## Advanced atom selection with VMD (experimental)

Similarly to the `@mda:` and `@mdt:` selectors above, you can use the two following selectors in order to
access to [VMD](https://www.ks.uiuc.edu/Research/vmd/) syntax for atoms selection:
- `@vmdexec:`: This selector launches an instance of VMD, so `vmd` executable should be in your execution path.
  Might be very slow or even crash your simulation. Notice that even if `vmd` executable is used,
  the implementation is still python based and so a working python interpreter should be provided.
- `@vmd:`: This selector tries to import the `vmd` python module. Notice that the best way to obtain this module
  is not within the standard VMD installer but rather by installing the python
  module that can be found at [this link](http://github.com/Eigenstate/vmd-python).
  The module is also available on [conda](https://anaconda.org/conda-forge/vmd-python).
  You should make sure the module is available in the python interpreter used by MOLINFO
  (check using the command `python -c "import vmd"`).

These two selectors are experimental and might be removed at some point.

*/
//+ENDPLUMEDOC


/*
This action is defined in core/ as it is used by other actions.
Anyway, it is registered here, so that excluding this module from
compilation will exclude it from plumed.
*/

typedef PLMD::GenericMolInfo MolInfo;

PLUMED_REGISTER_ACTION(MolInfo,"MOLINFO")

}
}
