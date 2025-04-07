/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
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
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/Group.h"
#include "AdjacencyMatrixBase.h"

//+PLUMEDOC MATRIX CONTACT_MATRIX
/*
Adjacency matrix in which two atoms are adjacent if they are within a certain cutoff.

A useful tool for developing complex collective variables is the notion of the
so called adjacency matrix. An adjacency matrix can be an $N \times N$ matrix in which the $i$th, $j$th element tells you whether
or not the $i$th and $j$th atoms/molecules from a set of $N$ atoms/molecules are adjacent or not. Alternatively, you can calculate
an adjacency matrix between a set of $N$ atoms and a second set of $M$ atoms.  For this type of matrix the $i$th, $j$th element tells you
whether the whether the $i$th atom in the first group and the $j$th atom in the second group are adjacent or not.  The adjacency matrix in this
case is thus $N \times M$.

The simplest type of adjacency matrix is a contact matrix.  The elements of a contact matrix are calculated using:

$$
a_{ij} = \sigma( r_{ij} )
$$

where $r_{ij}$ is the magnitude of the vector connecting atoms $i$ and $j$ and where $\sigma$ is a switching function. If you want to calculate a
contact matrix for one group of atoms the input would look something like this:

```plumed
c1: CONTACT_MATRIX GROUP=1-7 SWITCH={RATIONAL R_0=2.6 NN=6 MM=12}
```

Alternatively, if you want to calculate the contact matrix between two groups of atoms you would use an input like following:

```plumed
c2: CONTACT_MATRIX GROUPA=1-7 GROUPB=8-20 SWITCH={RATIONAL R_0=2.6 NN=6 MM=12}
```

Once you have calculated a contact matrix you can perform various matrix operations by using the tools in the matrixtools or clusters modules.

## Coordination numbers

If a contact matrix, $\mathbf{C}$, is multiplied from the back by a vector of ones the $i$th element of the resulting matrix tells you the number of atoms that are
within $r_c$ of atom $i$.  In other words, the coordination numbers of the atoms can be calculated from the contact matrix by doing matrix vector multiplication.

This realisation about the relationship between the contact map and the coordination number is heavily used in PLUMED.  For example, to calculate
and print the coordination numbers of the first 7 atoms in the system with themselves you would use an input something like this:

```plumed
c1: CONTACT_MATRIX GROUP=1-7 SWITCH={RATIONAL R_0=2.6 NN=6 MM=12}
ones: ONES SIZE=7
cc: MATRIX_VECTOR_PRODUCT ARG=c1,ones
PRINT ARG=cc FILE=colvar
```

Implmenting the coordination number this way is useful as there are many different ways to define whether two atoms/molecules and to construct a "contact" matrix based on
the result.  For example:

* You could say that two molecules are connected if they are within a certain distance of each other and if they have the same orientation (see [TORSIONS_MATRIX](TORSIONS_MATRIX.md)).
* You could say that two water molecules are connected if they are hydrogen bonded to each other (see [HBOND_MATRIX](HBOND_MATRIX.md)).
* You could say that two atoms are connected if they are within a certain distance of each other and if they have similar values for a CV (see [OUTER_PRODUCT](OUTER_PRODUCT.md)).

When the coordination numbers is implemented in the way described above (by doing the matrix-vector multiplication) you have the flexibility to define the contact matrix that
is used in the multiplication in whatever way you choose.  In other words, this implementation of the coordination number is much more flexible. For example, suppose you want
to calculate the number of atoms that have a coordination is greater than 3.0.  You can do this with PLUMED using the following input:

```plumed
# Calculate the contact matrix for the first seven atoms in the system
c1: CONTACT_MATRIX GROUP=1-7 SWITCH={RATIONAL R_0=2.6 NN=6 MM=12}
# Calculate the coordination numbers for the first seven atoms in the system
ones: ONES SIZE=7
cc: MATRIX_VECTOR_PRODUCT ARG=c1,ones
# Set the ith element of the vector mtc equal to one if the coordination number of atom i is greater than 3.
mtc: MORE_THAN ARG=cc SWITCH={RATIONAL D_0=3 R_0=1}
# Calculate the number of atoms with a coordination number greater than 3.
s: SUM ARG=mtc PERIODIC=NO
PRINT ARG=s FILE=colvar
```

Alternatively, consider the CV that was used to study perovskite nucleation in [this paper](https://pubs.acs.org/doi/abs/10.1021/acs.chemmater.9b04259).  The CV here measures the number of
methylamonium molecules that are attached to at least 5 other methylamoniusm molecules, at least 7 lead atoms and at least 11 ionide ions.  We can calculate something akin to this
CV and apply a restraint on the resulting quantity by using the following input file:

```plumed
# Lead ions
Pb: GROUP ATOMS=1-64
# Iodide atoms
I: GROUP ATOMS=65-256
# Methylamonium "atoms" -- in the real CV these are centers of mass rather than single atoms
cn: GROUP ATOMS=257-320

ones64: ONES SIZE=64
# Contact matrix that determines if methylamonium molecules are within 8 A of each other
cm_cncn: CONTACT_MATRIX GROUP=cn SWITCH={RATIONAL R_0=0.8}
# Coordination number of methylamounium with methylamonium
cc_cncn: MATRIX_VECTOR_PRODUCT ARG=cm_cncn,ones64
# Vector with elements that are one if coordiantion of methylamonium with methylamonium >5
mt_cncn: MORE_THAN ARG=cc_cncn SWITCH={RATIONAL R_0=5 NN=12 MM=24}

# Contact matrix that determines if methylamoinium moleulcule and Pb atom are within 7.5 A of each other
cm_cnpb: CONTACT_MATRIX GROUPA=cn GROUPB=Pb SWITCH={RATIONAL R_0=0.75}
# Coordination number of methylamonium with Pb
cc_cnpb: MATRIX_VECTOR_PRODUCT ARG=cm_cnpb,ones64
# Vector with elements that are one if coordination of methylamounium with lead is >7
mt_cnpb: MORE_THAN ARG=cc_cnpb SWITCH={RATIONAL R_0=7 NN=12 MM=24}

ones192: ONES SIZE=192
# Contact matrix that determines if methylamoinium moleulcule and I atom are within 6.5 A of each other
cm_cnI: CONTACT_MATRIX GROUPA=cn GROUPB=I SWITCH={RATIONAL R_0=0.65}
# Coordination number of methylamonium with I
cc_cnI: MATRIX_VECTOR_PRODUCT ARG=cm_cnI,ones192
# Vector with elements that are one if coordination of methylamounium with lead is >11
mt_cnI: MORE_THAN ARG=cc_cnI SWITCH={RATIONAL R_0=11 NN=12 MM=24}

# Element wise product of these three input vectors.
# mm[i]==1 if coordination number of corrsponding methylamounium with methylamonium is >5
# and if coordination of methylamounium with Pb is >7 and if coordination of methylamounium with I > 11
mm: CUSTOM ARG=mt_cncn,mt_cnpb,mt_cnI FUNC=x*y*z PERIODIC=NO

# Sum of coordination numbers and thus equal to number of methylamoniums with desired coordination numbers
ff: SUM ARG=mm PERIODIC=NO

rr: RESTRAINT ARG=ff AT=62 KAPPA=10
```

## COMPONENTS flag

If you add the flag COMPONENTS to the input as shown below:

```plumed
c4: CONTACT_MATRIX GROUP=1-7 COMPONENTS SWITCH={RATIONAL R_0=2.6 NN=6 MM=12}
```

then four matrices with the labels `c4.w`, `c4.x`, `c4.y` and `c4.z` are output by the action. The matrix with the label `c4.w` is the adjacency matrix
that would be output if you had not added the COMPONENTS flag. The $i,j$ component of the matrices `c4.x`, `c4.y` and `c4.z` contain the $x$, $y$ and $z$
components of the vector connecting atoms $i$ and $j$. Importantly, however, the components of these vectors are only stored in `c4.x`, `c4.y` and `c4.z`
if the elements of `c4.w` are non-zero.

## Optimisation details

Adjacency matrices are sparse.  Each atom is only be connected to a small number of neighbours and the vast majority of the elements of the contact matrix are thus zero.  To reduce
the amount of memory that PLUMED requires PLUMED uses sparse matrix storage.  Consequently, whenever you calculate and store a contact matrix only the elements of the matrix that are
non-zero are stored.  The same thing holds for the additional values that are created when you use the COMPONENTS flag. The components of the vectors connecting atoms are only stored
when the elements of `c4.w` are non-zero.

We can also use the sparsity of the adjacency matrix to make the time required to compute a contact matrix scale linearly rather than quadratically with the number of atoms. Element
$i,j$ of the contact matrix is only non-zero if two atoms are within a cutoff, $r_c$. We can determine that many pairs of atoms are further appart than $r_c$ without computing the
distance between these atoms by using divide and conquer strategies such as linked lists and neighbour lists.  __To turn on these features you need to set the `D_MAX` parameter in the
switching functions.__  The value you pass to the `D_MAX` keyword is used as the cutoff in the link cell algorithm.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace adjmat {

class ContactMatrixShortcut : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit ContactMatrixShortcut(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(ContactMatrixShortcut,"CONTACT_MATRIX")

void ContactMatrixShortcut::registerKeywords(Keywords& keys) {
  AdjacencyMatrixBase::registerKeywords( keys );
  keys.remove("GROUP");
  keys.add("numbered","GROUP","specifies the list of atoms that should be assumed indistinguishable");
  keys.add("compulsory","NN","6","The n parameter of the switching function ");
  keys.add("compulsory","MM","0","The m parameter of the switching function; 0 implies 2*NN");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
  keys.add("numbered","SWITCH","the input for the switching function that acts upon the distance between each pair of atoms");
  keys.linkActionInDocs("SWITCH","LESS_THAN");
  keys.addActionNameSuffix("_PROPER");
  keys.needsAction("TRANSPOSE");
  keys.needsAction("CONCATENATE");
}

ContactMatrixShortcut::ContactMatrixShortcut(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  std::vector<std::string> grp_str;
  std::string atomsstr="";
  std::vector<std::string> atomsvec;
  parseVector("ATOMS",atomsvec);
  if( atomsvec.size()>0 )  {
    for(unsigned i=0; i<atomsvec.size(); ++i) {
      Group* gg = plumed.getActionSet().selectWithLabel<Group*>( atomsvec[i] );
      if( gg ) {
        grp_str.push_back( atomsvec[i] );
      }
    }
    if( grp_str.size()!=atomsvec.size() ) {
      grp_str.resize(0);
      atomsstr = " ATOMS=" + atomsvec[0];
      for(unsigned i=1; i<atomsvec.size(); ++i) {
        atomsstr += "," + atomsvec[i];
      }
    }
  } else {
    std::string grp_inpt;
    for(unsigned i=1;; ++i) {
      if( !parseNumbered("GROUP",i,grp_inpt) ) {
        break;
      }
      grp_str.push_back( grp_inpt );
    }
  }
  if( grp_str.size()>9 ) {
    error("cannot handle more than 9 groups");
  }
  if( grp_str.size()==0 )  {
    readInputLine( getShortcutLabel() + ": CONTACT_MATRIX_PROPER " + atomsstr + " " + convertInputLineToString() );
    return;
  }

  for(unsigned i=0; i<grp_str.size(); ++i) {
    std::string sw_str, num;
    Tools::convert( i+1, num );
    parseNumbered("SWITCH", (i+1)*10 + 1 + i,  sw_str );
    if( sw_str.length()==0 ) {
      error("missing SWITCH" + num + num + " keyword");
    }
    readInputLine( getShortcutLabel() + num +  num + ": CONTACT_MATRIX_PROPER GROUP=" + grp_str[i] + " SWITCH={" + sw_str + "}" );
    for(unsigned j=0; j<i; ++j) {
      std::string sw_str2, jnum;
      Tools::convert( j+1, jnum );
      parseNumbered("SWITCH", (j+1)*10 + 1 + i, sw_str2);
      if( sw_str2.length()==0 ) {
        error("missing SWITCH" + jnum + num + " keyword");
      }
      readInputLine( getShortcutLabel() + jnum + num + ": CONTACT_MATRIX_PROPER GROUPA=" + grp_str[j] + " GROUPB=" + grp_str[i] + " SWITCH={" + sw_str2 +"}");
      readInputLine( getShortcutLabel() + num +  jnum + ": TRANSPOSE ARG=" + getShortcutLabel() + jnum + num );
    }
  }
  std::string join_matrices = getShortcutLabel() + ": CONCATENATE";
  for(unsigned i=0; i<grp_str.size(); ++i) {
    std::string inum;
    Tools::convert(i+1,inum);
    for(unsigned j=0; j<grp_str.size(); ++j) {
      std::string jnum;
      Tools::convert(j+1,jnum);
      if( i>j ) {
        join_matrices += " MATRIX" + inum + jnum + "=" + getShortcutLabel() + inum +  jnum;
      } else {
        join_matrices += " MATRIX" + inum + jnum + "=" + getShortcutLabel() + inum +  jnum;
      }
    }
  }
  readInputLine( join_matrices );
}

}
}
