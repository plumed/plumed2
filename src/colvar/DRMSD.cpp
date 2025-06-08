/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2024 The plumed team
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
#include "tools/PDB.h"

namespace PLMD {
namespace colvar {

//+PLUMEDOC DCOLVAR DRMSD
/*
Calculate the distance RMSD with respect to a reference structure.

To calculate the root-mean-square deviation between the atoms in two configurations
you must first superimpose the two structures in some ways.  Obviously, it is the internal vibrational
motions of the structure - i.e. not the translations and rotations - that are interesting. However,
aligning two structures by removing the translational and rotational motions is not easy.  Furthermore,
in some cases there can be alignment issues caused by so-called frame-fitting problems. It is thus
often easier to calculate the distances between all the pairs of atoms.  The distance
between the two structures, ${\bf X}^a$ and ${\bf X}^b$ can then be measured as:

$$
d({\bf X}^A, {\bf X}^B) = \sqrt{\frac{1}{N(N-1)} \sum_{i \ne j} [ d({\bf x}_i^a,{\bf x}_j^a) - d({\bf x}_i^b,{\bf x}_j^b) ]^2}
$$

where $N$ is the number of atoms and $d({\bf x}_i,{\bf x}_j)$ represents the distance between
atoms $i$ and $j$.  Clearly, this representation of the configuration is invariant to translation and rotation.
However, it can become expensive to calculate when the number of atoms is large.  This can be resolved
within the DRMSD colvar by setting `LOWER_CUTOFF` and `UPPER_CUTOFF`.  These keywords ensure that only
pairs of atoms that are within a certain range are incorporated into the above sum.

In PDB files the atomic coordinates and box lengths should be in Angstroms unless
you are working with natural units.  If you are working with natural units then the coordinates
should be in your natural length unit.  [Click here](http://www.wwpdb.org/docs.html) for more details on the PDB file format.

## Examples

The following tells plumed to calculate the distance RMSD between
the positions of the atoms in the reference file and their instantaneous
position. Only pairs of atoms whose distance in the reference structure is within
0.1 and 0.8 nm are considered.

```plumed
#SETTINGS INPUTFILES=regtest/basic/rt-drmsd/test1.pdb
d: DRMSD ...
  REFERENCE=regtest/basic/rt-drmsd/test1.pdb
  LOWER_CUTOFF=0.1 UPPER_CUTOFF=0.8
...
```

The following tells plumed to calculate the square of the DRMSD value for a pair of molecules.

```plumed
#SETTINGS INPUTFILES=regtest/basic/rt-drmsd/test0.pdb
d: DRMSD ...
  REFERENCE=regtest/basic/rt-drmsd/test0.pdb
  LOWER_CUTOFF=0.1 UPPER_CUTOFF=0.8 TYPE=INTER-DRMSD
  SQUARED
...
```

Notice that in the input reference file (which you can see by clicking on regtest/basic/rt-drmsd/test0.pdb )
the atoms in each of the two molecules are separated by a TER
command as shown below.

In this example the INTER-DRMSD type ensures that the set of distances from which the final
quantity is computed involve one atom from each of the two molecules.  If this is replaced
by INTRA-DRMSD then only those distances involving pairs of atoms that are both in the same
molecule are computed.

## Periodic boundary conditions

Notice that PLUMED does not use periodic boundary conditions when computing the reference distances.  However, the instantaneous
values of the distances are computed using periodic boundary conditions.  If you would like not to use periodic boundary conditions
when computing these instaneous distances you use the NOPBC flag as illustrated below:

```plumed
#SETTINGS INPUTFILES=regtest/basic/rt-drmsd/test1.pdb
d: DRMSD ...
  REFERENCE=regtest/basic/rt-drmsd/test1.pdb
  LOWER_CUTOFF=0.1 UPPER_CUTOFF=0.8 NOPBC
...
```

*/
//+ENDPLUMEDOC

class DRMSD : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit DRMSD(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(DRMSD,"DRMSD")

void DRMSD::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","REFERENCE","a file in pdb format containing the reference structure and the atoms involved in the CV.");
  keys.add("optional","LOWER_CUTOFF","only pairs of atoms further than LOWER_CUTOFF are considered in the calculation.");
  keys.add("optional","UPPER_CUTOFF","only pairs of atoms closer than UPPER_CUTOFF are considered in the calculation.");
  keys.add("compulsory","TYPE","DRMSD","what kind of DRMSD would you like to calculate.  You can use either the normal DRMSD involving all the distances between "
           "the atoms in your molecule.  Alternatively, if you have multiple molecules you can use the type INTER-DRMSD "
           "to compute DRMSD values involving only those distances between the atoms at least two molecules or the type INTRA-DRMSD "
           "to compute DRMSD values involving only those distances between atoms in the same molecule");
  keys.addFlag("SQUARED",false,"This should be setted if you want MSD instead of RMSD ");
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  // This is just ignored in reality which is probably bad
  keys.addDeprecatedFlag("NUMERICAL_DERIVATIVES","");
  keys.setValueDescription("scalar/vector","the DRMSD distance between the instantaneous structure and the reference structure");
  keys.needsAction("SUM");
  keys.needsAction("DISTANCE");
  keys.needsAction("CONSTANT");
  keys.needsAction("EUCLIDEAN_DISTANCE");
  keys.needsAction("CUSTOM");
}

DRMSD::DRMSD( const ActionOptions& ao ):
  Action(ao),
  ActionShortcut(ao) {
  // Read in the reference configuration
  std::string reference;
  parse("REFERENCE",reference);
  // First bit of input for the instantaneous distances
  bool numder;
  parseFlag("NUMERICAL_DERIVATIVES",numder);
  double fake_unit=0.1;
  FILE* fp2=fopen(reference.c_str(),"r");
  bool do_read=true;
  unsigned nframes=0;
  while( do_read ) {
    PDB mypdb;
    do_read=mypdb.readFromFilepointer(fp2,false,fake_unit);
    if( !do_read && nframes>0 ) {
      break ;
    }
    nframes++;
  }
  fclose(fp2);

  // Get cutoff information
  double lcut=0;
  parse("LOWER_CUTOFF",lcut);
  std::string drmsd_type;
  parse("TYPE",drmsd_type);
  double ucut=std::numeric_limits<double>::max();
  parse("UPPER_CUTOFF",ucut);
  bool nopbc;
  parseFlag("NOPBC",nopbc);
  std::string pbc_str;
  if(nopbc) {
    pbc_str="NOPBC";
  }
  // Open the pdb file
  FILE* fp=fopen(reference.c_str(),"r");
  do_read=true;
  if(!fp) {
    error("could not open reference file " + reference );
  }
  unsigned n=0;
  std::string allpairs="";
  std::vector<std::pair<unsigned,unsigned> > upairs;
  std::vector<std::string> refvals;
  while ( do_read ) {
    PDB mypdb;
    do_read=mypdb.readFromFilepointer(fp,false,fake_unit);
    if( !do_read && n>0 ) {
      break ;
    }
    std::vector<Vector> pos( mypdb.getPositions() );
    unsigned nn=1;
    if( pos.size()==0 ) {
      error("read no atoms from file named " + reference );
    }
    // This is what we do for the first frame
    if( n==0 ) {
      std::vector<AtomNumber> atoms( mypdb.getAtomNumbers() );
      if( drmsd_type=="DRMSD" ) {
        for(unsigned i=0; i<atoms.size()-1; ++i) {
          std::string istr;
          Tools::convert( atoms[i].serial(), istr );
          for(unsigned j=i+1; j<atoms.size(); ++j) {
            std::string jstr;
            Tools::convert( atoms[j].serial(), jstr );
            double distance = delta( pos[i], pos[j] ).modulo();
            if( distance < ucut && distance > lcut ) {
              std::string num;
              Tools::convert( nn, num );
              nn++;
              // Add this pair to list of pairs
              upairs.push_back( std::pair<unsigned,unsigned>(i,j) );
              // Add this distance to list of reference values
              std::string dstr;
              Tools::convert( distance, dstr );
              refvals.push_back( dstr );
              // Calculate this distance
              if( nframes==1 ) {
                allpairs += " ATOMS" + num + "=" + istr + "," + jstr;
              } else {
                readInputLine( getShortcutLabel() + "_d" + num + ": DISTANCE ATOMS=" + istr + "," + jstr + " " + pbc_str );
              }
            }
          }
        }
      } else {
        unsigned nblocks = mypdb.getNumberOfAtomBlocks();
        std::vector<unsigned> blocks( 1 + nblocks );
        if( nblocks==1 ) {
          blocks[0]=0;
          blocks[1]=atoms.size();
        } else {
          blocks[0]=0;
          for(unsigned i=0; i<nblocks; ++i) {
            blocks[i+1]=mypdb.getAtomBlockEnds()[i];
          }
        }
        if( drmsd_type=="INTRA-DRMSD" ) {
          for(unsigned i=0; i<nblocks; ++i) {
            for(unsigned iatom=blocks[i]+1; iatom<blocks[i+1]; ++iatom) {
              std::string istr;
              Tools::convert( atoms[iatom].serial(), istr );
              for(unsigned jatom=blocks[i]; jatom<iatom; ++jatom) {
                std::string jstr;
                Tools::convert( atoms[jatom].serial(), jstr );
                double distance = delta( pos[iatom], pos[jatom] ).modulo();
                if(distance < ucut && distance > lcut ) {
                  std::string num;
                  Tools::convert( nn, num );
                  nn++;
                  // Add this pair to list of pairs
                  upairs.push_back( std::pair<unsigned,unsigned>(iatom,jatom) );
                  // Add this distance to list of reference values
                  std::string dstr;
                  Tools::convert( distance, dstr );
                  refvals.push_back( dstr );
                  // Calculate this distance
                  if( nframes==1 ) {
                    allpairs += " ATOMS" + num + "=" + istr + "," + jstr;
                  } else {
                    readInputLine( getShortcutLabel() + "_d" + num + ": DISTANCE ATOMS=" + istr + "," + jstr + " " + pbc_str );
                  }
                }
              }
            }
          }
        } else if( drmsd_type=="INTER-DRMSD" ) {
          for(unsigned i=1; i<nblocks; ++i) {
            for(unsigned j=0; j<i; ++j) {
              for(unsigned iatom=blocks[i]; iatom<blocks[i+1]; ++iatom) {
                std::string istr;
                Tools::convert( atoms[iatom].serial(), istr );
                for(unsigned jatom=blocks[j]; jatom<blocks[j+1]; ++jatom) {
                  std::string jstr;
                  Tools::convert( atoms[jatom].serial(), jstr );
                  double distance = delta( pos[iatom], pos[jatom] ).modulo();
                  if(distance < ucut && distance > lcut ) {
                    std::string num;
                    Tools::convert( nn, num );
                    nn++;
                    // Add this pair to list of pairs
                    upairs.push_back( std::pair<unsigned,unsigned>(iatom,jatom) );
                    // Add this distance to list of reference values
                    std::string dstr;
                    Tools::convert( distance, dstr );
                    refvals.push_back( dstr );
                    // Calculate this distance
                    if( nframes==1 ) {
                      allpairs += " ATOMS" + num + "=" + istr + "," + jstr;
                    } else {
                      readInputLine( getShortcutLabel() + "_d" + num + ": DISTANCE ATOMS=" + istr + "," + jstr + " " + pbc_str );
                    }
                  }
                }
              }
            }
          }
        } else {
          plumed_merror( drmsd_type + " is not valid input to TYPE keyword");
        }
      }
      // This is for every subsequent frame
    } else {
      for(unsigned i=0; i<refvals.size(); ++i) {
        std::string dstr;
        Tools::convert( delta( pos[upairs[i].first], pos[upairs[i].second] ).modulo(), dstr );
        refvals[i] += "," + dstr;
      }
    }
    n++;
  }
  // Now create values that hold all the reference distances
  fclose(fp);

  if( nframes==1 ) {
    readInputLine( getShortcutLabel() + "_d: DISTANCE" + allpairs + " " + pbc_str );
    std::string refstr = refvals[0];
    for(unsigned i=1; i<refvals.size(); ++i) {
      refstr += "," + refvals[i];
    }
    readInputLine( getShortcutLabel() + "_ref: CONSTANT VALUES="  + refstr );
    readInputLine( getShortcutLabel() + "_diffs: CUSTOM ARG=" + getShortcutLabel() + "_d," + getShortcutLabel() + "_ref FUNC=(x-y)*(x-y) PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_u: SUM ARG=" + getShortcutLabel() + "_diffs PERIODIC=NO");
  } else {
    std::string arg_str1, arg_str2;
    for(unsigned i=0; i<refvals.size(); ++i ) {
      std::string inum;
      Tools::convert( i+1, inum );
      readInputLine( getShortcutLabel() + "_ref" + inum + ": CONSTANT VALUES=" + refvals[i] );
      if( i==0 ) {
        arg_str1 = getShortcutLabel() + "_d" + inum;
        arg_str2 = getShortcutLabel() + "_ref" + inum;
      } else {
        arg_str1 += "," + getShortcutLabel() + "_d" + inum;
        arg_str2 += "," + getShortcutLabel() + "_ref" + inum;
      }
    }
    // And calculate the euclidean distances between the true distances and the references
    readInputLine( getShortcutLabel() + "_u: EUCLIDEAN_DISTANCE SQUARED ARG2=" + arg_str1 + " ARG1=" + arg_str2 );
  }
  // And final value
  std::string nvals;
  Tools::convert( refvals.size(), nvals );
  bool squared;
  parseFlag("SQUARED",squared);
  if( squared ) {
    readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + getShortcutLabel() + "_u FUNC=x/" + nvals + " PERIODIC=NO");
  } else {
    readInputLine( getShortcutLabel() + "_2: CUSTOM ARG=" + getShortcutLabel() + "_u FUNC=(x/" + nvals + ") PERIODIC=NO");
    readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + getShortcutLabel() + "_2 FUNC=sqrt(x) PERIODIC=NO");
  }
}

}
}


