/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2020 The plumed team
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
#include "PathMSDBase.h"
#include "core/PlumedMain.h"

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR PROPERTYMAP
/*
Calculate generic property maps.

This Colvar calculates the property maps according to the work of Spiwok \cite Spiwok:2011ce.


Basically it calculates
\f{eqnarray*}{
X=\frac{\sum_i X_i*\exp(-\lambda D_i(x))}{\sum_i  \exp(-\lambda D_i(x))} \\
Y=\frac{\sum_i Y_i*\exp(-\lambda D_i(x))}{\sum_i  \exp(-\lambda D_i(x))} \\
\cdots\\
zzz=-\frac{1}{\lambda}\log(\sum_i  \exp(-\lambda D_i(x)))
\f}

where the parameters \f$X_i\f$  and  \f$Y_i\f$ are provided in the input pdb (allv.pdb in this case) and
 \f$D_i(x)\f$  is the mean squared displacement after optimal alignment calculated on the pdb frames you input (see Kearsley).


When running with periodic boundary conditions, the atoms should be
in the proper periodic image. This is done automatically since PLUMED 2.5,
by considering the ordered list of atoms and rebuilding molecules using a procedure
that is equivalent to that done in \ref WHOLEMOLECULES . Notice that
rebuilding is local to this action. This is different from \ref WHOLEMOLECULES
which actually modifies the coordinates stored in PLUMED.

In case you want to recover the old behavior you should use the NOPBC flag.
In that case you need to take care that atoms are in the correct
periodic image.

\par Examples

\plumedfile
p3: PROPERTYMAP REFERENCE=allv.pdb PROPERTY=X,Y LAMBDA=69087 NEIGH_SIZE=8 NEIGH_STRIDE=4
PRINT ARG=p3.X,p3.Y,p3.zzz STRIDE=1 FILE=colvar FMT=%8.4f
\endplumedfile

note that NEIGH_STRIDE=4 NEIGH_SIZE=8 control the neighbor list parameter (optional but
recommended for performance) and states that the neighbor list will be calculated every 4
steps and consider only the closest 8 member to the actual md snapshots.

In this case the input line instructs plumed to look for two properties X and Y with attached values in the REMARK
line of the reference pdb (Note: No spaces from X and = and 1 !!!!).
e.g.

\auxfile{allv.pdb}
REMARK X=1 Y=2
ATOM      1  CL  ALA     1      -3.171   0.295   2.045  1.00  1.00
ATOM      5  CLP ALA     1      -1.819  -0.143   1.679  1.00  1.00
END
REMARK X=2 Y=3
ATOM      1  CL  ALA     1      -3.175   0.365   2.024  1.00  1.00
ATOM      5  CLP ALA     1      -1.814  -0.106   1.685  1.00  1.00
END
\endauxfile

\note
The implementation of this collective variable and of \ref PATHMSD
is shared, as well as most input options.

*/
//+ENDPLUMEDOC

class PropertyMap : public PathMSDBase {
public:
  explicit PropertyMap(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(PropertyMap,"PROPERTYMAP")

void PropertyMap::registerKeywords(Keywords& keys) {
  PathMSDBase::registerKeywords(keys);
  keys.add("compulsory","PROPERTY","the property to be used in the indexing: this goes in the REMARK field of the reference");
  ActionWithValue::useCustomisableComponents(keys);
  keys.addOutputComponent("zzz","default","the minimum distance from the reference points");
}

PropertyMap::PropertyMap(const ActionOptions&ao):
  Action(ao),
  PathMSDBase(ao)
{
  // this is the only additional keyword needed
  parseVector("PROPERTY",labels);
  checkRead();
  log<<"  Bibliography "
     <<plumed.cite("Spiwok V, Kralova B  J. Chem. Phys. 135,  224504 (2011)")
     <<"\n";
  if(labels.size()==0) {
    char buf[500];
    sprintf(buf,"Need to specify PROPERTY with this action\n");
    plumed_merror(buf);
  } else {
    for(unsigned i=0; i<labels.size(); i++) {
      log<<" found custom propety to be found in the REMARK line: "<<labels[i].c_str()<<"\n";
      addComponentWithDerivatives(labels[i]); componentIsNotPeriodic(labels[i]);
    }
    // add distance anyhow
    addComponentWithDerivatives("zzz"); componentIsNotPeriodic("zzz");
    //reparse the REMARK field and pick the index
    for(unsigned i=0; i<pdbv.size(); i++) {
      // now look for X=1.34555 Y=5.6677
      std::vector<double> labelvals;
      for(unsigned j=0; j<labels.size(); j++) {
        double val;
        if( pdbv[i].getArgumentValue(labels[j],val) ) {labelvals.push_back(val);}
        else {
          char buf[500];
          sprintf(buf,"PROPERTY LABEL \" %s \" NOT FOUND IN REMARK FOR FRAME %u \n",labels[j].c_str(),i);
          plumed_merror(buf);
        };
      }
      indexvec.push_back(labelvals);
    }
  }
  requestAtoms(pdbv[0].getAtomNumbers());

}

}
}


