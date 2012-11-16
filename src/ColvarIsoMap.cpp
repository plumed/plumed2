/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#include <cmath>

#include "ColvarPathMSDBase.h"

#include "Colvar.h"
#include "ActionRegister.h"

#include <string>
#include <cstring>
#include <cassert>
#include <iostream>
#include "PDB.h"
#include "RMSD.h"
#include "Atoms.h"
#include "Tools.h"

using namespace std;

namespace PLMD{

//+PLUMEDOC COLVAR ISOMAP
/*
The implementation of this collective variable is based on the \ref PATHMSD.

PUT HERE DOCUMENTATION SPECIFIC FOR THIS VARIABLE

*/
//+ENDPLUMEDOC
   
class ColvarIsoMap : public ColvarPathMSDBase {
public:
  ColvarIsoMap(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(ColvarIsoMap,"ISOMAP")

void ColvarIsoMap::registerKeywords(Keywords& keys){
  ColvarPathMSDBase::registerKeywords(keys);
  keys.add("compulsory","PROPERTY","the property to be used in the indexing: this goes in the REMARK field of the reference");
}

ColvarIsoMap::ColvarIsoMap(const ActionOptions&ao):
Action(ao),
ColvarPathMSDBase(ao)
{
  // this is the only additional keyword needed 
  parseVector("PROPERTY",labels);
  checkRead();
  if(labels.size()==0){
	char buf[500];
        sprintf(buf,"Need to specify PROPERTY with this action\n");
        plumed_merror(buf);
        exit(0);
  }else{
      for(unsigned i=0;i<labels.size();i++){
	log<<" found custom propety to be found in the REMARK line: "<<labels[i].c_str()<<"\n";
        addComponentWithDerivatives(labels[i].c_str()); componentIsNotPeriodic(labels[i].c_str());
      }
      // add distance anyhow
      addComponentWithDerivatives("zzz"); componentIsNotPeriodic("zzz");
      //reparse the REMARK field and pick the index 
      for(unsigned i=0;i<pdbv.size();i++){
      	     vector<std::string> myv(pdbv[i].getRemark());	
              // now look for X=1.34555 Y=5.6677
              vector<double> labelvals; 
              for(unsigned j=0;j<labels.size();j++){
      	      double val;
                     if(Tools::parse(myv,labels[j],val)){labelvals.push_back(val);}
                     else{
      		   char buf[500];
      		   sprintf(buf,"PROPERTY LABEL \" %s \" NOT FOUND IN REMARK FOR FRAME %u \n",labels[j].c_str(),i);
      		   plumed_merror(buf);  
                     };
              }
              indexvec.push_back(labelvals);
      }
  }
  requestAtoms(pdbv[0].getAtomNumbers());  
 
};

}


