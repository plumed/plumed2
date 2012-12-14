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

#include "Colvar.h"
#include "ColvarPathMSDBase.h"
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

void ColvarPathMSDBase::registerKeywords(Keywords& keys){
  Colvar::registerKeywords(keys);
  keys.add("compulsory","LAMBDA","the lambda parameter is needed for smoothing, is in the units of plumed");
  keys.add("compulsory","REFERENCE","the pdb is needed to provide the various milestones");
  keys.add("optional","NEIGH_SIZE","size of the neighbor list");
  keys.add("optional","NEIGH_STRIDE","how often the neighbor list needs to be calculated in time units");
}

ColvarPathMSDBase::ColvarPathMSDBase(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
pbc(true),
neigh_size(-1),
neigh_stride(-1.),
nframes(0)
{
  parse("LAMBDA",lambda);
  parse("NEIGH_SIZE",neigh_size);
  parse("NEIGH_STRIDE",neigh_stride);
  parse("REFERENCE",reference);

  // open the file
  FILE* fp=fopen(reference.c_str(),"r");
  if (fp!=NULL)
  {
    log<<"Opening reference file "<<reference.c_str()<<"\n";
    bool do_read=true;
    while (do_read){
         PDB mypdb; 
         RMSD mymsd(log); 
         do_read=mypdb.readFromFilepointer(fp,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength());
         if(do_read){
            nframes++;
            log<<"Found PDB: "<<nframes<<" containing  "<<mypdb.getAtomNumbers().size()<<" atoms\n"; 
	    pdbv.push_back(mypdb); 
//            requestAtoms(mypdb.getAtomNumbers()); // is done in non base classes 
            derivs_s.resize(mypdb.getAtomNumbers().size());
            derivs_z.resize(mypdb.getAtomNumbers().size());
            mymsd.set(mypdb,"OPTIMAL");
            msdv.push_back(mymsd); // the vector that stores the frames
            //log<<mypdb; 
         }else{break ;}
    }
    fclose (fp);
    log<<"Found TOTAL "<<nframes<< " PDB in the file "<<reference.c_str()<<" \n"; 
  } 
  if(neigh_stride>0. || neigh_size>0){
           if(neigh_size>nframes){
           	log.printf(" List size required ( %d ) is too large: resizing to the maximum number of frames required: %d  \n",neigh_size,nframes);
 		neigh_size=nframes;
           }
           log.printf("  Neighbor list enabled: \n");
           log.printf("                size   :  %d elements\n",neigh_size);
           log.printf("                stride :  %f time \n",neigh_stride);
  }else{
           log.printf("  Neighbor list NOT enabled \n");
  }

}

void ColvarPathMSDBase::calculate(){

  //log.printf("NOW CALCULATE! \n");

  // clean vectors
  for(unsigned i=0;i< derivs_z.size();i++){derivs_z[i].zero();}

  // full list: resize and calculate the msd 
  if(imgVec.empty()){ // this is the signal that means: recalculate all 
      imgVec.resize(nframes);  
      derivs_v.clear();for(unsigned i=0;i<nframes;i++){derivs_v.push_back(derivs_z);}
      for(unsigned i=0;i<nframes;i++){
          imgVec[i].distder=derivs_z;
          imgVec[i].property=indexvec[i];
          imgVec[i].index=i;
          imgVec[i].distance=msdv[imgVec[i].index].calculate(getPositions(),imgVec[i].distder,true);
      }
  }else{// just recalculate rmsd for the things you have in the list and assume that property and msdv are ok 
	for(unsigned i=0;i<imgVec.size();i++){
		imgVec[i].distance=msdv[imgVec[i].index].calculate(getPositions(),imgVec[i].distder,true);
        }  
  }

  vector<Value*> val_s_path;
  if(labels.size()>0){
    for(unsigned i=0;i<labels.size();i++){ val_s_path.push_back(getPntrToComponent(labels[i].c_str()));}
  }else{
     val_s_path.push_back(getPntrToComponent("sss"));
  } 
  Value* val_z_path=getPntrToComponent("zzz");

  vector<double> s_path(val_s_path.size());for(unsigned i=0;i<s_path.size();i++)s_path[i]=0.;
  double partition=0.;
  double tmp;

  typedef  vector< class ImagePath  >::iterator imgiter;
  for(imgiter it=imgVec.begin();it!=imgVec.end();++it){ 
    (*it).similarity=exp(-lambda*((*it).distance));
    //log<<"DISTANCE "<<(*it).distance<<"\n";
    for(unsigned i=0;i<s_path.size();i++){
   	 s_path[i]+=((*it).property[i])*(*it).similarity;
    }
    partition+=(*it).similarity;
  }
  for(unsigned i=0;i<s_path.size();i++){ s_path[i]/=partition;  val_s_path[i]->set(s_path[i]) ;}
  val_z_path->set(-(1./lambda)*std::log(partition));
  for(unsigned j=0;j<s_path.size();j++){
    // clean up
    for(unsigned i=0;i< derivs_s.size();i++){derivs_s[i].zero();}
    // do the derivative 
    for(imgiter it=imgVec.begin();it!=imgVec.end();++it){ 
       double expval=(*it).similarity;
       tmp=lambda*expval*(s_path[j]-(*it).property[j])/partition;
       for(unsigned i=0;i< derivs_s.size();i++){ derivs_s[i]+=tmp*(*it).distder[i] ;} 
       if(j==0){for(unsigned i=0;i< derivs_z.size();i++){ derivs_z[i]+=(*it).distder[i]*expval/partition;}} 
    }
    for(unsigned i=0;i< derivs_s.size();i++){
          setAtomsDerivatives (val_s_path[j],i,derivs_s[i]); 
          if(j==0){setAtomsDerivatives (val_z_path,i,derivs_z[i]);} 
    }
  }
  for(unsigned i=0;i<val_s_path.size();++i) setBoxDerivativesNoPbc(val_s_path[i]);
  setBoxDerivativesNoPbc(val_z_path);
  //
  //  here set next round neighbors
  //
  if (neigh_size>0){
	if( int(getStep())%int(neigh_stride/getTimeStep())==0 ){
		// next round do it all:empty the vector	
		imgVec.clear();
        }
        // time to analyze the results: 
        if(imgVec.size()==nframes){
            //sort by msd
            sort(imgVec.begin(), imgVec.end(), imgOrderByDist()); 
            //resize
            imgVec.resize(neigh_size);
        } 
  }
  //log.printf("CALCULATION DONE! \n");
}

}

