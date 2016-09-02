/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2016 The plumed team
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
#include "Colvar.h"
#include "ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/PDB.h"
#include "tools/RMSD.h"
#include "tools/Tools.h"
#include <cmath>

using namespace std;

namespace PLMD{
namespace colvar{

void PathMSDBase::registerKeywords(Keywords& keys){
  Colvar::registerKeywords(keys);
  keys.remove("NOPBC");
  keys.add("compulsory","LAMBDA","the lambda parameter is needed for smoothing, is in the units of plumed");
  keys.add("compulsory","REFERENCE","the pdb is needed to provide the various milestones");
  keys.add("optional","NEIGH_SIZE","size of the neighbor list");
  keys.add("optional","NEIGH_STRIDE","how often the neighbor list needs to be calculated in time units");
}

PathMSDBase::PathMSDBase(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
neigh_size(-1),
neigh_stride(-1),
nframes(0)
{
  parse("LAMBDA",lambda);
  parse("NEIGH_SIZE",neigh_size);
  parse("NEIGH_STRIDE",neigh_stride);
  parse("REFERENCE",reference);

  // open the file
  FILE* fp=fopen(reference.c_str(),"r");
  std::vector<AtomNumber> aaa;
  if (fp!=NULL)
  {
    log<<"Opening reference file "<<reference.c_str()<<"\n";
    bool do_read=true;
    while (do_read){
         PDB mypdb; 
         RMSD mymsd; 
         do_read=mypdb.readFromFilepointer(fp,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength());
         if(do_read){
            nframes++;
            if(mypdb.getAtomNumbers().size()==0) error("number of atoms in a frame should be more than zero");
            unsigned nat=mypdb.getAtomNumbers().size();
            if(nat!=mypdb.getAtomNumbers().size()) error("frames should have the same number of atoms");
            if(aaa.empty()) aaa=mypdb.getAtomNumbers();
            if(aaa!=mypdb.getAtomNumbers()) error("frames should contain same atoms in same order");
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
    if(nframes==0) error("at least one frame expected");
  } 
  if(neigh_stride>0 || neigh_size>0){
           if(neigh_size>int(nframes)){
           	log.printf(" List size required ( %d ) is too large: resizing to the maximum number of frames required: %u  \n",neigh_size,nframes);
 		neigh_size=nframes;
           }
           log.printf("  Neighbor list enabled: \n");
           log.printf("                size   :  %d elements\n",neigh_size);
           log.printf("                stride :  %d timesteps \n",neigh_stride);
  }else{
           log.printf("  Neighbor list NOT enabled \n");
  }

}

void PathMSDBase::calculate(){

  if(neigh_size>0 && getExchangeStep()) error("Neighbor lists for this collective variable are not compatible with replica exchange, sorry for that!");

  //log.printf("NOW CALCULATE! \n");


  // resize the list to full
  if(imgVec.empty()){ // this is the signal that means: recalculate all 
      imgVec.resize(nframes);  
      for(unsigned i=0;i<nframes;i++){
          imgVec[i].property=indexvec[i];
          imgVec[i].index=i;
      }
  }

// THIS IS THE HEAVY PART (RMSD STUFF)
  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  unsigned nat=pdbv[0].size();
  plumed_assert(nat>0);
  plumed_assert(nframes>0);
  plumed_assert(imgVec.size()>0);

  std::vector<double> tmp_distances(imgVec.size(),0.0);
  std::vector<Vector> tmp_derivs;
// this array is a merge of all tmp_derivs, so as to allow a single comm.Sum below
  std::vector<Vector> tmp_derivs2(imgVec.size()*nat);

// if imgVec.size() is less than nframes, it means that only some msd will be calculated
  for(unsigned i=rank;i<imgVec.size();i+=stride){
// store temporary local results
    tmp_distances[i]=msdv[imgVec[i].index].calculate(getPositions(),tmp_derivs,true);
    plumed_assert(tmp_derivs.size()==nat);
    for(unsigned j=0;j<nat;j++) tmp_derivs2[i*nat+j]=tmp_derivs[j];
  }
// reduce over all processors
  comm.Sum(tmp_distances);
  comm.Sum(tmp_derivs2);
// assign imgVec[i].distance and imgVec[i].distder
  for(unsigned i=0;i<imgVec.size();i++){
    imgVec[i].distance=tmp_distances[i];
    imgVec[i].distder.assign(&tmp_derivs2[i*nat],nat+&tmp_derivs2[i*nat]);
  }

// END OF THE HEAVY PART

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

  // clean vector
  for(unsigned i=0;i< derivs_z.size();i++){derivs_z[i].zero();}

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
	//if( int(getStep())%int(neigh_stride/getTimeStep())==0 ){
	// enforce consistency: the stride is in time steps
	if( int(getStep())%int(neigh_stride)==0 ){

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

}
