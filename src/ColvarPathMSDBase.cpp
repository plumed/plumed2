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

class ColvarPathMSDBase : public Colvar {
/// this class is a general container for path stuff 
  class ImagePath {
     public:
        // cardinal indexing: needed to map over msd 
        unsigned index;
        // spiwok indexing
        vector<double> property;
        // distance
        double distance;
        // similarity (exp - lambda distance) or other
        double similarity;
        // derivatives of the distance
        vector<Vector> distder;
        // here one can add a pointer to a value (hypothetically providing a distance from a point) 
  };
  struct imgOrderByDist {
       bool operator ()(ImagePath const& a, ImagePath const& b) {
           return (a).distance < (b).distance;
       };
  };
  struct imgOrderBySimilarity {
       bool operator ()(ImagePath const& a, ImagePath const& b) {
           return (a).similarity > (b).similarity;
       };
  };

  double lambda;
  bool pbc;
  int neigh_size,nframes;
  unsigned propertypos; 
  double neigh_stride;
  vector<PDB> pdbv;
  vector<RMSD> msdv;
  vector<string> labels;
  string reference;
  vector< vector<double> > indexvec; // use double to allow isomaps
  vector<Vector> derivs_s;
  vector<Vector> derivs_z;
  vector< vector <Vector> > derivs_v;
  vector <ImagePath> imgVec; // this can be used for doing neighlist   
public:
  ColvarPathMSDBase(const ActionOptions&);
// active methods:
  virtual void calculate();
//  virtual void prepare();
  static void registerKeywords(Keywords& keys);
};

void ColvarPathMSDBase::registerKeywords(Keywords& keys){
  Colvar::registerKeywords(keys);
  keys.add("compulsory","LAMBDA","the lambda parameter is needed for smoothing ");
  keys.add("compulsory","REFERENCE","the pdb is needed to provide the various milestones");
  keys.add("optional","NEIGH_SIZE","size of the neighbor list");
  keys.add("optional","NEIGH_STRIDE","how often the neighbor list needs to be calculated");
  keys.add("optional","PROPERTY","the property to be used in the indexing: this goes in the REMARK field of the reference");
}

ColvarPathMSDBase::ColvarPathMSDBase(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
pbc(true),
neigh_size(-1),
neigh_stride(-1.),
nframes(0),
propertypos(0)
{
  parse("LAMBDA",lambda);
  parse("NEIGH_SIZE",neigh_size);
  parse("NEIGH_STRIDE",neigh_stride);
  parse("REFERENCE",reference);
  parseVector("PROPERTY",labels);
  //parse("PROPERTY",propertypos);
  checkRead();

  // note: first add the derivative and then request the atoms 

  if(labels.size()==0){
     addComponentWithDerivatives("sss"); componentIsNotPeriodic("sss");
     addComponentWithDerivatives("zzz"); componentIsNotPeriodic("zzz");
     //labels.push_back("s");labels.push_back("z")
  }else{
      for(unsigned i=0;i<labels.size();i++){
	log<<" found custom propety to be found in the REMARK line: "<<labels[i].c_str()<<"\n";
        addComponentWithDerivatives(labels[i].c_str()); componentIsNotPeriodic(labels[i].c_str());
      }
      // add distance anyhow
      addComponentWithDerivatives("zzz"); componentIsNotPeriodic("zzz");
  }

  // parse all the pdb
  // open the file
  FILE* fp=fopen(reference.c_str(),"r");
  if (fp!=NULL)
  {
    log<<"Opening reference file "<<reference.c_str()<<"\n";
    bool do_read=true;
    while (do_read){
         PDB mypdb; 
         RMSD mymsd(log); 
         do_read=mypdb.readFromFilepointer(fp,1.,1.);
         if(do_read){
            nframes++;
            log<<"Found PDB: "<<nframes<<" atoms: "<<mypdb.getAtomNumbers().size()<<"\n"; 
	    pdbv.push_back(mypdb); 
            requestAtoms(mypdb.getAtomNumbers()); 
            derivs_s.resize(getNumberOfAtoms());
            derivs_z.resize(getNumberOfAtoms());
            mymsd.set(mypdb,"OPTIMAL");
            msdv.push_back(mymsd); // the vector that stores the frames
            //log<<mypdb; 
         }else{break ;}
    }
    fclose (fp);
    log<<"Found TOTAL "<<nframes<< " PDB in the file "<<reference.c_str()<<" \n"; 
  } 
  // do some neighbor printout
  if(labels.size()!=0){
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
  }else{ // normal case: only one label
	  double i=1.;
	  for(unsigned it=0 ;it<nframes ;++it){
                        vector<double> v; v.push_back(i);
			indexvec.push_back(v);i+=1.; 
	  }
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

};

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
          if(j==0)setAtomsDerivatives (val_z_path,i,derivs_z[i]); 
    }
  }
  //
  //  here set next round neighbors
  //
  if (neigh_size>0){
	if( int(getStep())%int(neigh_stride/getTimeStep())==0 ){
		//log<<"CLEARING: NEXT STEP IS FULL\n ";
		// next round do it all:empty the vector	
		imgVec.clear();
        }
        // time to analyze the results: 
        if(imgVec.size()==nframes){
            //for(unsigned i=0;i<nframes;i++){
	    //    log<<"OLDERR "<<imgVec[i].distance<<"\n";
            //}
            //sort by msd
            sort(imgVec.begin(), imgVec.end(), imgOrderByDist()); 
            //resize
            //for(unsigned i=0;i<nframes;i++){
	    //    log<<"NEWERR "<<imgVec[i].distance<<"\n";
            //}
            imgVec.resize(neigh_size);
        } 
  }
  //log.printf("CALCULATION DONE! \n");
};

}

