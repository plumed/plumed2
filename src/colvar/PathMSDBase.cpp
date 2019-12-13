/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2019 The plumed team
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

namespace PLMD {
namespace colvar {

void PathMSDBase::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords(keys);
  keys.add("compulsory","LAMBDA","the lambda parameter is needed for smoothing, is in the units of plumed");
  keys.add("compulsory","REFERENCE","the pdb is needed to provide the various milestones");
  keys.add("optional","NEIGH_SIZE","size of the neighbor list");
  keys.add("optional","NEIGH_STRIDE","how often the neighbor list needs to be calculated in time units");
  keys.add("optional", "EPSILON", "(default=-1) the maximum distance between the close and the current structure, the positive value turn on the close structure method");
  keys.add("optional", "LOG_CLOSE", "(default=0) value 1 enables logging regarding the close structure");
  keys.add("optional", "DEBUG_CLOSE", "(default=0) value 1 enables extensive debugging info regarding the close structure, the simulation will run much slower");
}

PathMSDBase::PathMSDBase(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  nopbc(false),
  neigh_size(-1),
  neigh_stride(-1),
  epsilonClose(-1),
  debugClose(0),
  logClose(0),
  computeRefClose(false),
  nframes(0)
{
  parse("LAMBDA",lambda);
  parse("NEIGH_SIZE",neigh_size);
  parse("NEIGH_STRIDE",neigh_stride);
  parse("REFERENCE",reference);
  parse("EPSILON", epsilonClose);
  parse("LOG_CLOSE", logClose);
  parse("DEBUG_CLOSE", debugClose);
  parseFlag("NOPBC",nopbc);

  // open the file
  FILE* fp=fopen(reference.c_str(),"r");
  std::vector<AtomNumber> aaa;
  if (fp!=NULL)
  {
    log<<"Opening reference file "<<reference.c_str()<<"\n";
    bool do_read=true;
    while (do_read) {
      PDB mypdb;
      RMSD mymsd;
      do_read=mypdb.readFromFilepointer(fp,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength());
      if(do_read) {
        nframes++;
        if(mypdb.getAtomNumbers().size()==0) error("number of atoms in a frame should be more than zero");
        unsigned nat=mypdb.getAtomNumbers().size();
        if(nat!=mypdb.getAtomNumbers().size()) error("frames should have the same number of atoms");
        if(aaa.empty()) {
          aaa=mypdb.getAtomNumbers();
          log.printf("  found %z atoms in input \n",aaa.size());
          log.printf("  with indices : ");
          for(unsigned i=0; i<aaa.size(); ++i) {
            if(i%25==0) log<<"\n";
            log.printf("%d ",aaa[i].serial());
          }
          log.printf("\n");
        }
        if(aaa!=mypdb.getAtomNumbers()) error("frames should contain same atoms in same order");
        log<<"Found PDB: "<<nframes<<" containing  "<<mypdb.getAtomNumbers().size()<<" atoms\n";
        pdbv.push_back(mypdb);
        derivs_s.resize(mypdb.getAtomNumbers().size());
        derivs_z.resize(mypdb.getAtomNumbers().size());
        mymsd.set(mypdb,"OPTIMAL");
        msdv.push_back(mymsd); // the vector that stores the frames
      } else {break ;}
    }
    fclose (fp);
    log<<"Found TOTAL "<<nframes<< " PDB in the file "<<reference.c_str()<<" \n";
    if(nframes==0) error("at least one frame expected");
    //set up rmsdRefClose, initialize it to the first structure loaded from reference file
    rmsdPosClose.set(pdbv[0], "OPTIMAL");
    firstPosClose = true;
  }
  if(neigh_stride>0 || neigh_size>0) {
    if(neigh_size>int(nframes)) {
      log.printf(" List size required ( %d ) is too large: resizing to the maximum number of frames required: %u  \n",neigh_size,nframes);
      neigh_size=nframes;
    }
    log.printf("  Neighbor list enabled: \n");
    log.printf("                size   :  %d elements\n",neigh_size);
    log.printf("                stride :  %d timesteps \n",neigh_stride);
  } else {
    log.printf("  Neighbor list NOT enabled \n");
  }
  if (epsilonClose > 0) {
    log.printf(" Computing with the close structure, epsilon = %lf\n", epsilonClose);
    log << "  Bibliography " << plumed.cite("Pazurikova J, Krenek A, Spiwok V, Simkova M J. Chem. Phys. 146, 115101 (2017)") << "\n";
  }
  else {
    debugClose = 0;
    logClose = 0;
  }
  if (debugClose)
    log.printf(" Extensive debug info regarding close structure turned on\n");

  rotationRefClose.resize(nframes);
  savedIndices = vector<unsigned>(nframes);

  if(nopbc) log.printf("  without periodic boundary conditions\n");
  else      log.printf("  using periodic boundary conditions\n");

}

PathMSDBase::~PathMSDBase() {
}

void PathMSDBase::calculate() {

  if(neigh_size>0 && getExchangeStep()) error("Neighbor lists for this collective variable are not compatible with replica exchange, sorry for that!");

  //log.printf("NOW CALCULATE! \n");

  if(!nopbc) makeWhole();


  // resize the list to full
  if(imgVec.empty()) { // this is the signal that means: recalculate all
    imgVec.resize(nframes);
    #pragma omp simd
    for(unsigned i=0; i<nframes; i++) {
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

  std::vector<Tensor> tmp_rotationRefClose(nframes);

  if (epsilonClose > 0) {
    //compute rmsd between positions and close structure, save rotation matrix, drotation_drr01
    double posclose = rmsdPosClose.calc_Rot_DRotDRr01(getPositions(), rotationPosClose, drotationPosCloseDrr01, true);
    //if we compute for the first time or the existing close structure is too far from current structure
    if (firstPosClose || (posclose > epsilonClose)) {
      //set the current structure as close one for a few next steps
      if (logClose)
        log << "PLUMED_CLOSE: new close structure, rmsd pos close " << posclose << "\n";
      rmsdPosClose.clear();
      rmsdPosClose.setReference(getPositions());
      //as this is a new close structure, we need to save the rotation matrices fitted to the reference structures
      // and we need to accurately recalculate for all reference structures
      computeRefClose = true;
      imgVec.resize(nframes);
      for(unsigned i=0; i<nframes; i++) {
        imgVec[i].property=indexvec[i];
        imgVec[i].index=i;
      }
      firstPosClose = false;
    }
    else {
      //the current structure is pretty close to the close structure, so we use saved rotation matrices to decrease the complexity of rmsd comuptation
      if (debugClose)
        log << "PLUMED-CLOSE: old close structure, rmsd pos close " << posclose << "\n";
      computeRefClose = false;
    }
  }

  std::vector<double> tmp_distances(imgVec.size(),0.0);
  std::vector<Vector> tmp_derivs;
// this array is a merge of all tmp_derivs, so as to allow a single comm.Sum below
  std::vector<Vector> tmp_derivs2(imgVec.size()*nat);

// if imgVec.size() is less than nframes, it means that only some msd will be calculated
  if (epsilonClose > 0) {
    if (computeRefClose) {
      //recompute rotation matrices accurately
      for(unsigned i=rank; i<imgVec.size(); i+=stride) {
        tmp_distances[i] = msdv[imgVec[i].index].calc_Rot(getPositions(), tmp_derivs, tmp_rotationRefClose[imgVec[i].index], true);
        plumed_assert(tmp_derivs.size()==nat);
        #pragma omp simd
        for(unsigned j=0; j<nat; j++) tmp_derivs2[i*nat+j]=tmp_derivs[j];
      }
    }
    else {
      //approximate distance with saved rotation matrices
      for(unsigned i=rank; i<imgVec.size(); i+=stride) {
        tmp_distances[i] = msdv[imgVec[i].index].calculateWithCloseStructure(getPositions(), tmp_derivs, rotationPosClose, rotationRefClose[imgVec[i].index], drotationPosCloseDrr01, true);
        plumed_assert(tmp_derivs.size()==nat);
        #pragma omp simd
        for(unsigned j=0; j<nat; j++) tmp_derivs2[i*nat+j]=tmp_derivs[j];
        if (debugClose) {
          double withclose = tmp_distances[i];
          RMSD opt;
          opt.setType("OPTIMAL");
          opt.setReference(msdv[imgVec[i].index].getReference());
          vector<Vector> ders;
          double withoutclose = opt.calculate(getPositions(), ders, true);
          float difference = fabs(withoutclose-withclose);
          log<<"PLUMED-CLOSE: difference original "<<withoutclose;
          log<<" - with close "<<withclose<<" = "<<difference<<", step "<<getStep()<<", i "<<i<<" imgVec[i].index "<<imgVec[i].index<<"\n";
        }
      }
    }
  }
  else {
    // store temporary local results
    for(unsigned i=rank; i<imgVec.size(); i+=stride) {
      tmp_distances[i]=msdv[imgVec[i].index].calculate(getPositions(),tmp_derivs,true);
      plumed_assert(tmp_derivs.size()==nat);
      #pragma omp simd
      for(unsigned j=0; j<nat; j++) tmp_derivs2[i*nat+j]=tmp_derivs[j];
    }
  }

// reduce over all processors
  comm.Sum(tmp_distances);
  comm.Sum(tmp_derivs2);
  if (epsilonClose > 0 && computeRefClose) {
    comm.Sum(tmp_rotationRefClose);
    for (unsigned i=0; i<nframes; i++) {
      rotationRefClose[i] = tmp_rotationRefClose[i];
    }
  }
// assign imgVec[i].distance and imgVec[i].distder
  for(unsigned i=0; i<imgVec.size(); i++) {
    imgVec[i].distance=tmp_distances[i];
    imgVec[i].distder.assign(&tmp_derivs2[i*nat],nat+&tmp_derivs2[i*nat]);
  }

// END OF THE HEAVY PART

  vector<Value*> val_s_path;
  if(labels.size()>0) {
    for(unsigned i=0; i<labels.size(); i++) { val_s_path.push_back(getPntrToComponent(labels[i].c_str()));}
  } else {
    val_s_path.push_back(getPntrToComponent("sss"));
  }
  Value* val_z_path=getPntrToComponent("zzz");

  vector<double> s_path(val_s_path.size()); for(unsigned i=0; i<s_path.size(); i++)s_path[i]=0.;
  double partition=0.;
  double tmp;

  // clean vector
  for(unsigned i=0; i< derivs_z.size(); i++) {derivs_z[i].zero();}

  for(auto & it : imgVec) {
    it.similarity=exp(-lambda*(it.distance));
    //log<<"DISTANCE "<<(*it).distance<<"\n";
    for(unsigned i=0; i<s_path.size(); i++) {
      s_path[i]+=(it.property[i])*it.similarity;
    }
    partition+=it.similarity;
  }
  for(unsigned i=0; i<s_path.size(); i++) { s_path[i]/=partition;  val_s_path[i]->set(s_path[i]) ;}
  val_z_path->set(-(1./lambda)*std::log(partition));
  for(unsigned j=0; j<s_path.size(); j++) {
    // clean up
    #pragma omp simd
    for(unsigned i=0; i< derivs_s.size(); i++) {derivs_s[i].zero();}
    // do the derivative
    for(const auto & it : imgVec) {
      double expval=it.similarity;
      tmp=lambda*expval*(s_path[j]-it.property[j])/partition;
      #pragma omp simd
      for(unsigned i=0; i< derivs_s.size(); i++) { derivs_s[i]+=tmp*it.distder[i] ;}
      if(j==0) {
        #pragma omp simd
        for(unsigned i=0; i< derivs_z.size(); i++) { derivs_z[i]+=it.distder[i]*expval/partition;}
      }
    }
    for(unsigned i=0; i< derivs_s.size(); i++) {
      setAtomsDerivatives (val_s_path[j],i,derivs_s[i]);
      if(j==0) {setAtomsDerivatives (val_z_path,i,derivs_z[i]);}
    }
  }
  for(unsigned i=0; i<val_s_path.size(); ++i) setBoxDerivativesNoPbc(val_s_path[i]);
  setBoxDerivativesNoPbc(val_z_path);
  //
  //  here set next round neighbors
  //
  if (neigh_size>0) {
    //if( int(getStep())%int(neigh_stride/getTimeStep())==0 ){
    // enforce consistency: the stride is in time steps
    if( int(getStep())%int(neigh_stride)==0 ) {

      // next round do it all:empty the vector
      imgVec.clear();
    }
    // time to analyze the results:
    if(imgVec.size()==nframes) {
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
