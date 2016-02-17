/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016 The plumed team
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
#include "Colvar.h"
#include "ActionRegister.h"
#include "core/PlumedMain.h"
#include "tools/Matrix.h"
#include "core/SetupMolInfo.h"
#include "core/ActionSet.h"
#include "tools/File.h"

#include <string>
#include <cmath>
#include <unordered_map>
#include <numeric>

using namespace std;

namespace PLMD{
namespace colvar{

//+PLUMEDOC COLVAR OVERLAP
/*
Put Documentation here 

*/
//+ENDPLUMEDOC
   
class Overlap : public Colvar {

private:
 // model GMM
 vector<double> GMM_m_w_;
 vector< Matrix <double> > GMM_m_cov_;
 // data GMM
 vector<Vector> GMM_d_m_;
 vector<double> GMM_d_w_;
 vector< Matrix <double> > GMM_d_cov_;
 // temporary vectors
 vector<double> ovmd_;
 vector< vector< Vector > > ovmd_der_;
 
 void get_GMM_m(vector<AtomNumber> &atoms);
 void get_GMM_d(string gmm_file);
 void normalize_GMM(vector<double> &w);
 Vector get_overlap(Vector m_m, double m_w, Matrix <double> m_cov,
                    Vector d_m, double d_w, Matrix <double> d_cov,
                    double &ov, bool do_der=false);

public:
  static void registerKeywords( Keywords& keys );
  explicit Overlap(const ActionOptions&);
// active methods:
  virtual void calculate();
};

PLUMED_REGISTER_ACTION(Overlap,"OVERLAP")

void Overlap::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords( keys );
  keys.add("atoms","ATOMS","atoms for which we calculate the density map");
  keys.add("compulsory","GMM_FILE","file with the parameters of the GMM components");
  componentsAreNotOptional(keys);
  keys.addOutputComponent("ovmd", "COMPONENTS","overlap of the model with individual data GMM components");
  keys.addOutputComponent("ovdd", "COMPONENTS","overlap between individual data GMM components");
}

Overlap::Overlap(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao)

{
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  
  string GMM_file;
  parse("GMM_FILE",GMM_file);
  
  checkRead();

  log.printf("  atoms involved : ");
  for(unsigned i=0;i<atoms.size();++i) log.printf("%d ",atoms[i].serial());
  log.printf("\n");
  log.printf("  GMM data file : %s\n", GMM_file.c_str());

  // calculate model GMM constant parameters
  get_GMM_m(atoms);

  // read data GMM parameters
  get_GMM_d(GMM_file);
  log.printf("  number of GMM components : %d\n", GMM_d_m_.size());
  
  // add components
  for(unsigned i=0;i<GMM_d_m_.size();++i) {
      string num; Tools::convert(i,num);
      addComponentWithDerivatives("ovmd"+num); componentIsNotPeriodic("ovmd"+num);
      addComponent("ovdd"+num);                componentIsNotPeriodic("ovdd"+num);
  }
  
  // normalize GMMs
  normalize_GMM(GMM_m_w_);
  normalize_GMM(GMM_d_w_);
  
  // get overlap between data GMM components
  double ov;
  for(unsigned i=0;i<GMM_d_m_.size();++i) {
      get_overlap(GMM_d_m_[i], GMM_d_w_[i], GMM_d_cov_[i],
                  GMM_d_m_[i], GMM_d_w_[i], GMM_d_cov_[i],
                  ov);
      string num; Tools::convert(i,num);
      Value* value=getPntrToComponent("ovdd"+num);
      value->set(ov);
  }
  
  // initialize temporary vectors
  for(unsigned i=0;i<GMM_d_w_.size(); ++i){
     ovmd_.push_back(0.0);
     vector<Vector> der_tmp;
     for(unsigned j=0;j<getNumberOfAtoms();j++) der_tmp.push_back(Vector(0,0,0));
     ovmd_der_.push_back(der_tmp);
  }
     
  // request the atoms
  requestAtoms(atoms);
}

void Overlap::get_GMM_m(vector<AtomNumber> &atoms)
{
  vector<SetupMolInfo*> moldat=plumed.getActionSet().select<SetupMolInfo*>();
  Matrix <double> cov(3,3);

  // map of atom types to scattering radii - in nm
  unordered_map<string, double> radii_map;
  radii_map["C"] = 1.0;
  radii_map["O"] = 1.0;
  radii_map["N"] = 1.0;
  radii_map["H"] = 1.0;
  radii_map["S"] = 1.0;
  
  // check if MOLINFO line is present 
  if( moldat.size()==1 ){
    log<<"  MOLINFO DATA found, using proper atom names\n";
    for(unsigned i=0;i<atoms.size();++i){
      // get atom name
      string name = moldat[0]->getAtomName(atoms[i]);
      char type;
      // get atom type
      char first = name.at(0);
      // GOLDEN RULE: type is first letter, if not a number
      if (!isdigit(first)){
         type = first;
      // otherwise is the second
      } else {
         type = name.at(1);
      }
      // check if key in map
      std::string type_s = std::string(1,type);
      if(radii_map.find(type_s) != radii_map.end()){
        cov(0,0)=radii_map[type_s]; cov(0,1)=0.0; cov(0,2)=0.0;
        cov(1,0)=0.0; cov(1,1)=radii_map[type_s]; cov(1,2)=0.0;
        cov(2,0)=0.0; cov(2,1)=0.0; cov(2,2)=radii_map[type_s];
        GMM_m_cov_.push_back(cov);
        GMM_m_w_.push_back(1.0); 
      } else {
        error("Wrong atom type "+type_s+" from atom name "+name+"\n"); 
      }
    }
  } else {
    error("MOLINFO DATA not found\n");
  }
}

// read GMM data file in PLUMED format:
void Overlap::get_GMM_d(string GMM_file)
{
 int idcomp;
 double w,m0,m1,m2,c00,c01,c02,c10,c11,c12,c20,c21,c22;
 
 // open file
 IFile *ifile = new IFile();
 if(ifile->FileExist(GMM_file)){
    ifile->open(GMM_file);
    while(ifile->scanField("Id",idcomp)){
     ifile->scanField("Weight",w);
     ifile->scanField("Mean_0",m0);
     ifile->scanField("Mean_1",m1);
     ifile->scanField("Mean_2",m2);
     ifile->scanField("Cov_00",c00);
     ifile->scanField("Cov_01",c01);
     ifile->scanField("Cov_02",c02);
     ifile->scanField("Cov_10",c10);
     ifile->scanField("Cov_11",c11);
     ifile->scanField("Cov_12",c12);
     ifile->scanField("Cov_20",c20);
     ifile->scanField("Cov_21",c21);
     ifile->scanField("Cov_22",c22);
     // center of the Gaussian
     GMM_d_m_.push_back(Vector(m0,m1,m2));
     // covariance matrix
     Matrix <double> cov(3,3);
     cov(0,0)=c00; cov(0,1)=c01; cov(0,2)=c02;
     cov(1,0)=c10; cov(1,1)=c11; cov(1,2)=c12;
     cov(2,0)=c20; cov(2,1)=c21; cov(2,2)=c22;
     GMM_d_cov_.push_back(cov);
     // weight
     GMM_d_w_.push_back(w);
     // new line
     ifile->scanField();
    }
    ifile->close();
 } else {
    error("Cannot find GMM_FILE "+GMM_file+"\n"); 
 }
 delete ifile;

}

// normalize GMM to sum to 1
// since all the GMM components are individually normalized, we just need to 
// divide each weight for the sum of the weights
void Overlap::normalize_GMM(vector<double> &w)
 {
   double norm = accumulate(w.begin(), w.end(), 0.0);
   for(unsigned i=0; i<w.size(); ++i) w[i] /= norm;
 }

Vector Overlap::get_overlap(Vector m_m, double m_w, Matrix <double> m_cov,
                            Vector d_m, double d_w, Matrix <double> d_cov,
                            double &ov, bool do_der)
{
  Vector ov_der = Vector(0,0,0);
  
  return ov_der;
}

// calculator
void Overlap::calculate(){

  makeWhole();

  // clean temporary vectors  
  for(unsigned i=0;i<GMM_d_w_.size(); ++i){
    ovmd_[i]=0.0;
    for(unsigned j=0; j<getNumberOfAtoms(); ++j) ovmd_der_[i][j]=Vector(0,0,0);
  }
    
  // we have to cycle over all model and data GMM components - MPI 
  for(unsigned i=0;i<GMM_d_w_.size()*getNumberOfAtoms();++i) {
      // get indexes of data and model component
      unsigned id = 1;
      unsigned im = 1;
      // get overlap
      double ov = 0.0;
      Vector ov_der = get_overlap(getPosition(im), GMM_m_w_[im], GMM_m_cov_[im],
                                     GMM_d_m_[id], GMM_d_w_[id], GMM_d_cov_[id],
                                     ov, true);
      ovmd_[id] += ov;
      // derivatives
      ovmd_der_[id][im] += ov_der;
  }
  // MPI summation
  
  // put values and derivatives
  for(unsigned i=0;i<GMM_d_w_.size(); ++i) {
     string num; Tools::convert(i,num);
     Value* value=getPntrToComponent("ovmd"+num);
     value->set(ovmd_[i]);
     for(unsigned j=0;j<getNumberOfAtoms();j++){
       setAtomsDerivatives (value,j,ovmd_der_[i][j]);
     }
  }
  
}

}
}



