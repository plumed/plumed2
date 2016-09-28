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
#include <map>
#include <numeric>

using namespace std;

namespace PLMD{
namespace colvar{

//+PLUMEDOC COLVAR EM3DMAP
/*
Put Documentation here 

*/
//+ENDPLUMEDOC
   
class EM3Dmap : public Colvar {

private:

 // temperature in kbt
 double kbt_;
 // model GMM - weights and covariances
 vector<double>           GMM_m_w_;
 vector< Matrix<double> > GMM_m_cov_;
 // data GMM - means, weights, and covariances
 vector<Vector>           GMM_d_m_;
 vector<double>           GMM_d_w_;
 vector< Matrix<double> > GMM_d_cov_;
 // overlaps 
 vector<double> ovmd_;
 vector<double> ovdd_; 
 // and derivatives
 vector<Vector> ovmd_der_;
 vector<Vector> atom_der_;
 vector<double> ene_der_;
 // constant quantity;
 double cfact_;
 
 // list of prefactors for overlap between two components of model and data GMM
 // fact_md = w_m * w_d / (2pi)**1.5 / sqrt(det_md)
 vector< double > fact_md_;
 // inverse of the sum of model and data covariances matrices
 vector< Matrix<double> > inv_cov_md_;
 // neighbor list
 double   nl_cutoff_;
 unsigned nl_stride_;
 bool first_time_;
 vector < pair<unsigned, unsigned > > nl_;
 // parallel stuff
 bool serial_;
 unsigned size_;
 unsigned rank_;
 
 // calculate model GMM weights and covariances - these are constants
 void get_GMM_m(vector<AtomNumber> &atoms);
 // read data GMM file
 void get_GMM_d(string gmm_file);
 // normalize GMM
 void normalize_GMM(vector<double> &w);

 // get fact_md and inv_cov_md
 double get_prefactor_inverse (Matrix<double> &GMM_cov_0, Matrix<double> &GMM_cov_1,
        double &GMM_w_0, double &GMM_w_1, Matrix<double> &inv_sum);
 // calculate self overlaps between data GMM components - ovdd_
 double get_self_overlap(unsigned id);
 // calculate overlap between two components
 double get_overlap(Vector &m_m, Vector d_m, double &fact_md,
                    Matrix<double> &inv_cov_md, Vector &ov_der);
 double get_overlap(Vector &m_m, Vector d_m, double &fact_md,
                    Matrix<double> &inv_cov_md);
 // update the neighbor list
 void update_neighbor_list();
 // calculate overlap
 void calculate_overlap();

public:
  static void registerKeywords( Keywords& keys );
  explicit EM3Dmap(const ActionOptions&);
// active methods:
  virtual void calculate();
};

PLUMED_REGISTER_ACTION(EM3Dmap,"EM3DMAP")

void EM3Dmap::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords( keys );
  keys.add("atoms","ATOMS","atoms for which we calculate the density map");
  keys.add("compulsory","GMM_FILE","file with the parameters of the GMM components");
  keys.add("compulsory","TEMP","temperature in energy units");
  keys.addFlag("SERIAL",false,"perform the calculation in serial - for debug purpose");
  keys.add("optional","NL_CUTOFF","The cutoff in overlap for the neighbor list");
  keys.add("optional","NL_STRIDE","The frequency with which we are updating the neighbor list");
  componentsAreNotOptional(keys);
  keys.addOutputComponent("ovmd", "COMPONENTS","overlap of the model with individual data GMM components");
  keys.addOutputComponent("ovdd", "COMPONENTS","overlap between individual data GMM components");
}

EM3Dmap::EM3Dmap(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
nl_cutoff_(-1.0), nl_stride_(0),
first_time_(true), serial_(false)
{
  
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  
  string GMM_file;
  parse("GMM_FILE",GMM_file);
 
  parse("TEMP",kbt_);
 
  // neighbor list stuff
  parse("NL_CUTOFF",nl_cutoff_);
  if(nl_cutoff_<=0.0) error("NL_CUTOFF should be explicitly specified and positive");
  parse("NL_STRIDE",nl_stride_);
  if(nl_stride_<=0) error("NL_STRIDE should be explicitly specified and positive");
  
  // serial or parallel
  parseFlag("SERIAL",serial_);
  if(serial_){
    size_=1; rank_=0;
  } else {
    size_=comm.Get_size(); rank_=comm.Get_rank();
  }
  
  checkRead();

  log.printf("  atoms involved : ");
  for(unsigned i=0;i<atoms.size();++i) log.printf("%d ",atoms[i].serial());
  log.printf("\n");
  log.printf("  GMM data file : %s\n", GMM_file.c_str());
  if(serial_) log.printf("  serial calculation\n");
  log.printf("  neighbor list overlap cutoff : %lf\n", nl_cutoff_);
  log.printf("  neighbor list stride : %u\n",  nl_stride_);

  // set constant quantity before calculating stuff
  cfact_ = 1.0/pow( 2.0*pi, 1.5 );
  
  // calculate model GMM constant parameters
  get_GMM_m(atoms);

  // read data GMM parameters
  get_GMM_d(GMM_file);
  log.printf("  number of GMM components : %u\n", static_cast<unsigned>(GMM_d_m_.size()));
  
  // normalize GMMs
  normalize_GMM(GMM_m_w_);
  normalize_GMM(GMM_d_w_);
  
  // get self overlaps between data GMM components
  for(unsigned i=0;i<GMM_d_w_.size();++i) {
      double ov = get_self_overlap(i);
      ovdd_.push_back(ov);
  }

  // and prepare temporary vectors
  ovmd_.resize(GMM_d_w_.size());
  ene_der_.resize(GMM_d_w_.size());
  atom_der_.resize(GMM_m_w_.size());

  // request the atoms
  requestAtoms(atoms);

  // add value
  addValueWithDerivatives(); setNotPeriodic();
     
}

void EM3Dmap::get_GMM_m(vector<AtomNumber> &atoms)
{
  vector<SetupMolInfo*> moldat=plumed.getActionSet().select<SetupMolInfo*>();
  Matrix <double> cov(3,3);

  // map of atom types to A and B coefficients of scattering factor
  // f(s) = A * exp(-B*s**2)
  // B is in Angstrom squared
  // data from here
  // http://fg.oisin.rc-harwell.ac.uk/scm/loggerhead/cctbx/cctbx/view/head:/cctbx/eltbx/xray_scattering/n_gaussian_raw.cpp
  map<string, double> A_map, B_map;
  A_map["C"] = 5.96792806111; B_map["C"] = 14.8957682987;
  A_map["O"] = 7.9652690671;  B_map["O"] = 9.0526662027;
  A_map["N"] = 6.96715024214; B_map["N"] = 11.4372299305;
  A_map["S"] = 15.911119329;  B_map["S"] = 10.8469011094;
  
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
      if(A_map.find(type_s) != A_map.end()){
        // convert to sigma in nm
        // the Gaussian in density (real) space is the FT of scattering factor
        // f(r) = A * (pi/B)**1.5 * exp(-pi**2/B*r**2)
        double s = sqrt ( 0.5 * B_map[type_s] ) / pi * 0.1;
        // covariance matrix for spherical Gaussian
        cov(0,0)=s*s; cov(0,1)=0.0; cov(0,2)=0.0;
        cov(1,0)=0.0; cov(1,1)=s*s; cov(1,2)=0.0;
        cov(2,0)=0.0; cov(2,1)=0.0; cov(2,2)=s*s;
        GMM_m_cov_.push_back(cov);
        // this will be normalized to 1 in the final density
        GMM_m_w_.push_back(A_map[type_s]); 
      } else {
        error("Wrong atom type "+type_s+" from atom name "+name+"\n"); 
      }
    }
  } else {
    error("MOLINFO DATA not found\n");
  }
}

// read GMM data file in PLUMED format:
void EM3Dmap::get_GMM_d(string GMM_file)
{
 int idcomp;
 double w, m0, m1, m2;
 Matrix <double> cov(3,3);
 
 // open file
 IFile *ifile = new IFile();
 if(ifile->FileExist(GMM_file)){
    ifile->open(GMM_file);
    while(ifile->scanField("Id",idcomp)){
     ifile->scanField("Weight",w);
     ifile->scanField("Mean_0",m0);
     ifile->scanField("Mean_1",m1);
     ifile->scanField("Mean_2",m2);
     ifile->scanField("Cov_00",cov(0,0));
     ifile->scanField("Cov_01",cov(0,1));
     ifile->scanField("Cov_02",cov(0,2));
     ifile->scanField("Cov_10",cov(1,0));
     ifile->scanField("Cov_11",cov(1,1));
     ifile->scanField("Cov_12",cov(1,2));
     ifile->scanField("Cov_20",cov(2,0));
     ifile->scanField("Cov_21",cov(2,1));
     ifile->scanField("Cov_22",cov(2,2));
     // center of the Gaussian
     GMM_d_m_.push_back(Vector(m0,m1,m2));
     // covariance matrix
     GMM_d_cov_.push_back(cov);
     // weights
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
void EM3Dmap::normalize_GMM(vector<double> &w)
 {
   double norm = accumulate(w.begin(), w.end(), 0.0);
   for(unsigned i=0; i<w.size(); ++i) w[i] /= norm;
 }

// get prefactors
double EM3Dmap::get_prefactor_inverse
(Matrix<double> &GMM_cov_0, Matrix<double> &GMM_cov_1,
           double &GMM_w_0, double &GMM_w_1, Matrix<double> &inv_sum)
{
 // we need the sum of the covariance matrices
 Matrix<double> sum_0_1(3,3);
 for(unsigned k1=0; k1<3; ++k1){ 
  for(unsigned k2=0; k2<3; ++k2){ 
     sum_0_1[k1][k2] = GMM_cov_0[k1][k2] + GMM_cov_1[k1][k2];
  }
 }   
 // and to calculate its determinant
 double det = sum_0_1[0][0]*(sum_0_1[1][1]*sum_0_1[2][2]-sum_0_1[1][2]*sum_0_1[2][1]);
       det -= sum_0_1[0][1]*(sum_0_1[1][0]*sum_0_1[2][2]-sum_0_1[1][2]*sum_0_1[2][0]);
       det += sum_0_1[0][2]*(sum_0_1[1][0]*sum_0_1[2][1]-sum_0_1[1][1]*sum_0_1[2][0]);
 // the prefactor is 
 double pre_fact =  cfact_ / sqrt(det) * GMM_w_0 * GMM_w_1;
 // and its inverse
 inv_sum[0][0] = (sum_0_1[1][1]*sum_0_1[2][2] - sum_0_1[1][2]*sum_0_1[2][1])/det;
 inv_sum[0][1] = (sum_0_1[0][2]*sum_0_1[2][1] - sum_0_1[0][1]*sum_0_1[2][2])/det;
 inv_sum[0][2] = (sum_0_1[0][1]*sum_0_1[1][2] - sum_0_1[0][2]*sum_0_1[1][1])/det;
 inv_sum[1][0] = (sum_0_1[1][2]*sum_0_1[2][0] - sum_0_1[1][0]*sum_0_1[2][2])/det;
 inv_sum[1][1] = (sum_0_1[0][0]*sum_0_1[2][2] - sum_0_1[0][2]*sum_0_1[2][0])/det;
 inv_sum[1][2] = (sum_0_1[0][2]*sum_0_1[1][0] - sum_0_1[0][0]*sum_0_1[1][2])/det;
 inv_sum[2][0] = (sum_0_1[1][0]*sum_0_1[2][1] - sum_0_1[1][1]*sum_0_1[2][0])/det;
 inv_sum[2][1] = (sum_0_1[0][1]*sum_0_1[2][0] - sum_0_1[0][0]*sum_0_1[2][1])/det;
 inv_sum[2][2] = (sum_0_1[0][0]*sum_0_1[1][1] - sum_0_1[0][1]*sum_0_1[1][0])/det;
 // return pre-factor
 return pre_fact;
}


double EM3Dmap::get_self_overlap(unsigned id)
{

 double ov = 0.0;
 Matrix<double> inv_sum_id_i(3,3);
 // start loop
 for(unsigned i=0; i<GMM_d_w_.size(); ++i){
   // call auxiliary method
   double pre_fact = get_prefactor_inverse(GMM_d_cov_[id], GMM_d_cov_[i], 
                                             GMM_d_w_[id],   GMM_d_w_[i], inv_sum_id_i); 
   // calculate overlap
   double ov_tmp = get_overlap(GMM_d_m_[id], GMM_d_m_[i], pre_fact, inv_sum_id_i);
   // add to overlap
   ov += ov_tmp;
 }
 return ov;
}

double EM3Dmap::get_overlap(Vector &m_m, Vector d_m, double &fact_md,
                            Matrix<double> &inv_cov_md, Vector &ov_der)
{
  // calculate vector difference m_m-d_m
  Vector md = m_m - d_m;
  // calculate its transpose
  Matrix<double> md_t(1,3);
  for(unsigned i=0; i<3; ++i) md_t[0][i] = md[i];
  // calculate product of md_t and inv_cov_md 
  Matrix<double> prod(1,3);
  mult(md_t, inv_cov_md, prod);
  // calculate product of prod and md
  double ov = 0.0;
  for(unsigned i=0; i<3; ++i) ov += prod[0][i]*md[i];
  // final calculation
  ov = fact_md * exp(-0.5*ov);
  // derivatives
  double x = md[0]*inv_cov_md[0][0] + md[1]*inv_cov_md[0][1] + md[2]*inv_cov_md[0][2];
  double y = md[0]*inv_cov_md[0][1] + md[1]*inv_cov_md[1][1] + md[2]*inv_cov_md[1][2];
  double z = md[0]*inv_cov_md[0][2] + md[1]*inv_cov_md[1][2] + md[2]*inv_cov_md[2][2];
  ov_der = ov * Vector(-x, -y, -z); 
  return ov;
}

double EM3Dmap::get_overlap(Vector &m_m, Vector d_m, double &fact_md,
                            Matrix<double> &inv_cov_md)
{
  // calculate vector difference m_m-d_m
  Vector md = m_m - d_m;
  // calculate its transpose
  Matrix<double> md_t(1,3);
  for(unsigned i=0; i<3; ++i) md_t[0][i] = md[i];
  // calculate product of md_t and inv_cov_md 
  Matrix<double> prod(1,3);
  mult(md_t, inv_cov_md, prod);
  // calculate product of prod and md
  double ov = 0.0;
  for(unsigned i=0; i<3; ++i) ov += prod[0][i]*md[i];
  // final calculation
  ov = fact_md * exp(-0.5*ov);
  return ov;
}

void EM3Dmap::update_neighbor_list()
{
  // temporary stuff
  Matrix<double> inv_sum_i_j(3,3);
  // clear old neighbor list and auxiliary vectors
  nl_.clear(); fact_md_.clear(); inv_cov_md_.clear();
  // cycle on all overlaps
  for(unsigned i=0; i<GMM_d_w_.size(); ++i){
   for(unsigned j=0; j<GMM_m_w_.size(); ++j){
      // call auxiliary method 
      double pre_fact = get_prefactor_inverse(GMM_d_cov_[i], GMM_m_cov_[j], 
                                                GMM_d_w_[i],   GMM_m_w_[j], inv_sum_i_j);	
      // calculate overlap
      double ov = get_overlap(GMM_d_m_[i], getPosition(j), pre_fact, inv_sum_i_j);
      // fill the neighbor list and auxiliary vectors
      if(ov >= nl_cutoff_){
        nl_.push_back(make_pair(i,j));
        fact_md_.push_back(pre_fact);
        inv_cov_md_.push_back(inv_sum_i_j);
      }
   }
  }
  // now resize derivatives
  ovmd_der_.resize(nl_.size());
}

// overlap calculator
void EM3Dmap::calculate_overlap(){

  //makeWhole();
  
  // update neighbor list ?
  if(first_time_ || getStep()%nl_stride_==0){
     update_neighbor_list();
     first_time_=false;
  }
  
  // clean temporary vectors
  for(unsigned i=0; i<ovmd_.size(); ++i)     ovmd_[i] = 0.0;
  for(unsigned i=0; i<ovmd_der_.size(); ++i) ovmd_der_[i] = Vector(0,0,0);
  
  // we have to cycle over all model and data GMM components in the neighbor list
  for(unsigned i=rank_;i<nl_.size();i=i+size_) {
      // get indexes of data and model component
      unsigned id = nl_[i].first;
      unsigned im = nl_[i].second;
      // add overlap with im component of model GMM
      ovmd_[id] += get_overlap(GMM_d_m_[id], getPosition(im), fact_md_[i],
                               inv_cov_md_[i], ovmd_der_[i]);
  }
  // if parallel, communicate stuff
  if(!serial_){
   comm.Sum(&ovmd_[0], ovmd_.size());
   comm.Sum(&ovmd_der_[0][0], 3*ovmd_der_.size());
  }
}

void EM3Dmap::calculate(){

  // calculate CV 
  calculate_overlap();

  // calculate "restraint"
  double ene = 0.0;
  // count number of non-zero overlaps
  double ndata_zero = 0.0;
  for(unsigned i=0;i<ovmd_.size();++i){
    ene_der_[i] = 0.0;
    if(ovmd_[i] > 0.0){
     // individual term
     ene_der_[i] = std::log(ovmd_[i]/ovdd_[i]);
     // increment energy
     ene += ene_der_[i] * ene_der_[i];
     // increment counter
     ndata_zero += 1.0;
    }
  };
  
  // constant factor
  double fact = kbt_ * 0.5 * ndata_zero;

  // clear temporary vector
  for(unsigned i=0; i<atom_der_.size(); ++i) atom_der_[i] = Vector(0,0,0);

  // get derivatives of bias with respect to atoms
  for(unsigned i=rank_;i<nl_.size();i=i+size_) {
     // get indexes of data and model component
     unsigned id = nl_[i].first;
     unsigned im = nl_[i].second;
     // check for zero overlaps
     if(ovmd_[id] > 0.0 && ene > 0.0){
      double der = 2.0 * fact / ene * ene_der_[id] / ovmd_[id];
      // chain rule
      atom_der_[im] += der * ovmd_der_[i];
     }
  }
  // if parallel, communicate stuff
  if(!serial_) comm.Sum(&atom_der_[0][0], 3*atom_der_.size());
 
  // set derivative
  for(unsigned i=0;i<atom_der_.size();++i) setAtomsDerivatives(i, atom_der_[i]);

  // set value of the bias
  setValue(fact * std::log(ene));
}

}
}
