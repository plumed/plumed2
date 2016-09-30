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
 vector< VectorGeneric<9> > GMM_m_cov_;
 // data GMM - means, weights, and covariances
 vector<Vector>           GMM_d_m_;
 vector<double>           GMM_d_w_;
 vector< VectorGeneric<9> > GMM_d_cov_;
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
 vector< VectorGeneric<9> > inv_cov_md_;
 // neighbor list
 double   nl_cutoff_;
 unsigned nl_stride_;
 bool first_time_;
 vector < unsigned > nl_;
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
 double get_prefactor_inverse (const VectorGeneric<9> &GMM_cov_0, const VectorGeneric<9> &GMM_cov_1,
        double &GMM_w_0, double &GMM_w_1, 
        VectorGeneric<9> &sum, VectorGeneric<9> &inv_sum);
 // calculate self overlaps between data GMM components - ovdd_
 double get_self_overlap(unsigned id);
 // calculate overlap between two components
 double get_overlap(const Vector &m_m, const Vector &d_m, double &fact_md,
                    const VectorGeneric<9> &inv_cov_md, Vector &ov_der);
 double get_overlap(const Vector &m_m, const Vector &d_m, double &fact_md,
                    const VectorGeneric<9> &inv_cov_md);
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
  VectorGeneric<9> cov;

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
        cov[0]=s*s; cov[1]=0.0; cov[2]=0.0;
        cov[3]=0.0; cov[4]=s*s; cov[5]=0.0;
        cov[6]=0.0; cov[7]=0.0; cov[8]=s*s;
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
 VectorGeneric<9> cov;
 
 // open file
 IFile *ifile = new IFile();
 if(ifile->FileExist(GMM_file)){
    ifile->open(GMM_file);
    while(ifile->scanField("Id",idcomp)){
     ifile->scanField("Weight",w);
     ifile->scanField("Mean_0",m0);
     ifile->scanField("Mean_1",m1);
     ifile->scanField("Mean_2",m2);
     ifile->scanField("Cov_00",cov[0]);
     ifile->scanField("Cov_01",cov[1]);
     ifile->scanField("Cov_02",cov[2]);
     ifile->scanField("Cov_10",cov[3]);
     ifile->scanField("Cov_11",cov[4]);
     ifile->scanField("Cov_12",cov[5]);
     ifile->scanField("Cov_20",cov[6]);
     ifile->scanField("Cov_21",cov[7]);
     ifile->scanField("Cov_22",cov[8]);
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
(const VectorGeneric<9> &GMM_cov_0, const VectorGeneric<9> &GMM_cov_1,
 double &GMM_w_0, double &GMM_w_1, 
 VectorGeneric<9> &sum_0_1, VectorGeneric<9> &inv_sum)
{
 // we need the sum of the covariance matrices
 for(unsigned k=0; k<9; ++k) sum_0_1[k] = GMM_cov_0[k] + GMM_cov_1[k]; 
    
 // and to calculate its determinant
 double det = sum_0_1[0]*(sum_0_1[4]*sum_0_1[8]-sum_0_1[5]*sum_0_1[7]);
       det -= sum_0_1[1]*(sum_0_1[3]*sum_0_1[8]-sum_0_1[5]*sum_0_1[6]);
       det += sum_0_1[2]*(sum_0_1[3]*sum_0_1[7]-sum_0_1[4]*sum_0_1[6]);
 // the prefactor is 
 double pre_fact =  cfact_ / sqrt(det) * GMM_w_0 * GMM_w_1;
 // and its inverse
 inv_sum[0] = (sum_0_1[4]*sum_0_1[8] - sum_0_1[5]*sum_0_1[7])/det;
 inv_sum[1] = (sum_0_1[2]*sum_0_1[7] - sum_0_1[1]*sum_0_1[8])/det;
 inv_sum[2] = (sum_0_1[1]*sum_0_1[5] - sum_0_1[2]*sum_0_1[4])/det;
 inv_sum[3] = (sum_0_1[5]*sum_0_1[6] - sum_0_1[3]*sum_0_1[8])/det;
 inv_sum[4] = (sum_0_1[0]*sum_0_1[8] - sum_0_1[2]*sum_0_1[6])/det;
 inv_sum[5] = (sum_0_1[2]*sum_0_1[3] - sum_0_1[0]*sum_0_1[5])/det;
 inv_sum[6] = (sum_0_1[3]*sum_0_1[7] - sum_0_1[4]*sum_0_1[6])/det;
 inv_sum[7] = (sum_0_1[1]*sum_0_1[6] - sum_0_1[0]*sum_0_1[7])/det;
 inv_sum[8] = (sum_0_1[0]*sum_0_1[4] - sum_0_1[1]*sum_0_1[3])/det;
 // return pre-factor
 return pre_fact;
}


double EM3Dmap::get_self_overlap(unsigned id)
{
 double ov = 0.0;
 VectorGeneric<9> sum;
 VectorGeneric<9> inv_sum;
 // start loop
 for(unsigned i=0; i<GMM_d_w_.size(); ++i){
   // call auxiliary method
   double pre_fact = get_prefactor_inverse(GMM_d_cov_[id], GMM_d_cov_[i], 
                                             GMM_d_w_[id],   GMM_d_w_[i], sum, inv_sum); 
   // calculate overlap
   double ov_tmp = get_overlap(GMM_d_m_[id], GMM_d_m_[i], pre_fact, inv_sum);
   // add to overlap
   ov += ov_tmp;
 }
 return ov;
}

double EM3Dmap::get_overlap(const Vector &m_m, const Vector &d_m, double &fact_md,
                            const VectorGeneric<9> &inv_cov_md, Vector &ov_der)
{
  // calculate vector difference m_m-d_m
  double md_x = m_m[0] - d_m[0];
  double md_y = m_m[1] - d_m[1];
  double md_z = m_m[2] - d_m[2];
  // calculate product of transpose of md and inv_cov_md
  double p_x = md_x*inv_cov_md[0]+md_y*inv_cov_md[3]+md_z*inv_cov_md[6];
  double p_y = md_x*inv_cov_md[1]+md_y*inv_cov_md[4]+md_z*inv_cov_md[7];
  double p_z = md_x*inv_cov_md[2]+md_y*inv_cov_md[5]+md_z*inv_cov_md[8];
  // calculate product of prod and md
  double ov = md_x*p_x+md_y*p_y+md_z*p_z; 
  // final calculation
  ov = fact_md * exp(-0.5*ov);
  // derivatives
  double x = md_x*inv_cov_md[0] + md_y*inv_cov_md[1] + md_z*inv_cov_md[2];
  double y = md_x*inv_cov_md[1] + md_y*inv_cov_md[4] + md_z*inv_cov_md[5];
  double z = md_x*inv_cov_md[2] + md_y*inv_cov_md[5] + md_z*inv_cov_md[8];
  ov_der = ov * Vector(x, y, z); 
  return ov;
}

double EM3Dmap::get_overlap(const Vector &m_m, const Vector &d_m, double &fact_md, 
                            const VectorGeneric<9> &inv_cov_md)
                        
{
  // calculate vector difference m_m-d_m
  double md_x = m_m[0] - d_m[0];
  double md_y = m_m[1] - d_m[1];
  double md_z = m_m[2] - d_m[2];
  // calculate product of transpose of md and inv_cov_md
  double p_x = md_x*inv_cov_md[0]+md_y*inv_cov_md[3]+md_z*inv_cov_md[6];
  double p_y = md_x*inv_cov_md[1]+md_y*inv_cov_md[4]+md_z*inv_cov_md[7];
  double p_z = md_x*inv_cov_md[2]+md_y*inv_cov_md[5]+md_z*inv_cov_md[8];
  // calculate product of prod and md
  double ov = md_x*p_x+md_y*p_y+md_z*p_z; 
  // final calculation
  ov = fact_md * exp(-0.5*ov);
  return ov;
}

void EM3Dmap::update_neighbor_list()
{
  // temporary stuff
  VectorGeneric<9> sum;
  VectorGeneric<9> inv_sum;
  // clear old neighbor list and auxiliary vectors
  nl_.clear(); fact_md_.clear(); inv_cov_md_.clear();
  // local stuff
  vector< VectorGeneric<9> > inv_cov_md_l;
  vector < unsigned > nl_l;
  vector < double > fact_md_l;
  // cycle on all overlaps
  unsigned nover = GMM_d_w_.size() * GMM_m_w_.size();
  for(unsigned k=rank_; k<nover; k=k+size_){
      // get indexes
      unsigned i = k / GMM_m_w_.size();
      unsigned j = k % GMM_m_w_.size();
      // call auxiliary method 
      double pre_fact = get_prefactor_inverse(GMM_d_cov_[i], GMM_m_cov_[j], 
                                                GMM_d_w_[i],   GMM_m_w_[j], 
                                                sum, inv_sum);	
      // calculate overlap
      double ov = get_overlap(GMM_d_m_[i], getPosition(j), pre_fact, inv_sum);
      // fill the neighbor list and auxiliary vectors
      if(ov >= nl_cutoff_ * ovdd_[i]){
        nl_l.push_back(k);
        fact_md_l.push_back(pre_fact);
        inv_cov_md_l.push_back(inv_sum);
      }
  }
  // find total dimension of neighborlist
  vector <int> recvcounts(size_, 0);
  recvcounts[rank_] = nl_l.size();
  comm.Sum(&recvcounts[0], size_);
  int tot_size = accumulate(recvcounts.begin(), recvcounts.end(), 0);
  // resize stuff
  nl_.resize(tot_size); fact_md_.resize(tot_size); inv_cov_md_.resize(tot_size);
  // calculate vector of displacement
  vector<int> disp(size_);
  disp[0] = 0;
  int rank_size = 0;
  for(unsigned i=0; i<size_-1; ++i){
    rank_size += recvcounts[i];
    disp[i+1] = rank_size;
  }
  // Allgather
  comm.Allgatherv(&nl_l[0],      recvcounts[rank_], &nl_[0],      &recvcounts[0], &disp[0]);
  comm.Allgatherv(&fact_md_l[0], recvcounts[rank_], &fact_md_[0], &recvcounts[0], &disp[0]);
  // adapt for Vector generic
  for(unsigned i=0; i<size_; ++i){
   recvcounts[i] *= 9;
   disp[i] *= 9;
  }
  comm.Allgatherv(&inv_cov_md_l[0][0], recvcounts[rank_], &inv_cov_md_[0][0], &recvcounts[0], &disp[0]);
  // now resize derivatives
  ovmd_der_.resize(tot_size);
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
      unsigned id = nl_[i] / GMM_m_w_.size();
      unsigned im = nl_[i] % GMM_m_w_.size();
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
     unsigned id = nl_[i] / GMM_m_w_.size();
     unsigned im = nl_[i] % GMM_m_w_.size();
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
