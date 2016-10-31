/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#include "core/ActionSet.h"
#include "tools/File.h"

#include <string>
#include <cmath>
#include <map>
#include <numeric>
#include <ctime>

using namespace std;


namespace PLMD{
namespace colvar{

//+PLUMEDOC COLVAR QXL
/*
Put Documentation here 

*/
//+ENDPLUMEDOC

class QCrossLink : public Colvar
{
  // sigma parameter
  double sigma_;
  double sigma_min_;
  double sigma_max_;
  double Dsigma_;
  // beta parameter
  double beta_;
  double beta_min_;
  double beta_max_;
  double Dbeta_;
  // Rslope parameter
  double Rslope_;
  double Rslope_min_;
  double Rslope_max_;
  double DRslope_;
  // W parameter
  double W_;
  double W_min_;
  double W_max_;
  double DW_;
  // alpha parameter
  double alpha_;
  // temperature in kbt
  double kbt_;
  // Monte Carlo stuff
  int MCsteps_;
  int MCstride_;
  long int MCfirst_;
  unsigned int MCaccsig_;
  unsigned int MCaccbeta_;
  unsigned int MCaccRslope_;
  unsigned int MCaccW_;

  // parallel stuff
  unsigned rank_;
  unsigned rank_f_;
  bool is_state_A_;
  bool is_state_B_;
  // data map
  vector< pair <unsigned,unsigned> > pair_list_;
  vector< double > ratio_list_;
  // auxiliary stuff
  vector<Vector> dexp_;
  vector<Vector> dexp_f_;
  vector<Vector> atom_der_;
  
  void get_data(string data_file, vector<AtomNumber> atoms);
  double get_rho(double dist, double Rslope, double alpha);
  void doMonteCarlo();
  double getEnergy(double sigma, double w, double beta, double Rslope, double alpha);
  double proposeMove(double x, double xmin, double xmax, double dxmax);
  bool doAccept(double oldE, double newE);

public:
  static void registerKeywords( Keywords& keys );
  explicit QCrossLink(const ActionOptions&);
// active methods:
  virtual void calculate();
};


PLUMED_REGISTER_ACTION(QCrossLink,"QXL")

void QCrossLink::registerKeywords(Keywords& keys){
  Colvar::registerKeywords( keys );
  keys.add("atoms","ATOMS","atoms for which we calculate the density map");
  keys.add("compulsory","DATA_FILE","file with experimental data");
  keys.add("compulsory","TEMP","temperature");
  keys.add("compulsory","ALPHA","Value of the alpha parameter");
  keys.add("compulsory","RANK_A","MPI rank of state A");
  keys.add("compulsory","RANK_B","MPI rank of state B");
  componentsAreNotOptional(keys); 
  keys.addOutputComponent("sigma",    "default","uncertainty parameter");
  keys.addOutputComponent("accsig",   "default","MC acceptance sigma");
  keys.addOutputComponent("beta",     "default","beta parameter");
  keys.addOutputComponent("accbeta",  "default","MC acceptance beta");
  keys.addOutputComponent("Rslope",   "default","Rslope parameter");
  keys.addOutputComponent("accRslope","default","MC acceptance Rslope");
  keys.addOutputComponent("W",        "default","W parameter");
  keys.addOutputComponent("accW",     "default","MC acceptance W");
  keys.addOutputComponent("score",    "default","Bayesian score");
}

QCrossLink::QCrossLink(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
alpha_(1.0),
MCsteps_(1), MCstride_(1),MCfirst_(-1),
MCaccsig_(0), MCaccbeta_(0), MCaccRslope_(0), MCaccW_(0),
is_state_A_(false), is_state_B_(false)
{
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  // data file
  string data_file;
  parse("DATA_FILE",data_file);
  // sigma parameter
  parse("SIGMA0",   sigma_);
  parse("SIGMA_MIN",sigma_min_);
  parse("SIGMA_MAX",sigma_max_);
  parse("DSIGMA",   Dsigma_);
  // beta parameter
  parse("BETA0",    beta_);
  parse("BETA_MIN", beta_min_);
  parse("BETA_MAX", beta_max_);
  parse("DBETA",    Dbeta_);
  // Rslope parameter
  parse("RSLOPE0",   Rslope_);
  parse("RSLOPE_MIN",Rslope_min_);
  parse("RSLOPE_MAX",Rslope_max_);
  parse("DRSLOPE",   DRslope_);
  // W parameter
  parse("W0",   W_);
  parse("W_MIN",W_min_);
  parse("W_MAX",W_max_);
  parse("DW",   DW_);
  // alpha parameter (constant)
  parse("ALPHA", alpha_);
  // rank of state A
  unsigned rank_A;
  parse("RANK_A", rank_A);
  // rank of state B
  unsigned rank_B;
  parse("RANK_B", rank_B);
  // temperature
  double temp=0.0;
  parse("TEMP",temp);
  // convert temp to kbt
  if(temp>0.0) kbt_=plumed.getAtoms().getKBoltzmann()*temp;
  else kbt_=plumed.getAtoms().getKbT();

  checkRead();
  
  // get number of replicas
  if(comm.Get_rank()==0) {
    rank_ = multi_sim_comm.Get_rank();
  } else {
    rank_ = 0;
  }
  comm.Sum(&rank_,1);
  
  // find out if state A or B
  if(rank_==rank_A){
    is_state_A_ = true;
    rank_f_ = rank_B;
  }
  if(rank_==rank_B){
    is_state_B_ = true;
    rank_f_ = rank_A;
  }
  if(is_state_A_==is_state_B_) error("Each replica must be in either state A or B\n");
  
  log.printf("  atoms involved : ");
  for(unsigned i=0;i<atoms.size();++i) log.printf("%d ",atoms[i].serial());
  log.printf("\n");
  log.printf("  data file %s", data_file.c_str());
  if(is_state_A_) log.printf("  this replica is modelling state A\n");
  if(is_state_B_) log.printf("  this replica is modelling state B\n");
  log.printf("  initial value of uncertainty %f\n",sigma_);
  log.printf("  minimum value of uncertainty %f\n",sigma_min_);
  log.printf("  maximum value of uncertainty %f\n",sigma_max_);
  log.printf("  maximum MC move of the uncertainty parameter %f\n",Dsigma_);
  log.printf("  initial value of beta %f\n",beta_);
  log.printf("  minimum value of beta %f\n",beta_min_);
  log.printf("  maximum value of beta %f\n",beta_max_);
  log.printf("  maximum MC move of the beta parameter %f\n",Dbeta_);
  log.printf("  initial value of Rslope %f\n",Rslope_);
  log.printf("  minimum value of Rslope %f\n",Rslope_min_);
  log.printf("  maximum value of Rslope %f\n",Rslope_max_);
  log.printf("  maximum MC move of the Rslope parameter %f\n",DRslope_);
  log.printf("  initial value of W %f\n",W_);
  log.printf("  minimum value of W %f\n",W_min_);
  log.printf("  maximum value of W %f\n",W_max_);
  log.printf("  maximum MC move of the W parameter %f\n",DW_);
  log.printf("  alpha parameter %f\n", alpha_);
  log.printf("  temperature of the system in energy unit %f\n",kbt_);
  log.printf("  number of MC steps %d\n",MCsteps_);
  log.printf("  do MC every %d steps\n", MCstride_);

  addComponent("sigma");    componentIsNotPeriodic("sigma");
  addComponent("accsig");   componentIsNotPeriodic("accsig");
  addComponent("beta");     componentIsNotPeriodic("beta");
  addComponent("accbeta");  componentIsNotPeriodic("accbeta");
  addComponent("Rslope");   componentIsNotPeriodic("Rslope");
  addComponent("accRslope");componentIsNotPeriodic("accRslope");
  addComponent("W");        componentIsNotPeriodic("W");
  addComponent("accW");     componentIsNotPeriodic("accW");
  addComponentWithDerivatives("score"); componentIsNotPeriodic("score");
  

  // initialize random seed
  srand (time(NULL)+rank_);
  
  // read data file
  get_data(data_file, atoms);
  log.printf("  number of data points : %u\n", static_cast<unsigned>(ratio_list_.size()));
  
  // resize stuff
  dexp_.resize(ratio_list_.size());
  dexp_f_.resize(ratio_list_.size());
  atom_der_.resize(atoms.size());
  
  // request the atoms
  requestAtoms(atoms);
}

// read data file
void QCrossLink::get_data(string data_file, vector<AtomNumber> atoms)
{
 int id;
 int at0, at1;
 int id0, id1;
 double ratio;
 
 // open file
 IFile *ifile = new IFile();
 if(ifile->FileExist(data_file)){
    ifile->open(data_file);
    while(ifile->scanField("Id",id)){
     ifile->scanField("Atom0",at0);
     ifile->scanField("Atom1",at1);
     ifile->scanField("Ratio",ratio);
     // find indexes of at0 and at1 inside atoms
     for(unsigned i=0; i<atoms.size(); ++i){
       if(atoms[i].serial() == at0) id0 = i;
       if(atoms[i].serial() == at1) id1 = i;     
     }
     // add data to map
     pair_list_.push_back(make_pair(id0,id1));
     ratio_list_.push_back(ratio);
     // new line
     ifile->scanField();
    }
    ifile->close();
 } else {
    error("Cannot find data_FILE "+data_file+"\n"); 
 }
 delete ifile;
}

// get rho, geometric factor of XL rate coefficient
double QCrossLink::get_rho(double dist, double Rslope, double alpha)
{
 double rho = 1.0 - 1.0 / (1.0+exp(-alpha*(dist-Rslope)));
 return rho;
}

// used to update all the Bayesian parameters
double QCrossLink::getEnergy(double sigma, double w, double beta, double Rslope, double alpha)
{
  // calculate energy
  double rho_A, rho_B;
  double ene = 0.0;
  // cycle on arguments
  for(unsigned i=0;i<dexp_.size();++i){
    // calculate rho
    if(is_state_A_){ 
      rho_A = get_rho(dexp_[i].modulo(),   Rslope, alpha);
      rho_B = get_rho(dexp_f_[i].modulo(), Rslope, alpha);
    } else {
      rho_A = get_rho(dexp_f_[i].modulo(), Rslope, alpha);
      rho_B = get_rho(dexp_[i].modulo(),   Rslope, alpha);
    }
    // calculate forward model
    double fmod = w * (1.0-exp(-beta*rho_A)) / (1.0-exp(-beta*rho_B));
    // calculate argument
    double tmp0 = std::log(ratio_list_[i] / fmod);
    double tmp1 = tmp0 * tmp0 + 2.0 * sigma * sigma;
    // increment energy
    ene += kbt_ * std::log(tmp1);    
  }
  // add Jeffrey's priors + constant terms in sigma;
  ene += kbt_ * std::log(sigma) - kbt_ * static_cast<double>(dexp_.size()) * std::log(sigma);
  return ene;
}

double QCrossLink::proposeMove(double x, double xmin, double xmax, double dxmax)
{
 double r = static_cast<double>(rand()) / RAND_MAX;
 double dx = -dxmax + r * 2.0 * dxmax;
 double x_new = x + dx;
 // check boundaries
 if(x_new > xmax){x_new = 2.0 * xmax - x_new;}
 if(x_new < xmin){x_new = 2.0 * xmin - x_new;}
 return x_new;
}

bool QCrossLink::doAccept(double oldE, double newE){
  bool accept = false;
  // calculate delta energy 
  double delta = ( newE - oldE ) / kbt_;
  // if delta is negative always accept move
  if( delta < 0.0 ){ 
   accept = true;
  }else{
   // otherwise extract random number   
   double s = static_cast<double>(rand()) / RAND_MAX;
   if( s < exp(-delta) ) { accept = true; }
  }
  return accept;
}
    
void QCrossLink::doMonteCarlo()
{
 double oldE, newE;
 bool accept;
 // get step
 long int step = getStep();
 // this is needed when restarting simulations
 if(MCfirst_==-1) MCfirst_=step;
 // calculate acceptance
 double MCtrials = std::floor(static_cast<double>(step-MCfirst_) / static_cast<double>(MCstride_))+1.0;
 // store old energy
 oldE = getEnergy(sigma_, W_, beta_, Rslope_, alpha_);
 // cycle on MC steps 
 for(unsigned i=0;i<MCsteps_;++i){
  // propose move in sigma
  double new_sigma = proposeMove(sigma_,sigma_min_,sigma_max_,Dsigma_);
  // calculate new energy
  newE = getEnergy(new_sigma, W_, beta_, Rslope_, alpha_);
  // accept or reject
  accept = doAccept(oldE, newE);
  if(accept){
   sigma_ = new_sigma;
   oldE = newE;
   MCaccsig_++;
  }
  // propose move in beta
  double new_beta = proposeMove(beta_,beta_min_,beta_max_,Dbeta_);
  // calculate new energy
  newE = getEnergy(sigma_, W_, new_beta, Rslope_, alpha_);
  // accept or reject
  accept = doAccept(oldE, newE);
  if(accept){
   beta_ = new_beta;
   oldE = newE;
   MCaccbeta_++;
  }
  // propose move in W
  double new_W = proposeMove(W_,W_min_,W_max_,DW_);
  // calculate new energy
  newE = getEnergy(sigma_, new_W, beta_, Rslope_, alpha_);
  // accept or reject
  accept = doAccept(oldE, newE);
  if(accept){
   W_ = new_W;
   oldE = newE;
   MCaccW_++;
  }
  // propose move in Rslope
  double new_Rslope = proposeMove(Rslope_,Rslope_min_,Rslope_max_,DRslope_);
  // calculate new energy
  newE = getEnergy(sigma_, W_, beta_, new_Rslope, alpha_);
  // accept or reject
  accept = doAccept(oldE, newE);
  if(accept){
   Rslope_ = new_Rslope;
   oldE = newE;
   MCaccRslope_++;
  }
 }
 // set values
 // sigma
 getPntrToComponent("sigma")->set(sigma_);
 // sigma acceptance
 double accsig = static_cast<double>(MCaccsig_) / static_cast<double>(MCsteps_) / MCtrials;
 getPntrToComponent("accsig")->set(accsig);
 // beta
 getPntrToComponent("beta")->set(beta_);
 // beta acceptance
 double accbeta = static_cast<double>(MCaccbeta_) / static_cast<double>(MCsteps_) / MCtrials;
 getPntrToComponent("accbeta")->set(accbeta);
 // Rslope
 getPntrToComponent("Rslope")->set(Rslope_);
 // Rslope acceptance
 double accRslope = static_cast<double>(MCaccRslope_) / static_cast<double>(MCsteps_) / MCtrials;
 getPntrToComponent("accRslope")->set(accRslope);
 // W
 getPntrToComponent("W")->set(W_);
 // W acceptance
 double accW = static_cast<double>(MCaccW_) / static_cast<double>(MCsteps_) / MCtrials;
 getPntrToComponent("accW")->set(accW);
}

void QCrossLink::calculate()
{
  // calculate distances
  if(comm.Get_rank()==0){
   for(unsigned i=0; i<dexp_.size(); ++i){
     unsigned id0 = pair_list_[i].first;
     unsigned id1 = pair_list_[i].second;
     dexp_[i] = pbcDistance(getPosition(id0),getPosition(id1)); 
   }
   // send and receive dexp_ with partner
   multi_sim_comm.Isend(dexp_, rank_f_,123);
   multi_sim_comm.Recv(dexp_f_,rank_f_,123);
  } else {
   for(unsigned i=0; i<dexp_.size(); ++i){
     dexp_[i] = Vector(0,0,0);
     dexp_f_[i] = Vector(0,0,0);
   }
  }
  comm.Sum(&dexp_[0],   3*dexp_.size());
  comm.Sum(&dexp_f_[0], 3*dexp_f_.size());
  
  // get time step 
  long int step = getStep();
  // do MC stuff at the right time step
  if(step%MCstride_==0&&!getExchangeStep()) doMonteCarlo();

  // send values of parameters to state B
  if(comm.Get_rank()==0){
   if(is_state_A_){
    multi_sim_comm.Isend(sigma_, rank_f_,400);
    multi_sim_comm.Isend(beta_, rank_f_,401);
    multi_sim_comm.Isend(Rslope_, rank_f_,402);
    multi_sim_comm.Isend(W_, rank_f_,403);
   } else {
    multi_sim_comm.Recv(sigma_,rank_f_,400);
    multi_sim_comm.Recv(beta_,rank_f_,401);
    multi_sim_comm.Recv(Rslope_,rank_f_,402);
    multi_sim_comm.Recv(W_,rank_f_,403);
   }
  } else {
    sigma_  = 0.0;
    beta_   = 0.0;
    Rslope_ = 0.0;
    W_      = 0.0;
  }
  comm.Sum(&sigma_, 1);
  comm.Sum(&beta_, 1);
  comm.Sum(&Rslope_, 1);
  comm.Sum(&W_, 1);
  
  // clear derivatives
  for(unsigned i=0; i<atom_der_.size(); ++i) atom_der_[i] = Vector(0,0,0);
  
  // calculate energy
  double rho_A, rho_B;
  double ene = 0.0;
  // cycle on arguments
  for(unsigned i=0;i<dexp_.size();++i){
    // retrieve atom indexes
    unsigned id0 = pair_list_[i].first;
    unsigned id1 = pair_list_[i].second;
    // calculate rho
    if(is_state_A_){ 
      rho_A = get_rho(dexp_[i].modulo(),   Rslope_, alpha_);
      rho_B = get_rho(dexp_f_[i].modulo(), Rslope_, alpha_);
    } else {
      rho_A = get_rho(dexp_f_[i].modulo(), Rslope_, alpha_);
      rho_B = get_rho(dexp_[i].modulo(),   Rslope_, alpha_);
    }
    // calculate forward model
    double fmod_A = (1.0-exp(-beta_*rho_A));
    double fmod_B = (1.0-exp(-beta_*rho_B));
    double fmod  = W_ * fmod_A / fmod_B;
    // calculate argument
    double tmp0 = std::log(ratio_list_[i] / fmod);
    double tmp1 = tmp0 * tmp0 + 2.0 * sigma_ * sigma_;
    // increment energy
    ene += kbt_ * std::log(tmp1);
    // calculate derivative of score with respect to fmod
    double score_der_fmod = - kbt_ / tmp1 * 2.0 * tmp0 / fmod ;
    // derivative of fmod with respect to rho
    double fmod_der_rho; 
    if(is_state_A_){
      // calculate derivative of fmod with respect to rho_A
      fmod_der_rho = - W_ * (-fmod_A + 1.0 ) / fmod_B * beta_;
    } else {
      // calculate derivative of fmod with respect to rho_B
      fmod_der_rho =   W_ * fmod_A / fmod_B / fmod_B * (-fmod_B + 1.0) * beta_;
    }
    // calculate derivative of rho with respect to distance
    // rho = 1.0 - 1.0 / (1.0+exp(-alpha*(dist-Rslope)));
    double tmp_rho = exp(-alpha_*(dexp_[i].modulo()-Rslope_));
    double rho_der_dist = - 1.0 / ( 1.0 + tmp_rho ) / ( 1.0 + tmp_rho ) * tmp_rho * alpha_;
    // add to derivatives with respect to atom positions
    atom_der_[id0] -= score_der_fmod * fmod_der_rho * rho_der_dist * dexp_[i] / dexp_[i].modulo();
    atom_der_[id1] += score_der_fmod * fmod_der_rho * rho_der_dist * dexp_[i] / dexp_[i].modulo();
  }
  // add Jeffrey's priors + constant terms in sigma;
  ene += kbt_ * std::log(sigma_) - kbt_ * static_cast<double>(dexp_.size()) * std::log(sigma_);
  // set value of the bias
  getPntrToComponent("score")->set(ene);
  // set derivatives
  for(unsigned i=0; i<atom_der_.size(); ++i) setAtomsDerivatives(getPntrToComponent("score"), i, atom_der_[i]);
}


}
}


