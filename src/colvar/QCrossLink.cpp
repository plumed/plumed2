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
  // pbc
  bool nopbc_;
  // Monte Carlo stuff
  int MCsteps_;
  int MCstride_;
  long int MCfirst_;
  unsigned int MCaccsig_;
  unsigned int MCaccbeta_;
  unsigned int MCaccRslope_;
  unsigned int MCaccW_;
  // debug stuff
  bool debug_;
  string dfile_name_;
  OFile dfile_;

  // parallel stuff
  unsigned rank_;
  unsigned rank_f_;
  bool is_state_A_;
  bool is_state_B_;
  unsigned rank_local_;
  unsigned size_local_;
  bool serial_;
  
  // experimental data vectors
  vector< pair <unsigned,unsigned> > pair_list_;
  vector< double > ratio_list_;
  // auxiliary stuff
  vector<Vector> dexp_;
  vector<Vector> dexp_f_;
  vector<Vector> atom_der_;
  
  void get_data(string data_file, vector<AtomNumber> atoms);
  double get_rho(double dist, double Rslope, double alpha);
  void print_MC(long int MDstep, long int MCstep, double sigma, double beta, double Rslope, double W, double ene);
  void doMonteCarlo(long int step);
  double getEnergy(double sigma, double w, double beta, double Rslope, double alpha);
  double proposeMove(double x, double xmin, double xmax, double dxmax);
  bool doAccept(double oldE, double newE);

public:
  static void registerKeywords( Keywords& keys );
  explicit QCrossLink(const ActionOptions&);
  ~QCrossLink();
// active methods:
  virtual void calculate();
};


PLUMED_REGISTER_ACTION(QCrossLink,"QXL")

void QCrossLink::registerKeywords(Keywords& keys){
  Colvar::registerKeywords( keys );
  keys.add("atoms","ATOMS","atoms for which we have cross-linking data");
  keys.add("compulsory","DATA_FILE","file with experimental data");
  keys.add("compulsory","SIGMA0","initial value of the uncertainty parameter");
  keys.add("compulsory","SIGMA_MIN","minimum value of the uncertainty parameter");
  keys.add("compulsory","SIGMA_MAX","maximum value of the uncertainty parameter");
  keys.add("compulsory","DSIGMA","maximum MC move of the uncertainty parameter");
  keys.add("compulsory","BETA0","initial value of the beta parameter");
  keys.add("compulsory","BETA_MIN","minimum value of the beta parameter");
  keys.add("compulsory","BETA_MAX","maximum value of the beta parameter");
  keys.add("compulsory","DBETA","maximum MC move of the beta parameter");
  keys.add("compulsory","RSLOPE0","initial value of the Rslope parameter");
  keys.add("compulsory","RSLOPE_MIN","minimum value of the Rslope parameter");
  keys.add("compulsory","RSLOPE_MAX","maximum value of the Rslope parameter");
  keys.add("compulsory","DRSLOPE","maximum MC move of the Rslope parameter");
  keys.add("compulsory","W0","initial value of the W parameter");
  keys.add("compulsory","W_MIN","minimum value of the W parameter");
  keys.add("compulsory","W_MAX","maximum value of the W parameter");
  keys.add("compulsory","DW","maximum MC move of the W parameter");
  keys.add("compulsory","RANK_A","MPI rank of state A");
  keys.add("compulsory","RANK_B","MPI rank of state B");
  keys.add("compulsory","TEMP","temperature");
  keys.addFlag("SERIAL",false,"perform the calculation in serial - for debug purpose");
  keys.add("optional","ALPHA","Value of the alpha parameter");
  keys.add("optional","DEBUG_FILE","file with MC optimization details");
  keys.add("optional","MC_STEPS","number of MC steps");
  keys.add("optional","MC_STRIDE","MC stride");
  
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

QCrossLink::~QCrossLink()
{
  if(debug_) dfile_.close();
}

QCrossLink::QCrossLink(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao), alpha_(1.0), nopbc_(false),
MCsteps_(1), MCstride_(1),MCfirst_(-1),
MCaccsig_(0), MCaccbeta_(0), MCaccRslope_(0), MCaccW_(0), debug_(false),
is_state_A_(false), is_state_B_(false), serial_(false)
{
  // list of atoms
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
  // disable pbc
  parseFlag("NOPBC",nopbc_);
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
  // serial or parallel
  parseFlag("SERIAL",serial_);
  if(serial_){
    size_local_=1; rank_local_=0;
  } else {
    size_local_=comm.Get_size(); rank_local_=comm.Get_rank();
  }
  parse("MC_STEPS",MCsteps_);
  parse("MC_STRIDE",MCstride_);
  parse("DEBUG_FILE",dfile_name_);
  if(dfile_name_.length()>0) debug_ = true;

  checkRead();
  
  // get rank in replica parallelism
  if(comm.Get_rank()==0) {
    rank_ = multi_sim_comm.Get_rank();
  } else {
    rank_ = 0;
  }
  if(comm.Get_size()>1) comm.Sum(&rank_,1);
  
  // find out if this rank will model A or B
  if(rank_==rank_A){
    is_state_A_ = true;
    rank_f_ = rank_B;
  }
  if(rank_==rank_B){
    is_state_B_ = true;
    rank_f_ = rank_A;
  }
  // check that I am either state A or state B
  if(is_state_A_==is_state_B_) error("Each replica must be in either state A or B\n");
  
  log.printf("  atoms involved : ");
  for(unsigned i=0;i<atoms.size();++i) log.printf("%d ",atoms[i].serial());
  log.printf("\n");
  log.printf("  data file %s\n", data_file.c_str());
  if(is_state_A_) log.printf("  this replica is modelling state A\n");
  if(is_state_B_) log.printf("  this replica is modelling state B\n");
  if(serial_)     log.printf("  serial calculation\n");
  if(nopbc_)      log.printf("  distances are calculated without pbc\n");
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
  if(debug_) log.printf("  printing MC optimization details to %s\n", dfile_name_.c_str());

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
  
  // resize auxiliary vectors
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
     ifile->scanField("Atom1",at0);
     ifile->scanField("Atom2",at1);
     ifile->scanField("Ratio",ratio);
     // find indexes of at0 and at1 inside atoms
     for(unsigned i=0; i<atoms.size(); ++i){
       if(atoms[i].serial() == at0) id0 = i;
       if(atoms[i].serial() == at1) id1 = i;     
     }
     // add data to lists
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
  double rho_A, rho_B;
  // calculate energy
  double ene = 0.0;
  // cycle on arguments
  for(unsigned i=rank_local_;i<dexp_.size();i=i+size_local_){
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
  // if parallel, sum all the energies
  if(!serial_ && comm.Get_size()>1) comm.Sum(&ene,1);
  // add Jeffrey's prior + constant terms in sigma;
  ene += kbt_ * ( 1.0 - static_cast<double>(dexp_.size()) ) * std::log(sigma);
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

bool QCrossLink::doAccept(double oldE, double newE)
{
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

void QCrossLink::print_MC(long int MDstep, long int MCstep, double sigma,
  double beta, double Rslope, double W, double ene)
{
 // open the file if we are at the first MD and MC steps
 if(MDstep==0 && MCstep==0){
   dfile_.link(*this);
   dfile_.open(dfile_name_);
   dfile_.setHeavyFlush();
   dfile_.fmtField("%12.6f"); 
 }
 // write fields
 dfile_.printField("MD step", static_cast<int>(MDstep));
 dfile_.printField("MC step", static_cast<int>(MCstep));
 dfile_.printField("Sigma",   sigma);
 dfile_.printField("Beta",    beta);
 dfile_.printField("Rslope",  Rslope);
 dfile_.printField("W",       W);
 dfile_.printField("Energy",  ene);
 dfile_.printField();
}


void QCrossLink::doMonteCarlo(long int step)
{
 double oldE, newE;
 bool accept;
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
  // if debug, print out info about MC sampling
  if(debug_ && rank_==0 && comm.Get_rank()==0)
    print_MC(getStep(), i, sigma_, beta_, Rslope_, W_, newE);
 }
 // send values of parameters to state B, via buffer
 vector<double>       buff_d(4, 0.0), buff_d_r(4, 0.0);
 vector<unsigned int> buff_u(4, 0),   buff_u_r(4, 0);
  if(comm.Get_rank()==0){
   // prepare buffer double
   buff_d[0]=sigma_;    buff_d[1]=beta_;      buff_d[2]=Rslope_;      buff_d[3]=W_;
   // send and receive
   Communicator::Request req_d=multi_sim_comm.Isend(buff_d, rank_f_, 400);
   multi_sim_comm.Recv(buff_d_r, rank_f_, 400);
   req_d.wait();
   // prepare buffer unsigned 
   buff_u[0]=MCaccsig_; buff_u[1]=MCaccbeta_; buff_u[2]=MCaccRslope_; buff_u[3]=MCaccW_;
   // send and receive
   Communicator::Request req_u=multi_sim_comm.Isend(buff_u, rank_f_, 401);
   multi_sim_comm.Recv(buff_u_r, rank_f_, 401);
   req_u.wait();   
   // if state B write received buffer into variables
   if(is_state_B_){
    sigma_=buff_d_r[0];    beta_=buff_d_r[1];      Rslope_=buff_d_r[2];      W_=buff_d_r[3];
    MCaccsig_=buff_u_r[0]; MCaccbeta_=buff_u_r[1]; MCaccRslope_=buff_u_r[2]; MCaccW_=buff_u_r[3];
   }
  } else {
    sigma_  = 0.0;
    beta_   = 0.0;
    Rslope_ = 0.0;
    W_      = 0.0;
  }
  // wait for things to be done
  multi_sim_comm.Barrier();
  // local communication
  if(comm.Get_size()>1){
   comm.Sum(&sigma_, 1);
   comm.Sum(&beta_, 1);
   comm.Sum(&Rslope_, 1);
   comm.Sum(&W_, 1);
  }
 // set values of Bayesian parameters
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
     if(!nopbc_) dexp_[i] = pbcDistance(getPosition(id0),getPosition(id1));
     else        dexp_[i] = delta(getPosition(id0),getPosition(id1));
   }
   // send and receive dexp_ with partner
   Communicator::Request req=multi_sim_comm.Isend(dexp_, rank_f_, 123);
   multi_sim_comm.Recv(dexp_f_, rank_f_, 123);
   req.wait();
  } else {
   for(unsigned i=0; i<dexp_.size(); ++i){
     dexp_[i] = Vector(0,0,0);
     dexp_f_[i] = Vector(0,0,0);
   }
  }
  // wait for things to be done
  multi_sim_comm.Barrier();
  // local communication
  if(comm.Get_size()>1){
   comm.Sum(&dexp_[0],   3*dexp_.size());
   comm.Sum(&dexp_f_[0], 3*dexp_f_.size());
  }
  // get time step 
  long int step = getStep();
  // do MC stuff at the right time step
  if(step%MCstride_==0&&!getExchangeStep()) doMonteCarlo(step);
  
  // clear derivatives
  for(unsigned i=0; i<atom_der_.size(); ++i) atom_der_[i] = Vector(0,0,0);
  
  // calculate energy
  double rho_A, rho_B;
  double ene = 0.0;
  // cycle on data points, in parallel
  for(unsigned i=rank_local_;i<dexp_.size();i=i+size_local_){
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
    // calculate derivative of energy with respect to fmod
    double dene_dfmod = - kbt_ / tmp1 * 2.0 * tmp0 / fmod ;
    // derivative of fmod with respect to rho
    double dfmod_drho; 
    if(is_state_A_){
      // calculate derivative of fmod with respect to rho_A
      dfmod_drho =  - W_ / fmod_B * (fmod_A - 1.0) * beta_;
    } else {
      // calculate derivative of fmod with respect to rho_B
      dfmod_drho =  + W_ * fmod_A / fmod_B / fmod_B * (fmod_B - 1.0) * beta_;
    }
    // calculate derivative of rho with respect to distance
    // rho = 1.0 - 1.0 / (1.0+exp(-alpha*(dist-Rslope)));
    double tmp_rho = exp(-alpha_*(dexp_[i].modulo()-Rslope_));
    double drho_ddist = - 1.0 / ( 1.0 + tmp_rho ) / ( 1.0 + tmp_rho ) * tmp_rho * alpha_;
    // retrieve atom indexes
    unsigned id0 = pair_list_[i].first;
    unsigned id1 = pair_list_[i].second;
    // add to derivatives with respect to atom positions
    atom_der_[id0] -= dene_dfmod * dfmod_drho * drho_ddist * dexp_[i] / dexp_[i].modulo();
    atom_der_[id1] += dene_dfmod * dfmod_drho * drho_ddist * dexp_[i] / dexp_[i].modulo();
  }
  // gather contributions if parallel
  if(!serial_ && comm.Get_size()>1){
    comm.Sum(&ene, 1);
    comm.Sum(&atom_der_[0], 3*atom_der_.size());
  }
  // add Jeffrey's priors + constant terms in sigma;
  ene += kbt_ * ( 1.0 - static_cast<double>(dexp_.size()) ) * std::log(sigma_);
  // set value of the bias
  getPntrToComponent("score")->set(ene);
  // set derivatives
  for(unsigned i=0; i<atom_der_.size(); ++i) setAtomsDerivatives(getPntrToComponent("score"), i, atom_der_[i]);
}


}
}


