/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015,2016 The plumed team
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
/* 
 The sigma_mean optimisation method was developed by Thomas Loehr and
 Carlo Camilloni 
*/


#include "Bias.h"
#include "ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "core/Value.h"
#include <cmath>

using namespace std;

namespace PLMD{
namespace bias{

//+PLUMEDOC BIAS METAINFERENCE
/*
Calculate the Metainference Score for a set of back calculated experimental data.

*/
//+ENDPLUMEDOC

class Metainference : public Bias
{
  const double sqrt2_div_pi;
  // experimental values
  vector<double> parameters;
  // noise type
  unsigned noise_type_;
  enum { GAUSS, MGAUSS, OUTLIERS};
  // scale is data scaling factor
  bool   doscale_;
  double scale_;
  double scale_min_;
  double scale_max_;
  double Dscale_;
  // sigma is data uncertainty
  vector<double> sigma_;
  double sigma_min_;
  double sigma_max_;
  double Dsigma_;
  // sigma_mean is uncertainty in the mean estimate
  double sigma_mean_;
  double sigma_mean_min_;
  double sigma_mean_max_;
  double Dsigma_mean_;
  double max_plumed_force_;
  double max_md_force_;
  unsigned sum_nfmax_;
  // temperature in kbt
  double   kbt_;
  // number of data points
  unsigned ndata_;
  // Sigma mean minimisation
  bool optsigmamean_;
  // Monte Carlo stuff
  double   old_energy;
  unsigned MCsteps_;
  unsigned MCstride_;
  unsigned MCaccept_;
  long int MCfirst_;
  // output
  Value* valueBias;
  Value* valueScale;
  Value* valueAccept;
  vector<Value*> valueSigma;
  Value* valueSigmaMean;
  Value* valueSigmaMeanMin;
  Value* valueMaxForceMD;
  Value* valueMaxForcePLUMED;
  Value* valueFracViol;

  unsigned nrep_;
  unsigned replica_;

  Atoms& atoms;

  double getEnergySPE(const vector<double> &sigma, const double scale);
  double getEnergyGJE(const vector<double> &sigma, const double scale);
  void   doMonteCarlo();
  double getEnergyForceSPE();
  double getEnergyForceGJE();
  
public:
  explicit Metainference(const ActionOptions&);
  void calculate();
  void update();
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(Metainference,"METAINFERENCE")

void Metainference::registerKeywords(Keywords& keys){
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("optional","PARARG","reference values for the experimental data, these can be provided as arguments without derivatives"); 
  keys.add("optional","PARAMETERS","reference values for the experimental data");
  keys.add("compulsory","NOISETYPE","functional form of the noise (GAUSS,MGAUSS,OUTLIERS)");
  keys.addFlag("SCALEDATA",false,"Set to TRUE if you want to sample a scaling factor common to all values and replicas");  
  keys.addFlag("OPTSIGMAMEAN", false, "Set to minimize sigma_mean on the fly");
  keys.add("compulsory","SCALE0","initial value of the uncertainty parameter");
  keys.add("compulsory","SCALE_MIN","minimum value of the uncertainty parameter");
  keys.add("compulsory","SCALE_MAX","maximum value of the uncertainty parameter");
  keys.add("compulsory","DSCALE","maximum MC move of the uncertainty parameter");
  keys.add("compulsory","SIGMA0","initial value of the uncertainty parameter");
  keys.add("compulsory","SIGMA_MIN","minimum value of the uncertainty parameter");
  keys.add("compulsory","SIGMA_MAX","maximum value of the uncertainty parameter");
  keys.add("compulsory","DSIGMA","maximum MC move of the uncertainty parameter");
  keys.add("compulsory","SIGMA_MEAN0","starting value for the uncertainty in the mean estimate");
  keys.add("compulsory","SIGMA_MEAN_MIN","minimum uncertainty in the mean estimate");
  keys.add("compulsory","SIGMA_MEAN_MAX","maximum uncertainty in the mean estimate");
  keys.add("compulsory","DSIGMA_MEAN","maximum MC move of the uncertainty parameter");
  keys.add("compulsory","MAX_FORCE_MD","maximum allowable force");
  keys.add("compulsory","MAX_FORCE_PLUMED","maximum allowable force");
  keys.add("optional","TEMP","the system temperature - this is only needed if code doesnt' pass the temperature to plumed");
  keys.add("optional","MC_STEPS","number of MC steps");
  keys.add("optional","MC_STRIDE","MC stride");
  componentsAreNotOptional(keys);
  useCustomisableComponents(keys);
  keys.addOutputComponent("bias","default","the instantaneous value of the bias potential");
  keys.addOutputComponent("sigma", "default","uncertainty parameter");
  keys.addOutputComponent("scale", "default","scale parameter");
  keys.addOutputComponent("accept","default","MC acceptance");
  keys.addOutputComponent("sigmaMean","default","uncertainty in the mean estimate");
  keys.addOutputComponent("maxForceMD","default","max force on molecule");
  keys.addOutputComponent("maxForcePLUMED","default","max force on molecule");
  keys.addOutputComponent("fracViol","default","fraction of max force violations");
}

Metainference::Metainference(const ActionOptions&ao):
PLUMED_BIAS_INIT(ao), 
sqrt2_div_pi(0.45015815807855),
doscale_(false),
sum_nfmax_(0),
ndata_(getNumberOfArguments()),
optsigmamean_(false),
old_energy(0),
MCsteps_(1), 
MCstride_(1), 
MCaccept_(0), 
MCfirst_(-1),
atoms(plumed.getAtoms())
{
  parseVector("PARAMETERS",parameters);
  if(parameters.size()!=static_cast<unsigned>(getNumberOfArguments())&&!parameters.empty())
    error("Size of PARAMETERS array should be either 0 or the same as of the number of arguments in ARG1");

  vector<Value*> arg2;
  parseArgumentList("PARARG",arg2);

  if(!arg2.empty()) {
    if(parameters.size()>0) error("It is not possible to use PARARG and PARAMETERS together");
    if(arg2.size()!=getNumberOfArguments()) error("Size of PARARG array should be the same as number for arguments in ARG");
    for(unsigned i=0;i<arg2.size();i++){
      parameters.push_back(arg2[i]->get()); 
      if(arg2[i]->hasDerivatives()==true) error("PARARG can only accept arguments without derivatives");
    }
  }

  if(parameters.size()!=getNumberOfArguments()) 
    error("PARARG or PARAMETERS arrays should include the same number of elements as the arguments in ARG");

  string stringa_noise;
  parse("NOISETYPE",stringa_noise);
  if(stringa_noise=="GAUSS")       noise_type_ = GAUSS; 
  else if(stringa_noise=="MGAUSS") noise_type_ = MGAUSS;
  else if(stringa_noise=="OUTLIERS")  noise_type_ = OUTLIERS;
  else error("Unkwnow noise type"); 

  parseFlag("SCALEDATA", doscale_);
  if(doscale_) {
    parse("SCALE0",scale_);
    parse("SCALE_MIN",scale_min_);
    parse("SCALE_MAX",scale_max_);
    parse("DSCALE",Dscale_);
  } else {
    scale_=1.0;
  }

  parseFlag("OPTSIGMAMEAN", optsigmamean_);

  vector<double> readsigma;
  parseVector("SIGMA0",readsigma);
  if(noise_type_!=MGAUSS&&readsigma.size()>1) error("If you want to use more than one sigma you should use NOISETYPE=MGAUSS");
  if(noise_type_==MGAUSS) {
    if(readsigma.size()==getNumberOfArguments()) {
      sigma_.resize(getNumberOfArguments());
      sigma_=readsigma;
    } else if(readsigma.size()==1) {
      sigma_.resize(getNumberOfArguments(),readsigma[0]);
    } else {
      error("SIGMA0 can accept either one single value or as many values as the number of arguments (with NOISETYPE=MGAUSS)");
    } 
  } else sigma_.resize(1, readsigma[0]);

  parse("SIGMA_MIN",sigma_min_);
  parse("SIGMA_MAX",sigma_max_);
  parse("DSIGMA",Dsigma_);
  parse("MC_STEPS",MCsteps_);
  parse("MC_STRIDE",MCstride_);
  // get temperature
  double temp=0.0;
  parse("TEMP",temp);
  parse("SIGMA_MEAN0",sigma_mean_);

  // sigma_mean params
  if(optsigmamean_) {
    parse("SIGMA_MEAN_MIN", sigma_mean_min_);
    parse("SIGMA_MEAN_MAX", sigma_mean_max_);
    parse("DSIGMA_MEAN", Dsigma_mean_);
    parse("MAX_FORCE_PLUMED", max_plumed_force_);
    max_plumed_force_ *= max_plumed_force_;
    parse("MAX_FORCE_MD", max_md_force_);
    max_md_force_ *= max_md_force_; 
  }

  checkRead();

  if(temp>0.0) kbt_=plumed.getAtoms().getKBoltzmann()*temp;
  else kbt_=plumed.getAtoms().getKbT();

  // get number of replicas
  if(comm.Get_rank()==0) {
    nrep_ = multi_sim_comm.Get_size();
    replica_ = multi_sim_comm.Get_rank();
  } else {
    nrep_ = 0;
    replica_ = 0;
  }
  comm.Sum(&nrep_,1);
  comm.Sum(&replica_,1);

  // adjust for multiple-time steps
  MCstride_ *= getStride();

  switch(noise_type_) {
    case GAUSS:
      log.printf("  with gaussian noise and a single noise parameter for all the data\n");
      break;
    case MGAUSS:
      log.printf("  with gaussian noise and a noise parameter for each data point\n");
      break;
    case OUTLIERS:
      log.printf("  with long tailed gaussian noise and a single noise parameter for all the data\n");
      break;
  }

  if(doscale_) {
    log.printf("  sampling a common scaling factor with:\n");
    log.printf("    initial scale parameter %f\n",scale_);
    log.printf("    minimum scale parameter %f\n",scale_min_);
    log.printf("    maximum scale parameter %f\n",scale_max_);
    log.printf("    maximum MC move of scale parameter %f\n",Dscale_);
  }

  if(readsigma.size()==1) log.printf("  initial data uncertainty %f\n",sigma_[0]);
  else {
    log.printf("  initial data uncertainties");
    for(unsigned i=0;i<sigma_.size();++i) log.printf(" %f", sigma_[i]);
    log.printf("\n");
  }
  log.printf("  minimum data uncertainty %f\n",sigma_min_);
  log.printf("  maximum data uncertainty %f\n",sigma_max_);
  log.printf("  maximum MC move of data uncertainty %f\n",Dsigma_);
  log.printf("  uncertainty in the mean estimate %f\n",sigma_mean_);
  log.printf("  temperature of the system %f\n",kbt_);
  log.printf("  number of experimental data points %u\n",getNumberOfArguments());
  log.printf("  number of replicas %u\n",nrep_);
  log.printf("  MC steps %u\n",MCsteps_);
  log.printf("  MC stride %u\n",MCstride_);

  if(optsigmamean_) {
    log.printf("  minimum data uncertainty of the mean %f\n", sigma_mean_min_);
    log.printf("  maximum data uncertainty of the mean %f\n", sigma_mean_max_);
    log.printf("  maximum move of data uncertainty of the mean %f\n", Dsigma_mean_);
    log.printf("  maximum permissible force %f\n", sqrt(max_plumed_force_));
    log.printf("  maximum permissible force %f\n", sqrt(max_md_force_));
  }

  addComponent("bias");
  componentIsNotPeriodic("bias");
  valueBias=getPntrToComponent("bias");
  if(doscale_) { 
    addComponent("scale");  
    componentIsNotPeriodic("scale");
    valueScale=getPntrToComponent("scale");
  }
  addComponent("accept");
  componentIsNotPeriodic("accept");
  valueAccept=getPntrToComponent("accept");

  if(optsigmamean_) {
    addComponent("sigmaMean");
    componentIsNotPeriodic("sigmaMean");
    valueSigmaMean=getPntrToComponent("sigmaMean");
    addComponent("sigmaMeanMin");
    componentIsNotPeriodic("sigmaMeanMin");
    valueSigmaMeanMin=getPntrToComponent("sigmaMeanMin");
    addComponent("maxForceMD");
    componentIsNotPeriodic("maxForceMD");
    valueMaxForceMD=getPntrToComponent("maxForceMD");
    addComponent("maxForcePLUMED");
    componentIsNotPeriodic("maxForcePLUMED");
    valueMaxForcePLUMED=getPntrToComponent("maxForcePLUMED");
    addComponent("fracViol");
    componentIsNotPeriodic("fracViol");
    valueFracViol=getPntrToComponent("fracViol");
  }

  if(noise_type_==MGAUSS) {
    for(unsigned i=0;i<sigma_.size();++i){
      std::string num; Tools::convert(i,num);
      addComponent("sigma_"+num); componentIsNotPeriodic("sigma_"+num);
      valueSigma.push_back(getPntrToComponent("sigma_"+num));
    }
  } else {
    addComponent("sigma"); componentIsNotPeriodic("sigma");
    valueSigma.push_back(getPntrToComponent("sigma"));
  }
  // initialize random seed
  unsigned iseed;
  if(comm.Get_rank()==0) iseed = time(NULL)+replica_;
  else iseed = 0;     
  comm.Sum(&iseed, 1);
  srand(iseed);

  log<<"  Bibliography "<<plumed.cite("Bonomi, Camilloni, Cavalli, Vendruscolo, Sci. Adv. 2, e150117 (2016)");
}

double Metainference::getEnergySPE(const vector<double> &sigma, const double scale){
  // calculate effective sigma
  const double smean2 = sigma_mean_*sigma_mean_;
  const double s = sqrt( sigma[0]*sigma[0] + smean2 );
  // cycle on arguments
  double ene = 0.0;
  for(unsigned i=0;i<getNumberOfArguments();++i){
    const double dev = scale*getArgument(i)-parameters[i]; 
    // argument
    const double a2 = 0.5*dev*dev + s*s;
    // increment energy
    ene += std::log( 2.0 * a2 / ( 1.0 - exp(- a2 / smean2) ) );
  }
  // add normalization and Jeffrey's prior
  ene += std::log(s) - static_cast<double>(ndata_)*std::log(sqrt2_div_pi*s);
  return kbt_ * ene;
}

double Metainference::getEnergyGJE(const vector<double> &sigma, const double scale){
  // cycle on arguments
  double ene = 0.0;
  const double smean2 = sigma_mean_*sigma_mean_;
  double ss = sigma[0]*sigma[0] + smean2; 
  for(unsigned i=0;i<getNumberOfArguments();++i){
    if(noise_type_==MGAUSS){ 
      ss = sigma[i]*sigma[i] + smean2;
      // add Jeffrey's prior - one per sigma
      ene += 0.5*std::log(ss);
    }
    const double dev = scale*getArgument(i)-parameters[i]; 
    ene += 0.5*dev*dev/ss + 0.5*std::log(ss);
  }
  // add Jeffrey's prior in case one sigma for all data points
  if(noise_type_==GAUSS) ene += 0.5*std::log(ss);
  return kbt_ * ene;
}

void Metainference::doMonteCarlo(){
  // calculate old energy (first time only)
  if(MCfirst_==-1) {
    switch(noise_type_) {
      case GAUSS:
      case MGAUSS:
        old_energy = getEnergyGJE(sigma_,scale_);
        break;
      case OUTLIERS:
        old_energy = getEnergySPE(sigma_,scale_);
        break;
    }
  }
 
  // cycle on MC steps 
  for(unsigned i=0;i<MCsteps_;++i){
 
    // propose move for scale
    double new_scale = scale_;
    if(doscale_) {
      const double r1 = static_cast<double>(rand()) / RAND_MAX;
      const double ds1 = -Dscale_ + r1 * 2.0 * Dscale_;
      new_scale += ds1;
      // check boundaries
      if(new_scale > scale_max_){new_scale = 2.0 * scale_max_ - new_scale;}
      if(new_scale < scale_min_){new_scale = 2.0 * scale_min_ - new_scale;}
      // the scaling factor should be the same for all the replicas
      if(comm.Get_rank()==0) multi_sim_comm.Bcast(new_scale,0);
      comm.Bcast(new_scale,0);
    }
  
    // propose move for sigma
    vector<double> new_sigma(sigma_.size());
    for(unsigned j=0;j<sigma_.size();j++) {
      const double r2 = static_cast<double>(rand()) / RAND_MAX;
      const double ds2 = -Dsigma_ + r2 * 2.0 * Dsigma_;
      new_sigma[j] = sigma_[j] + ds2;
      // check boundaries
      if(new_sigma[j] > sigma_max_){new_sigma[j] = 2.0 * sigma_max_ - new_sigma[j];}
      if(new_sigma[j] < sigma_min_){new_sigma[j] = 2.0 * sigma_min_ - new_sigma[j];}
    }
 
    // calculate new energy
    double new_energy=0;
    switch(noise_type_) {
      case GAUSS:
      case MGAUSS:
        new_energy = getEnergyGJE(new_sigma,new_scale);
        break;
      case OUTLIERS:
        new_energy = getEnergySPE(new_sigma,new_scale);
        break;
    }
    // accept or reject
    const double delta = ( new_energy - old_energy ) / kbt_;
    // if delta is negative always accept move
    if( delta <= 0.0 ){
      old_energy = new_energy;
      scale_ = new_scale;
      sigma_ = new_sigma;
      MCaccept_++;
    // otherwise extract random number
    } else {
      const double s = static_cast<double>(rand()) / RAND_MAX;
      if( s < exp(-delta) ){
        old_energy = new_energy;
        scale_ = new_scale;
        sigma_ = new_sigma;
        MCaccept_++;
      }
    }
 
    if(doscale_) {
      // the scaling factor should be the same for all the replicas
      if(comm.Get_rank()==0) multi_sim_comm.Bcast(scale_,0);
      comm.Bcast(scale_,0);
    }
  }
  /* save the result of the sampling */
  if(doscale_) valueScale->set(scale_);
  for(unsigned i=0; i<sigma_.size(); i++) valueSigma[i]->set(sigma_[i]);
}

double Metainference::getEnergyForceSPE()
{
  double ene = 0.0;
  const unsigned narg=getNumberOfArguments();

  const double smean2 = sigma_mean_*sigma_mean_; 
  const double s = sqrt( sigma_[0]*sigma_[0] + smean2 );
  vector<double> f(narg,0);
  
  if(comm.Get_rank()==0){
   for(unsigned i=0;i<narg;++i){
     const double dev = scale_*getArgument(i)-parameters[i]; 
     const double a2 = 0.5*dev*dev + s*s;
     const double t = exp(-a2/smean2);
     const double dt = 1./t;
     const double it = 1./(1.-t);
     const double dit = 1./(1.-dt);
     ene += std::log(2.*a2*it);
     f[i] = -scale_*dev*(dit/smean2 + 1./a2);
   }
   // collect contribution to forces and energy from other replicas
   multi_sim_comm.Sum(&f[0],narg);
   multi_sim_comm.Sum(&ene,1);
   // add normalizations and priors of local replica
   ene += std::log(s) - static_cast<double>(ndata_)*std::log(sqrt2_div_pi*s);
  }
  // intra-replica summation
  comm.Sum(&f[0],narg);
  comm.Sum(&ene,1);

  for(unsigned i=0; i<narg; ++i) setOutputForce(i, kbt_ * f[i]);
  return ene;
}

double Metainference::getEnergyForceGJE()
{
  double ene = 0.0;

  const unsigned ssize = sigma_.size();
  vector<double> ss(ssize);
  vector<double> inv_s2(ssize, 0.);
  const double smean2 = sigma_mean_*sigma_mean_;

  for(unsigned i=0;i<ssize; ++i) {
    ss[i] = sigma_[i]*sigma_[i] + smean2;
    if(comm.Get_rank()==0) inv_s2[i] = 1.0/ss[i];
  }

  if(comm.Get_rank()==0) multi_sim_comm.Sum(&inv_s2[0],ssize); 
  comm.Sum(&inv_s2[0],ssize);  
  
  const unsigned narg=getNumberOfArguments();
  for(unsigned i=0;i<narg;++i){
    const double dev = scale_*getArgument(i)-parameters[i]; 
    unsigned sel_sigma=0;
    if(noise_type_==MGAUSS){
      sel_sigma=i;
      // add Jeffrey's prior - one per sigma
      ene += 0.5*std::log(ss[sel_sigma]);
    }
    ene += 0.5*dev*dev*inv_s2[sel_sigma] + 0.5*std::log(ss[sel_sigma]);
    setOutputForce(i, -kbt_*dev*scale_*inv_s2[sel_sigma]);
  }
  // add Jeffrey's prior in case one sigma for all data points
  if(noise_type_==GAUSS) ene += 0.5*std::log(ss[0]);
  return ene;
}

void Metainference::calculate(){
  /* MONTE CARLO */
  const long int step = getStep();
  if(step%MCstride_==0&&!getExchangeStep()) doMonteCarlo();
  // this is needed when restarting simulations
  if(MCfirst_==-1) MCfirst_=step;
  const double MCtrials = floor(static_cast<double>(step-MCfirst_) / static_cast<double>(MCstride_))+1.0;
  const double accept = static_cast<double>(MCaccept_) / static_cast<double>(MCsteps_) / MCtrials;
  valueAccept->set(accept);

  // calculate bias and forces
  double ene = 0; 
  switch(noise_type_) {
    case GAUSS:
    case MGAUSS:
      ene = getEnergyForceGJE();
      break;
    case OUTLIERS:
      ene = getEnergyForceSPE();
      break;
  }
  // set value of the bias
  valueBias->set(kbt_*ene);
}

void Metainference::update() {
  if(optsigmamean_) {
    const double EPS = 0.1;
    const long int step = getStep();
    // Get max force of whole system
    vector<Vector> md_forces;
    vector<Vector> plumed_forces;
    vector<double> masses;

    atoms.getLocalMDForces(md_forces);
    atoms.getLocalForces(plumed_forces);
    atoms.getLocalMasses(masses);

    vector<double> allforces_plumed;
    vector<double> allforces_md;
    for(unsigned i = 0; i < plumed_forces.size(); ++i) {
      const double pf2 = plumed_forces[i].modulo2();
      // we are only interested in forces plumed has an effect on
      if(pf2 > EPS && masses[i] > EPS ) {
        const double invm2 = 1./(masses[i]*masses[i]);
        allforces_md.push_back(md_forces[i].modulo2()*invm2);
        allforces_plumed.push_back(pf2*invm2);
      }
    }

    vector<double> fmax_tmp_md(comm.Get_size(),0);
    vector<double> fmax_tmp_plumed(comm.Get_size(),0);
    unsigned nfmax = 0;
    unsigned nfmax_md = 0;

    /* each local thread should look for the maximum force and number of violations */
    if(allforces_md.size()>0) {
      fmax_tmp_md[comm.Get_rank()]     = *max_element(allforces_md.begin(), allforces_md.end());
      fmax_tmp_plumed[comm.Get_rank()] = *max_element(allforces_plumed.begin(), allforces_plumed.end());
      for(unsigned i = 0; i < allforces_md.size(); ++i) {
        if(allforces_plumed[i] > max_plumed_force_) {
          nfmax += 1;
          if(allforces_md[i] > max_md_force_) {
            nfmax_md += 1;
          }
        }
      }
    }

    // the largest forces are shared among the local threads but not over the replicas 
    comm.Sum(fmax_tmp_md);
    comm.Sum(fmax_tmp_plumed);
    // these are the largest forces for a specific replica
    const double fmax_md     = *max_element(fmax_tmp_md.begin(), fmax_tmp_md.end());
    const double fmax_plumed = *max_element(fmax_tmp_plumed.begin(), fmax_tmp_plumed.end());

    // the number of violations is summed up over the local thread and over the replicas 
    comm.Sum(nfmax);
    comm.Sum(nfmax_md);
    if(comm.Get_rank() == 0) {
      multi_sim_comm.Sum(nfmax);
      multi_sim_comm.Sum(nfmax_md);
    }
    // and resent to all the threads
    comm.Bcast(nfmax,0);
    comm.Bcast(nfmax_md,0);

    // if no violations then we minimise sigma_mean_ (and also sigma_mean_min_ if needed)
    if(nfmax == 0) {
      sigma_mean_ -= Dsigma_mean_ * 0.01 * std::log(sigma_mean_/sigma_mean_min_);
      if((sigma_mean_-sigma_mean_min_)<Dsigma_mean_) {
        double sigma_mean_min_new = sigma_mean_min_ - Dsigma_mean_;
        if(sigma_mean_min_new > 0) sigma_mean_min_ = sigma_mean_min_new;
      }
    // otherwise we increase sigma_mean_ and also sigma_mean_min_ if needed
    } else {
      sum_nfmax_ += nfmax;
      sigma_mean_ += Dsigma_mean_ * nfmax;
      if(nfmax_md > 0) {
        sigma_mean_min_ += nfmax_md * Dsigma_mean_;
        sigma_mean_ += Dsigma_mean_ * nfmax_md;
      }
    }

    valueMaxForceMD->set(sqrt(fmax_md));
    valueMaxForcePLUMED->set(sqrt(fmax_plumed));
    valueSigmaMean->set(sigma_mean_);
    valueSigmaMeanMin->set(sigma_mean_min_);
    if(step > 0) { 
      valueFracViol->set(static_cast<double>(sum_nfmax_) / static_cast<double>(step));
    }
  }
}

}
}


