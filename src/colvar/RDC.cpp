/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014 The plumed team
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

#ifdef __PLUMED_HAS_GSL
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#endif

using namespace std;

namespace PLMD{

//+PLUMEDOC COLVAR RDC 
/*
Calculates the Residual Dipolar Coupling between two atoms. 

The RDC between two atomic nuclei depends on the \f$\theta\f$ angle between 
the inter-nuclear vector and the external magnetic field. In isotropic media RDCs average to zero because of the orientational 
averaging, but when the rotational symmetry is broken, either through the introduction of an alignment medium or for molecules 
with highly anisotropic paramagnetic susceptibility, RDCs become measurable.

\f[
D=D_{max}0.5(3\cos^2(\theta)-1)
\f]

where

\f[
D_{max}=-\mu_0\gamma_1\gamma_2h/(8\pi^3r^3)
\f]

that is the maximal value of the dipolar coupling for the two nuclear spins with gyromagnetic ratio \f$\gamma\f$. 
\f$\mu\f$ is the magnetic constant and h is the Planck constant. 

Common Gyromagnetic Ratios (C.G.S)
- H(1) 26.7513
- C(13) 6.7261
- N(15) -2.7116
- NH -72.5388
- CH 179.9319
- CN -18.2385
- CC 45.2404

This collective variable calculates the Residual Dipolar Coupling for a set of couple of atoms using the above definition. 
From the calculated RDCs and a set of experimental values it calculates either their correlation or the squared quality factor 

\f[
Q^2=\frac{\sum_i(D_i-D^{exp}_i)^2}{\sum_i(D^{exp}_i)^2}
\f]

RDCs report only on the fraction of molecules that is aligned, this means that comparing the RDCs from a single structure in
a MD simulation to the experimental values is not particularly meaningfull, from this point of view it is better to compare
their correlation. The fraction of aligned molecules can be obtained by maximising the correlation between the calculated and 
the experimental RDCs. This fraction can be used as a scaling factor in the calculation of the RDCs in order to compare their
values. The averaging of the RDCs calculated with the above definition from a standard MD should result to 0 because of
the rotational diffusion, but this variable can be used to break the rotational symmetry.

RDCs can also be calculated using a Single Value Decomposition approach, in this case the code rely on the
a set of function from the GNU Scientific Library (GSL). (With SVD forces are not currently implemented).

Replica-Averaged restrained simulations can be performed with this CV using the ENSEMBLE flag.

Additional material and examples can be also found in the tutorial \ref belfast-9

\par Examples
In the following example five N-H RDCs are defined and their correlation with
respect to a set of experimental data is calculated.  

\verbatim
RDC ...
GYROM=-72.5388
SCALE=1.0 
TYPE=CORRELATION
ATOMS1=20,21 COUPLING1=8.17
ATOMS2=37,38 COUPLING2=-8.271
ATOMS3=56,57 COUPLING3=-10.489
ATOMS4=76,77 COUPLING4=-9.871
ATOMS5=92,93 COUPLING5=-9.152
LABEL=nh
... RDC 

rdce: RESTRAINT ARG=nh KAPPA=0. SLOPE=-25000.0 AT=1.

PRINT ARG=nh,rdce.bias FILE=colvar
\endverbatim
(See also \ref PRINT, \ref RESTRAINT) 

*/
//+ENDPLUMEDOC

class RDC : public Colvar {
private:
  enum CV_TYPE   {QF, SDEV, CORR, COMP, SVD};
  const double   Const;
  vector<double> coupl;
  vector<double> mu_s;
  vector<double> scale;
  unsigned       ens_dim;
  unsigned       pperiod;
  unsigned       type;
  double         norm;
  bool           ensemble;
  bool           serial;
public:
  explicit RDC(const ActionOptions&);
  static void registerKeywords( Keywords& keys );
  virtual void calculate();
  virtual void qf_calc(bool qf, bool docomp);
  virtual void corr_calc();
  virtual void svd_calc();
};

PLUMED_REGISTER_ACTION(RDC,"RDC")

void RDC::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords( keys );
  componentsAreNotOptional(keys);
  keys.add("numbered","ATOMS","the couple of atoms involved in each of the bonds for which you wish to calculate the RDC. "
                             "Keywords like ATOMS1, ATOMS2, ATOMS3,... should be listed and one dipolar coupling will be "
                             "calculated for each ATOMS keyword you specify.");
  keys.reset_style("ATOMS","atoms");
  keys.add("compulsory","TYPE","QFACTOR","The type of calculation relative to the Residual Dipolar Coupling:"
                                         "QFACTOR, SQUAREDEV, CORRELATION, COMPONENTS, SVD.");
  keys.add("numbered","COUPLING","Add an experimental value for each coupling. ");
  keys.add("numbered","GYROM","Add the product of the gyromagnetic constants for each bond. ");
  keys.add("numbered","SCALE","Add a scaling factor to take into account concentration and other effects. ");
  keys.add("compulsory","WRITE_DC","0","Write the back-calculated dipolar couplings every # steps.");
  keys.addFlag("ENSEMBLE",false,"Set to TRUE if you want to average over multiple replicas.");  
  keys.addFlag("SERIAL",false,"Set to TRUE if you want to run the CV in serial.");  
  keys.addOutputComponent("bond_","COMPONENTS","the squared deviation of the RDC for bond #");
  keys.addOutputComponent("corr","CORRELATION","the correlation between calculated and experimental RDCs");
  keys.addOutputComponent("slope","CORRELATION","the slope between calculated and experimental RDCs");
  keys.addOutputComponent("inter","CORRELATION","the intercept between calculated and experimental RDCs");
}

RDC::RDC(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
Const(0.3356806),
type(QF),
norm(0.0)
{
  // Read in the atoms
  vector<AtomNumber> t, atoms;
  for(int i=1;;++i ){
     parseAtomList("ATOMS", i, t );
     if( t.empty() ) break;
     if( t.size()!=2 ){
         std::string ss; Tools::convert(i,ss);
         error("ATOMS" + ss + " keyword has the wrong number of atoms");
     }
     atoms.push_back(t[0]);
     atoms.push_back(t[1]);
     t.resize(0); 
  }

  std::string What;
  parse("TYPE",What);
  if(What=="QFACTOR")          type=QF;
  else if(What=="SQUAREDEV")   type=SDEV;
  else if(What=="CORRELATION") type=CORR;
  else if(What=="COMPONENTS")  type=COMP;
  else if(What=="SVD")         type=SVD;
  else error("Unrecognized calculation TYPE!\n");

#ifndef __PLUMED_HAS_GSL
  if(type==SVD) error("You CANNOT use SVD without GSL. Recompile PLUMED with GSL!\n");
#endif


  // define the number of bonds
  coupl.resize( atoms.size()/2 ); 
  unsigned ntarget=0;

  //if(type!=COMP) {
    // Read in RDC values
    for(unsigned i=0;i<coupl.size();++i){
       if( !parseNumbered( "COUPLING", i+1, coupl[i] ) ) break;
       ntarget++; 
    }
    if( ntarget!=coupl.size() ) error("found wrong number of COUPLING values");

    // Read in GYROMAGNETIC constants 
    mu_s.resize( coupl.size() ); 
    ntarget=0;
    for(unsigned i=0;i<coupl.size();++i){
       if( !parseNumbered( "GYROM", i+1, mu_s[i] ) ) break;
       ntarget++; 
    }
    if( ntarget==0 ){
        parse("GYROM",mu_s[0]);
        for(unsigned i=1;i<coupl.size();++i) mu_s[i]=mu_s[0];
    } else if( ntarget!=coupl.size() ) error("found wrong number of GYROM values");
  //}

  // Read in SCALING factors 
  scale.resize( coupl.size() );
  for(unsigned i=0;i<coupl.size();++i) scale[i]=1.0;
  ntarget=0;
  for(unsigned i=0;i<coupl.size();++i){
     if( !parseNumbered( "SCALE", i+1, scale[i] ) ) break;
     ntarget++; 
  }
  if( ntarget==0 ){
      parse("SCALE",scale[0]);
      for(unsigned i=1;i<coupl.size();++i) scale[i]=scale[0];
  } else if( ntarget!=coupl.size() ) error("found wrong number of SCALE values");

  ensemble=false;
  parseFlag("ENSEMBLE",ensemble);
  if(ensemble){
    if(comm.Get_rank()==0) {
      if(multi_sim_comm.Get_size()<2) error("You CANNOT run Replica-Averaged simulations without running multiple replicas!\n");
      ens_dim=multi_sim_comm.Get_size();
    } else ens_dim=0;
    comm.Sum(&ens_dim, 1);
  } else ens_dim=1;

  serial=false;
  parseFlag("SERIAL",serial);

  int w_period=0;
  parse("WRITE_DC", w_period);
  pperiod=w_period;

  // Ouput details of all contacts 
  for(unsigned i=0;i<coupl.size();++i){
    log.printf("  The %dth Bond Dipolar Coupling is calculated from atoms : %d %d.", i+1, atoms[2*i].serial(), atoms[2*i+1].serial()); 
    log.printf("  Dipolar Coupling is %f. Gyromagnetic moment is %f. Scaling factor is %f.\n",coupl[i],mu_s[i],scale[i]);
  }
  if(ensemble) { log.printf("  ENSEMBLE averaging over %u replicas\n", ens_dim); }

  checkRead();

  if(type==QF||type==SDEV) {
    addValueWithDerivatives();
    setNotPeriodic();
    if(ensemble) setEnsemble(ens_dim);
    else setNotEnsemble();

  } else if(type==CORR) {
      addComponentWithDerivatives("corr"); componentIsNotPeriodic("corr");
      addComponentWithDerivatives("slope"); componentIsNotPeriodic("slope");
      addComponentWithDerivatives("inter"); componentIsNotPeriodic("inter");
      if(ensemble) {
        componentIsEnsemble("corr", ens_dim);
        componentIsEnsemble("slope", ens_dim);
        componentIsEnsemble("inter", ens_dim);
      } else {
        componentIsNotEnsemble("corr");
        componentIsNotEnsemble("slope");
        componentIsNotEnsemble("inter");
      }

  } else if(type==COMP) {
    for(unsigned i=0;i<coupl.size();i++) {
      std::string num; Tools::convert(i,num);
      addComponentWithDerivatives("bond_"+num); componentIsNotPeriodic("bond_"+num);
      if(ensemble) {
        componentIsEnsemble("bond_"+num, ens_dim);
      } else {
        componentIsNotEnsemble("bond_"+num);
      }

    }
  } else if(type==SVD) {
    addValue();
    setNotPeriodic();
    if(ensemble) setEnsemble(ens_dim);
    else setNotEnsemble();
  }


  requestAtoms(atoms);

  norm = 0.; 
  for(unsigned i=0; i<coupl.size(); i++) norm += coupl[i]*coupl[i];

  log.printf("  DONE!\n"); log.flush();
}

void RDC::calculate()
{
  switch(type) {
    case QF:   qf_calc(true,false);  break;
    case SDEV: qf_calc(false,false); break;
    case CORR: corr_calc();          break;
    case COMP: qf_calc(false,true);  break;
    case SVD:  svd_calc();           break;
  }
}

void RDC::qf_calc(bool qf, bool docomp)
{
  double score=0.;
  vector<double> rdc( coupl.size() );
  unsigned N = getNumberOfAtoms();
  vector<Vector> dRDC(N);

  // internal parallelisation
  unsigned stride=2*comm.Get_size();
  unsigned rank=2*comm.Get_rank();
  if(serial){
    stride=2;
    rank=0;
  }

  /* RDC Calculations and forces */
  for(unsigned r=rank;r<N;r+=stride)
  {
    unsigned index=r/2;
    Vector distance;
    distance=pbcDistance(getPosition(r),getPosition(r+1));
    double d    = distance.modulo();
    double ind  = 1./d;
    double id3  = ind*ind*ind; 
    double id7  = id3*id3*ind;
    double id9  = id7*ind*ind;
    double max  = -Const*scale[index]*mu_s[index];
    double dmax = id3*max;
    double cos_theta = distance[2]*ind;
    rdc[index] = 0.5*dmax*(3.*cos_theta*cos_theta-1.);
    double x2=distance[0]*distance[0];
    double y2=distance[1]*distance[1];
    double z2=distance[2]*distance[2];
    double prod = -max*id7*(1.5*x2 +1.5*y2 -6.*z2);
    dRDC[r][0] = prod*distance[0];
    dRDC[r][1] = prod*distance[1];
    dRDC[r][2] = -max*id9*distance[2]*(4.5*x2*x2 + 4.5*y2*y2 + 1.5*y2*z2 - 3.*z2*z2 + x2*(9.*y2 + 1.5*z2));
    if(!ensemble) {
      double delta = rdc[index]-coupl[index];
      score += delta*delta;
      dRDC[r] *= delta;
    }
    dRDC[r+1] = -dRDC[r];
  }

  // share the calculated rdc
  if(!serial) comm.Sum(&rdc[0],rdc.size());

  // Printout
  bool printout=false;
  if(pperiod>0) printout = (!(getStep()%pperiod));
  if(printout) {
    // share the calculated rdc
    //if(!serial) comm.Sum(&rdc[0],rdc.size());
    // only the master replica is going to write
    if(comm.Get_rank()==0) {
      char tmp1[21]; sprintf(tmp1, "%ld", getStep()); 
      string dcfile = string("rdc-")+getLabel()+"-"+tmp1+string(".dat");
      FILE *outfile = fopen(dcfile.c_str(), "w");
      fprintf(outfile, "#index calc exp\n");
      double sum=0.;
      double Q=0.;
      for(unsigned r=0;r<coupl.size();r++) { 
        fprintf(outfile," %4u %10.6f %10.6f\n", r, rdc[r], coupl[r]);
        sum+=(rdc[r]-coupl[r])*(rdc[r]-coupl[r]);
        Q+=(coupl[r]*coupl[r]);
      }
      fprintf(outfile, "# Q factor = %6.4f\n", sqrt(sum/Q));
      fclose(outfile);
    }
  }

  // Ensemble averaging
  double fact=1.0;
  if(ensemble) {
    fact = 1./((double) ens_dim);
    // share the calculated rdc
    //if(!serial) comm.Sum(&rdc[0],rdc.size());
    // I am the master of my replica
    if(comm.Get_rank()==0) {
      // among replicas
      multi_sim_comm.Sum(&rdc[0], coupl.size() );
      for(unsigned i=0;i<coupl.size();i++) rdc[i] *= fact; 
    } else for(unsigned i=0;i<coupl.size();i++) rdc[i] = 0.;
    // inside each replica
    comm.Sum(&rdc[0], coupl.size() );
    score = 0.;
    for(unsigned r=rank;r<N;r+=stride) {
      unsigned index=r/2;
      double delta = rdc[index]-coupl[index];
      dRDC[r]   *= delta;
      dRDC[r+1] *= delta;
      score += delta*delta;
    } 
  }
 
  Tensor virial;
  virial.zero();
  vector<Vector> deriv(N);
  if(!serial) comm.Sum(score);

  double tmpder = 2.*fact;
  if(qf) {
    score  /= norm;
    tmpder /= norm;
  } 

  if(!docomp) {
    for(unsigned r=rank;r<N;r+=stride) {
      deriv[r]   = tmpder*dRDC[r];
      deriv[r+1] = tmpder*dRDC[r+1];
      virial=virial+(-1.*Tensor(getPosition(r),deriv[r]));
      virial=virial+(-1.*Tensor(getPosition(r+1),deriv[r+1]));
    }

    if(!serial){
      if(!deriv.empty()) comm.Sum(&deriv[0][0],3*deriv.size());
      comm.Sum(virial);
    }

    for(unsigned i=0;i<N;i++) setAtomsDerivatives(i,deriv[i]);
    setValue           (score);
    setBoxDerivatives  (virial);

  } else {

    for(unsigned r=rank;r<N;r+=stride) {
      deriv[r]   = tmpder*dRDC[r];
      deriv[r+1] = tmpder*dRDC[r+1];
    }

    if(!serial){
      if(!deriv.empty()) comm.Sum(&deriv[0][0],3*deriv.size());
    }

    unsigned index=0;
    for(unsigned r=0;r<N;r+=2) {
      Value* val=getPntrToComponent(index);
      virial= -Tensor(getPosition(r),deriv[r  ])-Tensor(getPosition(r+1),deriv[r+1]);
      setAtomsDerivatives( val, r  , deriv[r  ]);
      setAtomsDerivatives( val, r+1, deriv[r+1]); 
      setBoxDerivatives( val, virial );
      double delta = (rdc[index]-coupl[index])*(rdc[index]-coupl[index]);
      val->set(delta);
      index++;
    }

  }
}

void RDC::corr_calc()
{
  double scx=0., scx2=0., scy=0., scxy=0.;
  vector<double> rdc( coupl.size() );
  unsigned N = getNumberOfAtoms();
  vector<Vector> dRDC(N);

  // internal parallelisation
  unsigned stride=2*comm.Get_size();
  unsigned rank=2*comm.Get_rank();
  if(serial){
    stride=2;
    rank=0;
  }

  /* RDC Calculations and forces */
  for(unsigned r=rank;r<N;r+=stride)
  {
    unsigned index=r/2;
    Vector distance;
    distance=pbcDistance(getPosition(r),getPosition(r+1));
    double d    = distance.modulo();
    double ind  = 1./d;
    double id3  = ind*ind*ind; 
    double id7  = id3*id3*ind;
    double id9  = id7*ind*ind;
    double max  = -Const*scale[index]*mu_s[index];
    double dmax = id3*max;
    double cos_theta = distance[2]*ind;
    rdc[index] = 0.5*dmax*(3.*cos_theta*cos_theta-1.);
    double x2=distance[0]*distance[0];
    double y2=distance[1]*distance[1];
    double z2=distance[2]*distance[2];
    double prod = -max*id7*(1.5*x2 +1.5*y2 -6.*z2);
    dRDC[r][0] = prod*distance[0];
    dRDC[r][1] = prod*distance[1];
    dRDC[r][2] = -max*id9*distance[2]*(4.5*x2*x2 + 4.5*y2*y2 + 1.5*y2*z2 - 3.*z2*z2 + x2*(9.*y2 + 1.5*z2));
    if(!ensemble) {
      scx  += rdc[index];
      scx2 += rdc[index]*rdc[index];
      scy  += coupl[index];
      scxy += rdc[index]*coupl[index];
    }
    dRDC[r+1] = -dRDC[r];
  }

  // share the calculated rdc
  if(!serial) comm.Sum(&rdc[0],rdc.size());

  // Printout
  bool printout=false;
  if(pperiod>0) printout = (!(getStep()%pperiod));
  if(printout) {
    // share the calculated rdc
    //if(!serial) comm.Sum(&rdc[0],rdc.size());
    // only the master replica is going to write
    if(comm.Get_rank()==0) {
      char tmp1[21]; sprintf(tmp1, "%ld", getStep()); 
      string dcfile = string("rdc-")+getLabel()+"-"+tmp1+string(".dat");
      FILE *outfile = fopen(dcfile.c_str(), "w");
      fprintf(outfile, "#index calc exp\n");
      double sum=0.;
      double Q=0.;
      for(unsigned r=0;r<coupl.size();r++) { 
        fprintf(outfile," %4u %10.6f %10.6f\n", r, rdc[r], coupl[r]);
        sum+=(rdc[r]-coupl[r])*(rdc[r]-coupl[r]);
        Q+=(coupl[r]*coupl[r]);
      }
      fprintf(outfile, "# Q factor = %6.4f\n", sqrt(sum/Q));
      fclose(outfile);
    }
  }

  // Ensemble averaging
  double fact=1.0;
  if(ensemble) {
    fact = 1./((double) ens_dim);
    // share the calculated rdc
    //if(!serial) comm.Sum(&rdc[0],rdc.size());
    // I am the master of my replica
    if(comm.Get_rank()==0) {
      // among replicas
      multi_sim_comm.Sum(&rdc[0], coupl.size() );
      for(unsigned i=0;i<coupl.size();i++) rdc[i] *= fact; 
    } else for(unsigned i=0;i<coupl.size();i++) rdc[i] = 0.;
    // inside each replica
    comm.Sum(&rdc[0], coupl.size() );
    for(unsigned r=rank;r<N;r+=stride) {
      unsigned index=r/2;
      scx  += rdc[index];
      scx2 += rdc[index]*rdc[index];
      scy  += coupl[index];
      scxy += rdc[index]*coupl[index];
    } 
  }
 
  Tensor co_virial, sl_virial, int_virial;
  co_virial.zero();
  sl_virial.zero();
  int_virial.zero();
  vector<Vector> co_deriv(N);
  vector<Vector> sl_deriv(N);
  vector<Vector> int_deriv(N);

  if(!serial) { 
    comm.Sum(scxy);
    comm.Sum(scx);
    comm.Sum(scy);
    comm.Sum(scx2);
  }
  double ns = coupl.size();

  double num = ns*scxy - scx*scy;
  double scy2=norm;
  double idev2x = 1./(ns*scx2-scx*scx);
  double idevx = sqrt(idev2x);
  double idevy = 1./sqrt(ns*scy2-scy*scy);
  /* correlation */
  double correlation = num * idevx * idevy;
  /* slope and intercept */
  double slope = num * idev2x;
  double inter = (scy - slope * scx)/ns;

  /* derivatives */
  for(unsigned r=rank;r<N;r+=stride) {
    unsigned index=r/2;
    double common_d1 = (ns*coupl[index]-scy)*idevx;
    double common_d2 = num*(ns*rdc[index]-scx)*idev2x*idevx;
    double common_d3 = common_d1 - common_d2;
    /* correlation */
    double co_der = fact*common_d3*idevy;
    co_deriv[r]   = co_der*dRDC[r];
    co_deriv[r+1] = co_der*dRDC[r+1];
    co_virial=co_virial+(-1.*Tensor(getPosition(r),co_deriv[r]));
    co_virial=co_virial+(-1.*Tensor(getPosition(r+1),co_deriv[r+1]));
    /* end */

    /* slope */
    double common_sl = (common_d1-2.*common_d2)*idevx;
    double sl_der = fact*common_sl;
    sl_deriv[r]   = sl_der*dRDC[r];
    sl_deriv[r+1] = sl_der*dRDC[r+1];
    sl_virial=sl_virial+(-1.*Tensor(getPosition(r),sl_deriv[r]));
    sl_virial=sl_virial+(-1.*Tensor(getPosition(r+1),sl_deriv[r+1]));
    /* end */

    /* intercept */
    double int_der = -fact*(slope+ scx*common_sl)/ns;
    int_deriv[r]   = int_der*dRDC[r];
    int_deriv[r+1] = int_der*dRDC[r+1];
    int_virial=int_virial+(-1.*Tensor(getPosition(r),int_deriv[r]));
    int_virial=int_virial+(-1.*Tensor(getPosition(r+1),int_deriv[r+1]));
    /* end */
  }

  if(!serial){
    if(!co_deriv.empty())  comm.Sum(&co_deriv[0][0], 3*co_deriv.size());
    if(!sl_deriv.empty())  comm.Sum(&sl_deriv[0][0], 3*sl_deriv.size());
    if(!int_deriv.empty()) comm.Sum(&int_deriv[0][0],3*int_deriv.size());
    comm.Sum(co_virial);
    comm.Sum(sl_virial);
    comm.Sum(int_virial);
  }

  Value* valuea=getPntrToComponent("corr");
  Value* valueb=getPntrToComponent("slope");
  Value* valuec=getPntrToComponent("inter");

  for(unsigned i=0;i<N;i++) {
    setAtomsDerivatives(valuea,i,co_deriv[i]);
    setAtomsDerivatives(valueb,i,sl_deriv[i]);
    setAtomsDerivatives(valuec,i,int_deriv[i]);
  }
  valuea->set(correlation);
  valueb->set(slope);
  valuec->set(inter);
  setBoxDerivatives(valuea, co_virial);
  setBoxDerivatives(valueb, sl_virial);
  setBoxDerivatives(valuec, int_virial);
}

void RDC::svd_calc()
{
  double score=0.;
  vector<double> rdc( coupl.size() );
  unsigned N = getNumberOfAtoms();

#ifdef __PLUMED_HAS_GSL
  gsl_vector *rdc_vec, *S, *Stmp, *work, *bc;
  gsl_matrix *coef_mat, *A, *V;
  double Sxx,Syy,Szz,Sxy,Sxz,Syz;
  rdc_vec = gsl_vector_alloc(coupl.size());
  bc = gsl_vector_alloc(coupl.size());
  Stmp = gsl_vector_alloc(5);
  S = gsl_vector_alloc(5);
  work = gsl_vector_alloc(5);
  coef_mat = gsl_matrix_alloc(coupl.size(),5);
  A = gsl_matrix_alloc(coupl.size(),5);
  V = gsl_matrix_alloc(5,5);
  gsl_matrix_set_zero(coef_mat);
  gsl_vector_set_zero(bc);
  unsigned index=0;
  vector<double> dmax(coupl.size());
  for(unsigned r=0; r<N; r+=2) {
    Vector distance;
    distance=pbcDistance(getPosition(r),getPosition(r+1));
    double d    = distance.modulo();
    double d2   = d*d;
    double d3   = d2*d;
    double id3  = 1./d3; 
    double max  = -Const*mu_s[index]*scale[index];
    dmax[index] = id3*max;
    double mu_x = distance[0]/d;
    double mu_y = distance[1]/d;
    double mu_z = distance[2]/d;
    gsl_vector_set(rdc_vec,index,coupl[index]/dmax[index]);
    gsl_matrix_set(coef_mat,index,0,gsl_matrix_get(coef_mat,index,0)+(mu_x*mu_x-mu_z*mu_z));
    gsl_matrix_set(coef_mat,index,1,gsl_matrix_get(coef_mat,index,1)+(mu_y*mu_y-mu_z*mu_z));
    gsl_matrix_set(coef_mat,index,2,gsl_matrix_get(coef_mat,index,2)+(2.0*mu_x*mu_y));
    gsl_matrix_set(coef_mat,index,3,gsl_matrix_get(coef_mat,index,3)+(2.0*mu_x*mu_z));
    gsl_matrix_set(coef_mat,index,4,gsl_matrix_get(coef_mat,index,4)+(2.0*mu_y*mu_z));
    index++;
  }
  gsl_matrix_memcpy(A,coef_mat);
  gsl_linalg_SV_decomp(A, V, Stmp, work);
  gsl_linalg_SV_solve(A, V, Stmp, rdc_vec, S);
  Sxx = gsl_vector_get(S,0);
  Syy = gsl_vector_get(S,1);
  Szz = -Sxx-Syy;
  Sxy = gsl_vector_get(S,2);
  Sxz = gsl_vector_get(S,3);
  Syz = gsl_vector_get(S,4);
  gsl_blas_dgemv(CblasNoTrans, 1.0, coef_mat, S, 0., bc);
  for(index=0; index<coupl.size(); index++) {
    rdc[index] = gsl_vector_get(bc,index)*dmax[index];
    double delta = rdc[index]-coupl[index];
    score += delta*delta;
  }
  gsl_matrix_free(coef_mat);
  gsl_matrix_free(A);
  gsl_vector_free(rdc_vec);
  gsl_vector_free(bc);
  gsl_vector_free(Stmp);
  gsl_vector_free(S);
  gsl_vector_free(work);
#endif

  // Printout
  bool printout=false;
  if(pperiod>0) printout = (!(getStep()%pperiod));
  if(printout) {
    // only the master replica is going to write
    if(comm.Get_rank()==0) {
      char tmp1[21]; sprintf(tmp1, "%ld", getStep()); 
      string dcfile = string("rdc-")+getLabel()+"-"+tmp1+string(".dat");
      FILE *outfile = fopen(dcfile.c_str(), "w");
      fprintf(outfile, "#index calc exp\n");
      double sum=0.;
      double Q=0.;
      for(unsigned r=0;r<coupl.size();r++) { 
        fprintf(outfile," %4u %10.6f %10.6f\n", r, rdc[r], coupl[r]);
        sum+=(rdc[r]-coupl[r])*(rdc[r]-coupl[r]);
        Q+=(coupl[r]*coupl[r]);
      }
      fprintf(outfile, "# Q factor = %6.4f\n", sqrt(sum/Q));
      fclose(outfile);
    }
  }

  // Ensemble averaging
  double fact=1.0;
  if(ensemble) {
    fact = 1./((double) ens_dim);
    // I am the master of my replica
    if(comm.Get_rank()==0) {
      // among replicas
      multi_sim_comm.Sum(&rdc[0], coupl.size() );
      for(unsigned i=0;i<coupl.size();i++) rdc[i] *= fact; 
    } else for(unsigned i=0;i<coupl.size();i++) rdc[i] = 0.;
    // inside each replica
    comm.Sum(&rdc[0], coupl.size() );
    score = 0.;
    unsigned index=0;
    for(unsigned r=0;r<N;r+=2) {
      double delta = rdc[index]-coupl[index];
      score += delta*delta;
      index++;
    } 
  }
 
  score /= norm;
  setValue           (score);
}

}
