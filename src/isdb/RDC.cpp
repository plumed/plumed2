/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2019 The plumed team
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
#include "MetainferenceBase.h"
#include "core/ActionRegister.h"
#include "tools/Pbc.h"

#ifdef __PLUMED_HAS_GSL
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#endif

using namespace std;

namespace PLMD {
namespace isdb {

//+PLUMEDOC ISDB_COLVAR RDC
/*
Calculates the (Residual) Dipolar Coupling between two atoms.

The Dipolar Coupling between two nuclei depends on the \f$\theta\f$ angle between
the inter-nuclear vector and the external magnetic field.

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
and their products (this is what is given in input using the keyword GYROM)
- N-H -72.5388
- C-H 179.9319
- C-N -18.2385
- C-C 45.2404

In isotropic media DCs average to zero because of the rotational
averaging, but when the rotational symmetry is broken, either through the introduction of an alignment medium or for molecules
with highly anisotropic paramagnetic susceptibility, then the average of the DCs is not zero and it is possible to measure a Residual Dipolar Coupling (RDCs).

This collective variable calculates the Dipolar Coupling for a set of couple of atoms using
the above definition.

In a standard MD simulation the average over time of the DC should then be zero. If one wants to model the meaning of a set of measured RDCs it is possible to try to solve the following problem: "what is the distribution of structures and orientations that reproduce the measured RDCs".

This collective variable can then be use to break the rotational symmetry of a simulation by imposing that the average of the DCs over the conformational ensemble must be equal to the measured RDCs \cite Camilloni:2015ka . Since measured RDCs are also a function of the fraction of aligned molecules in the sample it is better to compare them modulo a constant or looking at the correlation.

Alternatively if the molecule is rigid it is possible to use the experimental data to calculate the alignment tensor and the use that to back calculate the RDCs, this is what is usually call the Single Value Decomposition approach. In this case the code rely on the
a set of function from the GNU Scientific Library (GSL). (With SVD forces are not currently implemented).

Replica-Averaged simulations can be performed using RDCs, \ref ENSEMBLE, \ref STATS and \ref RESTRAINT .
\ref METAINFERENCE can be activated using DOSCORE and the other relevant keywords.

Additional material and examples can be also found in the tutorial \ref belfast-9

\par Examples
In the following example five N-H RDCs are defined and averaged over multiple replicas,
their correlation is then calculated with respect to a set of experimental data and restrained.
In addition, and only for analysis purposes, the same RDCs each single conformation are calculated
using a Single Value Decomposition algorithm, then averaged and again compared with the experimental data.

\plumedfile
RDC ...
GYROM=-72.5388
SCALE=0.001
ATOMS1=20,21
ATOMS2=37,38
ATOMS3=56,57
ATOMS4=76,77
ATOMS5=92,93
LABEL=nh
... RDC

erdc: ENSEMBLE ARG=nh.*

st: STATS ARG=erdc.* PARAMETERS=8.17,-8.271,-10.489,-9.871,-9.152

rdce: RESTRAINT ARG=st.corr KAPPA=0. SLOPE=-25000.0 AT=1.

RDC ...
GYROM=-72.5388
SVD
ATOMS1=20,21 COUPLING1=8.17
ATOMS2=37,38 COUPLING2=-8.271
ATOMS3=56,57 COUPLING3=-10.489
ATOMS4=76,77 COUPLING4=-9.871
ATOMS5=92,93 COUPLING5=-9.152
LABEL=svd
... RDC

esvd: ENSEMBLE ARG=svd.*

st_svd: STATS ARG=esvd.* PARAMETERS=8.17,-8.271,-10.489,-9.871,-9.152


PRINT ARG=st.corr,st_svd.corr,rdce.bias FILE=colvar
\endplumedfile

*/
//+ENDPLUMEDOC

//+PLUMEDOC ISDB_COLVAR PCS
/*
Calculates the Pseudo-contact shift of a nucleus determined by the presence of a metal ion susceptible to anisotropic magnetization.

The PCS of an atomic nucleus depends on the \f$\theta\f$ angle between the vector from the spin-label to the nucleus
 and the external magnetic field and the module of the vector itself \cite Camilloni:2015jf . While in principle the averaging
resulting from the tumbling should remove the pseudo-contact shift, in presence of the NMR magnetic field the magnetically anisotropic molecule bound to system will break the rotational symmetry does resulting in measurable values for the PCS and RDC.

PCS values can also be calculated using a Single Value Decomposition approach, in this case the code rely on the
a set of function from the GNU Scientific Library (GSL). (With SVD forces are not currently implemented).

Replica-Averaged simulations can be performed using PCS values, \ref ENSEMBLE, \ref STATS and \ref RESTRAINT .
Metainference simulations can be performed with this CV and \ref METAINFERENCE .

\par Examples

In the following example five PCS values are defined and their correlation with
respect to a set of experimental data is calculated and restrained. In addition,
and only for analysis purposes, the same PCS values are calculated using a Single Value
Decomposition algorithm.

\plumedfile
PCS ...
ATOMS1=20,21
ATOMS2=20,38
ATOMS3=20,57
ATOMS4=20,77
ATOMS5=20,93
LABEL=nh
... PCS

enh: ENSEMBLE ARG=nh.*

st: STATS ARG=enh.* PARAMETERS=8.17,-8.271,-10.489,-9.871,-9.152

pcse: RESTRAINT ARG=st.corr KAPPA=0. SLOPE=-25000.0 AT=1.

PRINT ARG=st.corr,pcse.bias FILE=colvar
\endplumedfile

*/
//+ENDPLUMEDOC

class RDC :
  public MetainferenceBase
{
private:
  double         Const;
  double         mu_s;
  double         scale;
  vector<double> coupl;
  bool           svd;
  bool           pbc;

#ifdef __PLUMED_HAS_GSL
/// Auxiliary class to delete a gsl_vector.
/// If used somewhere else we can move it.
  struct gsl_vector_deleter {
    void operator()(gsl_vector* p) {
      gsl_vector_free(p);
    }
  };

/// unique_ptr to a gsl_vector.
/// Gets deleted when going out of scope.
  typedef std::unique_ptr<gsl_vector,gsl_vector_deleter> gsl_vector_unique_ptr;

/// Auxiliary class to delete a gsl_matrix.
/// If used somewhere else we can move it.
  struct gsl_matrix_deleter {
    void operator()(gsl_matrix* p) {
      gsl_matrix_free(p);
    }
  };

/// unique_ptr to a gsl_matrix.
/// Gets deleted when going out of scope.
  typedef std::unique_ptr<gsl_matrix,gsl_matrix_deleter> gsl_matrix_unique_ptr;
#endif


  void do_svd();
public:
  explicit RDC(const ActionOptions&);
  static void registerKeywords( Keywords& keys );
  virtual void calculate();
  void update();
};

PLUMED_REGISTER_ACTION(RDC,"RDC")
PLUMED_REGISTER_ACTION(RDC,"PCS")

void RDC::registerKeywords( Keywords& keys ) {
  componentsAreNotOptional(keys);
  useCustomisableComponents(keys);
  MetainferenceBase::registerKeywords(keys);
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  keys.add("numbered","ATOMS","the couple of atoms involved in each of the bonds for which you wish to calculate the RDC. "
           "Keywords like ATOMS1, ATOMS2, ATOMS3,... should be listed and one dipolar coupling will be "
           "calculated for each ATOMS keyword you specify.");
  keys.reset_style("ATOMS","atoms");
  keys.add("compulsory","GYROM","1.","Add the product of the gyromagnetic constants for the bond. ");
  keys.add("compulsory","SCALE","1.","Add the scaling factor to take into account concentration and other effects. ");
  keys.addFlag("SVD",false,"Set to TRUE if you want to back calculate using Single Value Decomposition (need GSL at compilation time).");
  keys.addFlag("ADDCOUPLINGS",false,"Set to TRUE if you want to have fixed components with the experimental values.");
  keys.add("numbered","COUPLING","Add an experimental value for each coupling (needed by SVD and useful for \\ref STATS).");
  keys.addOutputComponent("rdc","default","the calculated # RDC");
  keys.addOutputComponent("exp","SVD/ADDCOUPLINGS","the experimental # RDC");
}

RDC::RDC(const ActionOptions&ao):
  PLUMED_METAINF_INIT(ao),
  Const(1.),
  mu_s(1.),
  scale(1.),
  pbc(true)
{
  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  const double RDCConst = 0.3356806;
  const double PCSConst = 1.0;

  if( getName().find("RDC")!=std::string::npos) { Const *= RDCConst; }
  else if( getName().find("PCS")!=std::string::npos) { Const *= PCSConst; }

  // Read in the atoms
  vector<AtomNumber> t, atoms;
  for(int i=1;; ++i ) {
    parseAtomList("ATOMS", i, t );
    if( t.empty() ) break;
    if( t.size()!=2 ) {
      std::string ss; Tools::convert(i,ss);
      error("ATOMS" + ss + " keyword has the wrong number of atoms");
    }
    atoms.push_back(t[0]);
    atoms.push_back(t[1]);
    t.resize(0);
  }

  const unsigned ndata = atoms.size()/2;

  // Read in GYROMAGNETIC constants
  parse("GYROM", mu_s);
  if(mu_s==0.) error("GYROM cannot be 0");

  // Read in SCALING factors
  parse("SCALE", scale);
  if(scale==0.) error("SCALE cannot be 0");

  svd=false;
  parseFlag("SVD",svd);
#ifndef __PLUMED_HAS_GSL
  if(svd) error("You CANNOT use SVD without GSL. Recompile PLUMED with GSL!\n");
#endif
  if(svd&&getDoScore()) error("It is not possible to use SVD and METAINFERENCE together");

  bool addexp=false;
  parseFlag("ADDCOUPLINGS",addexp);
  if(getDoScore()||svd) addexp=true;

  if(addexp) {
    coupl.resize( ndata );
    unsigned ntarget=0;
    for(unsigned i=0; i<ndata; ++i) {
      if( !parseNumbered( "COUPLING", i+1, coupl[i] ) ) break;
      ntarget++;
    }
    if( ntarget!=ndata ) error("found wrong number of COUPLING values");
  }

  // Ouput details of all contacts
  log.printf("  Gyromagnetic moment is %f. Scaling factor is %f.",mu_s,scale);
  for(unsigned i=0; i<ndata; ++i) {
    log.printf("  The %uth Bond Dipolar Coupling is calculated from atoms : %d %d.", i+1, atoms[2*i].serial(), atoms[2*i+1].serial());
    if(addexp) log.printf(" Experimental coupling is %f.", coupl[i]);
    log.printf("\n");
  }

  log<<"  Bibliography "
     <<plumed.cite("Camilloni C, Vendruscolo M, J. Phys. Chem. B 119, 653 (2015)")
     <<plumed.cite("Camilloni C, Vendruscolo M, Biochemistry 54, 7470 (2015)");
  log<< plumed.cite("Bonomi, Camilloni, Bioinformatics, 33, 3999 (2017)") << "\n";


  if(!getDoScore()&&!svd) {
    for(unsigned i=0; i<ndata; i++) {
      std::string num; Tools::convert(i,num);
      addComponentWithDerivatives("rdc_"+num);
      componentIsNotPeriodic("rdc_"+num);
    }
    if(addexp) {
      for(unsigned i=0; i<ndata; i++) {
        std::string num; Tools::convert(i,num);
        addComponent("exp_"+num);
        componentIsNotPeriodic("exp_"+num);
        Value* comp=getPntrToComponent("exp_"+num);
        comp->set(coupl[i]);
      }
    }
  } else {
    for(unsigned i=0; i<ndata; i++) {
      std::string num; Tools::convert(i,num);
      addComponentWithDerivatives("rdc_"+num);
      componentIsNotPeriodic("rdc_"+num);
    }
    for(unsigned i=0; i<ndata; i++) {
      std::string num; Tools::convert(i,num);
      addComponent("exp_"+num);
      componentIsNotPeriodic("exp_"+num);
      Value* comp=getPntrToComponent("exp_"+num);
      comp->set(coupl[i]);
    }
  }

  if(svd) {
    addComponent("Sxx"); componentIsNotPeriodic("Sxx");
    addComponent("Syy"); componentIsNotPeriodic("Syy");
    addComponent("Szz"); componentIsNotPeriodic("Szz");
    addComponent("Sxy"); componentIsNotPeriodic("Sxy");
    addComponent("Sxz"); componentIsNotPeriodic("Sxz");
    addComponent("Syz"); componentIsNotPeriodic("Syz");
  }

  requestAtoms(atoms, false);
  if(getDoScore()) {
    setParameters(coupl);
    Initialise(coupl.size());
  }
  setDerivatives();
  checkRead();
}

void RDC::do_svd()
{
#ifdef __PLUMED_HAS_GSL
  gsl_vector_unique_ptr rdc_vec(gsl_vector_alloc(coupl.size())),
                        S(gsl_vector_alloc(5)),
                        Stmp(gsl_vector_alloc(5)),
                        work(gsl_vector_alloc(5)),
                        bc(gsl_vector_alloc(coupl.size()));

  gsl_matrix_unique_ptr coef_mat(gsl_matrix_alloc(coupl.size(),5)),
                        A(gsl_matrix_alloc(coupl.size(),5)),
                        V(gsl_matrix_alloc(5,5));

  gsl_matrix_set_zero(coef_mat.get());
  gsl_vector_set_zero(bc.get());

  unsigned index=0;
  vector<double> dmax(coupl.size());
  for(unsigned r=0; r<getNumberOfAtoms(); r+=2) {
    Vector  distance;
    if(pbc) distance = pbcDistance(getPosition(r),getPosition(r+1));
    else    distance = delta(getPosition(r),getPosition(r+1));
    double d    = distance.modulo();
    double d2   = d*d;
    double d3   = d2*d;
    double id3  = 1./d3;
    double max  = -Const*mu_s*scale;
    dmax[index] = id3*max;
    double mu_x = distance[0]/d;
    double mu_y = distance[1]/d;
    double mu_z = distance[2]/d;
    gsl_vector_set(rdc_vec.get(),index,coupl[index]/dmax[index]);
    gsl_matrix_set(coef_mat.get(),index,0,gsl_matrix_get(coef_mat.get(),index,0)+(mu_x*mu_x-mu_z*mu_z));
    gsl_matrix_set(coef_mat.get(),index,1,gsl_matrix_get(coef_mat.get(),index,1)+(mu_y*mu_y-mu_z*mu_z));
    gsl_matrix_set(coef_mat.get(),index,2,gsl_matrix_get(coef_mat.get(),index,2)+(2.0*mu_x*mu_y));
    gsl_matrix_set(coef_mat.get(),index,3,gsl_matrix_get(coef_mat.get(),index,3)+(2.0*mu_x*mu_z));
    gsl_matrix_set(coef_mat.get(),index,4,gsl_matrix_get(coef_mat.get(),index,4)+(2.0*mu_y*mu_z));
    index++;
  }
  gsl_matrix_memcpy(A.get(),coef_mat.get());
  gsl_linalg_SV_decomp(A.get(), V.get(), Stmp.get(), work.get());
  gsl_linalg_SV_solve(A.get(), V.get(), Stmp.get(), rdc_vec.get(), S.get());
  /* tensor */
  Value* tensor;
  tensor=getPntrToComponent("Sxx");
  double Sxx = gsl_vector_get(S.get(),0);
  tensor->set(Sxx);
  tensor=getPntrToComponent("Syy");
  double Syy = gsl_vector_get(S.get(),1);
  tensor->set(Syy);
  tensor=getPntrToComponent("Szz");
  double Szz = -Sxx-Syy;
  tensor->set(Szz);
  tensor=getPntrToComponent("Sxy");
  double Sxy = gsl_vector_get(S.get(),2);
  tensor->set(Sxy);
  tensor=getPntrToComponent("Sxz");
  double Sxz = gsl_vector_get(S.get(),3);
  tensor->set(Sxz);
  tensor=getPntrToComponent("Syz");
  double Syz = gsl_vector_get(S.get(),4);
  tensor->set(Syz);

  gsl_blas_dgemv(CblasNoTrans, 1.0, coef_mat.get(), S.get(), 0., bc.get());
  for(index=0; index<coupl.size(); index++) {
    double rdc = gsl_vector_get(bc.get(),index)*dmax[index];
    Value* val=getPntrToComponent(index);
    val->set(rdc);
  }
#endif
}

void RDC::calculate()
{
  if(svd) {
    do_svd();
    return;
  }

  const double max  = -Const*scale*mu_s;
  const unsigned N=getNumberOfAtoms();
  vector<Vector> dRDC(N/2, Vector{0.,0.,0.});

  /* RDC Calculations and forces */
  #pragma omp parallel num_threads(OpenMP::getNumThreads())
  {
    #pragma omp for
    for(unsigned r=0; r<N; r+=2)
    {
      const unsigned index=r/2;
      Vector  distance;
      if(pbc) distance   = pbcDistance(getPosition(r),getPosition(r+1));
      else    distance   = delta(getPosition(r),getPosition(r+1));
      const double d2    = distance.modulo2();
      const double ind   = 1./sqrt(d2);
      const double ind2  = 1./d2;
      const double ind3  = ind2*ind;
      const double x2    = distance[0]*distance[0]*ind2;
      const double y2    = distance[1]*distance[1]*ind2;
      const double z2    = distance[2]*distance[2]*ind2;
      const double dmax  = ind3*max;
      const double ddmax = dmax*ind2;

      const double rdc   = 0.5*dmax*(3.*z2-1.);
      const double prod_xy = (x2+y2-4.*z2);
      const double prod_z =  (3.*x2 + 3.*y2 - 2.*z2);

      dRDC[index] = -1.5*ddmax*distance;
      dRDC[index][0] *= prod_xy;
      dRDC[index][1] *= prod_xy;
      dRDC[index][2] *= prod_z;

      string num; Tools::convert(index,num);
      Value* val=getPntrToComponent("rdc_"+num);
      val->set(rdc);
      if(!getDoScore()) {
        setBoxDerivatives(val, Tensor(distance,dRDC[index]));
        setAtomsDerivatives(val, r,  dRDC[index]);
        setAtomsDerivatives(val, r+1, -dRDC[index]);
      } else setCalcData(index, rdc);
    }
  }

  if(getDoScore()) {
    /* Metainference */
    Tensor dervir;
    double score = getScore();
    setScore(score);

    /* calculate final derivatives */
    Value* val=getPntrToComponent("score");
    for(unsigned r=0; r<N; r+=2)
    {
      const unsigned index=r/2;
      Vector  distance;
      if(pbc) distance   = pbcDistance(getPosition(r),getPosition(r+1));
      else    distance   = delta(getPosition(r),getPosition(r+1));
      const Vector der = dRDC[index]*getMetaDer(index);
      dervir += Tensor(distance, der);
      setAtomsDerivatives(val, r,  der);
      setAtomsDerivatives(val, r+1, -der);
    }
    setBoxDerivatives(val, dervir);
  }
}

void RDC::update() {
  // write status file
  if(getWstride()>0&& (getStep()%getWstride()==0 || getCPT()) ) writeStatus();
}

}
}
