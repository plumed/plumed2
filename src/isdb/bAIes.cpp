/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2024 The plumed team
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

#include "colvar/Colvar.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "tools/Communicator.h"
#include "tools/Matrix.h"
#include "core/GenericMolInfo.h"
#include "core/ActionSet.h"
#include "tools/File.h"
#include "tools/OpenMP.h"
#include <string>
#include <cmath>
#include <map>
#include <ctime>
#include "tools/Random.h"

#ifndef M_PI
//why do not use PLMD::pi?
#define M_PI           3.14159265358979323846
#endif

namespace PLMD {
namespace isdb {

//+PLUMEDOC ISDB_COLVAR BAIES
/*
Bayesian refinement of AF models.

This action implements the Bayesian approach to refine AF models introduced here and here.
It can be used to generate conformational ensembles of IDPs or refine AF models prior to small-molecule virtual screening.

## Examples

Complete tutorials can be found <a href="https://github.com/COSBlab/bAIes-IDP">here</a> and here.

*/
//+ENDPLUMEDOC

class BAIES : public Colvar {

private:

// bool pbc
  bool pbc_;
// temperature in kbt
  double kbt_;
// number of threads
  unsigned nt_;
// positions
  std::vector<Vector> pos_;
// derivatives
  std::vector<Vector> atom_der_;
// AF2 parameters variables
  std::string mtype_;
  std::vector < std::pair<unsigned, unsigned> > atom_pairs_;
  std::vector<double> mus_;
  std::vector<double> sigmas_;
// constants
  double inv_sqrt2_, sqrt2_pi_, inv_pi2_;
// prior type
  unsigned prior_;
  enum { NONE, JEFFREYS, CAUCHY };
// private functions
  void read_data_file(std::string datafile, std::vector<AtomNumber> atoms, double smin);
// energie
  double getGauss();
  double getGaussJeffreys();
  double getGaussCauchy();
  double getLognorm();
  double getLognormJeffreys();
  double getLognormCauchy();

public:
  static void registerKeywords( Keywords& keys );
  explicit BAIES(const ActionOptions&);
// active methods:
  void calculate() override;
};

PLUMED_REGISTER_ACTION(BAIES,"BAIES")

void BAIES::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords( keys );
  keys.add("atoms","ATOMS","atoms used in the calculation of bAIes energy");
  keys.add("compulsory","DATA_FILE","file with AF2 fit parameters");
  keys.add("compulsory","PRIOR", "type of prior to use (NONE, JEFFREYS, CAUCHY");
  keys.add("optional","TEMP", "temperature in kBt units");
  keys.add("optional","SIGMA_MIN", "minimum value of sigma");
  keys.addOutputComponent("ene","default","scalar","Bayesian bAIes energy");
}

BAIES::BAIES(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc_(true) {
  // set constants
  inv_sqrt2_ = 1.0/sqrt(2.0);
  sqrt2_pi_  = sqrt(2.0 / M_PI);
  inv_pi2_   = 0.5 / M_PI / M_PI;

  // list of atoms
  std::vector<AtomNumber> atoms;
  parseAtomList("ATOMS", atoms);

  // file with AF2 fit parameters
  std::string datafile;
  parse("DATA_FILE", datafile);

  // prior type
  std::string prior;
  parse("PRIOR", prior);
  // check priors allowed
  if(prior=="NONE") {
    prior_ = NONE;
  } else if(prior=="JEFFREYS") {
    prior_ = JEFFREYS;
  } else if(prior=="CAUCHY") {
    prior_ = CAUCHY;
  } else {
    error("Unknown PRIOR type!");
  }

  bool nopbc=!pbc_;
  parseFlag("NOPBC",nopbc);
  pbc_=!nopbc;

  // set kbt
  //kbt_ = plumed.getAtoms().getKbT();
  kbt_ = getkBT();
  parse("TEMP", kbt_);

  // set sigma_min
  double sigma_min = -1.0;
  parse("SIGMA_MIN",sigma_min);

  // check read
  checkRead();

  // set parallel stuff
  unsigned mpisize=comm.Get_size();
  if(mpisize>1) {
    error("BAIES supports only OpenMP parallelization");
  }
  // set number of OpenMP threads
  nt_ = OpenMP::getNumThreads();

  // read the data file
  read_data_file(datafile, atoms, sigma_min);

  // print stuff to log
  log.printf("  number of atoms involved : %u\n", atoms.size());
  log.printf("  AF2 fit parameters file : %s\n", datafile.c_str());
  log.printf("  prior type : %s\n", prior.c_str());
  // print stuff from data file
  log.printf("  noise model type : %s\n", mtype_.c_str());
  log.printf("  number of pairs involved : %u\n", atom_pairs_.size());
  // other info
  log.printf("  temperature of the system in energy unit : %f\n", kbt_);
  log.printf("  minimum value of sigma : %f\n", sigma_min);
  if(pbc_) {
    log.printf("  using periodic boundary conditions\n");
  } else {
    log.printf("  without periodic boundary conditions\n");
  }

  // prepare other vectors: data and derivatives
  atom_der_.resize(atoms.size());

  // add components
  addComponentWithDerivatives("ene");
  componentIsNotPeriodic("ene");

  // request atoms
  requestAtoms(atoms);

  // print bibliography
  log<<"  Bibliography "<<plumed.cite("Bonomi, Camilloni, Bioinformatics, 33, 3999 (2017)");
  log<<plumed.cite("Schnapka, Morozova, Sen, Bonomi, bioRxiv (2025) doi: XXX");
  log<<"\n";
}

void BAIES::read_data_file(const std::string datafile, const std::vector<AtomNumber> atoms, double smin) {
  unsigned id, ai, aj;
  double mu, sigma;
// map serials to index in position array
  std::map< unsigned, unsigned > index;
  for(unsigned i=0; i<atoms.size(); ++i) {
    index[atoms[i].serial()]=i;
  }
// open file
  IFile *ifile = new IFile();
  if(ifile->FileExist(datafile)) {
    ifile->open(datafile);
    // read constant fields
    if( ifile->FieldExist("model") ) {
      ifile->scanField("model",mtype_);
      if( mtype_!="gaussian" && mtype_!="lognormal" ) {
        error("Unknown noise model type");
      }
    } else {
      error("Missing noise model type in DATA_FILE");
    }
    // read line-by-line
    while(ifile->scanField("Id",id)) {
      ifile->scanField("atom_i",ai);
      ifile->scanField("atom_j",aj);
      ifile->scanField("mu",mu);
      ifile->scanField("sigma",sigma);
      // list of pairs of atoms
      atom_pairs_.push_back(std::make_pair(index[ai],index[aj]));
      // list of mus
      mus_.push_back(mu);
      // list of sigmas
      sigmas_.push_back(std::max(smin,sigma));
      // new line
      ifile->scanField();
    }
    ifile->close();
  } else {
    error("Cannot find DATA_FILE "+datafile+"\n");
  }
  delete ifile;
}

// Energy models
double BAIES::getGauss() {
  double ene = 0.0;
  // loop over atoms pairs
  #pragma omp parallel num_threads(nt_) shared(ene)
  {
    std::vector<Vector> omp_der(atom_der_.size());
    #pragma omp for reduction( + : ene) nowait
    for(unsigned i=0; i<atom_pairs_.size(); ++i) {
      // get indexes
      unsigned ii = atom_pairs_[i].first;
      unsigned jj = atom_pairs_[i].second;
      // get distance
      Vector distance=delta(pos_[ii], pos_[jj]);
      const double value=distance.modulo();
      const double invvalue=1.0/value;
      // get energy
      double arg = ( value - mus_[i] ) / sigmas_[i];
      ene += arg * arg;
      // calculate derivatives
      double arg_der = kbt_ * arg / sigmas_[i];
      omp_der[ii] += - arg_der * invvalue * distance;
      omp_der[jj] +=   arg_der * invvalue * distance;
    }
    // gather from tasks
    #pragma omp critical
    for(unsigned i=0; i<atom_der_.size(); ++i) {
      atom_der_[i]+=omp_der[i];
    }
  }
  // missing Gaussian normalization constant
  return 0.5 * kbt_ * ene;
}

double BAIES::getGaussJeffreys() {
  double ene = 0.0;
  // loop over all atom pairs
  #pragma omp parallel num_threads(nt_) shared(ene)
  {
    std::vector<Vector> omp_der(atom_der_.size());
    #pragma omp for reduction( + : ene) nowait
    for(unsigned i=0; i<atom_pairs_.size(); ++i) {
      // get indexes
      unsigned ii = atom_pairs_[i].first;
      unsigned jj = atom_pairs_[i].second;
      // get distance
      Vector distance=delta(pos_[ii], pos_[jj]);
      const double value=distance.modulo();
      const double invvalue=1.0/value;
      // get energy
      double val_mu = value - mus_[i];
      double inv_val_mu = 1.0 / val_mu;
      double inv_sigma = 1.0 / sigmas_[i];
      double erfcontent = val_mu * inv_sigma * inv_sqrt2_;
      double erfval = std::erf(erfcontent);
      // increase energy
      ene += -std::log( erfval * inv_val_mu );
      // calculate derivatives
      double arg_der = - kbt_ * val_mu * ( (sqrt2_pi_ * std::exp(-erfcontent * erfcontent) * inv_val_mu * inv_sigma ) - ( erfval * inv_val_mu * inv_val_mu ) ) / erfval ;
      // increase derivatives
      omp_der[ii] += - arg_der * invvalue * distance;
      omp_der[jj] +=   arg_der * invvalue * distance;
    }
    // gather from tasks
    #pragma omp critical
    for(unsigned i=0; i<atom_der_.size(); ++i) {
      atom_der_[i]+=omp_der[i];
    }
  }
  // missing normalization constants
  return kbt_ * ene;
}

double BAIES::getGaussCauchy() {
  double ene = 0.0;
  // loop over all atom pairs
  #pragma omp parallel num_threads(nt_) shared(ene)
  {
    std::vector<Vector> omp_der(atom_der_.size());
    #pragma omp for reduction( + : ene) nowait
    for(unsigned i=0; i<atom_pairs_.size(); ++i) {
      // get indexes
      unsigned ii = atom_pairs_[i].first;
      unsigned jj = atom_pairs_[i].second;
      // get distance
      Vector distance=delta(pos_[ii], pos_[jj]);
      const double value=distance.modulo();
      const double invvalue=1.0/value;
      // get energy
      double val_mu = value - mus_[i];
      double lncontent = 0.5*val_mu*val_mu + sigmas_[i]*sigmas_[i];
      ene += std::log( lncontent );
      // calculate derivatives
      double arg_der =  kbt_ * val_mu / lncontent;
      omp_der[ii] += - arg_der * invvalue * distance;
      omp_der[jj] +=   arg_der * invvalue * distance;
    }
    // gather from tasks
    #pragma omp critical
    for(unsigned i=0; i<atom_der_.size(); ++i) {
      atom_der_[i]+=omp_der[i];
    }
  }
  return kbt_ * ene;
}

double BAIES::getLognorm() {
  double ene = 0.0;
  // loop over atoms pairs
  #pragma omp parallel num_threads(nt_) shared(ene)
  {
    std::vector<Vector> omp_der(atom_der_.size());
    #pragma omp for reduction( + : ene) nowait
    for(unsigned i=0; i<atom_pairs_.size(); ++i) {
      // get indexes
      unsigned ii = atom_pairs_[i].first;
      unsigned jj = atom_pairs_[i].second;
      // get distance
      Vector distance=delta(pos_[ii], pos_[jj]);
      const double value=distance.modulo();
      const double invvalue=1.0/value;
      // get energy
      double arg = ( std::log(value) - mus_[i] ) / sigmas_[i];
      ene += arg * arg;
      // calculate derivatives
      double arg_der = kbt_ * arg / (value * sigmas_[i]);
      omp_der[ii] += - arg_der * invvalue * distance;
      omp_der[jj] +=   arg_der * invvalue * distance;
    }
    // gather from tasks
    #pragma omp critical
    for(unsigned i=0; i<atom_der_.size(); ++i) {
      atom_der_[i]+=omp_der[i];
    }
  }
  return 0.5 * kbt_ * ene;
}

double BAIES::getLognormJeffreys() {
  double ene = 0.0;
  // loop over all atom pairs
  #pragma omp parallel num_threads(nt_) shared(ene)
  {
    std::vector<Vector> omp_der(atom_der_.size());
    #pragma omp for reduction( + : ene) nowait
    for(unsigned i=0; i<atom_pairs_.size(); ++i) {
      // get indexes
      unsigned ii = atom_pairs_[i].first;
      unsigned jj = atom_pairs_[i].second;
      // get distance
      Vector distance=delta(pos_[ii], pos_[jj]);
      const double value=distance.modulo();
      const double invvalue=1.0/value;
      // get energy
      double val_mu = std::log(value) - mus_[i];
      double inv_val_mu = 1 / val_mu;
      double inv_sigma = 1 / sigmas_[i];
      double erfcontent = val_mu * inv_sigma * inv_sqrt2_;
      double erfval = std::erf(erfcontent);
      // increase energy
      ene += -std::log( erfval * inv_val_mu );
      // calculate derivatives
      double arg_der = - kbt_ * val_mu * ( (sqrt2_pi_ * std::exp(-erfcontent * erfcontent) * inv_val_mu * inv_sigma ) - ( erfval * inv_val_mu * inv_val_mu ) ) / erfval ;
      // chain rule
      arg_der *= invvalue;
      // increase derivatives
      omp_der[ii] += - arg_der * invvalue * distance;
      omp_der[jj] +=   arg_der * invvalue * distance;
    }
    // gather from tasks
    #pragma omp critical
    for(unsigned i=0; i<atom_der_.size(); ++i) {
      atom_der_[i]+=omp_der[i];
    }
  }
  // missing constants
  return kbt_ * ene;
}

double BAIES::getLognormCauchy() {
  double ene = 0.0;
  // loop over all atom pairs
  #pragma omp parallel num_threads(nt_) shared(ene)
  {
    std::vector<Vector> omp_der(atom_der_.size());
    #pragma omp for reduction( + : ene) nowait
    for(unsigned i=0; i<atom_pairs_.size(); ++i) {
      // get indexes
      unsigned ii = atom_pairs_[i].first;
      unsigned jj = atom_pairs_[i].second;
      // get distance
      Vector distance=delta(pos_[ii], pos_[jj]);
      const double value=distance.modulo();
      const double invvalue=1.0/value;
      // get energy
      double val_mu = std::log(value) - mus_[i];
      double lncontent = 0.5*val_mu*val_mu + sigmas_[i]*sigmas_[i];
      ene += std::log( lncontent );
      // calculate derivatives
      double arg_der =  kbt_ * val_mu * invvalue / lncontent;
      omp_der[ii] += - arg_der * invvalue * distance;
      omp_der[jj] +=   arg_der * invvalue * distance;
    }
    // gather from tasks
    #pragma omp critical
    for(unsigned i=0; i<atom_der_.size(); ++i) {
      atom_der_[i]+=omp_der[i];
    }
  }
  // missing constants
  return kbt_ * ene;
}

void BAIES::calculate() {
  // fix PBCs
  if(pbc_) {
    makeWhole();
  }

  // get the positions at this step
  pos_ = getPositions();

  // reset derivatives
  #pragma omp parallel for num_threads(nt_)
  for(unsigned i=0; i<atom_der_.size(); ++i) {
    atom_der_[i] = Vector(0,0,0);
  }

  // calculate bAIes energy and derivatives
  double ene = 0.0;
  // Gaussian noise model
  if(mtype_=="gaussian") {
    switch(prior_) {
    case NONE:
      ene = getGauss();
      break;
    case JEFFREYS:
      ene = getGaussJeffreys();
      break;
    case CAUCHY:
      ene = getGaussCauchy();
      break;
    }
  } else {
    switch(prior_) {
    case NONE:
      ene = getLognorm();
      break;
    case JEFFREYS:
      ene = getLognormJeffreys();
      break;
    case CAUCHY:
      ene = getLognormCauchy();
      break;
    }
  }
  // set score
  Value* score = getPntrToComponent("ene");
  score->set(ene);
  // calculate virial
  Tensor virial;
  // declare omp reduction for Tensors
  #pragma omp declare reduction( sumTensor : Tensor : omp_out += omp_in )

  #pragma omp parallel for num_threads(nt_) reduction (sumTensor : virial)
  for(unsigned i=0; i<pos_.size(); ++i) {
    virial += Tensor(pos_[i], -atom_der_[i]);
  }
  // set virial
  setBoxDerivatives(score, virial);
  // set derivatives
  #pragma omp parallel for num_threads(nt_)
  for(unsigned i=0; i<atom_der_.size(); ++i) {
    setAtomsDerivatives(score, i, atom_der_[i]);
  }
}

}
}
