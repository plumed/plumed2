/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2026 The plumed team
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
#ifndef __PLUMED_isdb_MetainferenceCore_h
#define __PLUMED_isdb_MetainferenceCore_h

#include "core/Action.h"
#include "core/ActionWithValue.h"
#include "tools/File.h"
#include "tools/Random.h"
#include "tools/OpenMP.h"
#include <functional>
#include <string>
#include <vector>

namespace PLMD {

class Communicator;
class Keywords;
class Log;
class Value;

namespace isdb {

/// History of the estimates of the squared standard error of the mean of one datum:
/// a circular buffer of fixed capacity that keeps track of its maximum
class SigmaMeanHistory {
  std::vector<double> buf_;
  unsigned first_=0;
  unsigned n_=0;
  double max_=0.;
  void recomputeMax() {
    max_=get(0);
    for(unsigned k=1; k<n_; ++k) {
      if(max_<get(k)) {
        max_=get(k);
      }
    }
  }
public:
  explicit SigmaMeanHistory(const unsigned capacity=1): buf_(capacity>0?capacity:1) {}
  unsigned size() const {
    return n_;
  }
  /// the maximum of the stored values (the history must not be empty)
  double max() const {
    return max_;
  }
  /// the k-th value, from the oldest
  double get(const unsigned k) const {
    return buf_[(first_+k)%buf_.size()];
  }
  void clear() {
    first_=0;
    n_=0;
  }
  /// add a value, removing the oldest one if the history is full
  void push(const double v) {
    if(n_==buf_.size()) {
      const double old=buf_[first_];
      buf_[first_]=v;
      first_=(first_+1)%buf_.size();
      if(max_<v) {
        max_=v;
      } else if(old==max_) {
        recomputeMax();
      }
    } else {
      buf_[(first_+n_)%buf_.size()]=v;
      n_++;
      if(n_==1||max_<v) {
        max_=v;
      }
    }
  }
  /// replace the oldest value
  void setFirst(const double v) {
    buf_[first_]=v;
    recomputeMax();
  }
};

/**
The Metainference algorithm shared by the METAINFERENCE bias and by the
ISDB collective variables that inherit from MetainferenceBase.

This is not an Action: it is owned by an Action, which it uses for input
parsing, communication, log, restart and status files. The owner adds its
own value (the bias or the score) and calls setup() once the number of data
and the experimental values are known; setup() adds the components of
Metainference to the owner. At every step the owner fills data() with the
back-calculated values, calls calculate() and applies the derivatives of the
energy with respect to the data (getMetaDer()) and to the bias used for
reweighting (getBiasDer()).
*/
class MetainferenceCore {
public:
  /// get a pointer to a component of the owner action, given its name
  typedef std::function<Value*(const std::string&)> ComponentGetter;
private:
  Action& action_;
  Communicator& comm;
  Communicator& multi_sim_comm;
  Log& log;

  // true if Metainference is active (always for the bias, DOSCORE for the CVs)
  bool active_;
  unsigned write_stride_;
  // number of experimental data
  unsigned narg;
  // experimental data
  std::vector<double> parameters;
  // metainference derivatives
  std::vector<double> metader_;
  // vector of back-calculated experimental data
  std::vector<double> calc_data_;

  // noise type
  unsigned noise_type_;
  enum { GAUSS, MGAUSS, OUTLIERS, MOUTLIERS, GENERIC };
  unsigned gen_likelihood_;
  enum { LIKE_GAUSS, LIKE_LOGN };
  bool   doscale_;
  unsigned scale_prior_;
  enum { SC_GAUSS, SC_FLAT };
  double scale_;
  double scale_mu_;
  double scale_min_;
  double scale_max_;
  double Dscale_;
  // scale is data scaling factor
  // noise type
  unsigned offset_prior_;
  bool   dooffset_;
  double offset_;
  double offset_mu_;
  double offset_min_;
  double offset_max_;
  double Doffset_;
  // scale and offset regression
  bool doregres_zero_;
  int  nregres_zero_;
  // sigma is data uncertainty
  std::vector<double> sigma_;
  std::vector<double> sigma_min_;
  std::vector<double> sigma_max_;
  std::vector<double> Dsigma_;
  // sigma_mean is uncertainty in the mean estimate
  std::vector<double> sigma_mean2_;
  // this is the estimator of the mean value per replica for generic metainference
  std::vector<double> ftilde_;
  double Dftilde_;

  // temperature in kbt
  double   kbt_;

  // Monte Carlo stuff
  std::vector<Random> random;
  unsigned MCsteps_;
  long long unsigned MCaccept_;
  long long unsigned MCacceptScale_;
  long long unsigned MCacceptFT_;
  long long unsigned MCtrial_;
  unsigned MCchunksize_;

  // output
  Value*   valueScale=nullptr;
  Value*   valueOffset=nullptr;
  Value*   valueAccept=nullptr;
  Value*   valueAcceptScale=nullptr;
  Value*   valueAcceptFT=nullptr;
  Value*   valueBiasDer_=nullptr;
  Value*   valueWeight_=nullptr;
  Value*   valueNeff_=nullptr;
  std::vector<Value*> valueSigma;
  std::vector<Value*> valueSigmaMean;
  std::vector<Value*> valueFtilde;

  // restart
  std::string status_file_name_;
  OFile    sfile_;
  std::string fmt_;

  // others
  bool     firstTime;
  std::vector<bool> firstTimeW;
  bool     master;
  bool     do_reweight_;
  unsigned do_optsigmamean_;
  unsigned nrep_;
  unsigned replica_;

  // selector
  unsigned nsel_;
  std::string selector_;
  unsigned iselect;

  // optimize sigma mean
  std::vector< std::vector <SigmaMeanHistory> > sigma_mean2_last_;
  unsigned optsigmamean_stride_;
  // optimize sigma max
  unsigned N_optimized_step_;
  unsigned optimized_step_;
  bool sigmamax_opt_done_;
  std::vector<double> sigma_max_est_;
  // GENERIC: sigma mean (the minimum of sigma) has exceeded SIGMA_MAX
  bool warned_sigma_max_=false;

  // average weights
  double decay_w_;
  std::vector< std::vector <double> >  average_weights_;

  // derivative of the energy with respect to the bias used for reweighting
  double biasDer_;
  // value of the bias used for reweighting at this step
  double reweight_bias_;
  // the sigmas moved by the last MC move (MC_CHUNKSIZE)
  std::vector<unsigned> moved_;
  // the energy of the MC has been updated with differences
  bool energy_from_differences_=false;

  double getEnergyMIGEN(const std::vector<double> &mean, const std::vector<double> &ftilde, const std::vector<double> &sigma,
                        const double scale, const double offset);
  double getEnergySP(const std::vector<double> &mean, const std::vector<double> &sigma,
                     const double scale, const double offset);
  double getEnergySPE(const std::vector<double> &mean, const std::vector<double> &sigma,
                      const double scale, const double offset);
  double getEnergyGJ(const std::vector<double> &mean, const std::vector<double> &sigma,
                     const double scale, const double offset);
  double getEnergyGJE(const std::vector<double> &mean, const std::vector<double> &sigma,
                      const double scale, const double offset);
  double getEnergyTerm(const unsigned i, const double mean_i, const double sigma_i) const;
  void getEnergyForceSP(const std::vector<double> &mean, const std::vector<double> &dmean_x, const std::vector<double> &dmean_b);
  void getEnergyForceSPE(const std::vector<double> &mean, const std::vector<double> &dmean_x, const std::vector<double> &dmean_b);
  void getEnergyForceGJ(const std::vector<double> &mean, const std::vector<double> &dmean_x, const std::vector<double> &dmean_b);
  void getEnergyForceGJE(const std::vector<double> &mean, const std::vector<double> &dmean_x, const std::vector<double> &dmean_b);
  void getEnergyForceMIGEN(const std::vector<double> &mean, const std::vector<double> &dmean_x, const std::vector<double> &dmean_b);
  void get_weights(double &weight, double &norm, double &neff);
  void replica_averaging(const double weight, const double norm, std::vector<double> &mean, std::vector<double> &dmean_b);
  void get_sigma_mean(const double weight, const double norm, const double neff, const std::vector<double> &mean);
  void do_regression_zero(const std::vector<double> &mean);
  void moveTilde(const std::vector<double> &mean_, double &old_energy);
  void moveScaleOffset(const std::vector<double> &mean_, double &old_energy);
  void moveSigmas(const std::vector<double> &mean_, double &old_energy, const unsigned i, const std::vector<unsigned> &indices, bool &breaknow);
  double doMonteCarlo(const std::vector<double> &mean);

public:
  /// register the keywords and the components shared by all the Metainference actions
  static void registerKeywords( Keywords& keys );
  /// read the input of the owner action; active is false for a CV without DOSCORE
  MetainferenceCore(Action& action, bool active);
  ~MetainferenceCore();
  MetainferenceCore(const MetainferenceCore&) = delete;
  MetainferenceCore& operator=(const MetainferenceCore&) = delete;
  /// set the experimental values, read the status file at restart and add the components to the owner
  void setup(const std::vector<double>& params, ActionWithValue& av, const ComponentGetter& getComponent);
  /// the back-calculated data, to be set before calculate()
  std::vector<double>& data() {
    return calc_data_;
  }
  /// run Metainference with the current data and return the energy summed over the replicas
  /// isel is the selector index, reweight_bias the value of the bias used with REWEIGHT
  double calculate(unsigned isel, double reweight_bias);
  /// derivative of the energy with respect to the back-calculated datum i
  double getMetaDer(const unsigned i) const {
    return metader_[i];
  }
  /// derivative of the energy with respect to the bias used with REWEIGHT
  double getBiasDer() const {
    return biasDer_;
  }
  /// the selector index (0 without SELECTOR)
  unsigned getSelector() const;
  bool doReweight() const {
    return do_reweight_;
  }
  bool isActive() const {
    return active_;
  }
  unsigned getWriteStride() const {
    return write_stride_;
  }
  unsigned getNarg() const {
    return narg;
  }
  /// write the status file if this is a step at which it should be written
  void update();
  void writeStatus();
};

}
}

#endif
