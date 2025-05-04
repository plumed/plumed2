/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2021 The VES code team
   (see the PEOPLE-VES file at the root of this folder for a list of names)

   See http://www.ves-code.org for more information.

   This file is part of VES code module.

   The VES code module is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   The VES code module is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with the VES code module.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#include "TargetDistribution.h"
#include "GridIntegrationWeights.h"

#include "core/ActionRegister.h"
#include "tools/Grid.h"
#include "core/PlumedMain.h"



namespace PLMD {
namespace ves {

//+PLUMEDOC VES_TARGETDIST TD_WELLTEMPERED
/*
Well-tempered target distribution (dynamic).

Use as a target distribution the well-tempered distribution discussed in the first paper cited below,
which is given by

$$
p(\mathbf{s}) =
\frac{e^{-(\beta/\gamma) F(\mathbf{s})}}
{\int d\mathbf{s}\, e^{-(\beta/\gamma) F(\mathbf{s})}} =
\frac{[P_{0}(\mathbf{s})]^{1/\gamma}}
{\int d\mathbf{s}\, [P_{0}(\mathbf{s})]^{1/\gamma}}
$$

where $\gamma$ is a so-called bias factor and $P_{0}(\mathbf{s})$ is the
unbiased canonical distribution of the CVs. This target distribution thus
corresponds to a biased ensemble where, as compared to the unbiased one,
the probability peaks have been broaden and the fluctuations of the CVs are
enhanced.
The value of the bias factor $\gamma$ determines by how much the fluctuations
are enhanced.

The well-tempered distribution can be view as sampling on
an effective free energy surface $\tilde{F}(\mathbf{s}) = (1/\gamma) F(\mathbf{s})$
which has largely the same metastable states as the original $F(\mathbf{s})$
but with barriers that have been reduced by a factor of $\gamma$.
Generally one should use a value of $\gamma$ that results in
effective barriers on the order of few $k_{\mathrm{B}}T$
such that thermal fluctuations can easily induce transitions
between different metastable states.

At convergence the relationship between the bias potential and the free
energy surface is given by

$$
F(\mathbf{s}) = - \left(\frac{1}{1-\gamma^{-1}} \right) V(\mathbf{s})
$$

This target distribution depends directly on the free energy surface
$F(\mathbf{s})$ which is quantity that we do not know a-priori and
want to obtain. Therefore, this target distribution
is iteratively updated according to

$$
p^{(m+1)}(\mathbf{s}) =
\frac{e^{-(\beta/\gamma) F^{(m+1)}(\mathbf{s})}}
{\int d\mathbf{s}\, e^{-(\beta/\gamma) F^{(m+1)}(\mathbf{s})}}
$$

where $F^{(m+1)}(\mathbf{s})$ is the current best estimate of the
free energy surface obtained according to

$$
F^{(m+1)}(\mathbf{s}) =
- V^{(m+1)}(\mathbf{s}) - \frac{1}{\beta} \log p^{(m)}(\mathbf{s}) =
- V^{(m+1)}(\mathbf{s}) + \frac{1}{\gamma} F^{(m)}(\mathbf{s})
$$

The frequency of performing this update needs to be set in the
optimizer used in the calculation. Normally it is sufficient
to do it every 100-1000 bias update iterations.

## Examples

Employ a well-tempered target distribution with a bias factor of 10

```plumed
td_welltemp: TD_WELLTEMPERED BIASFACTOR=10
```

*/
//+ENDPLUMEDOC

class TD_WellTempered: public TargetDistribution {
private:
  double bias_factor_;
public:
  static void registerKeywords(Keywords&);
  explicit TD_WellTempered(const ActionOptions& ao);
  void updateGrid() override;
  double getValue(const std::vector<double>&) const override;
  ~TD_WellTempered() {}
};


PLUMED_REGISTER_ACTION(TD_WellTempered,"TD_WELLTEMPERED")


void TD_WellTempered::registerKeywords(Keywords& keys) {
  TargetDistribution::registerKeywords(keys);
  keys.add("compulsory","BIASFACTOR","The bias factor used for the well-tempered distribution.");
  keys.addDOI("10.1103/PhysRevLett.100.020603");
  keys.addDOI("10.1021/acs.jctc.5b00076");
}


TD_WellTempered::TD_WellTempered(const ActionOptions& ao):
  PLUMED_VES_TARGETDISTRIBUTION_INIT(ao),
  bias_factor_(0.0) {
  log.printf("  Well-tempered target distribution, see and cite ");
  log << plumed.cite("Valsson and Parrinello, J. Chem. Theory Comput. 11, 1996-2002 (2015)");
  log << plumed.cite("Barducci, Bussi, and Parrinello, Phys. Rev. Lett. 100, 020603 (2008)");
  log.printf("\n");
  parse("BIASFACTOR",bias_factor_);
  if(bias_factor_<=1.0) {
    plumed_merror("TD_WELLTEMPERED target distribution: the value of the bias factor doesn't make sense, it should be larger than 1.0");
  }
  setDynamic();
  setFesGridNeeded();
  checkRead();
}


double TD_WellTempered::getValue(const std::vector<double>& argument) const {
  plumed_merror("getValue not implemented for TD_WellTempered");
  return 0.0;
}


void TD_WellTempered::updateGrid() {
  double beta_prime = getBeta()/bias_factor_;
  plumed_massert(getFesGridPntr()!=NULL,"the FES grid has to be linked to use TD_WellTempered!");
  std::vector<double> integration_weights = GridIntegrationWeights::getIntegrationWeights(getTargetDistGridPntr());
  double norm = 0.0;
  for(Grid::index_t l=0; l<targetDistGrid().getSize(); l++) {
    double value = beta_prime * getFesGridPntr()->getValue(l);
    logTargetDistGrid().setValue(l,value);
    value = exp(-value);
    norm += integration_weights[l]*value;
    targetDistGrid().setValue(l,value);
  }
  targetDistGrid().scaleAllValuesAndDerivatives(1.0/norm);
  logTargetDistGrid().setMinToZero();
}


}
}
