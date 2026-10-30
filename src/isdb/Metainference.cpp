/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2023 The plumed team
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

#include "bias/Bias.h"
#include "core/ActionRegister.h"
#include "MetainferenceCore.h"

namespace PLMD {
namespace isdb {

//+PLUMEDOC ISDB_BIAS METAINFERENCE
/*
Calculates the Metainference energy for a set of experimental data.

The metainference method that is introduced in the first paper cited below is a Bayesian framework
to model heterogeneous systems by integrating prior information with noisy, ensemble-averaged data.
Metainference models a system and quantifies the level of noise in the data by considering a set of replicas of the system.

Calculated experimental data are given in input as ARG while reference experimental values
can be given either from fixed components of other actions using PARARG or as numbers using
PARAMETERS. The default behavior is that of averaging the data over the available replicas,
if this is not wanted the keyword NOENSEMBLE prevent this averaging.

Metadynamics Metainference (see second paper cited below) more in general biased Metainference requires the knowledge of
biasing potential in order to calculate the weighted average. In this case the value of the bias
can be provided as the last argument in ARG and adding the keyword REWEIGHT. To avoid the noise
resulting from the instantaneous value of the bias the weight of each replica can be averaged
over a give time using the keyword AVERAGING.

The data can be averaged by using multiple replicas and weighted for a bias if present.
The functional form of Metainference can be chosen among four variants selected
with NOISE=GAUSS,MGAUSS,OUTLIERS,MOUTLIERS,GENERIC which correspond to modelling the noise for
the arguments as a single gaussian common to all the data points, a gaussian per data
point, a single long-tailed gaussian common to all the data points, a log-tailed
 gaussian per data point or using two distinct noises as for the most general formulation of Metainference.
In this latter case the noise of the replica-averaging is gaussian (one per data point) and the noise for
the comparison with the experimental data can chosen using the keyword LIKELIHOOD
between gaussian or log-normal (one per data point), furthermore the evolution of the estimated average
over an infinite number of replicas is driven by DFTILDE.

As for Metainference theory there are two sigma values: SIGMA_MEAN0 represent the
error of calculating an average quantity using a finite set of replica and should
be set as small as possible following the guidelines for replica-averaged simulations
in the framework of the Maximum Entropy Principle. Alternatively, this can be obtained
automatically using the internal sigma mean optimization as introduced in the third paper cited below
(OPTSIGMAMEAN=SEM), in this second case sigma_mean is estimated from the maximum standard error
of the mean either over the simulation or over a defined time using the keyword AVERAGING.
SIGMA_BIAS is an uncertainty parameter, sampled by a MC algorithm in the bounded interval
defined by SIGMA_MIN and SIGMA_MAX. The initial value is set at SIGMA0. The MC move is a
random displacement of maximum value equal to DSIGMA. If the number of data point is
too large and the acceptance rate drops it is possible to make the MC move over mutually
exclusive, random subset of size MC_CHUNKSIZE and run more than one move setting MC_STEPS
in such a way that MC_CHUNKSIZE*MC_STEPS will cover all the data points.

Calculated and experimental data can be compared modulo a scaling factor and/or an offset
using SCALEDATA and/or ADDOFFSET, the sampling is obtained by a MC algorithm either using
a flat or a gaussian prior setting it with SCALE_PRIOR or OFFSET_PRIOR. With a flat prior the
moves are reflected into the interval given by SCALE_MIN and SCALE_MAX (OFFSET_MIN and OFFSET_MAX),
with a gaussian prior the prior is centered in SCALE0 (OFFSET0) and its width is DSCALE (DOFFSET),
which also sets the size of the moves.

The random numbers of the Monte Carlo are seeded from the clock, so that two runs are different;
SEED makes them reproducible. The status of Metainference (uncertainties, history of sigma mean,
scale, offset, averaged weights) is written every WRITE_STRIDE steps in the file MISTATUS followed
by the label of the action (STATUS_FILE changes MISTATUS) and is read back when restarting.

The same keywords can be used in the ISDB collective variables (e.g. [NOE](NOE.md), [RDC](RDC.md))
together with DOSCORE: in this case the score is calculated directly by the collective variable,
which is faster with many data.

## Examples

In the following example we calculate a set of [RDC](RDC.md), take the replica-average of
them and comparing them with a set of experimental values. RDCs are compared with
the experimental data but for a multiplication factor SCALE that is also sampled by
MC on-the-fly

```plumed
RDC ...
LABEL=rdc
SCALE=0.0001
GYROM=-72.5388
ATOMS1=22,23
ATOMS2=25,27
ATOMS3=29,31
ATOMS4=33,34
... RDC

METAINFERENCE ...
ARG=rdc.*
NOISETYPE=MGAUSS
PARAMETERS=1.9190,2.9190,3.9190,4.9190
SCALEDATA SCALE0=1 SCALE_MIN=0.1 SCALE_MAX=3 DSCALE=0.01
SIGMA0=0.01 SIGMA_MIN=0.00001 SIGMA_MAX=3 DSIGMA=0.01
SIGMA_MEAN0=0.001
LABEL=spe
... METAINFERENCE

PRINT ARG=spe.bias FILE=BIAS STRIDE=1
```

in the following example instead of using one uncertainty parameter per data point we use
a single uncertainty value in a long-tailed gaussian to take into account for outliers, furthermore
the data are weighted for the bias applied to other variables of the system.

```plumed
RDC ...
LABEL=rdc
SCALE=0.0001
GYROM=-72.5388
ATOMS1=22,23
ATOMS2=25,27
ATOMS3=29,31
ATOMS4=33,34
... RDC

cv1: TORSION ATOMS=1,2,3,4
cv2: TORSION ATOMS=2,3,4,5
mm: METAD ARG=cv1,cv2 HEIGHT=0.5 SIGMA=0.3,0.3 PACE=200 BIASFACTOR=8 WALKERS_MPI

METAINFERENCE ...
#SETTINGS NREPLICAS=2
ARG=rdc.*,mm.bias
REWEIGHT
NOISETYPE=OUTLIERS
PARAMETERS=1.9190,2.9190,3.9190,4.9190
SCALEDATA SCALE0=1 SCALE_MIN=0.1 SCALE_MAX=3 DSCALE=0.01
SIGMA0=0.01 SIGMA_MIN=0.00001 SIGMA_MAX=3 DSIGMA=0.01
SIGMA_MEAN0=0.001
LABEL=spe
... METAINFERENCE
```

(See also [RDC](RDC.md), [PBMETAD](PBMETAD.md)).

*/
//+ENDPLUMEDOC

class Metainference : public bias::Bias {
  // number of experimental data (the bias used with REWEIGHT is not counted)
  unsigned narg;
  // the Metainference algorithm
  MetainferenceCore core_;

public:
  explicit Metainference(const ActionOptions&);
  void calculate() override;
  void update() override;
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(Metainference,"METAINFERENCE")

void Metainference::registerKeywords(Keywords& keys) {
  Bias::registerKeywords(keys);
  keys.addInputKeyword("optional","PARARG","scalar","reference values for the experimental data, these can be provided as arguments without derivatives");
  keys.add("optional","PARAMETERS","reference values for the experimental data");
  MetainferenceCore::registerKeywords(keys);
  keys.addDOI("10.1126/sciadv.1501177");
  keys.addDOI("10.1038/srep31232");
  keys.addDOI("10.1063/1.4981211");
}

Metainference::Metainference(const ActionOptions&ao):
  PLUMED_BIAS_INIT(ao),
  narg(0),
  core_(*this,true) {
  narg = getNumberOfArguments();
  if(core_.doReweight()) {
    narg--;
  }

  std::vector<double> parameters;
  parseVector("PARAMETERS",parameters);
  if(parameters.size()!=static_cast<unsigned>(narg)&&!parameters.empty()) {
    error("Size of PARAMETERS array should be either 0 or the same as of the number of arguments in ARG1");
  }

  std::vector<Value*> arg2;
  parseArgumentList("PARARG",arg2);
  if(!arg2.empty()) {
    if(parameters.size()>0) {
      error("It is not possible to use PARARG and PARAMETERS together");
    }
    if(arg2.size()!=narg) {
      error("Size of PARARG array should be the same as number for arguments in ARG");
    }
    for(unsigned i=0; i<arg2.size(); i++) {
      parameters.push_back(arg2[i]->get());
      if(arg2[i]->hasDerivatives()==true) {
        error("PARARG can only accept arguments without derivatives");
      }
    }
  }

  if(parameters.size()!=narg) {
    error("PARARG or PARAMETERS arrays should include the same number of elements as the arguments in ARG");
  }

  checkRead();

  core_.setup(parameters, *this, [this](const std::string& name) {
    return getPntrToComponent(name);
  });
}

void Metainference::calculate() {
  std::vector<double>& data=core_.data();
  for(unsigned i=0; i<narg; ++i) {
    data[i]=getArgument(i);
  }

  const double ene = core_.calculate(core_.getSelector(), core_.doReweight() ? getArgument(narg) : 0.);

  for(unsigned i=0; i<narg; ++i) {
    setOutputForce(i, -core_.getMetaDer(i));
  }
  if(core_.doReweight()) {
    setOutputForce(narg, -core_.getBiasDer());
  }

  setBias(ene);
}

void Metainference::update() {
  core_.update();
}

}
}

