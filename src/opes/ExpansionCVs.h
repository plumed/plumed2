/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2020 of Michele Invernizzi.

The opes module is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The opes module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_opes_ExpansionCVs_h
#define __PLUMED_opes_ExpansionCVs_h

#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"

namespace PLMD {
namespace opes {

/*
\ingroup INHERIT
This is the abstract base class to use for implementing expansion CVs (ECVs).
ECVs should be used together with the OPESexpanded action.
They take as argument regular CVs, and output them as components without modification.
*/

class ExpansionCVs:
  public ActionWithValue,
  public ActionWithArguments
{
protected:
  bool isReady_; //true only after initECVs
  double kbt_;
  unsigned totNumECVs_;

//methods useful for linear expansions
  void setSteps(std::vector<double>&,const unsigned,const std::string&);
  unsigned estimateSteps(const double,const double,const std::vector<double>&,const std::string&) const;

public:
  explicit ExpansionCVs(const ActionOptions&);
  virtual ~ExpansionCVs() {};
  void apply() override;
  void calculate() override;
  static void registerKeywords(Keywords&);
  inline unsigned getNumberOfDerivatives() override {return 1;};

  inline double getKbT() const {return kbt_;};
  inline unsigned getTotNumECVs() const {plumed_massert(isReady_,"cannot ask for totNumECVs before ECV isReady"); return totNumECVs_;};
  virtual std::vector< std::vector<unsigned> > getIndex_k() const; //might need to override this

  virtual void calculateECVs(const double *) = 0;
  virtual const double * getPntrToECVs(unsigned) = 0;
  virtual const double * getPntrToDerECVs(unsigned) = 0;
  virtual std::vector<std::string> getLambdas() const = 0;
  virtual void initECVs_observ(const std::vector<double>&,const unsigned,const unsigned) = 0; //arg: all the observed CVs, the total numer of CVs, the first CV index referring to this ECV
  virtual void initECVs_restart(const std::vector<std::string>&) = 0; //arg: the lambdas read from DeltaF_name relative to this ECV
};

}
}

#endif

