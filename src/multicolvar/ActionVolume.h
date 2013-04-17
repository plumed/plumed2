/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#ifndef __PLUMED_multicolvar_ActionVolume_h
#define __PLUMED_multicolvar_ActionVolume_h

#include "core/ActionAtomistic.h"
#include "tools/HistogramBead.h"
#include "core/ActionWithValue.h"
#include "vesselbase/ActionWithVessel.h"
#include "vesselbase/BridgeVessel.h"
#include "MultiColvarBase.h"

namespace PLMD {
namespace multicolvar {

/**
\ingroup INHERIT
This is the abstract base class to use for implementing a new way of definining a particular region of the simulation
box. You can use this to calculate the number of atoms inside that part or the average value of a quantity like the 
coordination number inside that part of the cell. 
*/

class ActionVolume :
  public ActionAtomistic,
  public ActionWithValue,
  public vesselbase::ActionWithVessel
  {
friend class Region;
private:
/// The value of sigma
  double sigma;
/// Are we interested in the area outside the colvar
  bool not_in;
/// The bead for the histogram
  HistogramBead bead;
/// The action that is calculating the colvars of interest
  MultiColvarBase* mycolv;
/// The vessel that bridges
  vesselbase::BridgeVessel* myBridgeVessel;
/// Everything for controlling the updating of neighbor lists
  int updateFreq;
  unsigned lastUpdate;
protected:
  double getSigma() const ;
  MultiColvarBase* getPntrToMultiColvar();
public:
  static void registerKeywords( Keywords& keys );
  ActionVolume(const ActionOptions&);
/// Don't actually clear the derivatives when this is called from plumed main.  
/// They are calculated inside another action and clearing them would be bad  
  void clearDerivatives(){}
/// Get the number of derivatives for this action
  unsigned getNumberOfDerivatives();  // N.B. This is replacing the virtual function in ActionWithValue
/// Is the output quantity periodic
  bool isPeriodic();
/// Jobs to be done when the action is activated
  void prepare();
/// Do jobs required before tasks are undertaken
  void doJobsRequiredBeforeTaskList();
/// This calculates all the vessels and is called from within a bridge vessel
  void performTask(const unsigned& i );
/// Routines that have to be defined so as not to have problems with virtual methods 
  void deactivate_task();
  void calculate(){}
/// We need our own calculate numerical derivatives here
  void calculateNumericalDerivatives();
  virtual void setupRegion()=0;
  virtual bool derivativesOfFractionalCoordinates()=0;
  virtual double calculateNumberInside( const Vector& cpos, HistogramBead& bead, Vector& derivatives )=0;
  void apply(){};
};

inline
double ActionVolume::getSigma() const {
  return sigma;
}

inline
MultiColvarBase* ActionVolume::getPntrToMultiColvar(){
  return mycolv;
}

inline
unsigned ActionVolume::getNumberOfDerivatives(){
  return mycolv->getNumberOfDerivatives();
}

}
}
#endif
