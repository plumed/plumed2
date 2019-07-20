/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2017 The plumed team
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
#ifndef __PLUMED_core_AverageBase_h
#define __PLUMED_core_AverageBase_h
#include "ActionPilot.h"
#include "ActionAtomistic.h"
#include "ActionWithValue.h"
#include "ActionWithArguments.h"
#include "tools/RMSD.h"

namespace PLMD {

class AverageBase :
  public ActionPilot,
  public ActionAtomistic,
  public ActionWithValue,
  public ActionWithArguments {
private:
  bool clearnextstep, firststep;
  Tensor rot;
  PLMD::RMSD myrmsd;
  std::string rmsd_type;
  Matrix<std::vector<Vector> > DRotDPos;
  std::vector<double> data;
  std::vector<Vector> atom_pos, der, direction, centeredpos, centeredreference;
protected:
  bool clearnorm;
  unsigned clearstride;
  unsigned n_real_args;
  std::vector<double> align, displace;
/// This is used to setup the components for the actions that store data
  void setupComponents( const unsigned& nreplicas );
/// This is used to transfer the data in runFinalJobs for actions that collect data
  void transferCollectedDataToValue( const std::vector<std::vector<double> >& mydata, const std::vector<double>& myweights );
/// Get the number of atoms that we are averaging
  unsigned getNumberOfAtomsToAverage() const ;
/// Get the position of the ith atom in the reference configuration
  Vector getReferencePosition(const unsigned& i );
public:
  static void registerKeywords( Keywords& keys );
  explicit AverageBase( const ActionOptions& );
  void clearDerivatives( const bool& force=false ) {}
  virtual void resizeValues() {}
  unsigned getNumberOfDerivatives() const ;
  void getInfoForGridHeader( std::string& gtype, std::vector<std::string>& argn, std::vector<std::string>& min,
                             std::vector<std::string>& max, std::vector<unsigned>& nbin,
                             std::vector<double>& spacing, std::vector<bool>& pbc, const bool& dumpcube ) const ;
  void getGridPointIndicesAndCoordinates( const unsigned& ind, std::vector<unsigned>& indices, std::vector<double>& coords ) const ;
  void getGridPointAsCoordinate( const unsigned& ind, const bool& setlength, std::vector<double>& coords ) const ;
/// These are required because we inherit from both ActionAtomistic and ActionWithArguments
  void lockRequests();
  void unlockRequests();
  void calculateNumericalDerivatives( ActionWithValue* a=NULL ) { plumed_error(); }
  void calculate() {}
  void apply() {}
  void update();
  virtual void accumulateNorm( const double& cweight ) = 0 ;
  virtual void accumulateGrid( const double& cweight ){ plumed_error(); }
  virtual void accumulateValue( const double& cweight, const std::vector<double>& val ) = 0;
  virtual void setReferenceConfig();
  virtual void accumulateAtoms( const double& cweight, const std::vector<Vector>& dir ) = 0;
  std::string getStrideClearAndWeights() const ;
  std::string getAtomsData() const ;
  virtual AtomNumber getAtomNumber(const AtomNumber& num ) const { plumed_merror("No virtual atoms " + getLabel() ); }
};

inline
unsigned AverageBase::getNumberOfAtomsToAverage() const {
  return atom_pos.size();
}

inline
Vector AverageBase::getReferencePosition(const unsigned& i ) {
  return myrmsd.getReference()[i];
}

}
#endif
