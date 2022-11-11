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
#ifndef __PLUMED_core_CollectFrames_h
#define __PLUMED_core_CollectFrames_h
#include "AverageBase.h"

namespace PLMD {

class CollectFrames : public AverageBase {
private:
  bool save_all_bias;
  std::vector<unsigned> task_counts;
  unsigned nvals, ndata, task_start;
  std::vector<double> frame_weights, data, posdata;
  void retrieveDataPoint( const unsigned& itime, std::vector<double>& old_data );
  void computeCurrentBiasForData( const std::vector<double>& values, const bool& runserial, std::vector<double>& weights );
public:
  static void registerKeywords( Keywords& keys );
  explicit CollectFrames( const ActionOptions& );
  void getInfoForGridHeader( std::string& gtype, std::vector<std::string>& argn, std::vector<std::string>& min,
                             std::vector<std::string>& max, std::vector<unsigned>& nbin,
                             std::vector<double>& spacing, std::vector<bool>& pbc, const bool& dumpcube ) const ;
  void getGridPointIndicesAndCoordinates( const unsigned& ind, std::vector<unsigned>& indices, std::vector<double>& coords ) const ;
  void turnOnBiasHistory();
  unsigned getNumberOfColumns() const override ;
  void calculate();
  void finishComputations( const std::vector<double>& buf );
  void firstUpdate() override ;
  void accumulate( const std::vector<std::vector<Vector> >& dir ) override ;
  void setupCurrentTaskList() override ;
  void performTask( const unsigned& current, MultiValue& myvals ) const ;
};

}
#endif
