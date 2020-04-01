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
  std::vector<double> frame_weights, data, off_diag_bias, allweights, posdata;
  std::vector<std::vector<double> > alldata, allposdata;
  void retrieveDataPoint( const unsigned& itime, std::vector<double>& old_data );
  void computeCurrentBiasForData( const std::vector<double>& values, const bool& runserial, std::vector<double>& weights );
public:
  static void registerKeywords( Keywords& keys );
  explicit CollectFrames( const ActionOptions& );
  void turnOnBiasHistory();
  void calculate();
  void finishComputations( const std::vector<double>& buf );
  void accumulate( const std::vector<std::vector<Vector> >& dir );
  void buildCurrentTaskList( bool& forceAllTasks, std::vector<std::string>& actionsThatSelectTasks, std::vector<unsigned>& tflags );
  void performTask( const unsigned& current, MultiValue& myvals ) const ;
};

}
#endif
