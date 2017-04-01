/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2016 The plumed team
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
#ifndef __PLUMED_bias_ReweightWham_h
#define __PLUMED_bias_ReweightWham_h

#include "ReweightBase.h"

namespace PLMD {
namespace bias {

class ReweightWham : public ReweightBase {
private:
  double thresh;
  unsigned maxiter;
  std::vector<unsigned> wlists;
  std::vector<double> stored_biases;
  std::vector<double> final_weights;
public:
  static void registerKeywords(Keywords&);
  explicit ReweightWham(const ActionOptions&ao);
  void calculateWeights( const unsigned& nframes );
  void clearData();
  double getLogWeight();
  double getWeight( const unsigned& iweight ) const ;
};

inline
double ReweightWham::getWeight( const unsigned& iweight ) const {
  plumed_dbg_assert( calculatedWeights && iweight<final_weights.size() );
  return final_weights[iweight];
}

}
}
#endif
