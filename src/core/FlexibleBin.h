/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2020 The plumed team
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
#ifndef __PLUMED_core_FlexibleBin_h
#define __PLUMED_core_FlexibleBin_h

#include<vector>

namespace PLMD {

class ActionWithArguments;

class FlexibleBin {
private:
  int type;
  // this contains all the infos about the CVs including periodicity
  ActionWithArguments *paction;
  double sigma;
  // variance is the matrix that really matters
  std::vector<double> variance;
  // this is only there
  std::vector<double> average;
  // minimum and maximum values
  std::vector<double> sigmamin;
  std::vector<double> sigmamax;
  std::vector<bool>  limitmax;
  std::vector<bool>  limitmin;
public:
  /// a constructor that takes the pointer of the action that contains it
  FlexibleBin(int type,ActionWithArguments *paction, double const &d, std::vector<double> &sigmamin, std::vector<double> &sigmamax );
  /// a constructor for 1D FlexBin (for PBMETAD)
  FlexibleBin(int type, ActionWithArguments *paction, unsigned iarg, double const &d, std::vector<double> &sigmamin, std::vector<double> &sigmamax);
  /// update the average (always for diffusion) or calculate the geom covariance (  only when do_when_zero is zero)
  void update(bool nowAddAHill);
  void update(bool nowAddAHill, unsigned iarg);
  std::vector<double> getMatrix() const;
  std::vector<double> getInverseMatrix() const;
  std::vector<double> getInverseMatrix(unsigned iarg) const;
  enum AdaptiveHillsType { none, diffusion, geometry };
};



}
#endif
