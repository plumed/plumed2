/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2018 The VES code team
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
#ifndef __PLUMED_ves_VesTools_h
#define __PLUMED_ves_VesTools_h

#include <string>
#include <sstream>
#include <iomanip>
#include <limits>
#include <vector>

#include "core/ActionSet.h"


namespace PLMD {

class GridBase;

namespace ves {

class VesTools {
public:
  // Convert double into a string with more digits
  static void convertDbl2Str(const double value,std::string& str, unsigned int precision);
  static void convertDbl2Str(const double value,std::string& str);
  // copy grid values
  static void copyGridValues(GridBase* grid_pntr_orig, GridBase* grid_pntr_copy);
  static unsigned int getGridFileInfo(const std::string&, std::string&, std::vector<std::string>&, std::vector<std::string>&, std::vector<std::string>&, std::vector<bool>&, std::vector<unsigned int>&, bool&);
  //
  template<typename T> static std::vector<std::string> getLabelsOfAvailableActions(const ActionSet&);
  template<typename T> static T getPointerFromLabel(const std::string&, const ActionSet&, std::string&);
  template<typename T> static std::vector<T> getPointersFromLabels(const std::vector<std::string>&, const ActionSet&, std::string&);

};

inline
void VesTools::convertDbl2Str(const double value,std::string& str, unsigned int precision) {
  std::ostringstream ostr;
  ostr<<std::setprecision(precision)<<value;
  str=ostr.str();
}


inline
void VesTools::convertDbl2Str(const double value,std::string& str) {
  unsigned int precision = std::numeric_limits<double>::digits10 + 1;
  convertDbl2Str(value,str,precision);
}


template<typename T>
std::vector<std::string> VesTools::getLabelsOfAvailableActions(const ActionSet& actionset) {
  std::vector<std::string> avail_action_str(0);
  std::vector<T> avail_action_pntrs = actionset.select<T>();
  for(unsigned int i=0; i<avail_action_pntrs.size(); i++) {
    avail_action_str.push_back(avail_action_pntrs[i]->getLabel());
  }
  return avail_action_str;
}


template<typename T>
T VesTools::getPointerFromLabel(const std::string& action_label, const ActionSet& actionset, std::string& error_msg) {
  std::vector<std::string> action_labels(1);
  action_labels[0] = action_label;
  std::vector<T> action_pntrs = getPointersFromLabels<T>(action_labels,actionset,error_msg);
  return action_pntrs[0];
}


template<typename T>
std::vector<T> VesTools::getPointersFromLabels(const std::vector<std::string>& action_labels, const ActionSet& actionset, std::string& error_msg) {
  std::vector<T> action_pntrs(action_labels.size(),NULL);
  error_msg = "";
  std::vector<std::string> missing(0);
  for(unsigned int i=0; i<action_labels.size(); i++) {
    action_pntrs[i] = actionset.selectWithLabel<T>(action_labels[i]);
    if(action_pntrs[i]==NULL) {
      missing.push_back(action_labels[i]);
    }
  }
  // error handling
  if(missing.size()>0) {
    if(missing.size()==1) {
      error_msg = "label "+missing[0]+" does not exist\n";
    }
    else if(missing.size()>1) {
      std::string tmp="";
      for(unsigned int j=0; j<missing.size(); j++) {tmp +=missing[j]+" ";}
      error_msg = "labels "+tmp+"do not exist\n";
    }
    std::vector<T> avail_action_pntrs = actionset.select<T>();
    if(avail_action_pntrs.size()>0) {
      error_msg += "             Hint! the actions defined in the input file that can be used here are: \n";
      for(unsigned int i=0; i<avail_action_pntrs.size(); i++) {
        error_msg += "             " + avail_action_pntrs[i]->getName() + " with label " + avail_action_pntrs[i]->getLabel() + "\n";
      }
    }
    else {
      error_msg += "             Hint! no actions defined in the input file that can be used here, they should be defined before this actions\n";
    }
  }
  return action_pntrs;
}



}
}

#endif
