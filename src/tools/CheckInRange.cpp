/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#include "CheckInRange.h"
#include "Tools.h"

namespace PLMD {

bool CheckInRange::setBounds( const unsigned& n, const std::vector<std::string>& str_lower, const std::vector<std::string>& str_upper, std::string& errors ) {
  if( str_upper.size()!=n && str_upper.size()>0 ) {
    errors="wrong number of arguments for LESS_THAN_OR_EQUAL keyword";
    return false;
  }
  if( str_lower.size()!=n && str_lower.size()>0 ) {
    errors="wrong number of arguments for GREATER_THAN_OR_EQUAL keyword";
    return false;
  }
  if( str_upper.size()>0 && str_lower.size()>0 ) {
    lower.resize( str_lower.size() );
    upper.resize( str_upper.size() );
    for(unsigned i=0; i<upper.size(); ++i) {
      if( str_lower[i]=="none" ) {
        lower[i] = -std::numeric_limits<double>::max();
      } else {
        Tools::convert( str_lower[i], lower[i] );
      }
      if( str_upper[i]=="none" ) {
        upper[i] = std::numeric_limits<double>::max();
      } else {
        Tools::convert( str_upper[i], upper[i] );
      }
    }
  } else if( str_upper.size()>0 ) {
    upper.resize( str_upper.size() );
    for(unsigned i=0; i<upper.size(); ++i) {
      if( str_upper[i]=="none" ) {
        upper[i] = std::numeric_limits<double>::max();
      } else {
        Tools::convert( str_upper[i], upper[i] );
      }
    }
  } else if( str_lower.size()>0 ) {
    lower.resize( str_lower.size() );
    for(unsigned i=0; i<lower.size(); ++i) {
      if( str_lower[i]=="none" ) {
        lower[i] = -std::numeric_limits<double>::max();
      } else {
        Tools::convert( str_lower[i], lower[i] );
      }
    }
  }
  return true;
}

bool CheckInRange::wereSet() const {
  return lower.size()>0 || upper.size()>0;
}

std::string CheckInRange::report( const std::vector<std::string>& a ) const {
  if( upper.size()>0 && lower.size()>0 ) {
    std::string str_l, str_u;
    Tools::convert( upper[0], str_u );
    Tools::convert( lower[0], str_l );
    std::string out="only printing indices of atoms that have " + str_l + " <= " + a[0] + " <=" + str_u;
    for(unsigned i=1; i<upper.size(); ++i) {
      Tools::convert( upper[i], str_u );
      Tools::convert( lower[i], str_l );
      out += " and " + str_l + " <= " + a[i] + " <=" + str_u;
    }
    return out;
  }
  if( upper.size()>0 ) {
    std::string str_u;
    Tools::convert( upper[0], str_u );
    std::string out="only printing indices of atoms that have " + a[0] + " <= " + str_u;
    for(unsigned i=1; i<upper.size(); ++i) {
      Tools::convert( upper[i], str_u );
      out += " and " + a[i] + " <= " + str_u;
    }
    return out;
  }
  std::string str_l;
  Tools::convert( lower[0], str_l );
  std::string out="only printing indices of atoms that have " + str_l + " <= " + a[0];
  for(unsigned i=1; i<lower.size(); ++i) {
    Tools::convert( lower[i], str_l );
    out += " and " + str_l + " <= " + a[i];
  }
  return out;
}

bool CheckInRange::check( const std::vector<double>& vals ) const {
  for(unsigned j=0; j<vals.size(); ++j) {
    if( upper.size()>0 && vals[j]>upper[j] ) {
      return false;
    }
    if( lower.size()>0 && vals[j]<lower[j] ) {
      return false;
    }
  }
  return true;
}

}
