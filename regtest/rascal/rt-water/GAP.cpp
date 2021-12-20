#ifdef __PLUMED_HAS_RASCAL
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"

#include "rascal/utils/json_io.hh"

#include <chrono>
#include <cmath>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <fstream>
#include <list>
#include <string>
#include <cmath>

using namespace rascal;

namespace PLMD {
namespace rascal {

//+PLUMEDOC COLVAR GAP
/*
GAP potential implementation using RASCAL


*/
//+ENDPLUMEDOC

class GAP : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit GAP(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(GAP,"GAP")

void GAP::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords(keys); 
  keys.add("compulsory","FILE","the file that contains the details of the GAP potential");
  keys.add("numbered","SPECIES","the atoms in each species type"); keys.reset_style("SPECIES","atoms");
  keys.add("compulsory","NFEATURES","the number of SOAP features that are being computed");
}

GAP::GAP(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao)
{
  std::string filen; parse("FILE",filen);
  // Open the file
  std::ifstream ifile(filen); 
  // Read in the hyper parameters for the representation
  json mydata; ifile>>mydata;
  // Get the weights in a setup action so they can be used
  json weights=mydata["init_params"]["weights"]; std::string ddd, myvec; Tools::convert(weights[1][0][0],myvec); 
  for(unsigned i=1; i<weights[1].size(); ++i ) { Tools::convert(weights[1][i][0],ddd); myvec += "," + ddd; }
  readInputLine( getShortcutLabel() + "_weights: CONSTANT_VALUE VALUES=" + myvec );
  // Now get the values
  json vals=mydata["init_params"]["X_train"]["data"]["sparse_points"]["values"];

  // Read in the species that we are calculating this for
  std::string allspecies, spec; parse("SPECIES",spec); std::vector<double> mat_data;
  if( spec.length()>0 ) {
      allspecies = "SPECIES=" + spec;
      // Store all the matrix data
      for(unsigned i=0; i<vals[0][1][0][1].size(); ++i) mat_data.push_back( vals[0][1][0][1][i] );  
 
  } else {
   
      for(unsigned i=1;;++i) {
          if( !parseNumbered("SPECIES", i, spec ) ) break;
          std::string num; Tools::convert( i, num ); 
          allspecies += " SPECIES" + num + "=" + spec;
          // Store all the matrix data 
          for(unsigned j=0; j<vals[i-1][1][0][1].size(); ++j) mat_data.push_back( vals[i-1][1][0][1][j] );
      }
  }
  // Now print out the reference matrix
  OFile ofile; ofile.open( getShortcutLabel() + "_matrixfile"); 
  std::vector<unsigned> shape(2); shape[1]=weights[1].size(); shape[0] = mat_data.size() / shape[1];
  for(unsigned i=0; i<shape[0]; ++i) {
      for(unsigned j=0; j<shape[1]; ++j) ofile.printf( "%f ", mat_data[j*shape[0] + i] );
      ofile.printf("\n");
  }
  ofile.close(); 
  // Now read the matrix in
  readInputLine( getShortcutLabel() + "_matrix: CONSTANT_VALUE FILE=" + getShortcutLabel() + "_matrixfile" );
  // Create the SOAP object
  readInputLine( getShortcutLabel() + "_soap: SPHERICAL_INVARIANTS HYPERPARAMS=" + mydata["init_params"]["kernel"]["init_params"]["representation"]["init_params"].dump() + " " + allspecies );
  // Now do the matrix multiplication 
  readInputLine( getShortcutLabel() + "_prod: DOT ARG1=" + getShortcutLabel() + "_soap ARG2=" + getShortcutLabel() + "_matrix" ); 
  // And the vector matrix times the vector
  readInputLine( getShortcutLabel() + "_vprod: DOT ARG1=" + getShortcutLabel() + "_prod ARG2=" + getShortcutLabel() + "_weights");
  // And sum the vector matrix to get the total energy
  readInputLine( getShortcutLabel() + ": SUM ARG=" + getShortcutLabel() + "_vprod PERIODIC=NO");
  // And add a force on this value
  readInputLine( "BIASVALUE ARG=" + getShortcutLabel() );
}

}
}
#endif


