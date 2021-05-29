/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2018 The plumed team
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
#include "tools/IFile.h"

namespace PLMD {
namespace bias {

class ReadHills : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit ReadHills(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(ReadHills,"READ_HILLS")

void ReadHills::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","FILE","the file containing the hills to be read");
  keys.add("compulsory","GRID_MIN","auto","the lower bounds for the grid");
  keys.add("compulsory","GRID_MAX","auto","the upper bounds for the grid");
  keys.add("compulsory","GRID_BIN","the number of bins for the grid");
}

ReadHills::ReadHills(const ActionOptions& ao):
Action(ao),
ActionShortcut(ao)
{
  std::string hillsfile; parse("FILE",hillsfile); std::vector<std::string> values;
  IFile ifile; if( !ifile.FileExist(hillsfile) ) error("could not find file named " + hillsfile );
  ifile.allowIgnoredFields(); ifile.open(hillsfile); std::vector<std::string> fields; ifile.scanFieldList(fields);
  bool before_sigma=true;
  for(unsigned i=0; i<fields.size(); i++) {
    size_t founds=fields[i].find("sigma_");
    size_t foundm=fields[i].find("min_");
    size_t foundp=fields[i].find("max_");
    if ( founds!=std::string::npos || foundm!=std::string::npos ||  foundp!=std::string::npos )before_sigma=false;
    size_t found=fields[i].find("time");
    if( found==std::string::npos && before_sigma) {
        readInputLine( fields[i] + ": READ IGNORE_FORCES FILE=" + hillsfile + " VALUES=" + fields[i]  );
        values.push_back( fields[i] );
    }
  }
  // Do we have multivariate hills
  bool multivariate=false;
  if( ifile.FieldExist("multivariate")) {
      std::string sss; ifile.scanField("multivariate",sss);
      if(sss=="true") { multivariate=true;}
      else if(sss=="false") { multivariate=false;}
  }
  if( multivariate ) {
      std::string inum, jnum, col_string;
      readInputLine( getShortcutLabel() + "_zero: CONSTANT VALUES=0.0");
      for(unsigned i=0; i<values.size(); i++) {
          Tools::convert( i+1, inum ); 
          for(unsigned j=0;j<=i;++j) {
              Tools::convert( j+1, jnum );
              readInputLine( getShortcutLabel() + "_sigma_" + values[i] + "_" + values[j] + ": READ IGNORE_FORCES FILE=" + hillsfile +  
                             " VALUES=sigma_" + values[i] + "_" + values[j]  );
              if( i!=j ) col_string += " MATRIX" + inum + jnum + "=" + getShortcutLabel() + "_sigma_" + values[i] + "_" + values[j];
          }
          col_string += " MATRIX" + inum + inum + "=" + getShortcutLabel() + "_sigma_" + values[i] + "_" + values[i];
          for(unsigned j=i+1;j<values.size();++j) {
              Tools::convert( j+1, jnum ); col_string += " MATRIX" + inum + jnum + "=" + getShortcutLabel() + "_zero";
          }
      }
      // This is cholesky decomposition of matrix    
      readInputLine( getShortcutLabel() + "_chol: CONCATENATE " + col_string );
      // Transpose
      readInputLine( getShortcutLabel() + "_cholT: TRANSPOSE ARG=" + getShortcutLabel() + "_chol");
      // Recompose sigma matrix
      readInputLine( getShortcutLabel() + "_sigma: DOT ARG1=" + getShortcutLabel() + "_chol ARG2=" + getShortcutLabel() + "_cholT"); 
      // And compute final metric matrix
      readInputLine( getShortcutLabel() + "_icov: INVERT_MATRIX ARG=" + getShortcutLabel() + "_sigma");
  } else {
      std::string col_string = "ARG=" + getShortcutLabel() + "_icov_" + values[0];
      for(unsigned i=0; i<values.size(); i++) {
          readInputLine( getShortcutLabel() + "_sigma_" + values[i] + ": READ IGNORE_FORCES FILE=" + hillsfile + " VALUES=sigma_" + values[i]  );
          readInputLine( getShortcutLabel() + "_icov_" + values[i] + ": MATHEVAL ARG1=" + getShortcutLabel() + "_sigma_" + values[i] + " FUNC=1/(x*x) PERIODIC=NO" );
          if( i>0 ) col_string += "," + getShortcutLabel() + "_icov_" + values[i];
      }
      readInputLine( getShortcutLabel() + "_icov: CONCATENATE " + col_string );
  }
  readInputLine( getShortcutLabel() + "_height: READ IGNORE_FORCES FILE=" + hillsfile + " VALUES=height");
  // And sum the hills
  std::vector<std::string> gmin(values.size()), gmax(values.size()), grid_nbins(values.size());
  parseVector("GRID_MIN",gmin); parseVector("GRID_MAX",gmax); parseVector("GRID_BIN",grid_nbins);
  std::string input = getShortcutLabel() + "_kde: KDE_CALC METRIC=" + getShortcutLabel() + "_icov ARG1=" + values[0] + " HEIGHTS=" + getShortcutLabel() + "_height";
  std::string gminstr=" GRID_MIN=" + gmin[0]; std::string gmaxstr=" GRID_MAX=" + gmax[0]; std::string gbinstr=" GRID_BIN=" + grid_nbins[0];
  for(unsigned i=1;i<values.size();++i) { 
    std::string num; Tools::convert( i+1, num ); input += " ARG" + num + "=" + values[i];
    gminstr += "," + gmin[i]; gmaxstr += "," + gmax[i]; gbinstr += "," + grid_nbins[i];

  }
  readInputLine( input + " " + gminstr + " " + gmaxstr + " " + " " + gbinstr );
  readInputLine( getShortcutLabel() + ": AVERAGE ARG=" + getShortcutLabel() + "_kde NORMALIZATION=false STRIDE=1");
}

}
}
