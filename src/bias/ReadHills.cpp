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
#include "ReadHills.h"
#include "core/ActionRegister.h"
#include "tools/IFile.h"

namespace PLMD {
namespace bias {

PLUMED_REGISTER_ACTION(ReadHills,"READ_HILLS")

void ReadHills::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","FILE","HILLS","the file containing the hills to be read");
  keys.addFlag("IN_STARTUP",false,"read the full column of data to a constant to make all hills available during whole run");
  keys.needsAction("READ");
  keys.needsAction("READCOLUMN");
}

void ReadHills::registerGridKeywords(Keywords & keys) {
  keys.add("compulsory","GRID_MIN","auto","the lower bounds for the grid");
  keys.add("compulsory","GRID_MAX","auto","the upper bounds for the grid");
  keys.add("compulsory","GRID_BIN","the number of bins for the grid");
}

void ReadHills::createGrid( const std::string& argvals, const std::string& bandvals, const std::string heightvals, ActionShortcut* action ) {
  std::string gmin, gmax, grid_nbins;
  action->parse("GRID_MIN",gmin);
  action->parse("GRID_MAX",gmax);
  action->parse("GRID_BIN",grid_nbins);
  // Now lets create the KDE action
  action->readInputLine( action->getShortcutLabel() + "_bias: KDE ARG=" + argvals + " BANDWIDTH=" + bandvals + " HEIGHTS=" + heightvals + " GRID_MIN=" + gmin + " GRID_MAX=" + gmax + " GRID_BIN=" + grid_nbins );
}

std::string ReadHills::fixArgumentName( const std::string& argin ) {
  std::string argout = argin;
  std::size_t dot=argin.find(".");
  if( dot!=std::string::npos ) {
    argout = argin.substr(0,dot) + "_" + argin.substr(dot+1);
  }
  return argout;
}

void ReadHills::getHillsDetails( const std::string& hillsfilename, std::vector<std::string>& values, bool& multivariate, bool fixnames ) {
  IFile ifile;
  ifile.allowIgnoredFields();
  ifile.open(hillsfilename);
  std::vector<std::string> fields;
  ifile.scanFieldList(fields);
  bool before_sigma=true;
  for(unsigned i=0; i<fields.size(); i++) {
    std::size_t founds=fields[i].find("sigma_");
    std::size_t foundm=fields[i].find("min_");
    std::size_t foundp=fields[i].find("max_");
    if ( founds!=std::string::npos || foundm!=std::string::npos ||  foundp!=std::string::npos ) {
      before_sigma=false;
    }
    std::size_t found=fields[i].find("time");
    if( found==std::string::npos && before_sigma) {
      values.push_back( fields[i] );
    }
  }
  if( fixnames ) {
    for(unsigned i=0; i<values.size(); ++i) {
      values[i] = fixArgumentName(values[i]);
    }
  }
  multivariate=false;
  if( ifile.FieldExist("multivariate")) {
    std::string sss;
    ifile.scanField("multivariate",sss);
    if(sss=="true") {
      multivariate=true;
    } else if(sss=="false") {
      multivariate=false;
    }
  }
  ifile.scanField();
  ifile.close();
}

ReadHills::ReadHills(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  bool instartup;
  parseFlag("IN_STARTUP",instartup);
  std::string readcommand = "READ IGNORE_TIME IGNORE_COMPONENT VALUES=";
  if( instartup ) {
    readcommand = "READCOLUMN VALUE=";
  }
  bool multivariate;
  std::string hillsfile;
  parse("FILE",hillsfile);
  std::vector<std::string> values;
  IFile ifile;
  if( !ifile.FileExist(hillsfile) ) {
    error("could not find file named " + hillsfile );
  }
  getHillsDetails( hillsfile, values, multivariate, false );
  for(unsigned i=0; i<values.size(); ++i) {
    readInputLine( getShortcutLabel() + "_" + fixArgumentName(values[i]) + ": " + readcommand + values[i] + " FILE=" + hillsfile );
  }
  if( multivariate ) {
    error("multivariate gaussians not implemented");
  } else {
    // Read the widths
    for(unsigned i=0; i<values.size(); ++i) {
      readInputLine( getShortcutLabel() + "_sigma_" + fixArgumentName(values[i]) + ": " + readcommand + "sigma_" + values[i] + " FILE=" + hillsfile );
    }
    // Read the heights
    readInputLine( getShortcutLabel() + "_height: " + readcommand + "height FILE=" + hillsfile );
  }
}

}
}
