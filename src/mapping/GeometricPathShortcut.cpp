/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2018 The plumed team
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
#include "Path.h"
#include "tools/PDB.h"
#include "core/ActionRegister.h"

using namespace std;

namespace PLMD {
namespace mapping {

class GeometricPathShortcut : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit GeometricPathShortcut(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(GeometricPathShortcut,"GPATH")

void GeometricPathShortcut::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys ); Path::registerInputFileKeywords( keys );
  keys.add("optional","PROPERTY","read in path coordinates by finding option with this label in remark of pdb frames");
  keys.addFlag("UNFIX_FRAMES",false,"allows us to not fix frames for paths");
}

GeometricPathShortcut::GeometricPathShortcut( const ActionOptions& ao ):
  Action(ao),
  ActionShortcut(ao)
{
  std::vector<std::string> refactions; std::string mtype, refname;
  // Create list of reference configurations that PLUMED will use
  Path::readInputFrames( mtype, refname, true, this, refactions );
  // Now get coordinates on spath
  std::string pname, coord_str, ref_str; parse("PROPERTY",pname); 
  FILE* fp=std::fopen(refname.c_str(),"r"); bool do_read=true; double fake_unit=0.1; unsigned nfram = 0;
  while (do_read ) {
      PDB mypdb; do_read=mypdb.readFromFilepointer(fp,false,fake_unit);
      if( !do_read ) break ;

      if( nfram==0 ) { ref_str = " REFFRAMES=" + refactions[nfram]; } else { ref_str += "," + refactions[nfram]; }
      if( pname.length()>0 ) {
          double pval; std::string propstr;
          if( !mypdb.getArgumentValue(pname, pval) ) plumed_merror("could not find property named " + pname + " in input file " + refname );
          Tools::convert( pval, propstr ); 
          if( nfram==0 ) { coord_str = "COORDINATES=" + propstr; } else { coord_str += "," + propstr; }
      } else {
          std::string propstr; Tools::convert( nfram+1, propstr );
          if( nfram==0 ) { coord_str = "COORDINATES=" + propstr; } else { coord_str += "," + propstr; }
      }
      nfram++;
  }
  // Now setup action to compute distances between configurations
  std::string metric; unsigned nn = Path::getNumberOfFramesAndMetric( mtype, refname, metric );
  // Create action that computes the geometric path variables
  readInputLine( getShortcutLabel() + ": GEOMETRIC_PATH ARG=" + getShortcutLabel() + "_data " + coord_str + ref_str + metric ); 
}


}
}


