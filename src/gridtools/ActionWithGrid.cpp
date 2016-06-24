/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016 The plumed team
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
#include "ActionWithGrid.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

namespace PLMD {
namespace gridtools {

void ActionWithGrid::registerKeywords( Keywords& keys ){
  vesselbase::ActionWithAveraging::registerKeywords( keys );
  keys.add("compulsory","BANDWIDTH","the bandwidths for kernel density esimtation");
  keys.add("compulsory","KERNEL","gaussian","the kernel function you are using.  More details on  the kernels available "
                                            "in plumed plumed can be found in \\ref kernelfunctions.");
}

ActionWithGrid::ActionWithGrid( const ActionOptions& ao):
Action(ao),
ActionWithAveraging(ao),
mygrid(NULL)
{
}

void ActionWithGrid::createGrid( const std::string& type, const std::string& inputstr ){
  // Start creating the input for the grid
  std::string vstring = inputstr; 
  if( keywords.exists("KERNEL") ){
      std::string kstring; parse("KERNEL",kstring);
      if( kstring=="DISCRETE" ) vstring += " KERNEL=" + kstring;
      else vstring += " KERNEL=" + kstring + " " + getKeyword("BANDWIDTH");
  }

  vesselbase::VesselOptions da("mygrid","",-1,vstring,this);
  Keywords keys; gridtools::AverageOnGrid::registerKeywords( keys );
  vesselbase::VesselOptions dar( da, keys );
  if( type=="histogram" ){
     mygrid = new HistogramOnGrid(dar); 
  } else if( type=="average" ){
     mygrid = new AverageOnGrid(dar); 
  } else if( type=="grid" ){
     mygrid = new GridVessel(dar); 
  } else {
     plumed_merror("no way to create grid of type " + type );
  } 
}

void ActionWithGrid::performTask( const unsigned& task_index, const unsigned& current, MultiValue& myvals ) const {
  // Set the weight of this point
  myvals.setValue( 0, cweight ); compute( current, myvals );
}

}
}
