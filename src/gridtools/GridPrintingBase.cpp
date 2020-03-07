/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2020 The plumed team
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
#include "GridPrintingBase.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "vesselbase/ActionWithVessel.h"

namespace PLMD {
namespace gridtools {

void GridPrintingBase::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys ); ActionPilot::registerKeywords( keys );
  keys.add("compulsory","GRID","the action that creates the grid you would like to output");
  keys.add("compulsory","STRIDE","0","the frequency with which the grid should be output to the file.  The default "
           "value of 0 ensures that the grid is only output at the end of the trajectory");
  keys.add("compulsory","FILE","density","the file on which to write the grid.");
  keys.add("compulsory","REPLICA","0","the replica for which you would like to output this information");
  keys.add("optional","FMT","the format that should be used to output real numbers");
}

GridPrintingBase::GridPrintingBase(const ActionOptions&ao):
  Action(ao),
  ActionPilot(ao),
  fmt("%f"),
  output_for_all_replicas(false)
{
  std::string mlab; parse("GRID",mlab);
  vesselbase::ActionWithVessel* mves= plumed.getActionSet().selectWithLabel<vesselbase::ActionWithVessel*>(mlab);
  if(!mves) error("action labelled " +  mlab + " does not exist or does not have vessels");
  addDependency(mves);

  for(unsigned i=0; i<mves->getNumberOfVessels(); ++i) {
    ingrid=dynamic_cast<GridVessel*>( mves->getPntrToVessel(i) );
    if( ingrid ) break;
  }
  if( !ingrid ) error("input action does not calculate a grid");

  parse("FILE",filename);
  if(filename.length()==0) error("name out output file was not specified");
  log.printf("  outputting grid calculated by action %s to file named %s",mves->getLabel().c_str(), filename.c_str() );
  if( keywords.exists("FMT") ) {
    parse("FMT",fmt); log.printf(" with format %s \n", fmt.c_str() );
  } else {
    log.printf("\n");
  }
  std::string rep_data; parse("REPLICA",rep_data);
  if( rep_data=="all" ) output_for_all_replicas=true;
  else { preps.resize(1); Tools::convert( rep_data, preps[0] ); }
  if( output_for_all_replicas ) log.printf("  outputting files for all replicas \n");
  else {
    log.printf("  outputting data for replicas ");
    for(unsigned i=0; i<preps.size(); ++i) log.printf("%d ", preps[i] );
  }
}

void GridPrintingBase::update() {
  if( !output_for_all_replicas ) {
    bool found=false; unsigned myrep=plumed.multi_sim_comm.Get_rank();
    for(unsigned i=0; i<preps.size(); ++i) {
      if( myrep==preps[i] ) { found=true; break; }
    }
    if( !found ) return;
  }
  if( getStep()==0 || getStride()==0 ) return ;

  OFile ofile; ofile.link(*this);
  ofile.setBackupString("analysis");
  ofile.open( filename ); printGrid( ofile );
}

void GridPrintingBase::runFinalJobs() {
  if( !output_for_all_replicas ) {
    bool found=false; unsigned myrep=plumed.multi_sim_comm.Get_rank();
    for(unsigned i=0; i<preps.size(); ++i) {
      if( myrep==preps[i] ) { found=true; break; }
    }
    if( !found ) return;
  }
  if( getStride()>0 ) return;

  OFile ofile; ofile.link(*this);
  ofile.open( filename ); printGrid( ofile );
}

}
}
