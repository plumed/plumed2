/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2016 The plumed team
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
#include "core/ActionRegister.h"
#include "gridtools/ActionWithGrid.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "MultiColvarBase.h"

using namespace std;

namespace PLMD
{
namespace multicolvar {

//+PLUMEDOC ANALYSIS MULTIHISTOGRAM
/*
Evaluate the histogram for a particular multicolvar

\par Examples


*/
//+ENDPLUMEDOC

class MultiColvarHistogram : public gridtools::ActionWithGrid {
private:
  MultiColvarBase* mycolv; 
  vesselbase::StoreDataVessel* stash;
public:
  explicit MultiColvarHistogram(const ActionOptions&);
  static void registerKeywords( Keywords& keys );
  bool isPeriodic(){ return false; }
  void compute( const unsigned& current, MultiValue& myvals ) const ;
};

PLUMED_REGISTER_ACTION(MultiColvarHistogram,"MULTICOLVARHISTOGRAM")

void MultiColvarHistogram::registerKeywords( Keywords& keys ){
  gridtools::ActionWithGrid::registerKeywords( keys );
  keys.add("compulsory","GRID_MIN","the lower bounds for the grid");
  keys.add("compulsory","GRID_MAX","the upper bounds for the grid");
  keys.add("optional","GRID_BIN","the number of bins for the grid");
  keys.add("optional","GRID_SPACING","the approximate grid spacing (to be used as an alternative or together with GRID_BIN)");
}

MultiColvarHistogram::MultiColvarHistogram(const ActionOptions&ao):
  Action(ao),
  ActionWithGrid(ao)
{
  std::string mlab; parse("DATA",mlab);
  mycolv = plumed.getActionSet().selectWithLabel<MultiColvarBase*>(mlab);
  if(!mycolv) error("action labelled " +  mlab + " does not exist or is not a MultiColvar");
  stash = mycolv->buildDataStashes( NULL );

  // Read stuff for grid
  std::vector<std::string> gmin( getNumberOfArguments() ), gmax( getNumberOfArguments() );
  parseVector("GRID_MIN",gmin); parseVector("GRID_MAX",gmax);
  std::vector<unsigned> nbin; parseVector("GRID_BIN",nbin);
  std::vector<double> gspacing; parseVector("GRID_SPACING",gspacing);
  if( nbin.size()!=getNumberOfArguments() && gspacing.size()!=getNumberOfArguments() ){
      error("GRID_BIN or GRID_SPACING must be set");
  }
  
  // Input of name and labels
  std::string vstring="COMPONENTS=" + getLabel();
  vstring += " COORDINATES=" + getPntrToArgument(0)->getName();
  for(unsigned i=1;i<getNumberOfArguments();++i) vstring += "," + getPntrToArgument(i)->getName();
  // Input for PBC
  if( getPntrToArgument(0)->isPeriodic() ) vstring+=" PBC=T";
  else vstring+=" PBC=F";
  for(unsigned i=1;i<getNumberOfArguments();++i){
     if( getPntrToArgument(i)->isPeriodic() ) vstring+=",T";
     else vstring+=",F";
  }
  // And create the grid
  createGrid( "histogram", vstring );
  mygrid->setBounds( gmin, gmax, nbin, gspacing );
  // And finish the grid setup
  setAveragingAction( mygrid, false );

  // Create a task list
  for(unsigned i=0;i<mycolv->getFullNumberOfTasks();++i) addTaskToList(i);
  setAveragingAction( mygrid, false ); checkRead();
  log.printf("  for colvars calculated by action %s \n",mycolv->getLabel().c_str() );

  checkRead(); 
  // Add the dependency
  addDependency( mycolv );
}

void MultiColvarHistogram::compute( const unsigned& current, MultiValue& myvals ) const {
  std::vector<double> cvals( mycolv->getNumberOfQuantities() ); 
  stash->retrieveSequentialValue( current, false, cvals );
  myvals.setValue( 0, cvals[0] ); myvals.setValue( 1, cvals[1] );
}

}
}
