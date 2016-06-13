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

//+PLUMEDOC ANALYSIS MULTICOLVARHISTOGRAM
/*
Evaluate the histogram for a particular multicolvar

\par Examples


*/
//+ENDPLUMEDOC

class MultiColvarHistogram : public gridtools::ActionWithGrid {
private:
  double ww;
  MultiColvarBase* mycolv; 
  vesselbase::StoreDataVessel* stash;
public:
  explicit MultiColvarHistogram(const ActionOptions&);
  static void registerKeywords( Keywords& keys );
  unsigned getNumberOfQuantities() const ;
  bool isPeriodic(){ return false; }
  void prepareForAveraging();
  void compute( const unsigned& current, MultiValue& myvals ) const ;
};

PLUMED_REGISTER_ACTION(MultiColvarHistogram,"MULTICOLVARHISTOGRAM")

void MultiColvarHistogram::registerKeywords( Keywords& keys ){
  gridtools::ActionWithGrid::registerKeywords( keys );
  keys.add("compulsory","DATA","the multicolvar which you would like to calculate a histogram for");
  keys.add("compulsory","GRID_MIN","the lower bounds for the grid");
  keys.add("compulsory","GRID_MAX","the upper bounds for the grid");
  keys.add("optional","GRID_BIN","the number of bins for the grid");
  keys.add("optional","GRID_SPACING","the approximate grid spacing (to be used as an alternative or together with GRID_BIN)");
}

MultiColvarHistogram::MultiColvarHistogram(const ActionOptions&ao):
  Action(ao),
  ActionWithGrid(ao),
  ww(1.0)
{
  std::string mlab; parse("DATA",mlab);
  mycolv = plumed.getActionSet().selectWithLabel<MultiColvarBase*>(mlab);
  if(!mycolv) error("action labelled " +  mlab + " does not exist or is not a MultiColvar");
  stash = mycolv->buildDataStashes( NULL );

  // Read stuff for grid
  std::vector<std::string> gmin( 1 ), gmax( 1 );
  parseVector("GRID_MIN",gmin); parseVector("GRID_MAX",gmax);
  std::vector<unsigned> nbin; parseVector("GRID_BIN",nbin);
  std::vector<double> gspacing; parseVector("GRID_SPACING",gspacing);
  if( nbin.size()!=1 && gspacing.size()!=1 ){
      error("GRID_BIN or GRID_SPACING must be set");
  }
  
  // Input of name and labels
  std::string vstring="COMPONENTS=" + getLabel();
  vstring += " COORDINATES=" + mycolv->getLabel();
  // Input for PBC
  if( mycolv->isPeriodic() ) vstring+=" PBC=T";
  else vstring+=" PBC=F";
  // And create the grid
  createGrid( "histogram", vstring );
  mygrid->setBounds( gmin, gmax, nbin, gspacing );

  // Create a task list
  for(unsigned i=0;i<mycolv->getFullNumberOfTasks();++i) addTaskToList(i);
  // And finish the grid setup
  setAveragingAction( mygrid, true ); checkRead();
  log.printf("  for colvars calculated by action %s \n",mycolv->getLabel().c_str() );

  checkRead(); 
  // Add the dependency
  addDependency( mycolv );
}

void MultiColvarHistogram::prepareForAveraging(){
  deactivateAllTasks(); double norm=0; 
  std::vector<double> cvals( mycolv->getNumberOfQuantities() );
  for(unsigned i=0;i<stash->getNumberOfStoredValues();++i){
      taskFlags[i]=1; stash->retrieveSequentialValue(i, false, cvals ); 
      norm += cvals[0];
  }
  lockContributors(); ww = 1.0 / norm;
}

unsigned MultiColvarHistogram::getNumberOfQuantities() const {
  return 3;
}

void MultiColvarHistogram::compute( const unsigned& current, MultiValue& myvals ) const {
  std::vector<double> cvals( mycolv->getNumberOfQuantities() ); 
  stash->retrieveSequentialValue( current, false, cvals );
  myvals.setValue( 0, cvals[0] ); myvals.setValue( 1, cvals[1] ); myvals.setValue( 2, ww );
}

}
}
