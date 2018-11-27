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
#include "SetupReferenceBase.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/Atoms.h"

namespace PLMD {
namespace setup {

class CalculateReferenceValue : public SetupReferenceBase {
public:
  static void registerKeywords( Keywords& keys );
  explicit CalculateReferenceValue(const ActionOptions&ao); 
};

PLUMED_REGISTER_ACTION(CalculateReferenceValue,"CALCULATE_REFERENCE")

void CalculateReferenceValue::registerKeywords( Keywords& keys ) {
  SetupReferenceBase::registerKeywords( keys );
  keys.add("compulsory","CONFIG","the label of the READ_CONFIG command for which we are doing this calculation");
  keys.add("compulsory","INPUT","the file to use as input to PLUMED");
}

CalculateReferenceValue::CalculateReferenceValue(const ActionOptions&ao):
Action(ao),
SetupReferenceBase(ao)
{
   std::string lab; parse("CONFIG",lab); 
   SetupReferenceBase* as = plumed.getActionSet().selectWithLabel<SetupReferenceBase*>( lab );
   if( !as ) error("found no READ_CONFIG action with label " + lab );
   log.printf("  calculating reference values for positions read in by action with label %s \n", lab.c_str() );

   // Create a PlumedMain object to do the calculation
   PlumedMain p; int s=sizeof(double);
   p.cmd("setRealPrecision",&s);
   p.cmd("setNoVirial"); 
   p.cmd("setMDEngine","plumed");
   unsigned tatoms, targs; as->getNatomsAndNargs( tatoms, targs );
   int natoms = tatoms; p.cmd("setNatoms",&natoms);
   // Create copies of the values computed by the reference object
   for(unsigned i=0;i<as->getNumberOfComponents();++i) {
       std::vector<int> size(1+as->copyOutput(i)->getRank()); 
       size[0]=as->copyOutput(i)->getRank(); 
       for(unsigned j=0;j<size[0];++j) size[j+1]=as->copyOutput(i)->getShape()[j];
       p.cmd("createValue " + as->copyOutput(i)->getName(), &size[0] );
       if( !as->copyOutput(i)->isPeriodic() ) p.cmd("setValueNotPeriodic " + as->copyOutput(i)->getName());
   }
   double tstep=1.0; p.cmd("setTimestep",&tstep);
   // Now read the PLUMED command that we have to execute
   std::string inp; parse("INPUT",inp); const char* cinp=inp.c_str();
   std::vector<std::string> input=Tools::getWords(inp);
   if( input.size()==1 && !actionRegister().check(input[0]) ) {
       p.cmd("setPlumedDat",cinp); p.cmd("init");
   } else {
       p.cmd("init"); p.cmd("readInputLine",cinp);
   }
   // Setup the positions and masses using the indices from the read input
   int istep=0; p.cmd("setStep",&istep);
   std::vector<Vector> positions( natoms ), forces( natoms );
   std::vector<double> masses( natoms ), charges( natoms );
   as->getAtomsFromReference( 0, masses, charges, positions );
   p.cmd("setMasses",&masses[0]); if( atoms.chargesWereSet() ) p.cmd("setCharge",&charges[0]);
   p.cmd("setForces",&forces[0]); p.cmd("setPositions",&positions[0]);
   // Copy values from reference to PLUMED 
   for(unsigned i=0;i<as->getNumberOfComponents();++i) {
      unsigned nvals = as->copyOutput(i)->getSize();
      std::vector<double> valdata( nvals );
      for(unsigned j=0;j<nvals;++j) valdata[j] = as->copyOutput(i)->get(j); 
      p.cmd("setValue " + as->copyOutput(i)->getName(), &valdata[0] );
   }
   Tensor box( atoms.getPbc().getBox() ); p.cmd("setBox",&box[0][0]);
   // Now retrieve the final value
   ActionWithValue* fav = dynamic_cast<ActionWithValue*>( p.getActionSet()[p.getActionSet().size()-1].get() );
   if( !fav ) error("final value should calculate relevant value that you want as reference");
   if( fav->getNumberOfComponents()!=1 ) error("final action in input should have one component");
   std::string name = (fav->copyOutput(0))->getName();
   long rank; p.cmd("getDataRank " + name, &rank );
   if( rank==0 ) rank=1;
   std::vector<long> ishape( rank ); std::vector<unsigned> shape( rank ); 
   p.cmd("getDataShape " + name, &ishape[0] );
   unsigned nvals=1; for(unsigned i=0;i<shape.size();++i){ shape[i]=ishape[i]; nvals *= shape[i]; } 
   std::vector<double> data( nvals ); p.cmd("setMemoryForData " + name, &data[0] );
   // Do the calculation using the Plumed object
   p.cmd("calc"); addValue( shape ); getPntrToComponent(0)->buildDataStore( getLabel() );
   if( (fav->copyOutput(0))->isPeriodic() ) {
       std::string min, max; (fav->copyOutput(0))->getDomain( min, max ); 
       setPeriodic( min, max );
   } else setNotPeriodic();
   // And set it to what was calculated
   for(unsigned i=0;i<nvals;++i) getPntrToComponent(0)->set( i, data[i] ); 
   log.printf("  setup %d reference values \n", nvals);
   // Create a fake task list -- this ensures derivatives are looked at correctly
   for(unsigned i=0;i<shape[0];++i) addTaskToList( i );
}

}
}
