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
#include "ReadReferenceConfiguration.h"
#include "core/ActionWithValue.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/Atoms.h"

namespace PLMD {
namespace setup {

class CalculateReferenceValue : 
public ActionSetup,
public ActionWithValue,
public ActionAtomistic 
{
public:
  static void registerKeywords( Keywords& keys );
  explicit CalculateReferenceValue(const ActionOptions&ao); 
  void clearDerivatives( const bool& force=false ) {}
  unsigned getNumberOfDerivatives() const { return 0; }
};

PLUMED_REGISTER_ACTION(CalculateReferenceValue,"CALCULATE_REFERENCE")

void CalculateReferenceValue::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys ); ActionAtomistic::registerKeywords( keys ); ActionWithValue::registerKeywords( keys );
  keys.add("atoms","ATOMS","the label of the READ_ATOMS command for which we are doing this calculation");
  keys.add("compulsory","INPUT","the file to use as input to PLUMED");
}

CalculateReferenceValue::CalculateReferenceValue(const ActionOptions&ao):
Action(ao),
ActionSetup(ao),
ActionWithValue(ao),
ActionAtomistic(ao)
{
   std::string lab; std::vector<AtomNumber> atom_list; parseAtomList("ATOMS",atom_list); bool safe=true;
   for(unsigned i=0;i<atom_list.size();++i) {
       if( atoms.isVirtualAtom(atom_list[i]) ) {
           ReadReferenceConfiguration* as = dynamic_cast<ReadReferenceConfiguration*>( atoms.getVirtualAtomsAction(atom_list[i]) );
           if( !as ) safe=false;
           // Now check that we are only referencing one reference configuration
           if( i==0 ) lab=as->getLabel(); else if( lab!=as->getLabel() ) safe=false;
       } else safe=false;
   }
   if( !safe ) error("input must be to reference atom positions that are read in using READ_ATOMS"); 
   // And request the atoms
   requestAtoms( atom_list );
   // Create a PlumedMain object to do the calculation
   PlumedMain p; int s=sizeof(double);
   p.cmd("setRealPrecision",&s);
   p.cmd("setNoVirial"); 
   p.cmd("setMDEngine","plumed");
   ReadReferenceConfiguration* as = dynamic_cast<ReadReferenceConfiguration*>( atoms.getVirtualAtomsAction(atom_list[0]) ); 
   // Number of atoms is largest serial in input
   int natoms=as->myindices[0].serial();
   for(unsigned i=1;i<as->myindices.size();++i) {
       if( as->myindices[i].serial()>natoms ) natoms = as->myindices[i].serial();
   }
   p.cmd("setNatoms",&natoms);
   double tstep=1.0; p.cmd("setTimestep",&tstep);
   // Now read the PLUMED command that we have to execute
   std::string inp; parse("INPUT",inp); const char* cinp=inp.c_str();
   std::vector<std::string> input=Tools::getWords(inp);
   if( input.size()==1 ) {
       p.cmd("setPlumedDat",cinp); p.cmd("init");
   } else {
       p.cmd("init"); p.cmd("readInputLine",cinp);
   }
   // Setup the positions and masses using the indices from the read input
   int istep=0; p.cmd("setStep",&istep);
   retrieveAtoms(); std::vector<Vector> positions( natoms );
   std::vector<double> masses( natoms ), charges( natoms );
   for(unsigned i=0;i<as->myindices.size();++i) {
       masses[as->myindices[i].index()]=getMass(i); 
       if( atoms.chargesWereSet() ) charges[as->myindices[i].index()]=getCharge(i); 
       positions[as->myindices[i].index()]=getPosition(i);
   }
   p.cmd("setMasses",&masses[0]); 
   if( atoms.chargesWereSet() ) p.cmd("setCharge",&charges[0]);
   std::vector<Vector> forces( natoms ); p.cmd("setForces",&forces[0]);
   p.cmd("setPositions",&positions[0]);
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
