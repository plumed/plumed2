/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2015 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#include "SecondaryStructureRMSD.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/SetupMolInfo.h"
#include "core/Atoms.h"
#include "vesselbase/Vessel.h"
#include "reference/MetricRegister.h"
#include "reference/SingleDomainRMSD.h"

namespace PLMD {
namespace secondarystructure{

void SecondaryStructureRMSD::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  keys.add("residues","RESIDUES","this command is used to specify the set of residues that could conceivably form part of the secondary structure. "
                              "It is possible to use residues numbers as the various chains and residues should have been identified else using an instance of the "
                              "\\ref MOLINFO action. If you wish to use all the residues from all the chains in your system you can do so by "
                              "specifying all. Alternatively, if you wish to use a subset of the residues you can specify the particular residues "
                              "you are interested in as a list of numbers. Please be aware that to form secondary structure elements your chain "
                              "must contain at least N residues, where N is dependent on the particular secondary structure you are interested in. "
                              "As such if you define portions of the chain with fewer than N residues the code will crash."); 
  keys.add("compulsory","TYPE","DRMSD","the manner in which RMSD alignment is performed. Should be OPTIMAL, SIMPLE or DRMSD. "
                                       "For more details on the OPTIMAL and SIMPLE methods see \\ref RMSD. For more details on the "
                                       "DRMSD method see \\ref DRMSD.");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function.");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","NN","8","The n parameter of the switching function");
  keys.add("compulsory","MM","12","The m parameter of the switching function");
  keys.reserve("optional","STRANDS_CUTOFF","If in a segment of protein the two strands are further apart then the calculation "
                                       "of the actual RMSD is skipped as the structure is very far from being beta-sheet like. "
                                       "This keyword speeds up the calculation enormously when you are using the LESS_THAN option. "
                                       "However, if you are using some other option, then this cannot be used");
  keys.addFlag("VERBOSE",false,"write a more detailed output");
  keys.add("hidden","NL_STRIDE","the frequency with which the neighbor list should be updated. Between neighbour list update steps all quantities "
                                  "that contributed less than TOL at the previous neighbor list update step are ignored.");
  ActionWithVessel::registerKeywords( keys );
  keys.use("LESS_THAN"); keys.use("MIN"); keys.use("NL_TOL");
  keys.setComponentsIntroduction("By default this Action calculates the number of structural units that are within a certain "
                                 "distance of a idealised secondary structure element. This quantity can then be referenced "
                                 "elsewhere in the input by using the label of the action. However, thes Action can also be used to "
                                 "calculate the following quantities by using the keywords as described below.  The quantities then "
                                 "calculated can be referened using the label of the action followed by a dot and then the name "
                                 "from the table below.  Please note that you can use the LESS_THAN keyword more than once.  The resulting "
                                 "components will be labelled <em>label</em>.lessthan-1, <em>label</em>.lessthan-2 and so on unless you "
                                 "exploit the fact that these labels are customizable. In particular, by using the LABEL keyword in the "
                                 "description of you LESS_THAN function you can set name of the component that you are calculating"); 
}

SecondaryStructureRMSD::SecondaryStructureRMSD(const ActionOptions&ao):
Action(ao),
ActionAtomistic(ao),
ActionWithValue(ao),
ActionWithVessel(ao),
updateFreq(0),
align_strands(false),
s_cutoff(0),
align_atom_1(0),
align_atom_2(0)
{
  parse("TYPE",alignType);
  log.printf("  distances from secondary structure elements are calculated using %s algorithm\n",alignType.c_str() );
  log<<"  Bibliography "<<plumed.cite("Pietrucci and Laio, J. Chem. Theory Comput. 5, 2197 (2009)"); log<<"\n";

  parseFlag("VERBOSE",verbose_output);
  if( keywords.exists("NL_STRIDE") ) parse("NL_STRIDE",updateFreq);
  if(updateFreq>0){ 
     firsttime=true; 
     log.printf("  Updating contributors every %d steps.\n",updateFreq);
  } else {
     firsttime=false; contributorsAreUnlocked=true;   // This will lock during first prepare step methinks
     log.printf("  Updating contributors every step.\n");
  }

  if( keywords.exists("STRANDS_CUTOFF") ){
    parse("STRANDS_CUTOFF",s_cutoff); align_strands=true;
    if( s_cutoff>0) log.printf("  ignoring contributions from strands that are more than %f apart\n",s_cutoff);
  }
}

SecondaryStructureRMSD::~SecondaryStructureRMSD(){
  for(unsigned i=0;i<references.size();++i) delete references[i];
}

void SecondaryStructureRMSD::turnOnDerivatives(){
  ActionWithValue::turnOnDerivatives();
  needsDerivatives();
}

void SecondaryStructureRMSD::setAtomsFromStrands( const unsigned& atom1, const unsigned& atom2 ){
  align_atom_1=atom1; align_atom_2=atom2; 
}

void SecondaryStructureRMSD::readBackboneAtoms( const std::string& moltype, std::vector<unsigned>& chain_lengths ){
  std::vector<SetupMolInfo*> moldat=plumed.getActionSet().select<SetupMolInfo*>();
  if( moldat.size()==0 ) error("Unable to find MOLINFO in input");

  std::vector<std::string> resstrings; parseVector( "RESIDUES", resstrings );
  if( !verbose_output ){
      if(resstrings[0]=="all"){
         log.printf("  examining all possible secondary structure combinations\n");
      } else {
         log.printf("  examining secondary struture in residue poritions : %s \n",resstrings[0].c_str() );
         for(unsigned i=1;i<resstrings.size();++i) log.printf(", %s",resstrings[i].c_str() );
         log.printf("\n");
      }
  }
  std::vector< std::vector<AtomNumber> > backatoms;
  moldat[0]->getBackbone( resstrings, moltype, backatoms );

  chain_lengths.resize( backatoms.size() );
  for(unsigned i=0;i<backatoms.size();++i){
     chain_lengths[i]=backatoms[i].size();
     for(unsigned j=0;j<backatoms[i].size();++j) all_atoms.addIndexToList( backatoms[i][j] );
  }
}

void SecondaryStructureRMSD::addColvar( const std::vector<unsigned>& newatoms ){
  if( colvar_atoms.size()>0 ) plumed_assert( colvar_atoms[0].size()==newatoms.size() );
  if( verbose_output ){
     log.printf("  Secondary structure segment %u contains atoms : ", static_cast<unsigned>(colvar_atoms.size()+1));
     for(unsigned i=0;i<newatoms.size();++i) log.printf("%d ",all_atoms(newatoms[i]).serial() );
     log.printf("\n");
  }
  addTaskToList( colvar_atoms.size() );
  colvar_atoms.push_back( newatoms );
}

void SecondaryStructureRMSD::setSecondaryStructure( std::vector<Vector>& structure, double bondlength, double units ){
  // If we are in natural units get conversion factor from nm into natural length units
  if( plumed.getAtoms().usingNaturalUnits() ){
      error("cannot use this collective variable when using natural units");
  }
  plumed_massert( !(align_strands && align_atom_1==0 && align_atom_2==0), "you must use setAtomsFromStrands with strands cutoff"); 

  // Convert into correct units
  for(unsigned i=0;i<structure.size();++i){
     structure[i][0]*=units; structure[i][1]*=units; structure[i][2]*=units;
  }

  if( references.size()==0 ){
     pos.resize( structure.size() ); 
     finishTaskListUpdate();

     readVesselKeywords();
     if( getNumberOfVessels()==0 ){
         double r0; parse("R_0",r0); double d0; parse("D_0",d0);
         int nn; parse("NN",nn); int mm; parse("MM",mm);
         std::ostringstream ostr;
         ostr<<"RATIONAL R_0="<<r0<<" D_0="<<d0<<" NN="<<nn<<" MM="<<mm;
         std::string input=ostr.str(); addVessel( "LESS_THAN", input, -1 ); // -1 here means that this value will be named getLabel()
         readVesselKeywords();  // This makes sure resizing is done
     } 
  }

  // Set the reference structure
  references.push_back( metricRegister().create<SingleDomainRMSD>( alignType ) ); 
  unsigned nn=references.size()-1;
  std::vector<double> align( structure.size(), 1.0 ), displace( structure.size(), 1.0 );
  references[nn]->setBoundsOnDistances( true , bondlength );  // We always use pbc
  references[nn]->setReferenceAtoms( structure, align, displace );
  references[nn]->setNumberOfAtoms( structure.size() );
}

void SecondaryStructureRMSD::prepare(){
  if( contributorsAreUnlocked ) lockContributors();
  if( updateFreq>0 ){
      if( firsttime || getStep()%updateFreq==0 ){ firsttime=false; unlockContributors(); }
  }
}

void SecondaryStructureRMSD::finishTaskListUpdate(){
  all_atoms.deactivateAll();
  for(unsigned i=0;i<getCurrentNumberOfActiveTasks();++i){
      for(unsigned j=0;j<colvar_atoms[getActiveTask(i)].size();++j) all_atoms.activate( colvar_atoms[getActiveTask(i)][j] );
  }
  all_atoms.updateActiveMembers();
  ActionAtomistic::requestAtoms( all_atoms.retrieveActiveList() );
  forcesToApply.resize( getNumberOfDerivatives() );
}

void SecondaryStructureRMSD::calculate(){
  runAllTasks();
}

void SecondaryStructureRMSD::performTask(){
  // Retrieve the positions
  for(unsigned i=0;i<pos.size();++i) pos[i]=ActionAtomistic::getPosition( getAtomIndex(i) );

  // This does strands cutoff
  Vector distance=pbcDistance( pos[align_atom_1],pos[align_atom_2] ); 
  if( s_cutoff>0 ){
     if( distance.modulo()>s_cutoff ){
       setElementValue(1,0.0);
       return;
     }
  }

  // This aligns the two strands if this is required
  if( alignType!="DRMSD" && align_strands ){
     Vector origin_old, origin_new; origin_old=pos[align_atom_2];
     origin_new=pos[align_atom_1]+distance;
     for(unsigned i=15;i<30;++i){
         pos[i]+=( origin_new - origin_old );
     }
  } 

  // And now calculate the RMSD
  double r,nr; const Pbc& pbc=getPbc(); 
  closest=0; r=references[0]->calculate( pos, pbc, false );
  for(unsigned i=1;i<references.size();++i){
      nr=references[i]->calculate( pos, pbc, false );
      if( nr<r ){ closest=i; r=nr; }
  }
  setElementValue(1,1.0); 
  setElementValue(0,r);
  return;
}

void SecondaryStructureRMSD::mergeDerivatives( const unsigned& ider, const double& df ){
  plumed_dbg_assert( ider==0 );
  for(unsigned i=0;i<colvar_atoms[getCurrentTask()].size();++i){
     unsigned thisatom=getAtomIndex(i), thispos=3*thisatom; 
     Vector ader=references[closest]->getAtomDerivative(i);
     accumulateDerivative( thispos, df*ader[0] ); thispos++;
     accumulateDerivative( thispos, df*ader[1] ); thispos++;
     accumulateDerivative( thispos, df*ader[2] ); 
  }
  Tensor virial;
  if( !references[closest]->getVirial( virial ) ){ 
     virial.zero();
     for(unsigned i=0;i<colvar_atoms[getCurrentTask()].size();++i){
         virial+=(-1.0*Tensor( pos[i], references[closest]->getAtomDerivative(i) ));
     }
  } 

  // Easy to merge the virial
  unsigned outnat=3*getNumberOfAtoms();
  accumulateDerivative( outnat, df*virial(0,0) ); outnat++;
  accumulateDerivative( outnat, df*virial(0,1) ); outnat++;
  accumulateDerivative( outnat, df*virial(0,2) ); outnat++;
  accumulateDerivative( outnat, df*virial(1,0) ); outnat++;
  accumulateDerivative( outnat, df*virial(1,1) ); outnat++;
  accumulateDerivative( outnat, df*virial(1,2) ); outnat++;
  accumulateDerivative( outnat, df*virial(2,0) ); outnat++;
  accumulateDerivative( outnat, df*virial(2,1) ); outnat++;
  accumulateDerivative( outnat, df*virial(2,2) );
}

void SecondaryStructureRMSD::apply(){
  if( getForcesFromVessels( forcesToApply ) ) setForcesOnAtoms( forcesToApply );
}

void SecondaryStructureRMSD::clearDerivativesAfterTask( const unsigned& ival ){
  thisval_wasset[ival]=false; setElementValue( ival, 0.0 ); thisval_wasset[ival]=false;
}

}
}
