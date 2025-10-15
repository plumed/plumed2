/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
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
#ifndef __PLUMED_secondarystructure_SecondaryStructureBase_h
#define __PLUMED_secondarystructure_SecondaryStructureBase_h

#include "core/ActionWithVector.h"
#include "core/ActionShortcut.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/GenericMolInfo.h"
#include "core/ParallelTaskManager.h"
#include "tools/ColvarOutput.h"
#include <vector>

namespace PLMD {
namespace secondarystructure {

/// Base action for calculating things like AlphRMSD, AntibetaRMSD, etc
template <class T>
class SecondaryStructureBase: public ActionWithVector {
public:
  using input_type = T;
  using PTM = ParallelTaskManager<SecondaryStructureBase<T>>;
  static constexpr size_t virialSize=9;
  static constexpr unsigned customGatherStep=3;
  static constexpr unsigned customGatherStopBefore=virialSize;

private:
  PTM taskmanager;
public:
  static void registerKeywords( Keywords& keys );
  static void readBackboneAtoms( ActionShortcut* action, PlumedMain& plumed, const std::string& backnames, std::vector<unsigned>& chain_lengths, std::vector<std::string>& all_atoms );
  static bool readShortcutWords( std::string& ltmap, ActionShortcut* action );
  explicit SecondaryStructureBase(const ActionOptions&);
  unsigned getNumberOfDerivatives() override ;
  void calculate() override;
  void getInputData( std::vector<double>& inputdata ) const override ;
  void applyNonZeroRankForces( std::vector<double>& outforces ) override ;
  static void performTask( unsigned task_index, const T& actiondata, ParallelActionsInput& input, ParallelActionsOutput& output );
  static int getNumberOfValuesPerTask( std::size_t task_index, const T& actiondata );
  static void getForceIndices( std::size_t task_index, std::size_t colno, std::size_t ntotal_force, const T& actiondata, const ParallelActionsInput& input, ForceIndexHolder force_indices );
};

template <class T>
unsigned SecondaryStructureBase<T>::getNumberOfDerivatives() {
  return 3*getNumberOfAtoms()+virialSize;
}

template <class T>
bool SecondaryStructureBase<T>::readShortcutWords( std::string& ltmap, ActionShortcut* action ) {
  action->parse("LESS_THAN",ltmap);
  if( ltmap.length()==0 ) {
    std::string nn, mm, d_0, r_0;
    action->parse("R_0",r_0);
    if( r_0.length()==0 ) {
      r_0="0.08";
    }
    action->parse("NN",nn);
    action->parse("D_0",d_0);
    action->parse("MM",mm);
    ltmap = "RATIONAL R_0=" + r_0 + " D_0=" + d_0 + " NN=" + nn + " MM=" + mm;
    return false;
  }
  return true;
}

template <class T>
void SecondaryStructureBase<T>::registerKeywords( Keywords& keys ) {
  ActionWithVector::registerKeywords( keys );
  PTM::registerKeywords( keys );
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions");
  keys.addInputKeyword("optional","MASK","vector","a vector which is used to determine which elements of the secondary structure variable should be computed");
  keys.add("atoms","ATOMS","this is the full list of atoms that we are investigating");
  keys.add("numbered","SEGMENT","this is the lists of atoms in the segment that are being considered");
  if( keys.getDisplayName()=="SECONDARY_STRUCTURE_DRMSD" ) {
    keys.add("compulsory","BONDLENGTH","0.17","the length to use for bonds");
  }
  keys.add("numbered","STRUCTURE","the reference structure");
  if( keys.getDisplayName()=="SECONDARY_STRUCTURE_RMSD" ) {
    keys.add("compulsory","TYPE","OPTIMAL","the manner in which RMSD alignment is performed. Should be OPTIMAL or SIMPLE. "
             "For more details on the OPTIMAL and SIMPLE methods see \\ref RMSD.");
  } else if( keys.getDisplayName()!="SECONDARY_STRUCTURE_DRMSD" ) {
    keys.add("compulsory","TYPE","DRMSD","the manner in which RMSD alignment is performed. Should be OPTIMAL, SIMPLE or DRMSD. "
             "For more details on the OPTIMAL and SIMPLE methods see \\ref RMSD. For more details on the "
             "DRMSD method see \\ref DRMSD.");
  }
  keys.addFlag("VERBOSE",false,"write a more detailed output");
  keys.reset_style("VERBOSE","hidden");
  if( keys.getDisplayName()!="SECONDARY_STRUCTURE_DRMSD" && keys.getDisplayName()!="SECONDARY_STRUCTURE_RMSD" ) {
    keys.add("residues","RESIDUES","this command is used to specify the set of residues that could conceivably form part of the secondary structure. "
             "It is possible to use residues numbers as the various chains and residues should have been identified else using an instance of the "
             "\\ref MOLINFO action. If you wish to use all the residues from all the chains in your system you can do so by "
             "specifying all. Alternatively, if you wish to use a subset of the residues you can specify the particular residues "
             "you are interested in as a list of numbers. Please be aware that to form secondary structure elements your chain "
             "must contain at least N residues, where N is dependent on the particular secondary structure you are interested in. "
             "As such if you define portions of the chain with fewer than N residues the code will crash.");
    keys.add("optional","LESS_THAN","calculate the number of a residue segments that are within a certain target distance of this secondary structure type. "
             "This quantity is calculated using \\f$\\sum_i \\sigma(s_i)\\f$, where \\f$\\sigma(s)\\f$ is a \\ref switchingfunction.");
    keys.add("optional","R_0","The r_0 parameter of the switching function.");
    keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
    keys.add("compulsory","NN","8","The n parameter of the switching function");
    keys.add("compulsory","MM","12","The m parameter of the switching function");
  }
  keys.addFlag("ALIGN_STRANDS",false,"ensure that the two halves of a beta sheet are not broken by the periodic boundaries before doing alignment");
  keys.addOutputComponent("struct","default","scalar","the vectors containing the rmsd distances between the residues and each of the reference structures");
  keys.addOutputComponent("lessthan","default","scalar","the number blocks of residues that have an RMSD from the secondary structure that is less than the threshold");
  keys.needsAction("SECONDARY_STRUCTURE_RMSD");
  keys.needsAction("SECONDARY_STRUCTURE_DRMSD");
  keys.needsAction("LESS_THAN");
  keys.needsAction("SUM");
  keys.addDOI("10.1021/ct900202f");
}

template <class T>
void SecondaryStructureBase<T>::readBackboneAtoms( ActionShortcut* action, PlumedMain& plumed, const std::string& moltype, std::vector<unsigned>& chain_lengths, std::vector<std::string>& all_atoms ) {
  auto* moldat=plumed.getActionSet().selectLatest<GenericMolInfo*>(action);
  if( ! moldat ) {
    action->error("Unable to find MOLINFO in input");
  }

  std::vector<std::string> resstrings;
  action->parseVector( "RESIDUES", resstrings );
  if(resstrings.size()==0) {
    action->error("residues are not defined, check the keyword RESIDUES");
  } else if( Tools::caseInSensStringCompare(resstrings[0], "all") ) {
    resstrings[0]="all";
    action->log.printf("  examining all possible secondary structure combinations\n");
  } else {
    action->log.printf("  examining secondary structure in residue positions : %s ",resstrings[0].c_str() );
    for(unsigned i=1; i<resstrings.size(); ++i) {
      action->log.printf(", %s",resstrings[i].c_str() );
    }
    action->log.printf("\n");
  }
  std::vector< std::vector<AtomNumber> > backatoms;
  moldat->getBackbone( resstrings, moltype, backatoms );

  chain_lengths.resize( backatoms.size() );
  for(unsigned i=0; i<backatoms.size(); ++i) {
    chain_lengths[i]=backatoms[i].size();
    for(unsigned j=0; j<backatoms[i].size(); ++j) {
      std::string bat_str;
      Tools::convert( backatoms[i][j].serial(), bat_str );
      all_atoms.push_back( bat_str );
    }
  }
}

template <class T>
SecondaryStructureBase<T>::SecondaryStructureBase(const ActionOptions&ao):
  Action(ao),
  ActionWithVector(ao),
  taskmanager(this) {
  if( plumed.usingNaturalUnits() ) {
    error("cannot use this collective variable when using natural units");
  }

  input_type myinput;
  parseFlag("NOPBC",myinput.nopbc);
  std::string alignType="";
  if( getName()=="SECONDARY_STRUCTURE_RMSD" ) {
    parse("TYPE",alignType);
  }
  log.printf("  distances from secondary structure elements are calculated using %s algorithm\n", alignType.c_str() );
  log<<"  Bibliography "<<plumed.cite("Pietrucci and Laio, J. Chem. Theory Comput. 5, 2197 (2009)");
  log<<"\n";

  parseFlag("ALIGN_STRANDS",myinput.align_strands);
  bool verbose_output=false;
  parseFlag("VERBOSE",verbose_output);
  log.printf("  ensuring atoms 7 and 22 in each residue are not separated by pbc before doing alignment\n");

  // Read in the atoms
  std::vector<AtomNumber> all_atoms;
  parseAtomList("ATOMS",all_atoms);
  if( all_atoms.size()==0 ) {
    error("no atoms were specified -- use ATOMS");
  }
  requestAtoms( all_atoms );

  std::vector<std::vector<unsigned> > colvar_atoms;
  for(unsigned i=1;; ++i) {
    std::vector<unsigned> newatoms;
    if( !parseNumberedVector("SEGMENT",i,newatoms) ) {
      break;
    }
    if( verbose_output ) {
      log.printf("  Secondary structure segment %u contains atoms : ", static_cast<unsigned>(colvar_atoms.size()+1));
      for(unsigned ii=0; ii<newatoms.size(); ++ii) {
        log.printf("%d ",all_atoms[newatoms[ii]].serial() );
      }
      log.printf("\n");
    }
    colvar_atoms.push_back( newatoms );
  }
  if( colvar_atoms.size()==0 ) {
    error("did not find any SEGMENT keywords in input");
  }
  myinput.colvar_atoms.resize( colvar_atoms.size(), colvar_atoms[0].size() );
  for(unsigned i=0; i<colvar_atoms.size(); ++i) {
    for(unsigned j=0; j<colvar_atoms[i].size(); ++j) {
      myinput.colvar_atoms[i][j] = colvar_atoms[i][j];
    }
  }

  double bondlength=0.0;
  if( getName()=="SECONDARY_STRUCTURE_DRMSD" ) {
    parse("BONDLENGTH",bondlength);
    bondlength=bondlength/getUnits().getLength();
  }
  // Read in the reference structure
  for(unsigned ii=1;; ++ii) {
    std::vector<double> cstruct;
    if( !parseNumberedVector("STRUCTURE",ii,cstruct) ) {
      break ;
    }
    plumed_assert( cstruct.size()%3==0 && cstruct.size()/3==colvar_atoms[0].size() );
    std::vector<Vector> structure( cstruct.size()/3 );
    for(unsigned i=0; i<structure.size(); ++i) {
      for(unsigned j=0; j<3; ++j) {
        structure[i][j] = 0.1*cstruct[3*i+j]/getUnits().getLength();
      }
    }
    myinput.setReferenceStructure( alignType, bondlength, structure );
  }

  // And create values to hold everything
  plumed_assert( myinput.nstructures>0 );
  std::vector<std::size_t> shape(1);
  shape[0]=colvar_atoms.size();
  if( myinput.nstructures==1 ) {
    addValue( shape );
    setNotPeriodic();
  } else {
    std::string num;
    for(unsigned i=0; i<myinput.nstructures; ++i) {
      Tools::convert( i+1, num );
      addComponent( "struct-" + num, shape );
      componentIsNotPeriodic( "struct-" + num );
    }
  }
  for(unsigned i=0; i<getNumberOfComponents(); ++i) {
    getPntrToComponent(i)->setDerivativeIsZeroWhenValueIsZero();
  }
  myinput.nindices_per_task = colvar_atoms[0].size();
  taskmanager.setupParallelTaskManager( 3*colvar_atoms[0].size() + virialSize, 3*getNumberOfAtoms() + 9 );
  taskmanager.setActionInput( myinput );
}

template <class T>
void SecondaryStructureBase<T>::calculate() {
  taskmanager.runAllTasks();
}

template <class T>
void SecondaryStructureBase<T>::getInputData( std::vector<double>& inputdata ) const {
  if( inputdata.size()!=3*getNumberOfAtoms() ) {
    inputdata.resize( 3*getNumberOfAtoms() );
  }

  std::size_t k=0;
  for(unsigned i=0; i<getNumberOfAtoms(); ++i) {
    Vector mypos( getPosition(i) );
    inputdata[k] = mypos[0];
    k++;
    inputdata[k] = mypos[1];
    k++;
    inputdata[k] = mypos[2];
    k++;
  }
}

template <class T>
void SecondaryStructureBase<T>::performTask( unsigned task_index, const T& actiondata, ParallelActionsInput& input, ParallelActionsOutput& output ) {
  // std::vector<Vector> pos( actiondata.natoms );
  std::array<Vector,30> pos;

  for(unsigned i=0; i<actiondata.natoms; ++i) {
    const unsigned atno = actiondata.colvar_atoms(task_index,i);
    pos[i][0] = input.inputdata[3*atno+0];
    pos[i][1] = input.inputdata[3*atno+1];
    pos[i][2] = input.inputdata[3*atno+2];
  }
  // This aligns the two strands if this is required
  if( actiondata.align_strands ) {
    Vector distance=input.pbc->distance( pos[6],pos[21] );
    Vector origin_old, origin_new;
    origin_old=pos[21];
    origin_new=pos[6]+distance;
    for(unsigned i=15; i<30; ++i) {
      pos[i]+=( origin_new - origin_old );
    }
  } else if( !actiondata.nopbc ) {
    for(unsigned i=0; i<actiondata.natoms-1; ++i) {
      const Vector & first (pos[i]);
      Vector & second (pos[i+1]);
      second=first+input.pbc->distance(first,second);
    }
  }

  // Create a holder for the derivatives
  const unsigned rs = actiondata.nstructures;
  ColvarOutput rmsd_output( output.values,
                            3*pos.size()+virialSize,
                            output.derivatives.data() );
  // And now calculate the DRMSD
  for(unsigned i=0; i<rs; ++i) {
    T::calculateDistance( i, input.noderiv, actiondata, View{pos.data(),pos.size()}, rmsd_output );
  }
}

template <class T>
void SecondaryStructureBase<T>::applyNonZeroRankForces( std::vector<double>& outforces ) {
  taskmanager.applyForces( outforces );
}

template <class T>
int SecondaryStructureBase<T>::getNumberOfValuesPerTask( std::size_t task_index,
    const T& actiondata ) {
  return 1;
}

template <class T>
void SecondaryStructureBase<T>::getForceIndices( std::size_t task_index,
    std::size_t /* colno */,
    std::size_t ntotal_force,
    const T& actiondata,
    const ParallelActionsInput& input,
    ForceIndexHolder force_indices ) {
  for(unsigned i=0; i<input.ncomponents; ++i) {
    std::size_t m = 0;
    for(unsigned j=0; j<actiondata.nindices_per_task; ++j) {
      std::size_t base = 3*actiondata.colvar_atoms[task_index][j];
      force_indices.indices[i][m] = base + 0;
      ++m;
      force_indices.indices[i][m] = base + 1;
      ++m;
      force_indices.indices[i][m] = base + 2;
      ++m;
    }
    for(unsigned n=ntotal_force-virialSize; n<ntotal_force; ++n) {
      force_indices.indices[i][m] = n;
      ++m;
    }
    force_indices.threadsafe_derivatives_end[i] = 0;
    force_indices.tot_indices[i] = m;
  }
}

}
}

#endif
