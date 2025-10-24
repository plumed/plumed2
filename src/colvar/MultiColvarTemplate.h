/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#ifndef __PLUMED_colvar_MultiColvarTemplate_h
#define __PLUMED_colvar_MultiColvarTemplate_h

#include "core/ActionWithVector.h"
#include "core/ParallelTaskManager.h"
#include "tools/ColvarOutput.h"
#include "ColvarInput.h"

namespace PLMD {

class Colvar;

namespace colvar {

struct MultiColvarInput {
  bool usepbc;
  unsigned mode;
  unsigned nindices_per_task;
  //this is so simple that it better to use the Rule of 0
  //and the openacc helpers
  void toACCDevice()const {
#pragma acc enter data copyin(this[0:1], usepbc, mode, nindices_per_task)
  }
  void removeFromACCDevice() const  {
#pragma acc exit data delete(nindices_per_task, mode, usepbc, this[0:1])
  }
};

template <class T>
class MultiColvarTemplate : public ActionWithVector {
public:
  using input_type = MultiColvarInput;
  using PTM = ParallelTaskManager<MultiColvarTemplate<T>>;
  constexpr static size_t virialSize = 9;
private:
/// The parallel task manager
  PTM taskmanager;
/// An index that decides what we are calculating
  unsigned mode;
/// Are we using pbc to calculate the CVs
  bool usepbc;
/// Do we reassemble the molecule
  bool wholemolecules;
/// The number of atoms per task
  unsigned natoms_per_task;
public:
  static void registerKeywords(Keywords&);
  explicit MultiColvarTemplate(const ActionOptions&);
  unsigned getNumberOfDerivatives() override ;
  void addValueWithDerivatives( const std::vector<std::size_t>& shape=std::vector<std::size_t>() ) override ;
  void addComponentWithDerivatives( const std::string& name, const std::vector<std::size_t>& shape=std::vector<std::size_t>() ) override ;
  void getInputData( std::vector<double>& inputdata ) const override ;
  void calculate() override;
  void applyNonZeroRankForces( std::vector<double>& outforces ) override ;
  static void performTask( std::size_t task_index,
                           const MultiColvarInput& actiondata,
                           ParallelActionsInput& input,
                           ParallelActionsOutput& output );
  static int getNumberOfValuesPerTask( std::size_t task_index, const MultiColvarInput& actiondata );
  static void getForceIndices( std::size_t task_index,
                               std::size_t colno,
                               std::size_t ntotal_force,
                               const MultiColvarInput& actiondata,
                               const ParallelActionsInput& input,
                               ForceIndexHolder force_indices );
};

template <class T>
void MultiColvarTemplate<T>::registerKeywords(Keywords& keys ) {
  T::registerKeywords( keys );
  PTM::registerKeywords( keys );
  keys.addInputKeyword("optional","MASK","vector","the label for a sparse vector that should be used to determine which elements of the vector should be computed");
  unsigned nkeys = keys.size();
  for(unsigned i=0; i<nkeys; ++i) {
    if( keys.style( keys.getKeyword(i), "atoms" ) ) {
      keys.reset_style( keys.getKeyword(i), "numbered" );
    }
  }
  if( keys.outputComponentExists(".#!value") ) {
    keys.setValueDescription("vector","the " + keys.getDisplayName() + " for each set of specified atoms");
  }
}

template <class T>
MultiColvarTemplate<T>::MultiColvarTemplate(const ActionOptions&ao):
  Action(ao),
  ActionWithVector(ao),
  taskmanager(this),
  mode(0),
  usepbc(true),
  wholemolecules(false) {
  std::vector<AtomNumber> all_atoms;
  if( getName()=="POSITION_VECTOR" || getName()=="MASS_VECTOR" || getName()=="CHARGE_VECTOR" ) {
    parseAtomList( "ATOMS", all_atoms );
  }
  if( all_atoms.size()>0 ) {
    natoms_per_task=1;
  } else {
    std::vector<AtomNumber> t;
    for(int i=1;; ++i ) {
      T::parseAtomList( i, t, this );
      if( t.empty() ) {
        break;
      }

      if( i==1 ) {
        natoms_per_task=t.size();
      }
      if( t.size()!=natoms_per_task ) {
        std::string ss;
        Tools::convert(i,ss);
        error("ATOMS" + ss + " keyword has the wrong number of atoms");
      }
      for(unsigned j=0; j<natoms_per_task; ++j) {
        all_atoms.push_back( t[j] );
      }
      t.resize(0);
    }
  }
  if( all_atoms.size()==0 ) {
    error("No atoms have been specified");
  }
  requestAtoms(all_atoms);
  if( keywords.exists("NOPBC") ) {
    bool nopbc=!usepbc;
    parseFlag("NOPBC",nopbc);
    usepbc=!nopbc;
  }
  if( keywords.exists("WHOLEMOLECULES") ) {
    parseFlag("WHOLEMOLECULES",wholemolecules);
    if( wholemolecules ) {
      usepbc=false;
    }
  }
  if( usepbc ) {
    log.printf("  using periodic boundary conditions\n");
  } else {
    log.printf("  without periodic boundary conditions\n");
  }

  // Setup the values
  mode = T::getModeAndSetupValues( this );
  // This sets up an array in the parallel task manager to hold all the indices
  // Sets up the index list in the task manager
  taskmanager.setupParallelTaskManager( 3*natoms_per_task + virialSize, virialSize );
  taskmanager.setActionInput( MultiColvarInput{ usepbc, mode, natoms_per_task });
}

template <class T>
unsigned MultiColvarTemplate<T>::getNumberOfDerivatives() {
  return 3*getNumberOfAtoms()+9;
}

template <class T>
void MultiColvarTemplate<T>::calculate() {
  if( wholemolecules ) {
    makeWhole();
  }
  taskmanager.runAllTasks();
}

template <class T>
void MultiColvarTemplate<T>::applyNonZeroRankForces( std::vector<double>& outforces ) {
  taskmanager.applyForces( outforces );
}

template <class T>
void MultiColvarTemplate<T>::addValueWithDerivatives( const std::vector<std::size_t>& shape ) {
  std::vector<std::size_t> s(1);
  s[0]=getNumberOfAtoms() / natoms_per_task;
  addValue( s );
}

template <class T>
void MultiColvarTemplate<T>::addComponentWithDerivatives( const std::string& compName, const std::vector<std::size_t>& shape ) {
  std::vector<std::size_t> s(1);
  s[0]=getNumberOfAtoms() / natoms_per_task;
  addComponent( compName, s );
}

template <class T>
void MultiColvarTemplate<T>::getInputData( std::vector<double>& inputdata ) const {
  std::size_t ntasks = getConstPntrToComponent(0)->getNumberOfStoredValues();
  if( inputdata.size()!=5*natoms_per_task*ntasks ) {
    inputdata.resize( 5*natoms_per_task*ntasks );
  }

  std::size_t k=0;
  for(unsigned i=0; i<ntasks; ++i) {
    for(unsigned j=0; j<natoms_per_task; ++j) {
      Vector mypos( getPosition( natoms_per_task*i + j ) );
      inputdata[k] = mypos[0];
      k++;
      inputdata[k] = mypos[1];
      k++;
      inputdata[k] = mypos[2];
      k++;
    }
    for(unsigned j=0; j<natoms_per_task; ++j) {
      inputdata[k] = getMass( natoms_per_task*i + j );
      k++;
    }
    for(unsigned j=0; j<natoms_per_task; ++j) {
      inputdata[k] = getCharge( natoms_per_task*i + j );
      k++;
    }
  }
}

template <class T>
void MultiColvarTemplate<T>::performTask( std::size_t task_index,
    const MultiColvarInput& actiondata,
    ParallelActionsInput& input,
    ParallelActionsOutput& output ) {
  std::size_t pos_start = 5*actiondata.nindices_per_task*task_index;
  if( actiondata.usepbc ) {
    if( actiondata.nindices_per_task==1 ) {
      //this may be changed to input.pbc.apply() en mass or only on this one
      Vector fpos=input.pbc->distance(Vector(0.0,0.0,0.0),
                                      Vector(input.inputdata[pos_start],
                                             input.inputdata[pos_start+1],
                                             input.inputdata[pos_start+2]) );
      input.inputdata[pos_start]  =fpos[0];
      input.inputdata[pos_start+1]=fpos[1];
      input.inputdata[pos_start+2]=fpos[2];
    } else {
      //make whole?
      std::size_t apos_start = pos_start;
      //if accidentaly nindices_per_task is 0, this will work by looping on all possible unsigned integers!!!!
      for(unsigned j=0; j<actiondata.nindices_per_task-1; ++j) {
        Vector first(input.inputdata[apos_start],
                     input.inputdata[apos_start+1],
                     input.inputdata[apos_start+2]);
        Vector second(input.inputdata[apos_start+3],
                      input.inputdata[apos_start+4],
                      input.inputdata[apos_start+5]);
        //calling the pbc here gives problems
        second=first+input.pbc->distance(first,second);
        input.inputdata[apos_start+3]=second[0];
        input.inputdata[apos_start+4]=second[1];
        input.inputdata[apos_start+5]=second[2];
        apos_start += 3;
      }
    }
  } else if( actiondata.nindices_per_task==1 ) {
    //isn't this equivalent to x = x-0?
    //why this is needed?
    Vector fpos=delta(Vector(0.0,0.0,0.0),
                      Vector(input.inputdata[pos_start],
                             input.inputdata[pos_start+1],
                             input.inputdata[pos_start+2]));
    input.inputdata[pos_start]=fpos[0];
    input.inputdata[pos_start+1]=fpos[1];
    input.inputdata[pos_start+2]=fpos[2];
  }
  const size_t mass_start = pos_start + 3*actiondata.nindices_per_task;
  const size_t charge_start = mass_start + actiondata.nindices_per_task;
  const size_t local_ndev = 3*actiondata.nindices_per_task+virialSize;

  ColvarOutput cvout { output.values,
                       local_ndev,
                       output.derivatives.data() };
  T::calculateCV( ColvarInput{actiondata.mode,
                              actiondata.nindices_per_task,
                              input.inputdata+pos_start,
                              input.inputdata+mass_start,
                              input.inputdata+charge_start,
                              *input.pbc},
                  cvout );
}

template <class T>
int MultiColvarTemplate<T>::getNumberOfValuesPerTask( std::size_t task_index,
    const MultiColvarInput& actiondata ) {
  return 1;
}

template <class T>
void MultiColvarTemplate<T>::getForceIndices( std::size_t task_index,
    std::size_t colno,
    std::size_t ntotal_force,
    const MultiColvarInput& actiondata,
    const ParallelActionsInput& input,
    ForceIndexHolder force_indices ) {
  for(unsigned i=0; i<input.ncomponents; ++i) {
    std::size_t m=0;
    std::size_t base = 3*task_index*actiondata.nindices_per_task;
    for(unsigned j=0; j<actiondata.nindices_per_task; ++j) {
      force_indices.indices[i][m] = base + m;
      ++m;
      force_indices.indices[i][m] = base + m;
      ++m;
      force_indices.indices[i][m] = base + m;
      ++m;
    }
    force_indices.threadsafe_derivatives_end[i] = 3*actiondata.nindices_per_task;
    force_indices.indices[i][m+0] = ntotal_force - 9;
    force_indices.indices[i][m+1] = ntotal_force - 8;
    force_indices.indices[i][m+2] = ntotal_force - 7;
    force_indices.indices[i][m+3] = ntotal_force - 6;
    force_indices.indices[i][m+4] = ntotal_force - 5;
    force_indices.indices[i][m+5] = ntotal_force - 4;
    force_indices.indices[i][m+6] = ntotal_force - 3;
    force_indices.indices[i][m+7] = ntotal_force - 2;
    force_indices.indices[i][m+8] = ntotal_force - 1;
    force_indices.tot_indices[i] = 3*actiondata.nindices_per_task + virialSize;
  }
}

} // namespace colvar
} // namespace PLMD
#endif
