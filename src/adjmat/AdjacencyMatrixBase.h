/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2020 The plumed team
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
#ifndef __PLUMED_adjmat_AdjacencyMatrixBase_h
#define __PLUMED_adjmat_AdjacencyMatrixBase_h

#include <vector>
#include "core/ParallelTaskManager.h"
#include "core/ActionWithMatrix.h"
#include "core/PlumedMain.h"
#include "tools/LinkCells.h"

namespace PLMD {
namespace adjmat {

template <class T>
class AdjacencyMatrixData {
public:
  T matrixdata;
  bool usepbc{true};
  bool components{false};
  std::size_t nlists{0};
  unsigned natoms_per_list{0};
  std::vector<std::size_t> nlist;
  unsigned natoms_per_three_list{0};
  std::vector<std::size_t> nlist_three;
};

class AdjacencyMatrixInput {
public:
  bool noderiv{false};
  const Pbc* pbc;
  Vector pos;
  std::size_t natoms{0};
  std::vector<Vector> extra_positions;
  AdjacencyMatrixInput( bool n, const Pbc* b, Vector p ) : noderiv(n), pbc(b), pos(p) {}
};

class MatrixOutput {
public:
  View<double,1> val;
  View<double,helpers::dynamic_extent> deriv;
  MatrixOutput( std::size_t nd, double* vp, double* dp ) : val(vp), deriv(dp,nd) {}
};

template <class T>
class AdjacencyMatrixBase : public ActionWithMatrix {
public:
  using input_type = AdjacencyMatrixData<T>;
  using PTM = ParallelTaskManager<AdjacencyMatrixBase<T>>;
private:
  PTM taskmanager;
  bool nopbc, read_one_group;
  LinkCells linkcells, threecells;
  std::vector<unsigned> ablocks, threeblocks;
  double nl_cut, nl_cut2;
  unsigned nl_stride;
  void setupThirdAtomBlock( const std::vector<AtomNumber>& tc, std::vector<AtomNumber>& t );
public:
  static constexpr size_t virialSize = 9;
  static void registerKeywords( Keywords& keys );
  explicit AdjacencyMatrixBase(const ActionOptions&);
  unsigned getNumberOfDerivatives() override ;
  void calculate() override ;
  void applyNonZeroRankForces( std::vector<double>& outforces ) override ;
  void getInputData( std::vector<double>& inputdata ) const ;
  std::string writeInGraph() const override {
    if( getName()=="CONTACT_MATRIX_PROPER" ) {
      return "CONTACT_MATRIX";
    }
    return getName();
  }
  void setLinkCellCutoff( const bool& symmetric, const double& lcut, double tcut=-1.0 );
  static void performTask( std::size_t task_index,
                           const AdjacencyMatrixData<T>& actiondata,
                           ParallelActionsInput& input,
                           ParallelActionsOutput& output );
  static void gatherForces( std::size_t task_index,
                            const AdjacencyMatrixData<T>& actiondata,
                            const ParallelActionsInput& input,
                            const ForceInput& fdata,
                            ForceOutput forces );
};

template <class T>
void AdjacencyMatrixBase<T>::registerKeywords( Keywords& keys ) {
  ActionWithMatrix::registerKeywords( keys );
  keys.addInputKeyword("optional","MASK","vector","a vector that is used to used to determine which rows of the adjancency matrix to compute");
  keys.add("atoms","GROUP","the atoms for which you would like to calculate the adjacency matrix");
  keys.add("atoms","GROUPA","");
  keys.add("atoms","GROUPB","");
  keys.add("atoms-2","ATOMS","the atoms for which you would like to calculate the adjacency matrix. This is a depracated syntax that is equivalent to GROUP.  You are strongly recommened to use GROUP instead of ATOMS.");
  keys.reserve("atoms","GROUPC","");
  keys.addFlag("COMPONENTS",false,"also calculate the components of the vector connecting the atoms in the contact matrix");
  keys.addFlag("NOPBC",false,"don't use pbc");
  keys.add("compulsory","NL_CUTOFF","0.0","The cutoff for the neighbor list.  A value of 0 means we are not using a neighbor list");
  keys.add("compulsory","NL_STRIDE","1","The frequency with which we are updating the atoms in the neighbor list");
  T::registerKeywords( keys );
  PTM::registerKeywords( keys );
  keys.addOutputComponent("w","COMPONENTS","matrix","a matrix containing the weights for the bonds between each pair of atoms");
  keys.addOutputComponent("x","COMPONENTS","matrix","the projection of the bond on the x axis");
  keys.addOutputComponent("y","COMPONENTS","matrix","the projection of the bond on the y axis");
  keys.addOutputComponent("z","COMPONENTS","matrix","the projection of the bond on the z axis");
  keys.setValueDescription("matrix","a matrix containing the weights for the bonds between each pair of atoms");
}

template <class T>
AdjacencyMatrixBase<T>::AdjacencyMatrixBase(const ActionOptions& ao):
  Action(ao),
  ActionWithMatrix(ao),
  taskmanager(this),
  read_one_group(false),
  linkcells(comm),
  threecells(comm) {
  std::vector<std::size_t> shape(2);
  std::vector<AtomNumber> t;
  parseAtomList("GROUP", t );
  if( t.size()==0 ) {
    parseAtomList("ATOMS", t);
    if( t.size()>0 ) {
      warning("using depracated syntax for contact matrix.  You are strongly recommended to use GROUP instead of ATOMS");
    }
  }

  if( t.size()==0 ) {
    std::vector<AtomNumber> ta;
    parseAtomList("GROUPA",ta);
    if( ta.size()==0 && getName()=="HBOND_MATRIX") {
      parseAtomList("DONORS",ta);
    }
    std::vector<AtomNumber> tb;
    parseAtomList("GROUPB",tb);
    if( tb.size()==0 && getName()=="HBOND_MATRIX") {
      parseAtomList("ACCEPTORS",tb);
    }
    if( ta.size()==0 || tb.size()==0 ) {
      error("no atoms have been specified in input");
    }

    // Create list of tasks
    log.printf("  atoms in GROUPA ");
    for(unsigned i=0; i<ta.size(); ++i) {
      log.printf("%d ", ta[i].serial());
      t.push_back(ta[i]);
    }
    log.printf("\n");
    log.printf("  atoms in GROUPB ");
    ablocks.resize( tb.size() );
    unsigned n=0;
    for(unsigned i=0; i<tb.size(); ++i) {
      log.printf("%d ", tb[i].serial());
      ablocks[i]=ta.size()+n;
      t.push_back(tb[i]);
      n++;
    }
    log.printf("\n");
    shape[0]=ta.size();
    shape[1]=tb.size();
  } else {
    // Create list of tasks
    log.printf("  atoms in GROUP ");
    ablocks.resize( t.size() );
    for(unsigned i=0; i<t.size(); ++i) {
      log.printf("%d ", t[i].serial());
      ablocks[i]=i;
    }
    log.printf("\n");
    shape[0]=shape[1]=t.size();
    read_one_group=true;
  }
  if( keywords.exists("GROUPC") ) {
    std::vector<AtomNumber> tc;
    parseAtomList("GROUPC",tc);
    if( tc.size()==0 ) {
      error("no atoms in GROUPC specified");
    }
    log.printf("  atoms in GROUPC ");
    setupThirdAtomBlock( tc, t );
  } else if( keywords.exists("BRIDGING_ATOMS") ) {
    std::vector<AtomNumber> tc;
    parseAtomList("BRIDGING_ATOMS",tc);
    if( tc.size()==0 ) {
      error("no BRIDGING_ATOMS specified");
    }
    log.printf("  bridging atoms are ");
    setupThirdAtomBlock( tc, t );
  } else if( keywords.exists("HYDROGENS") ) {
    std::vector<AtomNumber> tc;
    parseAtomList("HYDROGENS",tc);
    if( tc.size()==0 ) {
      error("no HYDROGEN atoms specified");
    }
    log.printf("  hydrogen atoms are ");
    setupThirdAtomBlock( tc, t );
  } else if( keywords.exists("BACKGROUND_ATOMS") ) {
    std::vector<AtomNumber> tc;
    parseAtomList("BACKGROUND_ATOMS",tc);
    if( tc.size()==0 ) {
      error("no ATOMS atoms specified");
    }
    log.printf("  atoms for background density are ");
    setupThirdAtomBlock( tc, t );
  }
  // Request the atoms from the ActionAtomistic
  requestAtoms( t );
  bool components;
  parseFlag("COMPONENTS",components);
  parseFlag("NOPBC",nopbc);
  if( !components ) {
    addValue( shape );
    setNotPeriodic();
  } else {
    addComponent( "w", shape );
    componentIsNotPeriodic("w");
  }
  getPntrToComponent(0)->setDerivativeIsZeroWhenValueIsZero();
  // Stuff for neighbor list
  parse("NL_CUTOFF",nl_cut);
  nl_cut2=nl_cut*nl_cut;
  parse("NL_STRIDE",nl_stride);
  if( nl_cut==0 && nl_stride>1 ) {
    error("NL_CUTOFF must be set if NL_STRIDE is set greater than 1");
  }
  if( nl_cut>0 ) {
    log.printf("  using neighbor list with cutoff %f.  List is updated every %u steps.\n",nl_cut,nl_stride);
  }

  if( components ) {
    addComponent( "x", shape );
    componentIsNotPeriodic("x");
    addComponent( "y", shape );
    componentIsNotPeriodic("y");
    addComponent( "z", shape );
    componentIsNotPeriodic("z");
  }
  log<<"  Bibliography "<<plumed.cite("Tribello, Giberti, Sosso, Salvalaglio and Parrinello, J. Chem. Theory Comput. 13, 1317 (2017)")<<"\n";
  AdjacencyMatrixData<T> matdata;
  matdata.usepbc = !nopbc;
  matdata.components = components;
  matdata.nlists = getPntrToComponent(0)->getShape()[0];
  matdata.matrixdata.parseInput( this );
  taskmanager.setupParallelTaskManager( 1, 3 + virialSize, 3 + virialSize );
  taskmanager.setActionInput( matdata );
}

template <class T>
unsigned AdjacencyMatrixBase<T>::getNumberOfDerivatives() {
  return 3*getNumberOfAtoms() + 9;
}

template <class T>
void AdjacencyMatrixBase<T>::setupThirdAtomBlock( const std::vector<AtomNumber>& tc, std::vector<AtomNumber>& t ) {
  threeblocks.resize( tc.size() );
  unsigned base=t.size();
  for(unsigned i=0; i<tc.size(); ++i) {
    log.printf("%d ", tc[i].serial());
    t.push_back(tc[i]);
    threeblocks[i]=base+i;
  }
  log.printf("\n");
}

template <class T>
void AdjacencyMatrixBase<T>::setLinkCellCutoff( const bool& symmetric, const double& lcut, double tcut ) {
  if( read_one_group && symmetric ) {
    getPntrToComponent(0)->setSymmetric( true );
  }
  if( nl_cut>0 && lcut>nl_cut ) {
    error("D_MAX for switching functions should be shorter than neighbor list cutoff");
  }

  if( tcut<0 ) {
    tcut=lcut;
  }
  if( nl_cut>0 ) {
    linkcells.setCutoff( nl_cut );
  } else {
    linkcells.setCutoff( lcut );
  }
  if( linkcells.getCutoff()<std::numeric_limits<double>::max() ) {
    log.printf("  set link cell cutoff to %f \n", linkcells.getCutoff() );
  }
  threecells.setCutoff( tcut );
}

template <class T>
void AdjacencyMatrixBase<T>::getInputData( std::vector<double>& inputdata ) const {
  if( inputdata.size()!=3*getNumberOfAtoms() ) {
    inputdata.resize( 3*getNumberOfAtoms() );
  }

  std::size_t k=0;
  for(unsigned i=0; i<getNumberOfAtoms(); ++i) {
    Vector mypos( ActionAtomistic::getPosition(i) );
    inputdata[k] = mypos[0];
    k++;
    inputdata[k] = mypos[1];
    k++;
    inputdata[k] = mypos[2];
    k++;
  }
}

template <class T>
void AdjacencyMatrixBase<T>::calculate() {
  Value* myval = getPntrToComponent(0);
  // Retrieve the task list
  std::vector<unsigned> & pTaskList( getListOfActiveTasks(this) );
  // Get the number of tasks we have to deal with
  unsigned ntasks=myval->getShape()[0];
  if( nl_stride==1 ) {
    ntasks=pTaskList.size();
  } else {
    error("neighbour list non updates are not actually implemented or tested");
  }
  unsigned fbsize=0, lstart = getConstPntrToComponent(0)->getShape()[0];
  // This is triggered if you have GROUPA/GROUPB
  // in that case the second index for the matrix is recovered from the neighbour list
  // by subtracting the number of atoms in GROUPA as we do here.
  if( ablocks[0]>=lstart ) {
    fbsize = lstart;
  }

  // Get the atoms
  std::vector<Vector> ltmp_pos2( ntasks );
  for(unsigned i=0; i<ntasks; ++i) {
    ltmp_pos2[i] = ActionAtomistic::getPosition(pTaskList[i]);
  }
  AdjacencyMatrixData<T> & matdata = taskmanager.getActionInput();
  // Now update the neighbor lists
  if( getStep()%nl_stride==0 ) {
    // Build the link cells
    std::vector<Vector> ltmp_pos( ablocks.size() );
    for(unsigned i=0; i<ablocks.size(); ++i) {
      ltmp_pos[i]=ActionAtomistic::getPosition( ablocks[i] );
    }
    linkcells.createNeighborList( getConstPntrToComponent(0)->getShape()[0], ltmp_pos2, pTaskList, pTaskList, ltmp_pos, ablocks, getPbc(), matdata.natoms_per_list, matdata.nlist );
    if( threeblocks.size()>0 ) {
      std::vector<Vector> ltmp_pos3( threeblocks.size() );
      for(unsigned i=0; i<threeblocks.size(); ++i) {
        ltmp_pos3[i]=ActionAtomistic::getPosition( threeblocks[i] );
      }
      threecells.createNeighborList( getConstPntrToComponent(0)->getShape()[0], ltmp_pos2, pTaskList, pTaskList, ltmp_pos3, threeblocks, getPbc(), matdata.natoms_per_three_list, matdata.nlist_three );
    }
  } else {
    error("neighbour list non updates are not actually implemented or tested");
  }
  // And finally work out maximum number of columns to use
  unsigned maxcol = matdata.nlist[0];
  for(unsigned i=1; i<getConstPntrToComponent(0)->getShape()[0]; ++i) {
    if( matdata.nlist[i]>maxcol ) {
      maxcol = matdata.nlist[i];
    }
  }
  // The first element returned by the neighbour list is the central atom.  The
  // neighbour list thus always has one more atom in it than there are atoms in the
  // neighbourhood of the central atom
  maxcol = maxcol-1;

  // Reshape the matrix store if the number of columns has changed
  if( maxcol!=myval->getNumberOfColumns() ) {
    for(int i=0; i<getNumberOfComponents(); ++i) {
      getPntrToComponent(i)->reshapeMatrixStore( maxcol );
    }
  }
  // Clear the bookeeping array
  for(unsigned i=0; i<lstart; ++i) {
    myval->setMatrixBookeepingElement( i*(1+maxcol), 0 );
  }
  // Transfer neighbor list data to the bookeeping arrays
  for(unsigned i=0; i<ntasks; ++i) {
    unsigned bstart = pTaskList[i]*(1+maxcol);
    unsigned rstart = lstart + pTaskList[i]*(1+matdata.natoms_per_list);
    myval->setMatrixBookeepingElement( bstart, matdata.nlist[pTaskList[i]] - 1 );
    for(unsigned j=1; j<matdata.nlist[pTaskList[i]]; ++j) {
      myval->setMatrixBookeepingElement( bstart+j, matdata.nlist[rstart+j] - fbsize );
    }
  }
  for(unsigned i=1; i<getNumberOfComponents(); ++i) {
    getPntrToComponent(i)->copyBookeepingArrayFromArgument( getPntrToComponent(0) );
  }
  // We need to setup the task manager here because at this point we know the size of the sparse matrices
  unsigned nder = maxcol*(6 + 3*matdata.natoms_per_three_list + virialSize);
  if( matdata.components ) {
    taskmanager.setupParallelTaskManager( 0, nder, -1, 4*maxcol );
  } else {
    taskmanager.setupParallelTaskManager( 0, nder, -1, maxcol );
  }
  // And run all the tasks
  taskmanager.runAllTasks();
}

template <class T>
void AdjacencyMatrixBase<T>::performTask( std::size_t task_index,
    const AdjacencyMatrixData<T>& actiondata,
    ParallelActionsInput& input,
    ParallelActionsOutput& output ) {
  unsigned n3neigh=0, nneigh=actiondata.nlist[task_index];
  if( actiondata.natoms_per_three_list>0 ) {
    n3neigh = actiondata.nlist_three[task_index];
  }
  std::vector<Vector> atoms( nneigh + n3neigh );
  unsigned fstart = actiondata.nlists + task_index*(1+actiondata.natoms_per_list);
  Vector pos0( input.inputdata[3*task_index+0],input.inputdata[3*task_index+1],input.inputdata[3*task_index+2] );
  for(unsigned i=0; i<nneigh; ++i) {
    atoms[i][0] = input.inputdata[3*actiondata.nlist[fstart+i]+0] - pos0[0];
    atoms[i][1] = input.inputdata[3*actiondata.nlist[fstart+i]+1] - pos0[1];
    atoms[i][2] = input.inputdata[3*actiondata.nlist[fstart+i]+2] - pos0[2];
  }
  // Retrieve the set of third atoms
  unsigned fstart3 = actiondata.nlists + task_index*(1+actiondata.natoms_per_three_list);
  for(unsigned i=1; i<n3neigh; ++i) {
    atoms[nneigh+i-1][0] = input.inputdata[3*actiondata.nlist_three[fstart3+i]+0] - pos0[0];
    atoms[nneigh+i-1][1] = input.inputdata[3*actiondata.nlist_three[fstart3+i]+1] - pos0[1];
    atoms[nneigh+i-1][2] = input.inputdata[3*actiondata.nlist_three[fstart3+i]+2] - pos0[2];
  }
  // Apply periodic boundary conditions to all the atoms
  if( actiondata.usepbc ) {
    input.pbc->apply( atoms, atoms.size() );
  }
  AdjacencyMatrixInput adjinp( input.noderiv, input.pbc, atoms[0] );
  if( n3neigh>1 ) {
    adjinp.natoms = n3neigh-1;
    adjinp.extra_positions.resize( adjinp.natoms );
    for(unsigned i=0; i<adjinp.natoms; ++i) {
      adjinp.extra_positions[i] = atoms[nneigh+i];
    }
  }

  // And calculate this row of the matrices
  unsigned ncomponents = 1;
  std::size_t nderiv(6 + 3*adjinp.natoms + virialSize);
  if( actiondata.components ) {
    ncomponents = 4;
  }

  // Must clear the derivatives here as otherwise sparsity pattern
  // of derivatives does not match sparsity pattern for forces
  if( !input.noderiv ) {
    for(unsigned i=0; i<output.derivatives.size(); ++i) {
      output.derivatives[i] = 0.0;
    }
  }

  for(unsigned i=1; i<nneigh; ++i ) {
    adjinp.pos = atoms[i];
    std::size_t valpos = (i-1)*ncomponents;
    MatrixOutput adjout( nderiv, output.values.data() + valpos, output.derivatives.data() + valpos*nderiv );
    T::calculateWeight( actiondata.matrixdata, adjinp, adjout );
    if( !actiondata.components ) {
      continue ;
    }
    output.values[valpos+1] = atoms[i][0];
    output.values[valpos+2] = atoms[i][1];
    output.values[valpos+3] = atoms[i][2];
    if( input.noderiv ) {
      continue ;
    }
    // We store the distance here as we use it to calculate the virial contribution for the force on the x, y and z components of the bond vector
    output.derivatives[(valpos+1)*nderiv+0] = atoms[i][0];
    output.derivatives[(valpos+1)*nderiv+1] = atoms[i][1];
    output.derivatives[(valpos+1)*nderiv+2] = atoms[i][2];
  }
}

template <class T>
void AdjacencyMatrixBase<T>::applyNonZeroRankForces( std::vector<double>& outforces ) {
  taskmanager.applyForces( outforces );
}

template <class T>
void AdjacencyMatrixBase<T>::gatherForces( std::size_t task_index,
    const AdjacencyMatrixData<T>& actiondata,
    const ParallelActionsInput& input,
    const ForceInput& fdata,
    ForceOutput forces ) {
  unsigned n3neigh = 0;
  if( actiondata.natoms_per_three_list>0 ) {
    n3neigh = actiondata.nlist_three[task_index];
  }
  unsigned ncomponents = 1, nder = 6 + 3*actiondata.natoms_per_three_list + virialSize;
  if( actiondata.components ) {
    ncomponents = 4;
  }
  unsigned nneigh = actiondata.nlist[task_index];
  unsigned fstart = actiondata.nlists + task_index*(1+actiondata.natoms_per_list);
  unsigned fstart3 = actiondata.nlists + task_index*(1+actiondata.natoms_per_three_list);
  for(unsigned i=1; i<nneigh; ++i ) {
    std::size_t n = 0;
    unsigned fpos = (i-1)*ncomponents;
    double force = fdata.force[fpos];
    View<const double,helpers::dynamic_extent> deriv( fdata.deriv.data() + fpos*nder, nder );
    forces.thread_safe[ 3*task_index + 0 ] += force*deriv[n];
    n++;
    forces.thread_safe[ 3*task_index + 1 ] += force*deriv[n];
    n++;
    forces.thread_safe[ 3*task_index + 2 ] += force*deriv[n];
    n++;
    unsigned myatom = actiondata.nlist[fstart+i];
    forces.thread_safe[ 3*myatom + 0 ] += force*deriv[n];
    n++;
    forces.thread_safe[ 3*myatom + 1 ] += force*deriv[n];
    n++;
    forces.thread_safe[ 3*myatom + 2 ] += force*deriv[n];
    n++;
    for(unsigned j=1; j<n3neigh; ++j) {
      unsigned my3atom = actiondata.nlist_three[fstart3+j];
      forces.thread_safe[ 3*my3atom + 0 ] += force*deriv[n];
      n++;
      forces.thread_safe[ 3*my3atom + 1 ] += force*deriv[n];
      n++;
      forces.thread_safe[ 3*my3atom + 2 ] += force*deriv[n];
      n++;
    }
    unsigned virstart = forces.thread_safe.size()-virialSize;
    for(unsigned j=virstart; j<forces.thread_safe.size(); ++j) {
      forces.thread_safe[j] += force*deriv[n];
      n++;
    }
    if( !actiondata.components ) {
      continue ;
    }
    View<const double,3> dist( fdata.deriv.data() + (1+fpos)*nder );
    if( fabs(dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2])<epsilon ) {
      continue ;
    }

    forces.thread_safe[ 3*task_index + 0 ] += -fdata.force[fpos+1];
    forces.thread_safe[ 3*task_index + 1 ] += -fdata.force[fpos+2];
    forces.thread_safe[ 3*task_index + 2 ] += -fdata.force[fpos+3];
    forces.thread_safe[ 3*myatom + 0 ] += fdata.force[fpos+1];
    forces.thread_safe[ 3*myatom + 1 ] += fdata.force[fpos+2];
    forces.thread_safe[ 3*myatom + 2 ] += fdata.force[fpos+3];

    forces.thread_safe[virstart + 0] += -fdata.force[fpos+1]*dist[0];
    forces.thread_safe[virstart + 1] += -fdata.force[fpos+2]*dist[0];
    forces.thread_safe[virstart + 2] += -fdata.force[fpos+3]*dist[0];
    forces.thread_safe[virstart + 3] += -fdata.force[fpos+1]*dist[1];
    forces.thread_safe[virstart + 4] += -fdata.force[fpos+2]*dist[1];
    forces.thread_safe[virstart + 5] += -fdata.force[fpos+3]*dist[1];
    forces.thread_safe[virstart + 6] += -fdata.force[fpos+1]*dist[2];
    forces.thread_safe[virstart + 7] += -fdata.force[fpos+2]*dist[2];
    forces.thread_safe[virstart + 8] += -fdata.force[fpos+3]*dist[2];
  }
}

}
}

#endif
