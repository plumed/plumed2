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
struct AdjacencyMatrixData {
  T matrixdata;
  bool usepbc{true};
  bool components{false};
  std::size_t nlists{0};
  unsigned natoms_per_list{0};
  std::vector<std::size_t> nlist_v;
  std::size_t *nlist{nullptr};
  unsigned natoms_per_three_list{0};
  std::vector<std::size_t> nlist_three_v;
  std::size_t* nlist_three{nullptr};
  void update() {
    nlist=nlist_v.data();
    nlist_three=nlist_three_v.data();
  }
#ifdef __PLUMED_USE_OPENACC
  void toACCDevice() const {
#pragma acc enter data copyin(this[0:1],usepbc,components,nlists, \
                              natoms_per_list,nlist[0:nlist_v.size()], \
                              natoms_per_three_list, \
                              nlist_three[0:nlist_three_v.size()])
    matrixdata.toACCDevice();
  }
  void removeFromACCDevice() const {
    matrixdata.removeFromACCDevice();
#pragma acc exit data delete(nlist_three[0:nlist_three_v.size()], \
                             natoms_per_three_list, \
                             nlist[0:nlist_v.size()], natoms_per_list, \
                             nlists, components, usepbc, this[0:1])
  }
#endif //__PLUMED_USE_OPENACC
};

struct AdjacencyMatrixInput {
  bool noderiv{false};
  const Pbc* pbc;
  Vector pos;
  std::size_t natoms{0};
  VectorView extra_positions;
};

struct MatrixOutput {
  View<double,1> val;
  View<double> deriv;

  ///doing t= Tensor(v1,v2); deriv[x:x+9]=t with no extra memory allocation
  template <typename Iterable1, typename Iterable2>
  void assignOuterProduct(const std::size_t startingIndex,
                          const Iterable1& v1,
                          const Iterable2& v2 ) {
    deriv[startingIndex + 0] = v1[0] * v2[0];
    deriv[startingIndex + 1] = v1[0] * v2[1];
    deriv[startingIndex + 2] = v1[0] * v2[2];
    deriv[startingIndex + 3] = v1[1] * v2[0];
    deriv[startingIndex + 4] = v1[1] * v2[1];
    deriv[startingIndex + 5] = v1[1] * v2[2];
    deriv[startingIndex + 6] = v1[2] * v2[0];
    deriv[startingIndex + 7] = v1[2] * v2[1];
    deriv[startingIndex + 8] = v1[2] * v2[2];
  }
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
  void setupThirdAtomBlock( const std::vector<AtomNumber>& tc,
                            std::vector<AtomNumber>& t );
public:
  static constexpr size_t virialSize = 9;
  static void registerKeywords( Keywords& keys );
  explicit AdjacencyMatrixBase(const ActionOptions&);
  unsigned getNumberOfDerivatives() override ;
  void calculate() override ;
  void applyNonZeroRankForces( std::vector<double>& outforces ) override ;
  void getInputData( std::vector<double>& inputdata ) const override;
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
                           ParallelActionsOutput output );
  static int getNumberOfValuesPerTask( std::size_t task_index,
                                       const AdjacencyMatrixData<T>& actiondata );
  static void getForceIndices( std::size_t task_index,
                               std::size_t colno,
                               std::size_t ntotal_force,
                               const AdjacencyMatrixData<T>& actiondata,
                               const ParallelActionsInput& input,
                               ForceIndexHolder force_indices );
  void getMatrixColumnTitles( std::vector<std::string>& argnames ) const override ;
};

template <class T>
void AdjacencyMatrixBase<T>::registerKeywords( Keywords& keys ) {
  ActionWithMatrix::registerKeywords( keys );
  keys.addInputKeyword("optional","MASK","vector","a vector that is used to used to determine which rows of the adjancency matrix to compute");
  keys.add("atoms","GROUP","the atoms for which you would like to calculate the adjacency matrix");
  keys.add("atoms","GROUPA","when you are calculating the adjacency matrix between two sets of atoms this keyword is used to specify the atoms along with the keyword GROUPB");
  keys.add("atoms","GROUPB","when you are calculating the adjacency matrix between two sets of atoms this keyword is used to specify the atoms along with the keyword GROUPA");
  keys.addDeprecatedKeyword("ATOMS","GROUP");
  keys.reserve("atoms","GROUPC","a group of atoms that must be summed over when calculating each element of the adjacency matrix");
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
  keys.addDOI("10.1021/acs.jctc.6b01073");
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
  if( getName()!="HBOND_MATRIX" ) {
    if( t.size()==0 ) {
      parseAtomList("ATOMS", t);
      if( t.size()>0 ) {
        warning("using depracated syntax for contact matrix.  You are strongly recommended to use GROUP instead of ATOMS");
      }
    }
  } else if( t.size()>0 ) {
    warning("GROUP keyword has been deprecated for HBOND_MATRIX as it may lead users to wrongly assume that the matrices calculated by this action are symmetric.  We strongly recommend using DONORS/ACCEPTORS instead");
  }

  if( t.size()==0 ) {
    std::vector<AtomNumber> ta;
    if( getName()=="HBOND_MATRIX") {
      parseAtomList("DONORS",ta);
    } else {
      parseAtomList("GROUPA",ta);
    }
    std::vector<AtomNumber> tb;
    if( getName()=="HBOND_MATRIX") {
      parseAtomList("ACCEPTORS",tb);
    } else {
      parseAtomList("GROUPB",tb);
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
  requestAtoms( t, false );
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
  taskmanager.setActionInput( matdata );
}

template <class T>
unsigned AdjacencyMatrixBase<T>::getNumberOfDerivatives() {
  return 3*getNumberOfAtoms() + 9;
}

template <class T>
void AdjacencyMatrixBase<T>::setupThirdAtomBlock( const std::vector<AtomNumber>& tc,
    std::vector<AtomNumber>& t ) {
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
void AdjacencyMatrixBase<T>::setLinkCellCutoff( const bool& symmetric,
    const double& lcut,
    double tcut ) {
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
void AdjacencyMatrixBase<T>::getMatrixColumnTitles( std::vector<std::string>& argnames ) const {
  std::string num;
  for(unsigned i=0; i<getConstPntrToComponent(0)->getShape()[1]; ++i) {
    Tools::convert( i+1, num );
    argnames.push_back( num );
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
  unsigned fbsize=0;
  unsigned lstart = getConstPntrToComponent(0)->getShape()[0];
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
  auto & matdata = taskmanager.getActionInput();
  // Now update the neighbor lists
  if( getStep()%nl_stride==0 ) {
    // Build the link cells
    std::vector<Vector> ltmp_pos( ablocks.size() );
    for(unsigned i=0; i<ablocks.size(); ++i) {
      ltmp_pos[i]=ActionAtomistic::getPosition( ablocks[i] );
    }
    linkcells.createNeighborList( getConstPntrToComponent(0)->getShape()[0],
                                  make_const_view(ltmp_pos2),
                                  make_const_view(pTaskList),
                                  make_const_view(pTaskList),
                                  make_const_view(ltmp_pos),
                                  make_const_view(ablocks),
                                  getPbc(),
                                  matdata.natoms_per_list,
                                  matdata.nlist_v );
    if( threeblocks.size()>0 ) {
      std::vector<Vector> ltmp_pos3( threeblocks.size() );
      for(unsigned i=0; i<threeblocks.size(); ++i) {
        ltmp_pos3[i]=ActionAtomistic::getPosition( threeblocks[i] );
      }
      threecells.createNeighborList( getConstPntrToComponent(0)->getShape()[0],
                                     make_const_view(ltmp_pos2),
                                     make_const_view(pTaskList),
                                     make_const_view(pTaskList),
                                     make_const_view(ltmp_pos3),
                                     make_const_view(threeblocks),
                                     getPbc(),
                                     matdata.natoms_per_three_list,
                                     matdata.nlist_three_v);
    }
    matdata.update();
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
    for(unsigned i=0; i<getNumberOfComponents(); ++i) {
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
  unsigned nder = 6 + 3*matdata.natoms_per_three_list + virialSize;
  //there was an if on `matdata.components` but both branches were calling the next line
  taskmanager.setupParallelTaskManager( nder, getNumberOfDerivatives() );
  // Create the workspace that we use in performTask
  taskmanager.setWorkspaceSize(3*(maxcol + 2 + matdata.natoms_per_three_list) );
  // And run all the tasks
  taskmanager.runAllTasks();
}

template <class T>
void AdjacencyMatrixBase<T>::performTask( std::size_t task_index,
    const AdjacencyMatrixData<T>& actiondata,
    ParallelActionsInput& input,
    ParallelActionsOutput output ) {
  const unsigned nneigh=actiondata.nlist[task_index];
  const unsigned n3neigh=(actiondata.natoms_per_three_list>0) ?
                         actiondata.nlist_three[task_index]:0;
  VectorView atoms(output.buffer.data(), nneigh + n3neigh );

  const unsigned fstart = actiondata.nlists + task_index*(1+actiondata.natoms_per_list);
  const View<double,3> pos0( input.inputdata + 3 * task_index);
  for(unsigned i=0; i<nneigh; ++i) {
    atoms[i][0] = input.inputdata[3*actiondata.nlist[fstart+i]+0] - pos0[0];
    atoms[i][1] = input.inputdata[3*actiondata.nlist[fstart+i]+1] - pos0[1];
    atoms[i][2] = input.inputdata[3*actiondata.nlist[fstart+i]+2] - pos0[2];
  }

  if( actiondata.natoms_per_three_list>0 ) {
    // Retrieve the set of third atoms
    unsigned fstart3 = actiondata.nlists
                       + task_index*(1+actiondata.natoms_per_three_list);
    for(unsigned i=1; i<n3neigh; ++i) {
      atoms[nneigh+i-1][0] =
        input.inputdata[3*actiondata.nlist_three[fstart3+i]+0] - pos0[0];
      atoms[nneigh+i-1][1] =
        input.inputdata[3*actiondata.nlist_three[fstart3+i]+1] - pos0[1];
      atoms[nneigh+i-1][2] =
        input.inputdata[3*actiondata.nlist_three[fstart3+i]+2] - pos0[2];
    }
  }


  // Apply periodic boundary conditions to all the atoms
  if( actiondata.usepbc ) {
    input.pbc->apply( atoms, atoms.size() );
  }
  AdjacencyMatrixInput adjinp {input.noderiv,
                               input.pbc,
                               Vector{0.0,0.0,0.0},
                               0,
                               VectorView{ atoms.data() + 3*nneigh,n3neigh}};

  if( n3neigh>1 ) {
    adjinp.natoms = n3neigh-1;
  }

  // And calculate this row of the matrices
  std::size_t nderiv(6 + 3*adjinp.natoms + virialSize);
  const unsigned ncomponents =(actiondata.components) ? 4 : 1;

  // Must clear the derivatives here as otherwise sparsity pattern
  // of derivatives does not match sparsity pattern for forces
  if( !input.noderiv ) {
    for(unsigned i=0; i<output.derivatives.size(); ++i) {
      output.derivatives[i] = 0.0;
    }
  }

  for(unsigned i=1; i<nneigh; ++i ) {
    adjinp.pos = Vector(atoms[i][0],atoms[i][1],atoms[i][2]);
    const std::size_t valpos = (i-1)*ncomponents;
    MatrixOutput adjout{View<double,1>{&output.values[ valpos]},
                        View{output.derivatives.data() + valpos*nderiv,nderiv}};
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
    //sugar for not having to repeat [valpos*nderiv+something]
    for(int ii=1; ii<4; ++ii) {
      auto derivs  = output.derivatives.subview_n<5>(valpos*nderiv+ii*nderiv);
      derivs[0] = -1.0;
      derivs[1] =  1.0;
      derivs[2] = -atoms[i][0];
      derivs[3] = -atoms[i][1];
      derivs[4] = -atoms[i][2];
    }

  }
}

template <class T>
void AdjacencyMatrixBase<T>::applyNonZeroRankForces( std::vector<double>& outforces ) {
  taskmanager.applyForces( outforces );
}

template <class T>
int AdjacencyMatrixBase<T>::getNumberOfValuesPerTask( std::size_t task_index,
    const AdjacencyMatrixData<T>& actiondata ) {
  return actiondata.nlist[task_index] - 1;
}

template <class T>
void AdjacencyMatrixBase<T>::getForceIndices( const std::size_t task_index,
    const std::size_t colno,
    const std::size_t ntotal_force,
    const AdjacencyMatrixData<T>& actiondata,
    const ParallelActionsInput& /*input*/,
    ForceIndexHolder force_indices ) {

  const unsigned fstart = actiondata.nlists
                          + task_index*(1 + actiondata.natoms_per_list);

  const auto three_task_index= 3 * task_index;
  force_indices.indices[0][0] = three_task_index;
  force_indices.indices[0][1] = three_task_index + 1;
  force_indices.indices[0][2] = three_task_index + 2;

  const unsigned myatom = 3 * actiondata.nlist[fstart + 1 + colno];
  force_indices.indices[0][3] = myatom;
  force_indices.indices[0][4] = myatom + 1;
  force_indices.indices[0][5] = myatom + 2;
  if( actiondata.components ) {
    force_indices.indices[1][0] = three_task_index;
    force_indices.indices[1][1] = myatom;
    force_indices.indices[2][0] = three_task_index + 1;
    force_indices.indices[2][1] = myatom + 1;
    force_indices.indices[3][0] = three_task_index + 2;
    force_indices.indices[3][1] = myatom + 2;
  }

  force_indices.threadsafe_derivatives_end[0] = 0;
  unsigned n = 6;
  if( actiondata.natoms_per_three_list>0 ) {
    const unsigned n3neigh = actiondata.nlist_three[task_index];
    const unsigned fstart3 = actiondata.nlists
                             + task_index*(1 + actiondata.natoms_per_three_list);
    for(unsigned j=1; j<n3neigh; ++j) {
      unsigned my3atom = 3 * actiondata.nlist_three[fstart3 + j];
      force_indices.indices[0][n  ] = my3atom;
      force_indices.indices[0][n+1] = my3atom + 1;
      force_indices.indices[0][n+2] = my3atom + 2;
      n += 3;
    }
  }
  const unsigned virstart = ntotal_force - 9;
  force_indices.indices[0][n  ] = virstart + 0;
  force_indices.indices[0][n+1] = virstart + 1;
  force_indices.indices[0][n+2] = virstart + 2;
  force_indices.indices[0][n+3] = virstart + 3;
  force_indices.indices[0][n+4] = virstart + 4;
  force_indices.indices[0][n+5] = virstart + 5;
  force_indices.indices[0][n+6] = virstart + 6;
  force_indices.indices[0][n+7] = virstart + 7;
  force_indices.indices[0][n+8] = virstart + 8;

  force_indices.tot_indices[0] = n + 9;
  if( actiondata.components ) {
    force_indices.threadsafe_derivatives_end[1] = 0;
    force_indices.threadsafe_derivatives_end[2] = 0;
    force_indices.threadsafe_derivatives_end[3] = 0;
    force_indices.indices[1][2] = virstart;
    force_indices.indices[1][3] = virstart + 1;
    force_indices.indices[1][4] = virstart + 2;
    force_indices.indices[2][2] = virstart + 3;
    force_indices.indices[2][3] = virstart + 4;
    force_indices.indices[2][4] = virstart + 5;
    force_indices.indices[3][2] = virstart + 6;
    force_indices.indices[3][3] = virstart + 7;
    force_indices.indices[3][4] = virstart + 8;
    force_indices.tot_indices[1] = 5;
    force_indices.tot_indices[2] = 5;
    force_indices.tot_indices[3] = 5;
  }
}

} // namespace adjmat
} // namespace PLMD

#endif //__PLUMED_adjmat_AdjacencyMatrixBase_h
