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
#include "tools/TypeUtils.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <string>
#include <vector>
#include <type_traits>
#include <utility>

namespace PLMD {
namespace secondarystructure {

namespace helpers {
//this makes the variable needstype only needed to be declared explicitly
//needsType activate the reading of "TYPE" during the parsing
template <typename T, typename=void>
constexpr bool needsType=false;
//specialization to intercept needsType
template <typename T>
constexpr bool needsType<T,std::void_t<decltype(T::needsType)>> =T::needsType;

//this makes the variable needstype only needed to be declared explicitly
//needsBondLength activate the reading of "TYPE" during the parsing
template <typename T, typename=void>
constexpr bool needsBondLength=false;
//specialization to intercept needsBondLength
template <typename T>
constexpr bool needsBondLength<T,std::void_t<decltype(T::needsBondLength)>> =T::needsBondLength;
}

/// The largest number of atoms in a segment (the size of the position buffer in performTask).
/// ALPHARMSD, ANTIBETARMSD and PARABETARMSD all use segments of exactly this size
constexpr unsigned maxAtomsInSegment=30;

/// Base action for calculating things like AlphRMSD, AntibetaRMSD, etc
template <class T, typename myPTM=defaultPTM>
class SecondaryStructureBase: public ActionWithVector {
public:
  using input_type = T;
  using mytype=SecondaryStructureBase<T,myPTM>;
  using PTM =typename myPTM::template PTM<mytype>;
  typedef typename PTM::ParallelActionsInput ParallelActionsInput;
  typedef typename PTM::ParallelActionsOutput ParallelActionsOutput;
  using precision = cvprecision_t<T>;
  static constexpr size_t virialSize=9;
  static constexpr unsigned customGatherStep=3;
  static constexpr unsigned customGatherStopBefore=virialSize;
  static constexpr unsigned maxAtomsInSegment=secondarystructure::maxAtomsInSegment;

private:
  PTM taskmanager;
  static constexpr bool needsType = helpers::needsType<T>;
  static constexpr bool needsBondLength = helpers::needsBondLength<T>;
public:
  static void registerKeywords( Keywords& keys );
  static void readBackboneAtoms( ActionShortcut* action, PlumedMain& plumed, const std::string& backnames, std::vector<unsigned>& chain_lengths, std::vector<std::string>& all_atoms );
  /// Read LESS_THAN, or R_0, D_0, NN and MM, in a shortcut: ltmap is the switching function. The rational function defined
  /// by R_0, D_0, NN and MM gets the D_MAX that the other actions use for these keywords (where the function is 1e-5): it is
  /// then stretched to be zero at D_MAX. Returns true if LESS_THAN was used
  static bool readShortcutWords( std::string& ltmap, ActionShortcut* action );
  /// Read TYPE (returned in upper case), NOPBC and USEGPU in a shortcut: returns the name and the keywords of the action
  /// that computes the distances
  static std::string readDistanceAction( ActionShortcut* action, std::string& type, bool& nopbc, bool& usegpu );
  /// The pairs of strands of three residues that could form a parallel or antiparallel sheet with the given STYLE (all, inter
  /// or intra), in the order of the SEGMENT keywords of PARABETARMSD and ANTIBETARMSD. Each strand is given by the position
  /// of its first atom in the list of the backbone atoms (a strand is 15 consecutive atoms: N, CA, CB, C and O of three
  /// residues); chain_lengths contains the number of backbone atoms of each chain
  static void getSheetPairs( Action* action, const std::vector<unsigned>& chain_lengths, bool parallel, const std::string& style,
                             std::vector<std::pair<unsigned,unsigned> >& pairs );
  /// Append " SEGMENTnum=nlist[0],nlist[1],..." to line
  static void addSegment( std::string& line, unsigned num, const std::vector<unsigned>& nlist );
  /// Append " STRUCTUREnum=x0,y0,z0,x1,..." (in Angstrom) to line
  static void addStructure( std::string& line, unsigned num, const std::vector<Vector>& reference );
  /// Create the SUM of the switching functions in label_lt: the shortcut label if LESS_THAN was not used, label_lessthan otherwise
  static void addSum( ActionShortcut* action, bool uselessthan );
  explicit SecondaryStructureBase(const ActionOptions&);
  unsigned getNumberOfDerivatives() override ;
  void calculate() override;
  void getAtomicInputData( std::vector<double>& inputdata ) const override ;
  void getAtomicInputData( std::vector<float>& inputdata ) const override ;
  void applyNonZeroRankForces( std::vector<double>& outforces ) override ;
  static void performTask( unsigned task_index, const T& actiondata, ParallelActionsInput& input, ParallelActionsOutput& output );
  static int getNumberOfValuesPerTask( std::size_t task_index, const T& actiondata );
  static void getForceIndices( std::size_t task_index, std::size_t colno, std::size_t ntotal_force, const T& actiondata, const ParallelActionsInput& input, ForceIndexHolder force_indices );
};

//this is used to get some of the static non template dependant functionalities for Secondary structure
using SecondaryStructureInput=SecondaryStructureBase<Empty>;

template <class T,typename myPTM>
unsigned SecondaryStructureBase<T,myPTM>::getNumberOfDerivatives() {
  return 3*getNumberOfAtoms()+virialSize;
}

template <class T,typename myPTM>
bool SecondaryStructureBase<T,myPTM>::readShortcutWords( std::string& ltmap, ActionShortcut* action ) {
  action->parse("LESS_THAN",ltmap);
  if( ltmap.length()==0 ) {
    std::string nn, mm, d_0, r_0;
    action->parse("R_0",r_0);
    if( r_0.length()==0 ) {
      // The default is 0.08 nm. The distances are in the units that were set with UNITS, so it has to be converted
      Tools::convert( 0.08/action->plumed.getUnits().getLength(), r_0 );
    }
    action->parse("NN",nn);
    action->parse("D_0",d_0);
    action->parse("MM",mm);
    ltmap = "RATIONAL R_0=" + r_0 + " D_0=" + d_0 + " NN=" + nn + " MM=" + mm;
    // as SwitchingFunction::set(nn,mm,r0,d0)
    double r0, d0;
    int n, m;
    Tools::convert( r_0, r0 );
    Tools::convert( d_0, d0 );
    Tools::convert( nn, n );
    Tools::convert( mm, m );
    if( m==0 ) {
      m=2*n;
    }
    if( n==m ) {
      action->error("NN and MM must be different");
    }
    std::string dmax;
    Tools::convert( d0+r0*std::pow( 0.00001, 1.0/(n-m) ), dmax );
    ltmap += " D_MAX=" + dmax;
    return false;
  }
  return true;
}

template <class T,typename myPTM>
std::string SecondaryStructureBase<T,myPTM>::readDistanceAction( ActionShortcut* action, std::string& type, bool& nopbc, bool& usegpu ) {
  action->parse("TYPE",type);
  // RMSD only accepts upper case names
  std::transform( type.begin(), type.end(), type.begin(), [](unsigned char c) {
    return std::toupper(c);
  } );
  // BONDLENGTH is in nm, see the constructor
  std::string line = type=="DRMSD" ? "SECONDARY_STRUCTURE_DRMSD BONDLENGTH=0.17" : "SECONDARY_STRUCTURE_RMSD TYPE=" + type;
  action->parseFlag("NOPBC",nopbc);
  if( nopbc ) {
    line += " NOPBC";
  }
  action->parseFlag("USEGPU",usegpu);
  if( usegpu ) {
    line += " USEGPU";
  }
  return line;
}

template <class T,typename myPTM>
void SecondaryStructureBase<T,myPTM>::getSheetPairs( Action* action, const std::vector<unsigned>& chains, const bool parallel,
    const std::string& style, std::vector<std::pair<unsigned,unsigned> >& pairs ) {
  bool intra_chain(false), inter_chain(false);
  if( Tools::caseInSensStringCompare(style, "all") ) {
    intra_chain=true;
    inter_chain=true;
  } else if( Tools::caseInSensStringCompare(style, "inter") ) {
    inter_chain=true;
  } else if( Tools::caseInSensStringCompare(style, "intra") ) {
    intra_chain=true;
  } else {
    action->error( style + " is not a valid directive for the STYLE keyword");
  }
  pairs.resize(0);
  if( intra_chain ) {
    unsigned nprevious=0;
    for(unsigned i=0; i<chains.size(); ++i) {
      if( chains[i]<40 ) {
        action->error( std::string("segment of backbone is not long enough to form ") + (parallel ? "a parallel beta sheet" : "an antiparallel beta hairpin")
                       + ". Each backbone fragment must contain a minimum of 8 residues");
      }
      const unsigned nres=chains[i]/5;
      if( chains[i]%5!=0 ) {
        action->error("backbone segment received does not contain a multiple of five residues");
      }
      if( parallel ) {
        // the two strands are separated by at least three residues
        for(unsigned ires=0; ires+8<nres; ires++) {
          for(unsigned jres=ires+6; jres+2<nres; jres++) {
            pairs.emplace_back( nprevious+5*ires, nprevious+5*jres );
          }
        }
      } else {
        // the strands are i,i+1,i+2 and j-2,j-1,j, separated by at least two residues
        for(unsigned ires=0; ires+7<nres; ires++) {
          for(unsigned jres=ires+7; jres<nres; jres++) {
            pairs.emplace_back( nprevious+5*ires, nprevious+5*(jres-2) );
          }
        }
      }
      nprevious+=chains[i];
    }
  }
  if( inter_chain ) {
    if( chains.size()==1 && !Tools::caseInSensStringCompare(style, "all") ) {
      action->error("there is only one chain defined so cannot use inter_chain option");
    }
    for(unsigned ichain=1; ichain<chains.size(); ++ichain) {
      unsigned iprev=0;
      for(unsigned i=0; i<ichain; ++i) {
        iprev+=chains[i];
      }
      const unsigned inres=chains[ichain]/5;
      if( chains[ichain]%5!=0 ) {
        action->error("backbone segment received does not contain a multiple of five residues");
      }
      if( inres<3 ) {
        // A strand contains three residues. Without this the unsigned subtraction in the loop below wraps around
        continue;
      }
      for(unsigned ires=0; ires<inres-2; ++ires) {
        for(unsigned jchain=0; jchain<ichain; ++jchain) {
          unsigned jprev=0;
          for(unsigned i=0; i<jchain; ++i) {
            jprev+=chains[i];
          }
          const unsigned jnres=chains[jchain]/5;
          if( chains[jchain]%5!=0 ) {
            action->error("backbone segment received does not contain a multiple of five residues");
          }
          if( jnres<3 ) {
            continue;
          }
          for(unsigned jres=0; jres<jnres-2; ++jres) {
            pairs.emplace_back( iprev+5*ires, jprev+5*jres );
          }
        }
      }
    }
  }
  if( pairs.empty() ) {
    action->error("it is not possible to build any segment of backbone that could form a sheet with the residues and the STYLE that were given");
  }
}

template <class T,typename myPTM>
void SecondaryStructureBase<T,myPTM>::addSegment( std::string& line, unsigned num, const std::vector<unsigned>& nlist ) {
  // std::to_string is much faster than Tools::convert, and there can be hundreds of thousands of segments
  line += " SEGMENT" + std::to_string(num) + "=" + std::to_string(nlist[0]);
  for(unsigned i=1; i<nlist.size(); ++i) {
    line += "," + std::to_string(nlist[i]);
  }
}

template <class T,typename myPTM>
void SecondaryStructureBase<T,myPTM>::addStructure( std::string& line, unsigned num, const std::vector<Vector>& reference ) {
  line += " STRUCTURE" + std::to_string(num) + "=";
  std::string coord;
  for(unsigned i=0; i<reference.size(); ++i) {
    for(unsigned k=0; k<3; ++k) {
      Tools::convert( reference[i][k], coord );
      line += (i==0 && k==0 ? "" : ",") + coord;
    }
  }
}

template <class T,typename myPTM>
void SecondaryStructureBase<T,myPTM>::addSum( ActionShortcut* action, bool uselessthan ) {
  const std::string lab = action->getShortcutLabel();
  action->readInputLine( lab + (uselessthan ? "_lessthan" : "") + ": SUM ARG=" + lab + "_lt PERIODIC=NO");
}

template <class T,typename myPTM>
void SecondaryStructureBase<T,myPTM>::registerKeywords( Keywords& keys ) {
  ActionWithVector::registerKeywords( keys );
  T::registerKeywords( keys );
  PTM::registerKeywords( keys );
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions");
  keys.addInputKeyword("optional","MASK","vector","a vector which is used to determine which elements of the secondary structure variable should be computed");
  keys.add("atoms","ATOMS","this is the full list of atoms that we are investigating");
  keys.add("numbered","SEGMENT","the indices in the list given with ATOMS (starting from zero) of the atoms of a segment. All the segments must contain the same number of atoms (between 2 and 30)");
  if constexpr ( needsBondLength /*keys.getDisplayName()=="SECONDARY_STRUCTURE_DRMSD"*/ ) {
    keys.add("compulsory","BONDLENGTH","0.17","the pairs of atoms that are closer than this distance (in nm) in the reference structure are not used to compute the DRMSD");
  }
  keys.add("numbered","STRUCTURE","the coordinates (in Angstrom, three for each atom of a segment) of a reference structure. Use STRUCTURE1, STRUCTURE2, ... for multiple references");
  if constexpr ( needsType /*keys.getDisplayName()=="SECONDARY_STRUCTURE_RMSD"*/ ) {
    keys.add("compulsory","TYPE","OPTIMAL","the manner in which RMSD alignment is performed. Should be OPTIMAL or SIMPLE. "
             "For more details on the OPTIMAL and SIMPLE methods see RMSD.");
  } else if constexpr ( !needsBondLength /*keys.getDisplayName()!="SECONDARY_STRUCTURE_DRMSD"*/ ) {
    keys.add("compulsory","TYPE","DRMSD","the manner in which RMSD alignment is performed. Should be OPTIMAL, SIMPLE or DRMSD. "
             "For more details on the OPTIMAL and SIMPLE methods see RMSD. For more details on the "
             "DRMSD method see DRMSD.");
  }
  keys.addFlag("VERBOSE",false,"write a more detailed output");
  keys.reset_style("VERBOSE","hidden");
  if constexpr ( !needsBondLength && !needsType ) {
    keys.add("residues","RESIDUES","this command is used to specify the set of residues that could conceivably form part of the secondary structure. "
             "It is possible to use residues numbers as the various chains and residues should have been identified else using an instance of the "
             "MOLINFO action. If you wish to use all the residues from all the chains in your system you can do so by "
             "specifying all. Alternatively, if you wish to use a subset of the residues you can specify the particular residues "
             "you are interested in as a list of numbers. Please be aware that to form secondary structure elements your chain "
             "must contain at least N residues, where N is dependent on the particular secondary structure you are interested in. "
             "As such if you define portions of the chain with fewer than N residues PLUMED stops with an error.");
    keys.add("optional","LESS_THAN","calculate the number of segments of residues that are within a certain target distance of this secondary structure type. "
             "This quantity is calculated using $\\sum_i \\sigma(s_i)$, where $\\sigma(s)$ is the switching function given here. "
             "With this keyword the value of the action is the vector of the distances and the number of segments is the lessthan component. "
             "Give a D_MAX in this switching function (e.g. LESS_THAN={RATIONAL R_0=0.08 NN=8 MM=12 D_MAX=1.42}) to make it zero beyond D_MAX, "
             "as the one defined by R_0, NN, MM and D_0");
    keys.add("optional","R_0","The r_0 parameter of the switching function, in the length units of UNITS. The default is 0.08 nm. "
             "The switching function is stretched to be zero at the distance where it is 1e-5 (its D_MAX in the expanded input)");
    keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
    keys.add("compulsory","NN","8","The n parameter of the switching function");
    keys.add("compulsory","MM","12","The m parameter of the switching function");
  }
  keys.addFlag("ALIGN_STRANDS",false,"the segments contain two strands of 15 atoms each: make the two strands whole separately and move the second one to the periodic image where its 7th atom (the 22nd of the segment) is closest to the 7th atom of the first strand");
  if constexpr ( !needsBondLength && !needsType ) {
    // ALPHARMSD, ANTIBETARMSD and PARABETARMSD
    keys.addOutputComponent("lessthan","LESS_THAN","scalar","the number of segments of residues that have a distance from the secondary structure that is less than the threshold given with LESS_THAN");
  } else {
    keys.setValueDescription("vector","the distances between each segment and the reference structure (if only one STRUCTURE is given)");
    keys.addOutputComponent("struct","STRUCTURE","vector","struct-n is the vector of distances between each segment and the n-th reference structure (if more than one STRUCTURE is given)");
  }
  keys.needsAction("SECONDARY_STRUCTURE_RMSD");
  keys.needsAction("SECONDARY_STRUCTURE_DRMSD");
  keys.needsAction("LESS_THAN");
  keys.needsAction("SUM");
  keys.addDOI("10.1021/ct900202f");
}

template <class T,typename myPTM>
void SecondaryStructureBase<T,myPTM>::readBackboneAtoms( ActionShortcut* action, PlumedMain& plumed, const std::string& moltype, std::vector<unsigned>& chain_lengths, std::vector<std::string>& all_atoms ) {
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
  if( all_atoms.empty() ) {
    action->error("no backbone atoms were found for the residues specified with RESIDUES");
  }
}

template <class T,typename myPTM>
SecondaryStructureBase<T,myPTM>::SecondaryStructureBase(const ActionOptions&ao):
  Action(ao),
  ActionWithVector(ao),
  taskmanager(this) {
  if( plumed.usingNaturalUnits() ) {
    error("cannot use this collective variable when using natural units");
  }

  input_type myinput;
  parseFlag("NOPBC",myinput.nopbc);
  std::string alignType="";
  if constexpr ( needsType ) {
    parse("TYPE",alignType);
  }
  log.printf("  distances from secondary structure elements are calculated using %s algorithm\n", needsType ? alignType.c_str() : "DRMSD" );
  log<<"  Bibliography "<<plumed.cite("Pietrucci and Laio, J. Chem. Theory Comput. 5, 2197 (2009)");
  log<<"\n";

  parseFlag("ALIGN_STRANDS",myinput.align_strands);
  bool verbose_output=false;
  parseFlag("VERBOSE",verbose_output);
  if( myinput.nopbc ) {
    log.printf("  without periodic boundary conditions\n");
  } else if( myinput.align_strands ) {
    log.printf("  the two strands of each segment are made whole separately, and atom 22 is moved to the periodic image that is closest to atom 7\n");
  } else {
    log.printf("  each segment is made whole before computing the distance\n");
  }

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
    // The segments are copied into a rectangular matrix and the positions of their atoms into a
    // buffer of fixed size, so everything is checked here rather than in performTask
    if( newatoms.size()<2 || newatoms.size()>maxAtomsInSegment ) {
      error("SEGMENT" + std::to_string(i) + " contains " + std::to_string(newatoms.size()) + " atoms. A segment must contain between 2 and "
            + std::to_string(maxAtomsInSegment) + " atoms");
    }
    if( !colvar_atoms.empty() && newatoms.size()!=colvar_atoms[0].size() ) {
      error("all the SEGMENT keywords must contain the same number of atoms but SEGMENT1 contains "
            + std::to_string(colvar_atoms[0].size()) + " atoms and SEGMENT" + std::to_string(i) + " contains " + std::to_string(newatoms.size()) );
    }
    for(unsigned ii=0; ii<newatoms.size(); ++ii) {
      if( newatoms[ii]>=all_atoms.size() ) {
        error("SEGMENT" + std::to_string(i) + " contains index " + std::to_string(newatoms[ii]) + " but only " + std::to_string(all_atoms.size())
              + " atoms were specified with ATOMS (indices start from zero)");
      }
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
  if constexpr ( needsBondLength ) {
    parse("BONDLENGTH",bondlength);
    bondlength=bondlength/getUnits().getLength();
  }
  // Read in the reference structure
  for(unsigned ii=1;; ++ii) {
    std::vector<double> cstruct;
    if( !parseNumberedVector("STRUCTURE",ii,cstruct) ) {
      break ;
    }
    if( cstruct.size()%3!=0 || cstruct.size()/3!=colvar_atoms[0].size() ) {
      error("STRUCTURE" + std::to_string(ii) + " contains " + std::to_string(cstruct.size()) + " numbers. It must contain the three coordinates of each of the "
            + std::to_string(colvar_atoms[0].size()) + " atoms in a segment");
    }
    std::vector<Vector> structure( cstruct.size()/3 );
    for(unsigned i=0; i<structure.size(); ++i) {
      for(unsigned j=0; j<3; ++j) {
        structure[i][j] = 0.1*cstruct[3*i+j]/getUnits().getLength();
      }
    }
    myinput.setReferenceStructure( alignType, bondlength, structure );
  }

  // And create values to hold everything
  if( myinput.nstructures==0 ) {
    error("no reference structure was found in input -- use STRUCTURE");
  }
  if( myinput.align_strands && colvar_atoms[0].size()!=maxAtomsInSegment ) {
    error("ALIGN_STRANDS assumes that the two strands are the first and second halves of a segment containing " + std::to_string(maxAtomsInSegment)
          + " atoms, but the segments contain " + std::to_string(colvar_atoms[0].size()) + " atoms");
  }
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

template <class T,typename myPTM>
void SecondaryStructureBase<T,myPTM>::calculate() {
  taskmanager.runAllTasks();
}

template <class T,typename myPTM>
void SecondaryStructureBase<T,myPTM>::getAtomicInputData( std::vector<double>& inputdata ) const {
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

template <class T,typename myPTM>
void SecondaryStructureBase<T,myPTM>::getAtomicInputData( std::vector<float>& inputdata ) const {
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

template <class T,typename myPTM>
void SecondaryStructureBase<T,myPTM>::performTask( unsigned task_index, const T& actiondata, ParallelActionsInput& input, ParallelActionsOutput& output ) {
  // Empty is only used to reach the static shortcut helpers (see the
  // SecondaryStructureInput alias); this body is never called for it. MSVC
  // instantiates member bodies eagerly, hence the guard.
  if constexpr ( !std::is_same_v<T,Empty> ) {
    // The constructor checks that natoms<=maxAtomsInSegment
    std::array<Vector,maxAtomsInSegment> pos;

    for(unsigned i=0; i<actiondata.natoms; ++i) {
      const unsigned atno = actiondata.colvar_atoms(task_index,i);
      pos[i][0] = input.inputdata[3*atno+0];
      pos[i][1] = input.inputdata[3*atno+1];
      pos[i][2] = input.inputdata[3*atno+2];
    }
    if( !actiondata.nopbc ) {
      // Reconstruct the segment by making each pair of consecutive atoms whole. The two strands of a sheet
      // (atoms 0-14 and atoms 15-29) are made whole separately, because the last atom of the first strand
      // and the first atom of the second strand are not bonded
      const unsigned nstrand = actiondata.align_strands ? actiondata.natoms/2 : actiondata.natoms;
      for(unsigned start=0; start<actiondata.natoms; start+=nstrand) {
        for(unsigned i=start; i+1<start+nstrand; ++i) {
          const Vector & first (pos[i]);
          Vector & second (pos[i+1]);
          second=first+input.pbc->distance(first,second);
        }
      }
      // This aligns the two strands if this is required: atom 22 is moved to the periodic image that is
      // closest to atom 7 and the rest of the second strand follows
      if( actiondata.align_strands ) {
        Vector distance=input.pbc->distance( pos[6],pos[21] );
        distance += pos[6] - pos[21];
        for(unsigned i=actiondata.natoms/2; i<actiondata.natoms; ++i) {
          pos[i]+=distance;
        }
      }
    }

    // Create a holder for the derivatives. The number of derivatives per component has to be the one that
    // was given to setupParallelTaskManager (3*natoms+virialSize), and the view must contain just the atoms of
    // the segment: the size of the buffer is only an upper limit.
    const unsigned rs = actiondata.nstructures;
    ColvarOutput rmsd_output( output.values,
                              3*actiondata.natoms+virialSize,
                              output.derivatives.data() );
    // And now calculate the DRMSD
    for(unsigned i=0; i<rs; ++i) {
      T::calculateDistance( i, input.noderiv, actiondata, View{pos.data(),actiondata.natoms}, rmsd_output );
    }
  }
}

template <class T,typename myPTM>
void SecondaryStructureBase<T,myPTM>::applyNonZeroRankForces( std::vector<double>& outforces ) {
  taskmanager.applyForces( outforces );
}

template <class T,typename myPTM>
int SecondaryStructureBase<T,myPTM>::getNumberOfValuesPerTask( std::size_t task_index,
    const T& actiondata ) {
  return 1;
}

template <class T,typename myPTM>
void SecondaryStructureBase<T,myPTM>::getForceIndices( std::size_t task_index,
    std::size_t /* colno */,
    std::size_t ntotal_force,
    const T& actiondata,
    const ParallelActionsInput& input,
    ForceIndexHolder force_indices ) {
  // Guarded for the same reason as performTask().
  if constexpr ( !std::is_same_v<T,Empty> ) {
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
}

#endif
