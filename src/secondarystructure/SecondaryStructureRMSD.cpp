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
#include "SecondaryStructureRMSD.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/GenericMolInfo.h"
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"

//+PLUMEDOC MCOLVAR SECONDARY_STRUCTURE_RMSD
/*
Calclulate the distance between segments of a protein and a reference structure of interest

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace secondarystructure {

PLUMED_REGISTER_ACTION(SecondaryStructureRMSD,"SECONDARY_STRUCTURE_RMSD");

bool SecondaryStructureRMSD::readShortcutWords( std::string& ltmap, ActionShortcut* action ) {
  action->parse("LESS_THAN",ltmap);
  if( ltmap.length()==0 ) {
    std::string nn, mm, d_0, r_0; action->parse("R_0",r_0);
    if( r_0.length()==0 ) r_0="0.08";
    action->parse("NN",nn); action->parse("D_0",d_0); action->parse("MM",mm);
    ltmap = "RATIONAL R_0=" + r_0 + " D_0=" + d_0 + " NN=" + nn + " MM=" + mm;
    return false;
  }
  return true;
}

void SecondaryStructureRMSD::expandShortcut( const bool& uselessthan, const std::string& labout, const std::string& labin, const std::string& ltmap, ActionShortcut* action ) {
  action->readInputLine( labout + "_lt: LESS_THAN ARG=" + labin + " SWITCH={" + ltmap  +"}");
  if( uselessthan ) action->readInputLine( labout + "_lessthan: SUM ARG=" + labout + "_lt PERIODIC=NO");
  else action->readInputLine( labout + ": SUM ARG=" + labout + "_lt PERIODIC=NO");
}

void SecondaryStructureRMSD::registerKeywords( Keywords& keys ) {
  ActionWithVector::registerKeywords( keys );
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions");
  keys.add("residues","RESIDUES","this command is used to specify the set of residues that could conceivably form part of the secondary structure. "
           "It is possible to use residues numbers as the various chains and residues should have been identified else using an instance of the "
           "\\ref MOLINFO action. If you wish to use all the residues from all the chains in your system you can do so by "
           "specifying all. Alternatively, if you wish to use a subset of the residues you can specify the particular residues "
           "you are interested in as a list of numbers. Please be aware that to form secondary structure elements your chain "
           "must contain at least N residues, where N is dependent on the particular secondary structure you are interested in. "
           "As such if you define portions of the chain with fewer than N residues the code will crash.");
  keys.add("atoms","ATOMS","this is the full list of atoms that we are investigating");
  keys.add("numbered","SEGMENT","this is the lists of atoms in the segment that are being considered");
  keys.add("compulsory","BONDLENGTH","the length to use for bonds");
  keys.add("numbered","STRUCTURE","the reference structure");
  keys.add("compulsory","TYPE","DRMSD","the manner in which RMSD alignment is performed. Should be OPTIMAL, SIMPLE or DRMSD. "
           "For more details on the OPTIMAL and SIMPLE methods see \\ref RMSD. For more details on the "
           "DRMSD method see \\ref DRMSD.");
  keys.add("optional","STRANDS_CUTOFF","If in a segment of protein the two strands are further apart then the calculation "
           "of the actual RMSD is skipped as the structure is very far from being beta-sheet like. "
           "This keyword speeds up the calculation enormously when you are using the LESS_THAN option. "
           "However, if you are using some other option, then this cannot be used");
  keys.add("optional","CUTOFF_ATOMS","the pair of atoms that are used to calculate the strand cutoff");
  keys.addFlag("VERBOSE",false,"write a more detailed output");
  keys.add("optional","LESS_THAN","calculate the number of a residue segments that are within a certain target distance of this secondary structure type. "
           "This quantity is calculated using \\f$\\sum_i \\sigma(s_i)\\f$, where \\f$\\sigma(s)\\f$ is a \\ref switchingfunction.");
  keys.add("optional","R_0","The r_0 parameter of the switching function.");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","NN","8","The n parameter of the switching function");
  keys.add("compulsory","MM","12","The m parameter of the switching function");
  keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
  keys.addOutputComponent("struct","default","the vectors containing the rmsd distances between the residues and each of the reference structures");
  keys.addOutputComponent("lessthan","default","the number blocks of residues that have an RMSD from the secondary structure that is less than the threshold");
  keys.needsAction("SECONDARY_STRUCTURE_RMSD"); keys.needsAction("LESS_THAN"); keys.needsAction("SUM");
}

void SecondaryStructureRMSD::readBackboneAtoms( ActionShortcut* action, PlumedMain& plumed, const std::string& moltype, std::vector<unsigned>& chain_lengths, std::string& all_atoms ) {
  auto* moldat=plumed.getActionSet().selectLatest<GenericMolInfo*>(action);
  if( ! moldat ) action->error("Unable to find MOLINFO in input");

  std::vector<std::string> resstrings; action->parseVector( "RESIDUES", resstrings );
  if(resstrings.size()==0) action->error("residues are not defined, check the keyword RESIDUES");
  else if( Tools::caseInSensStringCompare(resstrings[0], "all") ) {
    resstrings[0]="all";
    action->log.printf("  examining all possible secondary structure combinations\n");
  } else {
    action->log.printf("  examining secondary structure in residue positions : %s ",resstrings[0].c_str() );
    for(unsigned i=1; i<resstrings.size(); ++i) action->log.printf(", %s",resstrings[i].c_str() );
    action->log.printf("\n");
  }
  std::vector< std::vector<AtomNumber> > backatoms;
  moldat->getBackbone( resstrings, moltype, backatoms );

  chain_lengths.resize( backatoms.size() );
  for(unsigned i=0; i<backatoms.size(); ++i) {
    chain_lengths[i]=backatoms[i].size();
    for(unsigned j=0; j<backatoms[i].size(); ++j) {
      std::string bat_str; Tools::convert( backatoms[i][j].serial(), bat_str );
      if( i==0 && j==0 ) all_atoms = "ATOMS=" + bat_str;
      else all_atoms += "," + bat_str;
    }
  }
}


SecondaryStructureRMSD::SecondaryStructureRMSD(const ActionOptions&ao):
  Action(ao),
  ActionWithVector(ao),
  nopbc(false),
  align_strands(false),
  s_cutoff2(0),
  align_atom_1(0),
  align_atom_2(0)
{
  if( plumed.usingNaturalUnits() ) error("cannot use this collective variable when using natural units");

  parse("TYPE",alignType); parseFlag("NOPBC",nopbc);
  log.printf("  distances from secondary structure elements are calculated using %s algorithm\n",alignType.c_str() );
  log<<"  Bibliography "<<plumed.cite("Pietrucci and Laio, J. Chem. Theory Comput. 5, 2197 (2009)"); log<<"\n";

  parseFlag("VERBOSE",verbose_output);

  if( keywords.exists("STRANDS_CUTOFF") ) {
    double s_cutoff = 0;
    parse("STRANDS_CUTOFF",s_cutoff); align_strands=true;
    if( s_cutoff>0) {
      log.printf("  ignoring contributions from strands that are more than %f apart\n",s_cutoff);
      std::vector<unsigned> cutatoms; parseVector("CUTOFF_ATOMS",cutatoms);
      if( cutatoms.size()==2 ) {
        align_atom_1=cutatoms[0]; align_atom_2=cutatoms[1];
      } else error("did not find CUTOFF_ATOMS in input");
    }
    s_cutoff2=s_cutoff*s_cutoff;
  }

  // Read in the atoms
  std::vector<AtomNumber> all_atoms; parseAtomList("ATOMS",all_atoms); requestAtoms( all_atoms );

  for(unsigned i=1;; ++i) {
    std::vector<unsigned> newatoms;
    if( !parseNumberedVector("SEGMENT",i,newatoms) ) break;
    if( verbose_output ) {
      log.printf("  Secondary structure segment %u contains atoms : ", static_cast<unsigned>(colvar_atoms.size()+1));
      for(unsigned i=0; i<newatoms.size(); ++i) log.printf("%d ",all_atoms[newatoms[i]].serial() );
      log.printf("\n");
    }
    colvar_atoms.push_back( newatoms );
  }

  double bondlength; parse("BONDLENGTH",bondlength); bondlength=bondlength/getUnits().getLength();

  // Read in the reference structure
  for(unsigned ii=1;; ++ii) {
    std::vector<double> cstruct;
    if( !parseNumberedVector("STRUCTURE",ii,cstruct) ) break ;
    plumed_assert( cstruct.size()%3==0 && cstruct.size()/3==colvar_atoms[0].size() );
    std::vector<Vector> structure( cstruct.size()/3 );
    for(unsigned i=0; i<structure.size(); ++i) {
      for(unsigned j=0; j<3; ++j) structure[i][j] = 0.1*cstruct[3*i+j]/getUnits().getLength();
    }
    if( alignType=="DRMSD" ) {
      std::map<std::pair<unsigned,unsigned>, double> targets;
      for(unsigned i=0; i<structure.size()-1; ++i) {
        for(unsigned j=i+1; j<structure.size(); ++j) {
          double distance = delta( structure[i], structure[j] ).modulo();
          if(distance > bondlength) targets[std::make_pair(i,j)] = distance;
        }
      }
      drmsd_targets.push_back( targets );
    } else {
      Vector center; std::vector<double> align( structure.size(), 1.0 ), displace( structure.size(), 1.0 );
      for(unsigned i=0; i<structure.size(); ++i) center+=structure[i]*align[i];
      for(unsigned i=0; i<structure.size(); ++i) structure[i] -= center;
      RMSD newrmsd; newrmsd.clear();
      newrmsd.set(align,displace,structure,alignType,true,true);
      myrmsd.push_back( newrmsd );
    }
  }

  // And create values to hold everything
  unsigned nref = myrmsd.size(); if( alignType=="DRMSD" ) nref=drmsd_targets.size();
  plumed_assert( nref>0 );
  std::vector<unsigned> shape(1); shape[0]=colvar_atoms.size();
  if( nref==1 ) { addValue( shape ); setNotPeriodic(); }
  else {
    std::string num;
    for(unsigned i=0; i<nref; ++i) {
      Tools::convert( i+1, num ); addComponent( "struct-" + num, shape );
      componentIsNotPeriodic( "struct-" + num );
    }
  }
}

void SecondaryStructureRMSD::areAllTasksRequired( std::vector<ActionWithVector*>& task_reducing_actions ) {
  if( s_cutoff2>0 ) task_reducing_actions.push_back(this);
}

int SecondaryStructureRMSD::checkTaskStatus( const unsigned& taskno, int& flag ) const {
  if( s_cutoff2>0 ) {
    Vector distance=pbcDistance( ActionAtomistic::getPosition( getAtomIndex(taskno,align_atom_1) ),
                                 ActionAtomistic::getPosition( getAtomIndex(taskno,align_atom_2) ) );
    if( distance.modulo2()<s_cutoff2 ) return 1;
    return 0;
  } return flag;
}

void SecondaryStructureRMSD::calculate() {
  runAllTasks();
}

void SecondaryStructureRMSD::performTask( const unsigned& current, MultiValue& myvals ) const {
  // Resize the derivatives if need be
  unsigned nderi = 3*getNumberOfAtoms()+9;
  if( myvals.getNumberOfDerivatives()!=nderi ) myvals.resize( myvals.getNumberOfValues(), nderi, 0, 0 );
  // Retrieve the positions
  const unsigned natoms = colvar_atoms[current].size();
  std::vector<Vector> pos( natoms ), deriv( natoms );
  for(unsigned i=0; i<natoms; ++i) pos[i]=ActionAtomistic::getPosition( getAtomIndex(current,i) );

  // This aligns the two strands if this is required
  Vector distance=pbcDistance( pos[align_atom_1],pos[align_atom_2] );
  if( align_strands ) {
    Vector origin_old, origin_new; origin_old=pos[align_atom_2];
    origin_new=pos[align_atom_1]+distance;
    for(unsigned i=15; i<30; ++i) {
      pos[i]+=( origin_new - origin_old );
    }
  } else if( !nopbc ) {
    for(unsigned i=0; i<natoms-1; ++i) {
      const Vector & first (pos[i]);
      Vector & second (pos[i+1]);
      second=first+pbcDistance(first,second);
    }
  }
  // Create a holder for the derivatives
  if( alignType=="DRMSD" ) {
    // And now calculate the DRMSD
    const unsigned rs = drmsd_targets.size();
    for(unsigned i=0; i<rs; ++i) {
      double drmsd=0; Vector distance; Tensor vir; vir.zero();
      for(unsigned j=0; j<natoms; ++j) deriv[j].zero();
      for(const auto & it : drmsd_targets[i] ) {
        const unsigned k=it.first.first;
        const unsigned j=it.first.second;

        distance=delta( pos[k], pos[j] );
        const double len = distance.modulo();
        const double diff = len - it.second;
        const double der = diff / len;
        drmsd += diff*diff;

        if( !doNotCalculateDerivatives() ) {
          deriv[k] += -der*distance; deriv[j] += der*distance;
          vir += -der*Tensor(distance,distance);
        }
      }

      const double inpairs = 1./static_cast<double>(drmsd_targets[i].size());
      unsigned ostrn = getConstPntrToComponent(i)->getPositionInStream();
      drmsd = sqrt(inpairs*drmsd); myvals.setValue( ostrn, drmsd );

      if( !doNotCalculateDerivatives() ) {
        double scalef = inpairs / drmsd;
        for(unsigned j=0; j<natoms; ++j) {
          const unsigned ja = getAtomIndex( current, j );
          myvals.addDerivative( ostrn, 3*ja + 0, scalef*deriv[j][0] ); myvals.updateIndex( ostrn, 3*ja+0 );
          myvals.addDerivative( ostrn, 3*ja + 1, scalef*deriv[j][1] ); myvals.updateIndex( ostrn, 3*ja+1 );
          myvals.addDerivative( ostrn, 3*ja + 2, scalef*deriv[j][2] ); myvals.updateIndex( ostrn, 3*ja+2 );
        }
        unsigned nbase = myvals.getNumberOfDerivatives() - 9;
        for(unsigned k=0; k<3; ++k) {
          for(unsigned j=0; j<3; ++j) {
            myvals.addDerivative( ostrn, nbase + 3*k + j, scalef*vir(k,j) );
            myvals.updateIndex( ostrn, nbase + 3*k + j );
          }
        }
      }
    }
  } else {
    const unsigned rs = myrmsd.size();
    for(unsigned i=0; i<rs; ++i) {
      double nr = myrmsd[i].calculate( pos, deriv, false );
      unsigned ostrn = getConstPntrToComponent(i)->getPositionInStream();
      myvals.setValue( ostrn, nr );

      if( !doNotCalculateDerivatives() ) {
        Tensor vir; vir.zero();
        for(unsigned j=0; j<natoms; ++j) {
          const unsigned ja = getAtomIndex( current, j );
          myvals.addDerivative( ostrn, 3*ja + 0, deriv[j][0] ); myvals.updateIndex( ostrn, 3*colvar_atoms[current][j]+0 );
          myvals.addDerivative( ostrn, 3*ja + 1, deriv[j][1] ); myvals.updateIndex( ostrn, 3*colvar_atoms[current][j]+1 );
          myvals.addDerivative( ostrn, 3*ja + 2, deriv[j][2] ); myvals.updateIndex( ostrn, 3*colvar_atoms[current][j]+2 );
          vir+=(-1.0*Tensor( pos[j], deriv[j] ));
        }
        unsigned nbase = myvals.getNumberOfDerivatives() - 9;
        for(unsigned k=0; k<3; ++k) {
          for(unsigned j=0; j<3; ++j) {
            myvals.addDerivative( ostrn, nbase + 3*k + j, vir(k,j) );
            myvals.updateIndex( ostrn, nbase + 3*k + j );
          }
        }
      }
    }
  }
  return;
}

}
}
