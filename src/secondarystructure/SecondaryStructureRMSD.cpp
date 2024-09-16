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

void SecondaryStructureRMSD::registerKeywords( Keywords& keys ) {
  ActionWithVector::registerKeywords( keys ); keys.use("MASK");
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
  keys.addFlag("VERBOSE",false,"write a more detailed output");
  keys.add("optional","LESS_THAN","calculate the number of a residue segments that are within a certain target distance of this secondary structure type. "
           "This quantity is calculated using \\f$\\sum_i \\sigma(s_i)\\f$, where \\f$\\sigma(s)\\f$ is a \\ref switchingfunction.");
  keys.add("optional","R_0","The r_0 parameter of the switching function.");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","NN","8","The n parameter of the switching function");
  keys.add("compulsory","MM","12","The m parameter of the switching function");
  keys.addFlag("ALIGN_STRANDS",false,"ensure that the two halves of a beta sheet are not broken by the periodic boundaries before doing alignment");
  keys.addOutputComponent("struct","default","the vectors containing the rmsd distances between the residues and each of the reference structures");
  keys.addOutputComponent("lessthan","default","the number blocks of residues that have an RMSD from the secondary structure that is less than the threshold");
  keys.needsAction("SECONDARY_STRUCTURE_RMSD"); keys.needsAction("LESS_THAN"); keys.needsAction("SUM");
}

void SecondaryStructureRMSD::readBackboneAtoms( ActionShortcut* action, PlumedMain& plumed, const std::string& moltype, std::vector<unsigned>& chain_lengths, std::vector<std::string>& all_atoms ) {
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
      all_atoms.push_back( bat_str );
    }
  }
}


SecondaryStructureRMSD::SecondaryStructureRMSD(const ActionOptions&ao):
  Action(ao),
  ActionWithVector(ao),
  nopbc(false)
{
  if( plumed.usingNaturalUnits() ) error("cannot use this collective variable when using natural units");

  parse("TYPE",alignType); parseFlag("NOPBC",nopbc);
  log.printf("  distances from secondary structure elements are calculated using %s algorithm\n",alignType.c_str() );
  log<<"  Bibliography "<<plumed.cite("Pietrucci and Laio, J. Chem. Theory Comput. 5, 2197 (2009)"); log<<"\n";

  parseFlag("VERBOSE",verbose_output); parseFlag("ALIGN_STRANDS",align_strands);
  log.printf("  ensuring atoms 7 and 22 in each residue are not separated by pbc before doing alignment\n");

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
  for(unsigned i=0; i<getNumberOfComponents(); ++i) getPntrToComponent(i)->setDerivativeIsZeroWhenValueIsZero();
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
  if( align_strands ) {
    Vector distance=pbcDistance( pos[6],pos[21] );
    Vector origin_old, origin_new; origin_old=pos[21];
    origin_new=pos[6]+distance;
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
      drmsd = sqrt(inpairs*drmsd); myvals.setValue( i, drmsd );

      if( !doNotCalculateDerivatives() ) {
        double scalef = inpairs / drmsd;
        for(unsigned j=0; j<natoms; ++j) {
          const unsigned ja = getAtomIndex( current, j );
          myvals.addDerivative( i, 3*ja + 0, scalef*deriv[j][0] ); myvals.updateIndex( i, 3*ja+0 );
          myvals.addDerivative( i, 3*ja + 1, scalef*deriv[j][1] ); myvals.updateIndex( i, 3*ja+1 );
          myvals.addDerivative( i, 3*ja + 2, scalef*deriv[j][2] ); myvals.updateIndex( i, 3*ja+2 );
        }
        unsigned nbase = myvals.getNumberOfDerivatives() - 9;
        for(unsigned k=0; k<3; ++k) {
          for(unsigned j=0; j<3; ++j) {
            myvals.addDerivative( i, nbase + 3*k + j, scalef*vir(k,j) );
            myvals.updateIndex( i, nbase + 3*k + j );
          }
        }
      }
    }
  } else {
    const unsigned rs = myrmsd.size();
    for(unsigned i=0; i<rs; ++i) {
      double nr = myrmsd[i].calculate( pos, deriv, false );
      myvals.setValue( i, nr );

      if( !doNotCalculateDerivatives() ) {
        Tensor vir; vir.zero();
        for(unsigned j=0; j<natoms; ++j) {
          const unsigned ja = getAtomIndex( current, j );
          myvals.addDerivative( i, 3*ja + 0, deriv[j][0] ); myvals.updateIndex( i, 3*colvar_atoms[current][j]+0 );
          myvals.addDerivative( i, 3*ja + 1, deriv[j][1] ); myvals.updateIndex( i, 3*colvar_atoms[current][j]+1 );
          myvals.addDerivative( i, 3*ja + 2, deriv[j][2] ); myvals.updateIndex( i, 3*colvar_atoms[current][j]+2 );
          vir+=(-1.0*Tensor( pos[j], deriv[j] ));
        }
        unsigned nbase = myvals.getNumberOfDerivatives() - 9;
        for(unsigned k=0; k<3; ++k) {
          for(unsigned j=0; j<3; ++j) {
            myvals.addDerivative( i, nbase + 3*k + j, vir(k,j) );
            myvals.updateIndex( i, nbase + 3*k + j );
          }
        }
      }
    }
  }
  return;
}

}
}
