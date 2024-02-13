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
#include "RMSDShortcut.h"
#include "core/ActionRegister.h"
#include "core/ActionWithArguments.h"
#include "tools/Pdb.h"
#include "core/PlumedMain.h"

namespace PLMD {
namespace colvar {

PLUMED_REGISTER_ACTION(RMSDShortcut,"RMSD")

void RMSDShortcut::registerRMSD(Keywords& keys) {
  keys.add("compulsory","TYPE","SIMPLE","the manner in which RMSD alignment is performed.  Should be OPTIMAL or SIMPLE.");
  keys.addFlag("SQUARED",false," This should be setted if you want MSD instead of RMSD ");
}

void RMSDShortcut::registerKeywords(Keywords& keys){
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","REFERENCE","a file in pdb format containing the reference structure and the atoms involved in the CV");
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  keys.addFlag("NUMERICAL_DERIVATIVES", false, "calculate the derivatives for these quantities numerically");
  keys.addFlag("DISPLACEMENT",false,"Calculate the vector of displacements instead of the length of this vector");
  keys.add("compulsory","NUMBER","0","if there are multiple structures in the pdb file you can specify that you want the RMSD from a specific structure by specifying its place in the file here. If NUMBER=0 then the RMSD from all structures are computed");
  registerRMSD( keys );
}

void RMSDShortcut::createPosVector( const std::string& lab, const PDB& pdb ) {
  bool nopbc; parseFlag("NOPBC",nopbc); std::vector<AtomNumber> anum( pdb.getAtomNumbers() );
  if( !nopbc ) {
      std::string num; Tools::convert( anum[0].serial(), num ); std::string wm_line = "WHOLEMOLECULES ENTITY0=" + num;
      for(unsigned i=1; i<anum.size(); ++i) { Tools::convert( anum[i].serial(), num ); wm_line += "," + num; }
      readInputLine( wm_line );
  }
  std::string num; Tools::convert( anum[0].serial(), num ); std::string pos_line = lab + "_pos: POSITION NOPBC ATOMS=" + num;
  for(unsigned i=1; i<anum.size(); ++i) { Tools::convert( anum[i].serial(), num ); pos_line += "," + num; }
  readInputLine( pos_line );
  // Concatenate the three positions together
  readInputLine( lab + ": CONCATENATE ARG=" + lab + "_pos.x," + lab + "_pos.y," + lab + "_pos.z");
}

void RMSDShortcut::readAlignAndDisplace( ActionWithArguments* action, const bool& norm_weights, std::vector<double>& align, std::vector<double>& displace, std::vector<double>& sqrtdisplace ) {
  unsigned natoms = (action->getPntrToArgument(0))->getShape()[0] / 3;
  align.resize( natoms ); action->parseVector("ALIGN",align);
  displace.resize( natoms ); action->parseVector("DISPLACE",displace);
  double wa=0, wd=0; sqrtdisplace.resize( displace.size() );
  for(unsigned i=0; i<align.size(); ++i) { wa+=align[i]; wd+=displace[i]; }
      
  if( wa>epsilon ) {
      double iwa = 1. / wa;
      for(unsigned i=0; i<align.size(); ++i) align[i] *= iwa;
  } else {
      double iwa = 1. / natoms;
      for(unsigned i=0; i<align.size(); ++i) align[i] = iwa;
  }   
  if( wd>epsilon ) {
      if( !norm_weights ) { wd = 1; } double iwd = 1. / wd;
      for(unsigned i=0; i<align.size(); ++i) displace[i] *= iwd;
  } else {
      double iwd = 1. / natoms;
      for(unsigned i=0; i<align.size(); ++i) displace[i] = iwd;
  }
  for(unsigned i=0; i<align.size(); ++i) sqrtdisplace[i] = sqrt(displace[i]);
}

void RMSDShortcut::setReferenceConfiguration( const unsigned& num, const Value* refarg, const std::vector<double>& align, const std::vector<double>& displace, const std::string& type, const bool& norm_weights, PLMD::RMSD& myrmsd ) {
  unsigned natoms = refarg->getShape()[0] / 3; if( refarg->getRank()==2 ) natoms = refarg->getShape()[1] / 3;
  Vector center; std::vector<Vector> pos( natoms );
  for(unsigned i=0; i<pos.size(); ++i) {
      for(unsigned j=0; j<3; ++j) pos[i][j] = refarg->get( (3*num+j)*pos.size() + i );
      center+=pos[i]*align[i];
  }
  for(unsigned i=0; i<pos.size(); ++i) pos[i] -= center;
  myrmsd.clear(); myrmsd.set(align,displace,pos,type,true,norm_weights);
}

RMSDShortcut::RMSDShortcut(const ActionOptions& ao):
Action(ao),
ActionShortcut(ao)
{
  bool disp; parseFlag("DISPLACEMENT",disp);
  std::string reference; parse("REFERENCE",reference);
  // Read the reference pdb file
  PDB pdb; if( !pdb.read(reference,plumed.usingNaturalUnits(),0.1/plumed.getUnits().getLength()) ) plumed_merror("missing file " + reference );
  unsigned frame; parse("NUMBER",frame); unsigned nf=1;
  if( frame==0 ) { 
      FILE* fp=std::fopen(reference.c_str(),"r"); bool do_read=true; nf=0;
      while ( do_read ) {
          PDB mypdb; do_read=mypdb.readFromFilepointer(fp,plumed.usingNaturalUnits(),0.1/plumed.getUnits().getLength());
          if( !do_read && nf>0 ) break ;
          nf++;
      }
  } 
  // Now create the RMSD object
  std::string num, rmsd_line = getShortcutLabel() + ": ";
  if( nf==1 && !disp ) {
      rmsd_line += "RMSD_SCALAR REFERENCE=" + reference;
  } else { 
      std::string ffnum; Tools::convert( frame, ffnum );
      readInputLine( getShortcutLabel() + "_ref: PDB2CONSTANT REFERENCE=" + reference + " NUMBER=" + ffnum );
      createPosVector( getShortcutLabel() + "_pos", pdb );
      if( nf==1 && disp ) rmsd_line += "RMSD_DISPLACEMENT_VECTOR ARG=" + getShortcutLabel() + "_pos," + getShortcutLabel() + "_ref";
      // Now align
      std::vector<double> align( pdb.getOccupancy() ); Tools::convert( align[0], num ); rmsd_line += " ALIGN=" + num;
      for(unsigned i=1; i<align.size(); ++i) { Tools::convert( align[i], num ); rmsd_line += "," + num; }
      // And displace
      std::vector<double> displace( pdb.getBeta() ); Tools::convert( displace[0], num ); rmsd_line += " DISPLACE=" + num;
      for(unsigned i=1; i<displace.size(); ++i) { Tools::convert( displace[i], num ); rmsd_line += "," + num; }
  }
  // And create the RMSD object
  bool numder; parseFlag("NUMERICAL_DERIVATIVES",numder); 
  if(numder && nf==1 && !disp ) rmsd_line += " NUMERICAL_DERIVATIVES"; else if( numder ) error("can only use NUMERICAL_DERIVATIVES flag when RMSD is calculating a single scalar value");
  bool squared; parseFlag("SQUARED",squared); if(squared) rmsd_line += " SQUARED";
  bool nopbc; parseFlag("NOPBC",nopbc); if(nopbc) rmsd_line += " NOPBC";
  std::string tt; parse("TYPE",tt); readInputLine( rmsd_line + " TYPE=" + tt );
}

double RMSDShortcut::calculateDisplacement( const std::string& type, const std::vector<double>& align, const std::vector<double>& displace, const std::vector<double>& sqrtdisplace, 
                                            const std::vector<Vector>& pos, PLMD::RMSD& myrmsd, std::vector<Vector>& direction, std::vector<Vector>& der, const bool& squared ) {
  if(type=="SIMPLE") return myrmsd.simpleAlignment( align, displace, pos, myrmsd.getReference(), der, direction, squared );

  unsigned natoms = pos.size();
  Tensor rot; Matrix<std::vector<Vector> > DRotDPos(3,3); std::vector<Vector> centeredpos( natoms ), centeredreference( natoms );
  double r = myrmsd.calc_PCAelements( pos, der, rot, DRotDPos, direction, centeredpos, centeredreference, squared ); std::vector<Vector> ref( myrmsd.getReference() );
  for(unsigned i=0;i<direction.size();++i) direction[i] = sqrtdisplace[i]*( direction[i] - ref[i] );

  return r;
}

void RMSDShortcut::addDisplacementForces( const std::string& type, const std::vector<double>& align, const std::vector<double>& displace, const std::vector<double>& sqrtdisplace,
                                          const std::vector<Vector>& pos, PLMD::RMSD& myrmsd, std::vector<Vector>& direction, std::vector<Vector>& der, Value* myval, const bool& squared ) {
  unsigned natoms = pos.size();
  if( type=="SIMPLE" ) {
      double i = myrmsd.simpleAlignment( align, displace, pos, myrmsd.getReference(), der, direction, squared );
      Vector comforce; comforce.zero();
      for(unsigned i=0; i<natoms; i++) {
          for(unsigned k=0; k<3; ++k) comforce[k] += align[i]*myval->getForce( k*natoms + i);
      }
      for(unsigned i=0; i<natoms; i++) {
          for(unsigned k=0; k<3; ++k) myval->addForce( k*natoms + i, -comforce[k] );
      }
      return; 
  }
  Tensor rot; Matrix<std::vector<Vector> > DRotDPos(3,3); std::vector<Vector> centeredpos( natoms ), centeredreference( natoms );
  double r = myrmsd.calc_PCAelements( pos, der, rot, DRotDPos, direction, centeredpos, centeredreference, squared ); std::vector<Vector> ref( myrmsd.getReference() );
  Tensor trot=rot.transpose(); double prefactor = 1 / static_cast<double>( natoms ); Vector v1; v1.zero();
  for(unsigned n=0; n<natoms; n++) { 
       Vector ff; for(unsigned k=0; k<3; ++k ) ff[k] = myval->getForce( k*natoms + n );
       v1+=prefactor*matmul(trot,ff);
  }
  // Notice that we use centreredreference here to accumulate the true forces
  for(unsigned n=0; n<natoms; n++) { 
       Vector ff; for(unsigned k=0; k<3; ++k ) ff[k] = myval->getForce( k*natoms + n );
       centeredreference[n] = sqrtdisplace[n]*( matmul(trot,ff) - v1 );
  }
  for(unsigned a=0; a<3; a++) {
      for(unsigned b=0; b<3; b++) {
          double tmp1=0.; for(unsigned m=0; m<natoms; m++) tmp1+=centeredpos[m][b]*myval->getForce( a*natoms + m );
          for(unsigned i=0; i<natoms; i++) centeredreference[i] += sqrtdisplace[i]*tmp1*DRotDPos[a][b][i];
      }
  }
  // Now subtract the current force and add on the true force
  for(unsigned n=0; n<natoms; n++) {
      for(unsigned k=0; k<3; ++k) myval->addForce( k*natoms + n, centeredreference[n][k]-myval->getForce( k*natoms + n ) );
  }
}

}
}
