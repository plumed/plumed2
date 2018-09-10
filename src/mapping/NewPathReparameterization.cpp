/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2018 The plumed team
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
#include "core/ActionPilot.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionRegister.h"
#include "setup/SetupReferenceBase.h"

namespace PLMD {
namespace mapping {

class NPathReparameterization : public ActionPilot {
private:
/// The object that calculates the distances between the frames
  PlumedMain metric; 
/// Number of cycles of the optimization algorithm to run
  unsigned maxcycles;
/// The points on the path to fix
  unsigned ifix1, ifix2;
/// Tolerance for the minimization algorithm
  double TOL;
/// Tempory masses and charges for the atoms and the displacement
  std::vector<double> masses, charges, data;
/// Positions and forces
  std::vector<Vector> positions, forces;
/// Used to store current spacing between frames in path
  std::vector<double> len, sumlen, sfrac;
/// The list of setup actions that contain the details of the path
  std::vector<setup::SetupReferenceBase*> mypath;
///
  bool loopEnd( const int& index, const int& end, const int& inc ) const ;
///
  void getDisplaceVector( const unsigned& ifrom, const unsigned& ito );
///
  double computeSpacing( const unsigned& ifrom, const unsigned& ito );
///
  void calcCurrentPathSpacings( const int& istart, const int& iend );
///
  void reparameterizePart( const int& istart, const int& iend, const double& target );
public:
  static void registerKeywords( Keywords& keys );
  NPathReparameterization(const ActionOptions&);
  void calculate(){}
  void apply(){} 
  void update();
};

PLUMED_REGISTER_ACTION(NPathReparameterization,"REPARAMETERIZE_PATH")

void NPathReparameterization::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys ); ActionPilot::registerKeywords( keys );
  keys.add("numbered","FRAME","list of frames that are to be adjusted and recentered on the path");
  keys.add("compulsory","METRIC","the way to compute the distances between the frames");
  keys.add("compulsory","FIXED","0","the frames in the path to fix");
  keys.add("compulsory","MAXCYLES","100","number of cycles of the algorithm to run");
  keys.add("compulsory","TOL","1E-4","the tolerance to use for the path reparameterization algorithm");
}

NPathReparameterization::NPathReparameterization(const ActionOptions&ao):
Action(ao),
ActionPilot(ao)
{
  parse("MAXCYLES",maxcycles); parse("TOL",TOL); unsigned natoms, nargs;
  log.printf("  running till change is less than %f or until there have been %d optimization cycles \n", TOL, maxcycles);
  // Setup the reference frames
  for(unsigned i=1;;++i) {
     std::string fname;
     if( !parseNumbered("FRAME",i,fname) ) break;

     // Now get a reference object to readjust
     setup::SetupReferenceBase* myset = plumed.getActionSet().selectWithLabel<setup::SetupReferenceBase*>( fname );
     if( !myset ) error("action with label " + fname + " is not a reference object");
     unsigned tatoms, targs; myset->getNatomsAndNargs( tatoms, targs );
     if( i==1 ) { natoms=tatoms; nargs=targs; }
     else if( natoms!=tatoms || nargs!=targs ) error("mismatched reference configurations");

     log.printf("  %dth reference configuration is in action %s \n", i, fname.c_str() );
     mypath.push_back( myset );
  }
  len.resize( mypath.size() ); sumlen.resize( mypath.size() ); sfrac.resize( mypath.size() );
  std::vector<unsigned> fixed; parseVector("FIXED",fixed);
  if( fixed.size()==1 ) {
      if( fixed[0]!=0 ) error("input to FIXED should be two integers");
      ifix1=0; ifix2=mypath.size()-1;
  } else if( fixed.size()==2 ) {
      if( fixed[0]<0 || fixed[1]<0 || fixed[0]>(mypath.size()-1) || fixed[1]>(mypath.size()-1) ) {
        error("input to FIXED should be two numbers between 0 and the number of frames-1");
      }
      ifix1=fixed[0]; ifix2=fixed[1];
  } else error("input to FIXED should be two integers");
  log.printf("  fixing frames %d and %d when reparameterizing \n", ifix1, ifix2 );
  
  // Now set up the metric object 
  int s=sizeof(double);
  metric.cmd("setRealPrecision",&s);
  metric.cmd("setNoVirial");
  metric.cmd("setMDEngine","plumed");
  int nat=2*natoms; metric.cmd("setNatoms",&nat);
  positions.resize(nat); masses.resize(nat); forces.resize(nat); charges.resize(nat);
  if( nargs>0 ) {
      std::vector<int> size(2); size[0]=1; size[1]=nargs;
      metric.cmd("createValue arg1",&size[0]); metric.cmd("setValueNotPeriodic arg1");
      metric.cmd("createValue arg2",&size[0]); metric.cmd("setValueNotPeriodic arg2");
  }
  std::string mtype; parse("METRIC",mtype); const char* cinp=mtype.c_str();
  std::vector<std::string> input=Tools::getWords(mtype);
  if( mtype.size()==1 && !actionRegister().check(input[0]) ) {
      metric.cmd("setPlumedDat",cinp); metric.cmd("init");
  } else {
      metric.cmd("init"); metric.cmd("readInputLine",cinp);
  }
  // Now setup stuff to retrieve the final displacement
  ActionWithValue* fav = dynamic_cast<ActionWithValue*>( metric.getActionSet()[metric.getActionSet().size()-1].get() );
  if( !fav ) error("final value should calculate relevant value that you want as reference");
  std::string name = (fav->copyOutput(0))->getName(); long rank; metric.cmd("getDataRank " + name, &rank );
  if( rank==0 ) rank=1;
  std::vector<long> ishape( rank ); metric.cmd("getDataShape " + name, &ishape[0] );
  unsigned nvals=1; for(unsigned i=0;i<ishape.size();++i) nvals *= ishape[i];
  data.resize( nvals ); metric.cmd("setMemoryForData " + name, &data[0] );
}

bool NPathReparameterization::loopEnd( const int& index, const int& end, const int& inc ) const {
  if( inc>0 && index<end ) return false;
  else if( inc<0 && index>end ) return false;
  return true;
}

void NPathReparameterization::getDisplaceVector( const unsigned& ifrom, const unsigned& ito ) {
  int step = getStep(); metric.cmd("setStep",&step);
  mypath[ifrom]->transferDataToPlumed( 0, masses, charges, positions, "arg1", metric );
  mypath[ito]->transferDataToPlumed( positions.size()/2, masses, charges, positions, "arg2", metric );
  metric.cmd("setMasses",&masses[0]); metric.cmd("setCharges",&charges[0]);
  metric.cmd("setPositions",&positions[0]); metric.cmd("setForces",&forces[0]);
  Tensor box( plumed.getAtoms().getPbc().getBox() ); metric.cmd("setBox",&box[0][0]);
  metric.cmd("calc");  
}

double NPathReparameterization::computeSpacing( const unsigned& ifrom, const unsigned& ito ) {
  getDisplaceVector( ifrom, ito ); double len=0; for(unsigned i=0;i<data.size();++i) len += data[i]*data[i]; 
  return sqrt( len );
}

void NPathReparameterization::calcCurrentPathSpacings( const int& istart, const int& iend ) {
  plumed_dbg_assert( istart<len.size() && iend<len.size() );
  len[istart] = sumlen[istart]=0;
  //printf("HELLO PATH SPACINGS ARE CURRENTLY \n");

  // Get the spacings given we can go forward and backwards
  int incr=1; if( istart>iend ) { incr=-1; }

  for(int i=istart+incr; loopEnd(i,iend+incr,incr)==false; i+=incr) {
    len[i] = computeSpacing( i-incr, i ); sumlen[i] = sumlen[i-incr] + len[i];
    //printf("FRAME %d TO FRAME %d EQUALS %f : %f \n",i-incr,i,len[i],sumlen[i] );
  }
}

void NPathReparameterization::reparameterizePart( const int& istart, const int& iend, const double& target ) {
  calcCurrentPathSpacings( istart, iend ); unsigned cfin;
  // If a target separation is set we fix where we want the nodes
  int incr=1; if( istart>iend ) { incr=-1; }

  if( target>0 ) {
    if( iend>istart ) {
      for(unsigned i=istart; i<iend+1; ++i) sfrac[i] = target*(i-istart);
    } else {
      for(int i=istart-1; i>iend-1; --i) sfrac[i]=target*(istart-i);
    }
    cfin = iend+incr;
  } else {
    cfin = iend;
  }

  double prevsum=0.; Matrix<double> newmatrix( mypath.size(), data.size() );
  for(unsigned iter=0; iter<maxcycles; ++iter) {
    if( fabs(sumlen[iend] - prevsum)<=TOL ) break ;
    prevsum = sumlen[iend];
    // If no target is set we redistribute length
    if( target<0 ) {
      plumed_assert( istart<iend );
      double dr = sumlen[iend] / static_cast<double>( iend - istart );
      for(unsigned i=istart; i<iend; ++i) sfrac[i] = dr*(i-istart);
    }

    // Now compute positions of new nodes in path
    for(int i=istart+incr; loopEnd(i,cfin,incr)==false; i+=incr) {
      int k = istart;
      while( !((sumlen[k] < sfrac[i]) && (sumlen[k+incr]>=sfrac[i])) ) {
        k+=incr;
        if( cfin==iend && k>= iend+1 ) plumed_merror("path reparameterization error");
        else if( cfin==(iend+1) && k>=iend ) { k=iend-1; break; }
        else if( cfin==(iend-1) && k<=iend ) { k=iend+1; break; }
      }
      double dr = (sfrac[i]-sumlen[k])/len[k+incr];
      // Copy the reference configuration to the row of a matrix
      mypath[k]->getReferenceConfiguration( data );
      for(unsigned j=0;j<data.size();++j) newmatrix(i,j) = data[j]; 
      getDisplaceVector( k, k+incr );
      // Shift the reference configuration by this ammount
      for(unsigned j=0;j<data.size();++j) newmatrix(i,j) += dr*data[j];
    }

    // Copy the positions of the new path to the new paths
    for(int i=istart+incr; loopEnd(i,cfin,incr)==false; i+=incr) {
      for(unsigned j=0;j<data.size();++j) data[j] = newmatrix(i,j);
      mypath[i]->setReferenceConfiguration( data );
    }

    // Recompute the separations between frames
    calcCurrentPathSpacings( istart, iend );
  }
}

void NPathReparameterization::update() {
  // First reparameterize the part between the fixed frames
  reparameterizePart( ifix1, ifix2, -1.0 );

  // Get the separation between frames which we will use to set the remaining frames
  double target = sumlen[ifix2] / ( ifix2 - ifix1 );

  // And reparameterize the begining and end of the path
  if( ifix1>0 ) reparameterizePart( ifix1, 0, target );
  if( ifix2<(mypath.size()-1) ) reparameterizePart( ifix2, mypath.size()-1, target );
}

}
}
