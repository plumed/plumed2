/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2020 The plumed team
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
#include "PathReparameterization.h"

namespace PLMD {
namespace mapping {

PathReparameterization::PathReparameterization( const Pbc& ipbc, const std::vector<Value*>& iargs, std::vector<std::unique_ptr<ReferenceConfiguration>>& pp ):
  mydpack( 1, pp[0]->getNumberOfReferenceArguments() + 3*pp[0]->getNumberOfReferencePositions() + 9 ),
  mypack( pp[0]->getNumberOfReferenceArguments(), pp[0]->getNumberOfReferencePositions(), mydpack ),
  mydir(ReferenceConfigurationOptions("DIRECTION")),
  pbc(ipbc),
  args(iargs),
  mypath(pp),
  len(pp.size()),
  sumlen(pp.size()),
  sfrac(pp.size()),
  MAXCYCLES(100)
{
  mypdb.setAtomNumbers(  pp[0]->getAbsoluteIndexes() ); mypdb.addBlockEnd( pp[0]->getAbsoluteIndexes().size() );
  if( pp[0]->getArgumentNames().size()>0 ) mypdb.setArgumentNames( pp[0]->getArgumentNames() );
  mydir.read( mypdb ); mydir.zeroDirection(); pp[0]->setupPCAStorage( mypack );
}

bool PathReparameterization::loopEnd( const int& index, const int& end, const int& inc ) const {
  if( inc>0 && index<end ) return false;
  else if( inc<0 && index>end ) return false;
  return true;
}

void PathReparameterization::calcCurrentPathSpacings( const int& istart, const int& iend ) {
  plumed_dbg_assert( istart<len.size() && iend<len.size() );
  len[istart] = sumlen[istart]=0;
  //printf("HELLO PATH SPACINGS ARE CURRENTLY \n");

  // Get the spacings given we can go forward and backwards
  int incr=1; if( istart>iend ) { incr=-1; }

  for(int i=istart+incr; loopEnd(i,iend+incr,incr)==false; i+=incr) {
    len[i] = mypath[i-incr]->calc( mypath[i]->getReferencePositions(), pbc, args, mypath[i]->getReferenceArguments(), mypack, false );
    sumlen[i] = sumlen[i-incr] + len[i];
    //printf("FRAME %d TO FRAME %d EQUALS %f : %f \n",i-incr,i,len[i],sumlen[i] );
  }
}

void PathReparameterization::reparameterizePart( const int& istart, const int& iend, const double& target, const double& TOL ) {
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

  std::vector<Direction> newpath;
  for(unsigned i=0; i<mypath.size(); ++i) {
    newpath.push_back( Direction(ReferenceConfigurationOptions("DIRECTION")) ); newpath[i].read( mypdb );
  }

  double prevsum=0.;
  for(unsigned iter=0; iter<MAXCYCLES; ++iter) {
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
      // Calculate the displacement between the appropriate points
      // double dd = mypath[k]->calc( mypath[k+incr]->getReferencePositions(), pbc, args, mypath[k+incr]->getReferenceArguments(), mypack, true );
      // Copy the reference configuration from the configuration to a tempory direction
      newpath[i].setDirection( mypath[k]->getReferencePositions(), mypath[k]->getReferenceArguments() );
      // Get the displacement of the path
      mypath[k]->extractDisplacementVector( mypath[k+incr]->getReferencePositions(), args, mypath[k+incr]->getReferenceArguments(), false, mydir );
      // Set our direction equal to the displacement
      // mydir.setDirection( mypack );
      // Shift the reference configuration by this amount
      newpath[i].displaceReferenceConfiguration( dr, mydir );
    }

    // Copy the positions of the new path to the new paths
    for(int i=istart+incr; loopEnd(i,cfin,incr)==false; i+=incr) {
      mypdb.setAtomPositions( newpath[i].getReferencePositions() );
      for(unsigned j=0; j<newpath[i].getNumberOfReferenceArguments(); ++j) mypdb.setArgumentValue( mypath[i]->getArgumentNames()[j], newpath[i].getReferenceArgument(j) );
      mypath[i]->read( mypdb );
    }

    // Recompute the separations between frames
    calcCurrentPathSpacings( istart, iend );
  }
}

void PathReparameterization::reparameterize( const int& ifix1, const int& ifix2, const double& TOL ) {
  plumed_dbg_assert( ifix1<ifix2 );
  // First reparameterize the part between the fixed frames
  reparameterizePart( ifix1, ifix2, -1.0, TOL );

  // Get the separation between frames which we will use to set the remaining frames
  double target = sumlen[ifix2] / ( ifix2 - ifix1 );

  // And reparameterize the beginning and end of the path
  if( ifix1>0 ) reparameterizePart( ifix1, 0, target, TOL );
  if( ifix2<(mypath.size()-1) ) reparameterizePart( ifix2, mypath.size()-1, target, TOL );
//  calcCurrentPathSpacings( 0, mypath.size()-1 );
}

}
}
