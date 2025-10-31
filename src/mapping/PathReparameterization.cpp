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
#include "core/ActionWithValue.h"
#include "core/ActionPilot.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "tools/Pbc.h"
#include "tools/Matrix.h"
#include "PathProjectionCalculator.h"

//+PLUMEDOC ANALYSIS REPARAMETERIZE_PATH
/*
Take an input path with frames that are not equally spaced and make the frames equally spaced

This action is used by [ADAPTIVE_PATH](ADAPTIVE_PATH.md) and [pathtools](pathtools.md). The algorithm in this action takes
a set of trajectory frames along a [PATH](PATH.md) in input and adjusts them until they are all are equally spaced.  To call
this algorithm from [pathtools](pathtools.md) you would use an input like this one:

```plumed
plumed pathtools --path in_path.pdb --metric EUCLIDEAN --out final_path.pdb
```

If you are using this action directly and not through [pathtools](pathtools.md) or [ADAPTIVE_PATH](ADAPTIVE_PATH.md) you will
probably be using it in an input something like this if your path is defined using atomic coordinates:

```plumed
#SETTINGS INPUTFILES=regtest/trajectories/path_msd/all.pdb
rmsd: RMSD REFERENCE=regtest/trajectories/path_msd/all.pdb DISPLACEMENT TYPE=OPTIMAL
# Accumulate the average displacement between the reference path and the trajectories that have sampled the transition
disp: AVERAGE_PATH_DISPLACEMENT ...
  ARG=rmsd.disp STRIDE=1
  METRIC={RMSD DISPLACEMENT TYPE=OPTIMAL ALIGN=1,1,1,1,1,1,1,1,1,1,1,1,1 DISPLACE=1,1,1,1,1,1,1,1,1,1,1,1,1}
  METRIC_COMPONENT=disp REFERENCE=rmsd_ref
...
# Now displace the original path by the accumulated displacement and reparameterize so that all frames are equally spaced
REPARAMETERIZE_PATH ...
  DISPLACE_FRAMES=disp FIXED=1,42
  METRIC={RMSD DISPLACEMENT TYPE=OPTIMAL ALIGN=1,1,1,1,1,1,1,1,1,1,1,1,1 DISPLACE=1,1,1,1,1,1,1,1,1,1,1,1,1}
  METRIC_COMPONENT=disp REFERENCE=rmsd_ref
  MAXCYLES=100 TOL=1E-4 STRIDE=0
...
# And output the final reparameterized path at the end of the simulation
DUMPPDB DESCRIPTION=PATH STRIDE=0 FILE=outpatb.pdb ATOMS=rmsd_ref ATOM_INDICES=1-13
```

The method for reparameterizing paths that is implemented in this input is discussed in the example documentation for
[AVERAGE_PATH_DISPLACEMENT](AVERAGE_PATH_DISPLACEMENT.md).

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace mapping {

class PathReparameterization : public ActionPilot {
private:
/// Number of cycles of the optimization algorithm to run
  unsigned maxcycles;
/// The points on the path to fix
  unsigned ifix1, ifix2;
/// Tolerance for the minimization algorithm
  double TOL;
/// Value containing ammount to displace each reference configuration by
  Value* displace_value;
/// The action for calculating the distances between the frames
  PathProjectionCalculator path_projector;
/// Used to store current spacing between frames in path
  std::vector<double> data, len, sumlen, sfrac;
///
  bool loopEnd( int index, int end, int inc ) const;
///
  double computeSpacing( unsigned ifrom, unsigned ito );
///
  void calcCurrentPathSpacings( int istart, int iend );
///
  void reparameterizePart( int istart, int iend, double target );
public:
  static void registerKeywords( Keywords& keys );
  PathReparameterization(const ActionOptions&);
  void calculate() {}
  void apply() {}
  void update();
};

PLUMED_REGISTER_ACTION(PathReparameterization,"REPARAMETERIZE_PATH")

void PathReparameterization::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  PathProjectionCalculator::registerKeywords(keys);
  keys.remove("ARG");
  keys.add("compulsory","STRIDE","1","the frequency with which to reparameterize the path");
  keys.add("compulsory","FIXED","0","the frames in the path to fix");
  keys.add("compulsory","MAXCYLES","100","number of cycles of the algorithm to run");
  keys.add("compulsory","TOL","1E-4","the tolerance to use for the path reparameterization algorithm");
  keys.add("optional","DISPLACE_FRAMES","label of an action that tells us how to displace the frames.  These displacements are applied before "
           "running the reparameterization algorith");
}

PathReparameterization::PathReparameterization(const ActionOptions&ao):
  Action(ao),
  ActionPilot(ao),
  displace_value(NULL),
  path_projector(this) {
  parse("MAXCYLES",maxcycles);
  parse("TOL",TOL);
  log.printf("  running till change is less than %f or until there have been %d optimization cycles \n", TOL, maxcycles);
  len.resize( path_projector.getNumberOfFrames()  );
  sumlen.resize( path_projector.getNumberOfFrames() );
  sfrac.resize( path_projector.getNumberOfFrames() );
  std::vector<unsigned> fixed;
  parseVector("FIXED",fixed);
  if( fixed.size()==1 ) {
    if( fixed[0]!=0 ) {
      error("input to FIXED should be two integers");
    }
    ifix1=0;
    ifix2=path_projector.getNumberOfFrames()-1;
  } else if( fixed.size()==2 ) {
    if( fixed[0]<1 || fixed[1]<1 || fixed[0]>path_projector.getNumberOfFrames() || fixed[1]>path_projector.getNumberOfFrames() ) {
      error("input to FIXED should be two numbers between 1 and the number of frames");
    }
    ifix1=fixed[0]-1;
    ifix2=fixed[1]-1;
  } else {
    error("input to FIXED should be two integers");
  }
  log.printf("  fixing frames %d and %d when reparameterizing \n", ifix1, ifix2 );
  std::string dframe;
  parse("DISPLACE_FRAMES",dframe);
  if( dframe.length()>0 ) {
    ActionWithValue* av = plumed.getActionSet().selectWithLabel<ActionWithValue*>( dframe );
    if( !av ) {
      error("could not find action with label " + dframe + " specified to DISPLACE_FRAMES keyword in input file");
    }
    if( av->getName()!="AVERAGE_PATH_DISPLACEMENT" ) {
      error("displace object is not of correct type");
    }
    displace_value = av->copyOutput(0);
  }
}

bool PathReparameterization::loopEnd( const int index, const int end, const int inc ) const {
  if( inc>0 && index<end ) {
    return false;
  } else if( inc<0 && index>end ) {
    return false;
  }
  return true;
}

double PathReparameterization::computeSpacing( const unsigned ifrom, const unsigned ito ) {
  path_projector.getDisplaceVector( ifrom, ito, data );
  double length=0;
  for(unsigned i=0; i<data.size(); ++i) {
    length += data[i]*data[i];
  }
  return sqrt( length );
}

void PathReparameterization::calcCurrentPathSpacings( const int istart, const int iend ) {
  plumed_dbg_assert( static_cast<unsigned>(istart)<len.size() && static_cast<unsigned>(iend)<len.size() );
  len[istart] = sumlen[istart]=0;
  //printf("HELLO PATH SPACINGS ARE CURRENTLY \n");

  // Get the spacings given we can go forward and backwards
  int incr=1;
  if( istart>iend ) {
    incr=-1;
  }

  for(int i=istart+incr; loopEnd(i,iend+incr,incr)==false; i+=incr) {
    len[i] = computeSpacing( i-incr, i );
    sumlen[i] = sumlen[i-incr] + len[i];
    //printf("FRAME %d TO FRAME %d EQUALS %f : %f \n",i-incr,i,len[i],sumlen[i] );
  }
}

void PathReparameterization::reparameterizePart( const int istart, const int iend, const double target ) {
  calcCurrentPathSpacings( istart, iend );
  int cfin;
  // If a target separation is set we fix where we want the nodes
  int incr=1;
  if( istart>iend ) {
    incr=-1;
  }

  if( target>0 ) {
    if( iend>istart ) {
      for(int i=istart; i<iend+1; ++i) {
        sfrac[i] = target*(i-istart);
      }
    } else {
      for(int i=istart-1; i>iend-1; --i) {
        sfrac[i]=target*(istart-i);
      }
    }
    cfin = iend+incr;
  } else {
    cfin = iend;
  }

  double prevsum=0.;
  Matrix<double> newmatrix( path_projector.getNumberOfFrames(), data.size() );
  for(unsigned iter=0; iter<maxcycles; ++iter) {
    if( fabs(sumlen[iend] - prevsum)<=TOL ) {
      break ;
    }
    prevsum = sumlen[iend];
    // If no target is set we redistribute length
    if( target<0 ) {
      plumed_assert( istart<iend );
      double dr = sumlen[iend] / static_cast<double>( iend - istart );
      for(int i=istart; i<iend; ++i) {
        sfrac[i] = dr*(i-istart);
      }
    }

    // Now compute positions of new nodes in path
    for(int i=istart+incr; loopEnd(i,cfin,incr)==false; i+=incr) {
      int k = istart;
      while( !((sumlen[k] < sfrac[i]) && (sumlen[k+incr]>=sfrac[i])) ) {
        k+=incr;
        if( cfin==iend && k>= iend+1 ) {
          plumed_merror("path reparameterization error");
        } else if( cfin==(iend+1) && k>=iend ) {
          k=iend-1;
          break;
        } else if( cfin==(iend-1) && k<=iend ) {
          k=iend+1;
          break;
        }
      }
      double dr = (sfrac[i]-sumlen[k])/len[k+incr];
      // Copy the reference configuration to the row of a matrix
      path_projector.getReferenceConfiguration( k, data );
      for(unsigned j=0; j<data.size(); ++j) {
        newmatrix(i,j) = data[j];
      }
      path_projector.getDisplaceVector( k, k+incr, data );
      // Shift the reference configuration by this ammount
      for(unsigned j=0; j<data.size(); ++j) {
        newmatrix(i,j) += dr*data[j];
      }
    }

    // Copy the positions of the new path to the new paths
    for(int i=istart+incr; loopEnd(i,cfin,incr)==false; i+=incr) {
      for(unsigned j=0; j<data.size(); ++j) {
        data[j] = newmatrix(i,j);
      }
      path_projector.setReferenceConfiguration( i, data );
    }

    // Recompute the separations between frames
    calcCurrentPathSpacings( istart, iend );
  }
}

void PathReparameterization::update() {
  // We never run this on the first step
  if( getStep()==0 ) {
    return ;
  }

  // Shift the frames using the displacements
  if( displace_value ) {
    for(unsigned i=0; i<path_projector.getNumberOfFrames(); ++i) {
      if( i==ifix1 || i==ifix2 ) {
        continue ;
      }
      // Retrieve the current position of the frame
      path_projector.getReferenceConfiguration( i, data );
      // Shift using the averages accumulated in the action that accumulates the displacements
      unsigned kstart = i*data.size();
      for(unsigned j=0; j<data.size(); ++j) {
        data[j] += displace_value->get( kstart + j );
      }
      // And now set the new position of the refernce frame
      path_projector.setReferenceConfiguration( i, data );
    }
  }

  // First reparameterize the part between the fixed frames
  reparameterizePart( ifix1, ifix2, -1.0 );

  // Get the separation between frames which we will use to set the remaining frames
  double target = sumlen[ifix2] / ( ifix2 - ifix1 );

  // And reparameterize the begining and end of the path
  if( ifix1>0 ) {
    reparameterizePart( ifix1, 0, target );
  }
  if( ifix2<(path_projector.getNumberOfFrames()-1) ) {
    reparameterizePart( ifix2, path_projector.getNumberOfFrames()-1, target );
  }
  // And update any RMSD objects that depend on input values
  path_projector.updateDepedentRMSDObjects();
}

}
}
