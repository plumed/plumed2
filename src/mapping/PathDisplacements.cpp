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
#include "core/ActionWithValue.h"
#include "core/ActionPilot.h"
#include "core/ActionWithArguments.h"
#include "core/ActionRegister.h"
#include "PathProjectionCalculator.h"

//+PLUMEDOC COLVAR ADAPTIVE_PATH
/*
Compute path collective variables that adapt to the lowest free energy path connecting states A and B.

The Path Collective Variables developed by Branduardi and co-workers \cite brand07 allow one
to compute the progress along a high-dimensional path and the distance from the high-dimensional
path.  The progress along the path (s) is computed using:

\f[
s = i_2 + \textrm{sign}(i_2-i_1) \frac{ \sqrt{( \mathbf{v}_1\cdot\mathbf{v}_2 )^2 - |\mathbf{v}_3|^2(|\mathbf{v}_1|^2 - |\mathbf{v}_2|^2) } }{2|\mathbf{v}_3|^2} - \frac{\mathbf{v}_1\cdot\mathbf{v}_3 - |\mathbf{v}_3|^2}{2|\mathbf{v}_3|^2}
\f]

In this expression \f$\mathbf{v}_1\f$ and \f$\mathbf{v}_3\f$ are the vectors connecting the current position to the closest and second closest node of the path,
respectfully and \f$i_1\f$ and \f$i_2\f$ are the projections of the closest and second closest frames of the path. \f$\mathbf{v}_2\f$, meanwhile, is the
vector connecting the closest frame to the second closest frame.  The distance from the path, \f$z\f$ is calculated using:

\f[
z = \sqrt{ \left[ |\mathbf{v}_1|^2 - |\mathbf{v}_2| \left( \frac{ \sqrt{( \mathbf{v}_1\cdot\mathbf{v}_2 )^2 - |\mathbf{v}_3|^2(|\mathbf{v}_1|^2 - |\mathbf{v}_2|^2) } }{2|\mathbf{v}_3|^2} - \frac{\mathbf{v}_1\cdot\mathbf{v}_3 - |\mathbf{v}_3|^2}{2|\mathbf{v}_3|^2} \right) \right]^2 }
\f]

Notice that these are the definitions of \f$s\f$ and \f$z\f$ that are used by \ref PATH when the GPATH option is employed.  The reason for this is that
the adaptive path method implemented in this action was inspired by the work of Diaz and Ensing in which these formula were used \cite BerndAdaptivePath.
To learn more about how the path is adapted we strongly recommend reading this paper.

\par Examples

The input below provides an example of how the adaptive path works in practise. The path is updated every 50 steps of
MD based on the data accumulated during the preceding 50 time steps.

\plumedfile
d1: DISTANCE ATOMS=1,2 COMPONENTS
pp: ADAPTIVE_PATH TYPE=EUCLIDEAN FIXED=5,15 UPDATE=50 WFILE=out-path.pdb WSTRIDE=50 REFERENCE=mypath.pdb
PRINT ARG=d1.x,d1.y,pp.* FILE=colvar
\endplumedfile

In the case above the distance between frames is calculated based on the \f$x\f$ and \f$y\f$ components of the vector connecting
atoms 1 and 2.  As such an extract from the input reference path (mypath.pdb) would look as follows:

\verbatim
REMARK ARG=d1.x,d1.y d1.x=1.12 d1.y=-.60
END
REMARK ARG=d1.x,d1.y d1.x=.99 d1.y=-.45
END
\endverbatim

Notice that one can also use RMSD frames in place of arguments like those above.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace mapping {

class PathDisplacements : public ActionWithValue, public ActionPilot, public ActionWithArguments {
private:
  bool clearnextstep;
  unsigned clearstride;
  double fadefact;
  std::vector<double> wsum, displace_v;
  Matrix<double> displacements;
  PathProjectionCalculator path_projector;
public:
  static void registerKeywords( Keywords& keys );
  explicit PathDisplacements(const ActionOptions&);
  unsigned getNumberOfDerivatives() const ;
  void clearDerivatives( const bool& force=false ) {}
  void calculate(){}
  void apply(){}
  void update();
};

PLUMED_REGISTER_ACTION(PathDisplacements,"AVERAGE_PATH_DISPLACEMENT")

void PathDisplacements::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys ); ActionWithValue::registerKeywords( keys ); ActionPilot::registerKeywords( keys ); 
  ActionWithArguments::registerKeywords( keys ); keys.use("ARG"); PathProjectionCalculator::registerKeywords( keys );
  keys.add("compulsory","STRIDE","1","the frequency with which the average displacements should be collected and added to the average displacements");
  keys.add("compulsory","HALFLIFE","-1","the number of MD steps after which a previously measured path distance weighs only 50% in the average. This option may increase convergence by allowing to \"forget\" the memory of a bad initial guess path. The default is to set this to infinity");
  keys.add("compulsory","CLEAR","0","the frequency with which to clear all the accumulated data.  The default value "
                                    "of 0 implies that all the data will be used and that the grid will never be cleared");
}

PathDisplacements::PathDisplacements(const ActionOptions& ao):
  Action(ao),
  ActionWithValue(ao),
  ActionPilot(ao),
  ActionWithArguments(ao),
  clearnextstep(false),
  path_projector(this)
{
  plumed_assert( !actionInChain() );
  // Read in clear instructions
  parse("CLEAR",clearstride);
  if( clearstride>0 ) {
    if( clearstride%getStride()!=0 ) error("CLEAR parameter must be a multiple of STRIDE");
    log.printf("  clearing average every %u steps \n",clearstride);
  }
  if( arg_ends.size()>0 ) error("makes no sense to use ARG1, ARG2... with this action use single ARG keyword");
  double halflife; parse("HALFLIFE",halflife);
  log.printf("  weight of contribution to frame halves every %f steps \n",halflife);
  if( halflife<0 ) fadefact=1.0;
  else fadefact = exp( -0.693147180559945 / static_cast<double>(halflife) );
  // Now create the weights vector and displacements matrix 
  unsigned nrows = getPntrToArgument(0)->getShape()[0];
  unsigned ncols = getPntrToArgument(0)->getShape()[1];
  wsum.resize( nrows ); displacements.resize( nrows, ncols );
  for(unsigned i=0;i<nrows;++i) {
      wsum[i]=0; for(unsigned j=0;j<ncols;++j) displacements(i,j)=0;
  }
  // Add bibliography
  log<<"  Bibliography "<<plumed.cite("Diaz Leines and Ensing, Phys. Rev. Lett. 109, 020601 (2012)")<<"\n";
  // And create a value to hold the displacements
  std::vector<unsigned> shape(2); shape[0]=nrows; shape[1]=ncols;
  addValue( shape ); setNotPeriodic(); 
  getPntrToOutput(0)->alwaysStoreValues(); 
  getPntrToOutput(0)->setShape( shape ); 
  getPntrToOutput(0)->clearDerivatives();
}

unsigned PathDisplacements::getNumberOfDerivatives() const {
  return 0;
}

void PathDisplacements::update() {
  unsigned nrows = getPntrToArgument(0)->getShape()[0];
  unsigned ncols = getPntrToArgument(0)->getShape()[1];

  if( clearnextstep ) {
      unsigned k=0;
      for(unsigned i=0;i<nrows;++i) {
          for(unsigned j=0;j<ncols;++j) { displacements(i,j)=0; getPntrToOutput(0)->set(k,0); k++; }
      }
      clearnextstep=false;
  }

  unsigned k=0, iclose1, iclose2; double v1v1, v3v3;
  for(unsigned i=0;i<nrows;++i) {
      double dist = 0;
      for(unsigned j=0;j<ncols;++j) {
          double tmp = getPntrToArgument(0)->get(k);
          dist += tmp*tmp; k++;
      }
      if( i==0 ) { v1v1 = dist; iclose1 = 0; }
      else if( dist<v1v1 ) { v3v3=v1v1; v1v1=dist; iclose2=iclose1; iclose1=i; }
      else if( i==1 ) { v3v3=dist; iclose2=1; }
      else if( dist<v3v3 ) { v3v3=dist; iclose2=i; }
  }
  // And find third closest point
  int isign = iclose1 - iclose2;
  if( isign>1 ) isign=1; else if( isign<-1 ) isign=-1;
  int iclose3 = iclose1 + isign; 
  unsigned ifrom=iclose1, ito=iclose3; if( iclose3<0 || iclose3>=nrows ) { ifrom=iclose2; ito=iclose1; } 

  // Calculate the dot product of v1 with v2
  Tensor box( plumed.getAtoms().getPbc().getBox() );
  path_projector.getDisplaceVector( ifrom, ito, box, displace_v );
  double v2v2=0, v1v2=0; unsigned kclose1 = iclose1*ncols;
  for(unsigned i=0;i<displace_v.size();++i) { v2v2 += displace_v[i]*displace_v[i]; v1v2 += displace_v[i]*getPntrToArgument(0)->get(kclose1+i); }

  double root = sqrt( v1v2*v1v2 - v2v2 * ( v1v1 - v3v3) );
  double dx = 0.5 * ( (root + v1v2) / v2v2 - 1.);
  double weight2 = -1.* dx; double weight1 = 1.0 + dx;
  if( weight1>1.0 ) {
    weight1=1.0; weight2=0.0;
  } else if( weight2>1.0 ) {
    weight1=0.0; weight2=1.0;
  }

  // Accumulate displacements for path
  for(unsigned i=0;i<ncols;++i) {
      double displace = getPntrToArgument(0)->get(kclose1+i) - dx*displace_v[i];
      displacements(iclose1,i) += weight1 * displace; displacements(iclose2,i) += weight2 * displace; 
  } 

  // Update weight accumulators
  wsum[iclose1] *= fadefact; wsum[iclose2] *= fadefact;
  wsum[iclose1] += weight1; wsum[iclose2] += weight2;

  // Update numbers in values
  if( wsum[iclose1] > epsilon ) {
      for(unsigned i=0;i<ncols;++i) getPntrToOutput(0)->set( kclose1+i, displacements(iclose1,i) / wsum[iclose1] );
  } 
  if( wsum[iclose2] > epsilon ) { 
      unsigned kclose2 = iclose2*ncols;
      for(unsigned i=0;i<ncols;++i) getPntrToOutput(0)->set( kclose2+i, displacements(iclose2,i) / wsum[iclose2] );
  }

  // Clear if required
  if( (getStep()>0 && clearstride>0 && getStep()%clearstride==0) ) clearnextstep=true;
}

}
}
