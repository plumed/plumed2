/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016,2017 The plumed team
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
#include "core/ActionWithArguments.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/Atoms.h"
#include "tools/Pbc.h"

namespace PLMD {
namespace mapping {

class GeometricPath : public ActionWithValue, public ActionWithArguments {
private:
  PlumedMain metric;
  std::vector<double> masses;
  std::vector<double> charges;
  std::vector<Vector> positions;
  std::vector<Vector> forces;
  std::vector<double> pcoords;
  std::vector<double> data;
  std::vector<std::vector<AtomNumber> > atind_for_frames;
  double getProjectionOnPath( const unsigned& ifrom, const unsigned& ito, const unsigned& closest, double& len ); 
public:
  static void registerKeywords(Keywords& keys);
  explicit GeometricPath(const ActionOptions&);
  void calculate();
  unsigned getNumberOfDerivatives() const ;
  void apply() {}
};

PLUMED_REGISTER_ACTION(GeometricPath,"GEOMETRIC_PATH")

void GeometricPath::registerKeywords(Keywords& keys) {
  Action::registerKeywords(keys); ActionWithValue::registerKeywords(keys); ActionWithArguments::registerKeywords(keys); keys.use("ARG");
  keys.add("compulsory","METRIC","the method to use for computing the displacement vectors between the reference frames");
  keys.add("compulsory","REFFRAMES","labels for actions that contain reference coordinates for each point on the path");
  keys.add("compulsory","COORDINATES","a vector of coordinates describing the position of each point along the path.  The default "
           "is to place these coordinates at 1, 2, 3, ...");
  componentsAreNotOptional(keys);
  keys.addOutputComponent("s","default","the position on the path");
  keys.addOutputComponent("z","default","the distance from the path");
}

GeometricPath::GeometricPath(const ActionOptions&ao):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao)
{
  // Ensure that values are stored in base calculation and that PLUMED doesn't try to calculate this in the stream
  plumed_assert( !actionInChain() ); getPntrToArgument(0)->buildDataStore( getLabel() );
  if( arg_ends.size()>0 ) error("makes no sense to use ARG1, ARG2... with this action use single ARG keyword");
  // Check that we have only one argument as input
  if( getNumberOfArguments()!=1 ) error("should only have one argument to this function");
  // Check that the input is a matrix
  if( getPntrToArgument(0)->getRank()!=2 ) error("the input to this action should be a matrix");
  // Check none of the arguments are periodic
  if( getPntrToArgument(0)->isPeriodic() ) error("cannot use this function on periodic functions");
  // Get the labels for the reference points
  std::vector<std::string> reflabs( getPntrToArgument(0)->getShape()[0] ); parseVector("REFFRAMES", reflabs );
  if( plumed.getAtoms().getAllGroups().count(reflabs[0]) ) {
      for(unsigned i=0;i<reflabs.size();++i) {
          if( !plumed.getAtoms().getAllGroups().count(reflabs[i]) ) error("count not find group with name " + reflabs[i] );
          atind_for_frames.push_back( plumed.getAtoms().getAllGroups().find(reflabs[i])->second );
          if( i>0 && atind_for_frames[i].size()!=atind_for_frames[0].size() ) error("mismatched numbers of atoms in reference frames");
      }
  }
  // Get the coordinates in the low dimensional space
  pcoords.resize( getPntrToArgument(0)->getShape()[0] ); parseVector("COORDINATES", pcoords );
  for(unsigned i=0;i<reflabs.size();++i) log.printf("  projecting frame read in by action %s at %f \n", reflabs[i].c_str(), pcoords[i] );
  // Create the values to store the output
  addComponentWithDerivatives("s"); componentIsNotPeriodic("s");
  addComponentWithDerivatives("z"); componentIsNotPeriodic("z");
  // Create a plumed main object to compute distances between reference configurations
  int s=sizeof(double);
  metric.cmd("setRealPrecision",&s);
  metric.cmd("setNoVirial"); 
  metric.cmd("setMDEngine","plumed");
  int natoms=0; if( atind_for_frames.size()>0 ) natoms = 2*atind_for_frames[0].size();
  metric.cmd("setNatoms",&natoms);
  positions.resize(natoms); masses.resize(natoms); forces.resize(natoms); charges.resize(natoms);
  double tstep=1.0; metric.cmd("setTimestep",&tstep);
  std::string inp; parse("METRIC",inp); const char* cinp=inp.c_str();
  std::vector<std::string> input=Tools::getWords(inp);
  if( input.size()==1 && !actionRegister().check(input[0]) ) {
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

unsigned GeometricPath::getNumberOfDerivatives() const {
  return 0;
}

double GeometricPath::getProjectionOnPath( const unsigned& ifrom, const unsigned& ito, const unsigned& closest, double& len ) {
  if( atind_for_frames.size()>0 ) {
      unsigned nref_at = atind_for_frames[0].size();
      Atoms& atoms( plumed.getAtoms() ); int istep=0; metric.cmd("setStep",&istep);
      // Set the atomic positions of the coordinates
      for(unsigned i=0;i<nref_at;++i) {
          // masses[i] = m[atind_for_frames[ifrom][i].index()]; masses[nref_at + i] = m[atind_for_frames[ito][i].index()]; 
          // charges[i] = c[atind_for_frames[ifrom][i].index()]; charges[nref_at + i] = c[atind_for_frames[ito][i].index()]; 
          positions[i] = atoms.getVatomPosition(atind_for_frames[ifrom][i]); positions[nref_at + i] = atoms.getVatomPosition(atind_for_frames[ito][i]);
      }
      metric.cmd("setMasses",&masses[0]);
      metric.cmd("setCharges",&charges[0]);
      metric.cmd("setPositions",&positions[0]);
      metric.cmd("setForces",&forces[0]);
      Tensor box( plumed.getAtoms().getPbc().getBox() ); metric.cmd("setBox",&box[0][0]);
  } else {
      plumed_error();
  }
  metric.cmd("calc");  
  double fval=0; len=0; Value* arg=getPntrToArgument(0); unsigned k=arg->getShape()[1]*closest;
  for(unsigned i=0;i<data.size();++i) { len += data[i]*data[i]; fval += data[i]*arg->get(k+i); }
  return fval;
}

void GeometricPath::calculate() {
  unsigned k=0, iclose1, iclose2; double v1v1, v3v3;
  unsigned nrows = getPntrToArgument(0)->getShape()[0];
  unsigned ncols = getPntrToArgument(0)->getShape()[1];
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

  // And calculate projection of vector connecting current point to closest frame on vector connecting nearest two frames
  double v2v2, v1v2 = getProjectionOnPath( ifrom, ito, iclose1, v2v2 ); 

  // This computes s value
  double spacing = pcoords[iclose1] - pcoords[iclose2];
  double root = sqrt( v1v2*v1v2 - v2v2 * ( v1v1 - v3v3) );
  double dx = 0.5 * ( (root + v1v2) / v2v2 - 1.);
  double path_s = pcoords[iclose1] + spacing * dx;
  double fact = 0.25*spacing / v2v2; Value* sp = getPntrToComponent(0); sp->set( path_s );

  // This computes z value
  double v4v4, proj = getProjectionOnPath( iclose2, iclose1, iclose1, v4v4 ); 
  double path_z = v1v1 + dx*dx*v4v4 - 2*dx*proj; path_z = sqrt(path_z);
  Value* zp = getPntrToComponent(1); zp->set( path_z );
}

}
}
