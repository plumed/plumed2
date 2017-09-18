/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2017 The plumed team
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
#include "SymmetryFunctionBase.h"
#include "multicolvar/MultiColvarBase.h"
#include "core/ActionRegister.h"
#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace symfunc {

//+PLUMEDOC MCOLVARF LOCAL_AVERAGE
/*
Calculate averages over spherical regions centered on atoms

As is explained in <a href="http://www.youtube.com/watch?v=iDvZmbWE5ps"> this video </a> certain multicolvars
calculate one scalar quantity or one vector for each of the atoms in the system.  For example
\ref COORDINATIONNUMBER measures the coordination number of each of the atoms in the system and \ref Q4 measures
the 4th order Steinhardt parameter for each of the atoms in the system.  These quantities provide tell us something about
the disposition of the atoms in the first coordination sphere of each of the atoms of interest.  Lechner and Dellago \cite dellago-q6
have suggested that one can probe local order in a system by taking the average value of such symmetry functions over
the atoms within a spherical cutoff of each of these atoms in the systems.  When this is done with Steinhardt parameters
they claim this gives a coordinate that is better able to distinguish solid and liquid configurations of Lennard-Jones atoms.

You can calculate such locally averaged quantities within plumed by using the LOCAL_AVERAGE command.  This command calculates
the following atom-centered quantities:

\f[
s_i = \frac{ c_i + \sum_j \sigma(r_{ij})c_j }{ 1 + \sum_j \sigma(r_{ij}) }
\f]

where the \f$c_i\f$ and \f$c_j\f$ values can be for any one of the symmetry functions that can be calculated using plumed
multicolvars.  The function \f$\sigma( r_{ij} )\f$ is a \ref switchingfunction that acts on the distance between
atoms \f$i\f$ and \f$j\f$.  Lechner and Dellago suggest that the parameters of this function should be set so that it the function is equal to one
when atom \f$j\f$ is in the first coordination sphere of atom \f$i\f$ and is zero otherwise.

The \f$s_i\f$ quantities calculated using the above command can be again thought of as atom-centred symmetry functions.  They
thus operate much like multicolvars.  You can thus calculate properties of the distribution of \f$s_i\f$ values using MEAN, LESS_THAN, HISTOGRAM
and so on.  You can also probe the value of these averaged variables in regions of the box by using the command in tandem with the
\ref AROUND command.

\par Examples

This example input calculates the coordination numbers for all the atoms in the system.  These coordination numbers are then averaged over
spherical regions.  The number of averaged coordination numbers that are greater than 4 is then output to a file.

\plumedfile
COORDINATIONNUMBER SPECIES=1-64 D_0=1.3 R_0=0.2 LABEL=d1
LOCAL_AVERAGE ARG=d1 SWITCH={RATIONAL D_0=1.3 R_0=0.2} MORE_THAN={RATIONAL R_0=4} LABEL=la
PRINT ARG=la.* FILE=colvar
\endplumedfile

This example input calculates the \f$q_4\f$ (see \ref Q4) vectors for each of the atoms in the system.  These vectors are then averaged
component by component over a spherical region.  The average value for this quantity is then outputeed to a file.  This calculates the
quantities that were used in the paper by Lechner and Dellago \cite dellago-q6

\plumedfile
Q4 SPECIES=1-64 SWITCH={RATIONAL D_0=1.3 R_0=0.2} LABEL=q4
LOCAL_AVERAGE ARG=q4 SWITCH={RATIONAL D_0=1.3 R_0=0.2} MEAN LABEL=la
PRINT ARG=la.* FILE=colvar
\endplumedfile

*/
//+ENDPLUMEDOC

class MatrixTimesVector : public SymmetryFunctionBase {
public:
  static void shortcutKeywords( Keywords& keys );
  static void expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                              const std::map<std::string,std::string>& keys,
                              std::vector<std::vector<std::string> >& actions );
  static void registerKeywords( Keywords& keys );
  explicit MatrixTimesVector(const ActionOptions&);
  void compute( const double& weight, const Vector& vec, MultiValue& myvals ) const ;
};

PLUMED_REGISTER_ACTION(MatrixTimesVector,"MATRIX_VECTOR_PRODUCT")
PLUMED_REGISTER_SHORTCUT(MatrixTimesVector,"LOCAL_AVERAGE_Q6")

void MatrixTimesVector::shortcutKeywords( Keywords& keys ) {
  SymmetryFunctionBase::shortcutKeywords( keys );
} 

void MatrixTimesVector::expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                                        const std::map<std::string,std::string>& keys,
                                        std::vector<std::vector<std::string> >& actions ) {
  SymmetryFunctionBase::expandMatrix( false, lab, words, keys, actions );
  std::vector<std::string> coord; coord.push_back( lab + "_coord:"); coord.push_back("COORDINATIONNUMBER");
  coord.push_back("WEIGHT=" + lab + "_mat.w"); actions.push_back( coord );
  if( words[0].find("LOCAL_AVERAGE_Q")!=std::string::npos ) {
      int l; Tools::convert( words[0].substr(15), l );
      for(int i=-l;i<=l;++i) {
          std::string num; Tools::convert( i, num );
          std::vector<std::string> realdata; realdata.push_back( lab + "_prod-rmn-[" + num + "]:");
          realdata.push_back("MATRIX_VECTOR_PRODUCT"); realdata.push_back("WEIGHT=" + lab + "_mat.w");
          realdata.push_back("VECTOR=" + keys.find("SPECIES")->second + "_rmn-[" + num + "]"); 
          actions.push_back( realdata );
          std::vector<std::string> realmath; realmath.push_back( lab + "_av-rmn-[" + num + "]:");
          realmath.push_back("MATHEVAL"); realmath.push_back("ARG1=" + lab + "_prod-rmn-[" + num + "]");
          realmath.push_back("ARG2=" + keys.find("SPECIES")->second + "_rmn-[" + num + "]");
          realmath.push_back("ARG3=" + lab + "_coord"); realmath.push_back("FUNC=(x+y)/(1+z)"); 
          realmath.push_back("PERIODIC=NO"); actions.push_back( realmath );
          std::vector<std::string> imagdata; imagdata.push_back( lab + "_prod-imn-[" + num + "]:");
          imagdata.push_back("MATRIX_VECTOR_PRODUCT"); imagdata.push_back("WEIGHT=" + lab + "_mat.w");
          imagdata.push_back("VECTOR=" + keys.find("SPECIES")->second + "_imn-[" + num + "]");
          actions.push_back( imagdata );
          std::vector<std::string> imagmath; imagmath.push_back( lab + "_av-imn-[" + num + "]:");
          imagmath.push_back("MATHEVAL"); imagmath.push_back("ARG1=" + lab + "_prod-imn-[" + num + "]");
          imagmath.push_back("ARG2=" + keys.find("SPECIES")->second + "_imn-[" + num + "]");
          imagmath.push_back("ARG3=" + lab + "_coord"); imagmath.push_back("FUNC=(x+y)/(1+z)"); 
          imagmath.push_back("PERIODIC=NO"); actions.push_back( imagmath );
      }
      std::vector<std::string> comb; comb.push_back( lab + "_2:" ); comb.push_back("COMBINE");
      unsigned k=0; std::string powstr="POWERS=2";
      for(int i=-l;i<=l;++i) {
          std::string num, anum; Tools::convert(i,num);
          k++; Tools::convert(k,anum); comb.push_back("ARG" + anum + "=" + lab + "_av-rmn-[" + num + "]" );
          k++; Tools::convert(k,anum); comb.push_back("ARG" + anum + "=" + lab + "_av-imn-[" + num + "]" );
          if( i==-l ) powstr += ",2"; else powstr += ",2,2";
      }
      comb.push_back(powstr); comb.push_back("PERIODIC=NO"); actions.push_back( comb );
      std::vector<std::string> mcomb; mcomb.push_back( lab + ":"); mcomb.push_back("MATHEVAL");
      mcomb.push_back("ARG1=" + lab + "_2"); mcomb.push_back("FUNC=sqrt(x)"); mcomb.push_back("PERIODIC=NO");
      actions.push_back( mcomb );  
  }
  multicolvar::MultiColvarBase::expandFunctions( lab, lab, "", words, keys, actions );
}

void MatrixTimesVector::registerKeywords( Keywords& keys ) {
  SymmetryFunctionBase::registerKeywords( keys ); keys.remove("VECTORS");
  keys.add("compulsory","VECTOR","the vector that you would like to multiply by the input matrix");
}

MatrixTimesVector::MatrixTimesVector(const ActionOptions&ao):
  Action(ao),
  SymmetryFunctionBase(ao)
{
  std::vector<Value*> vecs; parseArgumentList("VECTOR",vecs);
  if( vecs.size()!=1 ) error("keyword VECTOR shoudl only be provided with the label of a singl action");
  if( vecs[0]->getShape()[0]!=getPntrToArgument(0)->getShape()[1] ) error("shape of input vector should match second dimension of input WEIGHT matrix");
  vecs[0]->buildDataStore(); log.printf("  calculating product of input weight matrix with vector of weights labelled %s \n",vecs[0]->getName().c_str() );
  std::vector<Value*> args( getArguments() ); args.push_back( vecs[0] ); 
  requestArguments( args, true ); addValueWithDerivatives();
}

void MatrixTimesVector::compute( const double& val, const Vector& dir, MultiValue& myvals ) const {
  unsigned tindex = myvals.getSecondTaskIndex(); if( tindex>=getFullNumberOfTasks() ) tindex -= getFullNumberOfTasks();
  double func = getPntrToArgument(1)->get(tindex); addToValue( 0, val*func, myvals ); 
  if( doNotCalculateDerivatives() ) return;  
  addWeightDerivative( 0, func, myvals ); 
  myvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), arg_deriv_starts[1] + tindex, func );
}

}
}
