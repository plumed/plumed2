/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2020 The plumed team
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
#include "AnalysisBase.h"
#include "tools/PDB.h"
#include "core/ActionRegister.h"
#include "reference/MetricRegister.h"
#include "reference/ReferenceConfiguration.h"

//+PLUMEDOC ANALYSIS EUCLIDEAN_DISSIMILARITIES
/*
Calculate the matrix of dissimilarities between a trajectory of atomic configurations.

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace analysis {

class EuclideanDissimilarityMatrix : public AnalysisBase {
private:
  PDB mypdb;
  std::string mtype;
  Matrix<double> dissimilarities;
public:
  static void registerKeywords( Keywords& keys );
  explicit EuclideanDissimilarityMatrix( const ActionOptions& ao );
/// Do the analysis
  void performAnalysis() override;
/// This ensures that classes that use this data know that dissimilarities were set
  bool dissimilaritiesWereSet() const override { return true; }
/// Get information on how to calculate dissimilarities
  std::string getDissimilarityInstruction() const override;
/// Get the squared dissimilarity between two reference configurations
  double getDissimilarity( const unsigned& i, const unsigned& j ) override;
/// This is just to deal with ActionWithVessel
  void performTask( const unsigned&, const unsigned&, MultiValue& ) const override { plumed_error(); }
};

PLUMED_REGISTER_ACTION(EuclideanDissimilarityMatrix,"EUCLIDEAN_DISSIMILARITIES")

void EuclideanDissimilarityMatrix::registerKeywords( Keywords& keys ) {
  AnalysisBase::registerKeywords( keys ); keys.use("ARG"); keys.reset_style("ARG","optional");
  keys.add("compulsory","METRIC","EUCLIDEAN","the method that you are going to use to measure the distances between points");
  keys.add("atoms","ATOMS","the list of atoms that you are going to use in the measure of distance that you are using");
}

EuclideanDissimilarityMatrix::EuclideanDissimilarityMatrix( const ActionOptions& ao ):
  Action(ao),
  AnalysisBase(ao)
{
  parse("METRIC",mtype); std::vector<AtomNumber> atoms;
  if( my_input_data->getNumberOfAtoms()>0 ) {
    parseAtomList("ATOMS",atoms);
    if( atoms.size()!=0 ) {
      mypdb.setAtomNumbers( atoms );
      for(unsigned i=0; i<atoms.size(); ++i) {
        bool found=false;
        for(unsigned j=0; j<my_input_data->getAtomIndexes().size(); ++j) {
          if( my_input_data->getAtomIndexes()[j]==atoms[i] ) { found=true; break; }
        }
        if( !found ) {
          std::string num; Tools::convert( atoms[i].serial(), num );
          error("atom number " + num + " is not stored in any action that has been input");
        }
      }
      mypdb.addBlockEnd( atoms.size() );
    } else if( getNumberOfArguments()==0 ) {
      mypdb.setAtomNumbers( my_input_data->getAtomIndexes() );
      mypdb.addBlockEnd( my_input_data->getAtomIndexes().size() );
      if( mtype=="EUCLIDEAN" ) mtype="OPTIMAL";
    }
  }
  log.printf("  measuring distances using %s metric \n",mtype.c_str() );
  if( my_input_data->getArgumentNames().size()>0 ) {
    if( getNumberOfArguments()==0 && atoms.size()==0 ) {
      std::vector<std::string> argnames( my_input_data->getArgumentNames() );
      mypdb.setArgumentNames( argnames ); requestArguments( my_input_data->getArgumentList() );
    } else {
      std::vector<Value*> myargs( getArguments() );
      std::vector<std::string> inargnames( my_input_data->getArgumentNames() );
      std::vector<std::string> argnames( myargs.size() );
      for(unsigned i=0; i<myargs.size(); ++i) {
        argnames[i]=myargs[i]->getName();
        bool found=false;
        for(unsigned j=0; j<inargnames.size(); ++j) {
          if( argnames[i]==inargnames[j] ) { found=true; break; }
        }
        if( !found ) error("input named " + my_input_data->getLabel() + " does not store/calculate quantity named " + argnames[i] );
      }
      mypdb.setArgumentNames( argnames ); requestArguments( myargs );
    }
  }
}

void EuclideanDissimilarityMatrix::performAnalysis() {
  // Resize dissimilarities matrix and set all elements to zero
  if( !usingLowMem() ) {
    dissimilarities.resize( getNumberOfDataPoints(), getNumberOfDataPoints() ); dissimilarities=0;
  }
}

std::string EuclideanDissimilarityMatrix::getDissimilarityInstruction() const {
  return "TYPE=" + mtype;
}

double EuclideanDissimilarityMatrix::getDissimilarity( const unsigned& iframe, const unsigned& jframe ) {
  plumed_dbg_assert( iframe<getNumberOfDataPoints() && jframe<getNumberOfDataPoints() );
  if( !usingLowMem() ) {
    if( dissimilarities(iframe,jframe)>0. ) { return dissimilarities(iframe,jframe); }
  }
  if( iframe!=jframe ) {
    double dd;
    getStoredData( iframe, true ).transferDataToPDB( mypdb );
    auto myref1=metricRegister().create<ReferenceConfiguration>(mtype, mypdb);
    getStoredData( jframe, true ).transferDataToPDB( mypdb );
    auto myref2=metricRegister().create<ReferenceConfiguration>(mtype, mypdb);
    if( !usingLowMem() ) dd=dissimilarities(iframe,jframe) = dissimilarities(jframe,iframe) = distance( getPbc(), getArguments(), myref1.get(), myref2.get(), true );
    else dd=distance( getPbc(), getArguments(), myref1.get(), myref2.get(), true );
    return dd;
  }
  return 0.0;
}

}
}
