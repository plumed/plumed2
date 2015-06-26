/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013,2014 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#include "tools/IFile.h"
#include "DissimilarityMatrixBase.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace analysis {

// PLUMED_REGISTER_ACTION(DissimilarityMatrixBase,"READ_DISSIMILARITIES")

void DissimilarityMatrixBase::registerKeywords( Keywords& keys ){
  Analysis::registerKeywords( keys );
  keys.add("compulsory","DATA","the action that stores the reference configurations which we are calculating the distance between");
//   keys.add("compulsory","FILE","an input file containing the matrix of dissimilarities");
//   keys.add("compulsory","NPOINTS","number of rows/columns in dissimilarity matrix");
//   keys.add("compulsory","WFILE","input file containing weights of points");
}

DissimilarityMatrixBase::DissimilarityMatrixBase( const ActionOptions& ao ):
Action(ao),
Analysis(ao)
{
//   if( keywords.exists("FILE") ) parse("FILE",fname);
//   if( fname.length()>0 ){
//       parse("NPOINTS",nnodes);
//       parse("WFILE",wfile);     // Don't quite know how to deal with this yet  
//   } else {
//       std::string aname; parse("DATA",aname);
//       myanalysis_obj = plumed.getActionSet().selectWithLabel<Analysis*>(aname);
//       log.printf("  calculating dissimilarity matrix based on data stored by action %s \n",aname.c_str() );
//       if(!myanalysis_obj) error("action labelled " +  aname + " does not exist or is not of Analysis type");
//       addDependency( myanalysis_obj );
//       setStride( myanalysis_obj->getRunFrequency() );
//   }
}

void DissimilarityMatrixBase::performAnalysis(){
  mydissimilarities.resize( getNumberOfDataPoints(), getNumberOfDataPoints() ); mydissimilarities=0.0;
}

// void DissimilarityMatrixBase::runFinalJobs(){
// abort();
//   if( myanalysis_obj ){
//       mydissimilarities.resize( myanalysis_obj->getNumberOfDataPoints(), myanalysis_obj->getNumberOfDataPoints() );
//       mydissimilarities=0.0; setupDissimilarityCalc();
//   } else {
//       plumed_assert( fname.length()>0 ); mydissimilarities.resize( nnodes, nnodes );
//       IFile mfile; mfile.open(fname); std::vector<std::string> words;
//       for(unsigned i=0;i<nnodes;++i){
//           Tools::getParsedLine( mfile, words );
//           if( words.size()!=nnodes ) error("bad formatting in matrix file");
//           for(unsigned j=0;j<nnodes;++j) Tools::convert( words[j], mydissimilarities(i,j) ); 
//       }
//       mfile.close();
//   }
// }

}
}
