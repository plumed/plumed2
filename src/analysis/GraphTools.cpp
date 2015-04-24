/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2014 The plumed team
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
#include "core/CLTool.h"
#include "core/PlumedMain.h"
#include "cltools/CLToolRegister.h"
#include "tools/Tools.h"
#include "tools/Matrix.h"
#include "tools/IFile.h"
#include "core/ActionRegister.h"
#include "Analysis.h"
#include <cstdio>
#include <string>
#include <vector>
#include <iostream>

using namespace std;

namespace PLMD {
namespace analysis{

//+PLUMEDOC TOOLS graphtools
/*
graphtools is a tool that allows you to analyse read in matrices of pairwise dissimilarities

\par Examples
 

*/
//+ENDPLUMEDOC

class GraphTools:
public CLTool
{
public:
  static void registerKeywords( Keywords& keys );
  GraphTools(const CLToolOptions& co );
  int main(FILE* in, FILE*out,Communicator& pc);
  string description()const{
    return "use plumed's analysis tools to analyse a matrix of pairwise disimilarities";
  }
};

PLUMED_REGISTER_CLTOOL(GraphTools,"graphtools")

void GraphTools::registerKeywords( Keywords& keys ){
  CLTool::registerKeywords( keys );
  keys.add("compulsory","--plumed","plumed.dat","the plumed input file");
  keys.add("compulsory","--nnodes","number of nodes in your graph");
  keys.add("compulsory","--matrix","the file containing the matrix of dissimilarities");
  keys.add("optional","--weights","the file containing the vector of weights.  If absent all weights are assumed equal to one.");
  keys.addFlag("--squared-dissimilarities",false,"are the disimilarities in your matrix file squared distances");
}

GraphTools::GraphTools(const CLToolOptions& co ):
CLTool(co)
{
  inputdata=commandline;
}

int GraphTools::main(FILE* in, FILE*out,Communicator& pc){

 unsigned nnodes; parse("--nnodes",nnodes);

 std::string plumedin; parse("--plumed",plumedin); 

 PlumedMain fakemain; fakemain.init( false );
 std::vector<Analysis*> analysis_set; PLMD::IFile ifile; 
 ifile.open(plumedin); std::vector<std::string> words; 
 while( Tools::getParsedLine(ifile,words) ){
    if( words.empty() ) continue;
    else if( words[0]=="ENDPLUMED" ) break;
    else {
        std::vector<std::string> interpreted(words);
        Tools::interpretLabel(interpreted);

        bool found_copy=false;
        for(unsigned i=0;i<interpreted.size();++i){
            if( interpreted[i]=="USE_ALL_DATA" ) error("do not set USE_ALL_DATA in input"); 
            std::size_t dot=interpreted[i].find_first_of("=");
            if( interpreted[i].substr(0,dot)=="RUN" ) error("do not set RUN in input");
            if( interpreted[i].substr(0,dot)=="STRIDE" ) error("do not set STRIDE in input");
            if( interpreted[i].substr(0,dot)=="REUSE_DATA_FROM" || interpreted[i].substr(0,dot)=="USE_DIMRED_DATA_FROM" ) found_copy=true;
        }
        if( analysis_set.size()==0 ){
            if( found_copy ) error("first action should use data from read in matrix");
            interpreted.push_back("USE_ALL_DATA");
        } else if( !found_copy ){
             error("later actions should refer to data in first action");  
        }

        Analysis* analysis=dynamic_cast<Analysis*>( actionRegister().create(ActionOptions(fakemain,interpreted)) );
        if( !analysis ){
            std::string errstring;
            errstring = "ERROR:\n  I cannot understand line:";
            for(unsigned i=0;i<interpreted.size();++i) errstring += " " + interpreted[i];
            errstring += "\n";
            error(errstring);
        }
        analysis->checkRead();
        analysis_set.push_back( analysis );
    }
 } 
 ifile.close(); fakemain.finish_init_print();

 std::string matin; parse("--matrix",matin); 
 bool sq; parseFlag("--squared-dissimilarities",sq);
 Matrix<double> mymatrix(nnodes, nnodes);
 IFile mfile; mfile.open(matin);
 for(unsigned i=0;i<nnodes;++i){
     Tools::getParsedLine( mfile, words );
     if( words.size()!=nnodes ) error("bad formatting in matrix file");
     for(unsigned j=0;j<nnodes;++j) Tools::convert( words[j], mymatrix(i,j) );
 }
 mfile.close();

 std::string win; std::vector<double> myweights(nnodes); parse("--weights",win);
 if( win.length()==0 ) myweights.assign( myweights.size(), 1.0 );

 analysis_set[0]->setPairwiseDisimilarityMatrix( mymatrix, myweights, sq ); 
 for(unsigned i=0;i<analysis_set.size();++i) analysis_set[i]->runAnalysis();
 for(unsigned i=0;i<analysis_set.size();++i) delete analysis_set[i]; 

 return 0;
}

} // End of namespace
}
