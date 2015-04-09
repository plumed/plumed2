/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#include "AnalysisWithLandmarks.h"
#include "reference/ReferenceAtoms.h"
#include "reference/ReferenceArguments.h"
#include "reference/MultiReferenceBase.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"

//+PLUMEDOC ANALYSIS OUTPUT_LANDMARKS
/*
Select a set of the stored configurations using a landmark selection and output them to a file

For methods such as \ref CLASSICAL_MDS it can be expensive to run the analysis calculation 
with a large number of configurations.  What might be required is to run the analysis on a 
subset of frames.  One may then use the results from the performed analysis to do some further
analysis on the stored trajectory.  When running \ref CLASSICAL_MDS for example one may subsquently
project the remainder of the trajectory using some form of out of sample extension.  There are 
various ways of selecting the subset of frames on which to perform the analysis.  
These various methods are described below:

<table align=center frame=void width=95%% cellpadding=5%%>
<tr> 
<td> TYPE </td> <td> DESCRIPTION </td> <td> EXAMPLE INPUT </td>
</tr>
<tr>
<td> ALL </td> <td> use all the stored frames </td> <td> LANDMARKS={ALL} </td> </tr> <tr>
<td> STRIDE </td> <td> pick frames from the trajectory using a uniform stride  </td> <td> LANDMARKS={STRIDE N=\f$n\f$} </td> </tr> <tr>
<td> FPS </td> <td> pick \f$n\f$ frames using farthest point sampling </td> <td> LANMARKS={FPS N=\f$n\f$} </td> </tr> <tr>
</tr>
</table>

This action simply performs the required form of landmark selection and outputs the selected frames to a file.
The above landmark selection algorithms can be invoked in other analysis methods (e.g. \ref CLASSICAL_MDS or \ref SKETCHMAP).
When this is done the expensive part of the analysis will only be done on the selected landmark frames.

Weights are ascribed to each of the the points by doing a Voronoi analysis over all the fraems in the trajectory
unless this features is explicitally turned off using the keyword NOVORONOI.  As such a landmarks point with 
20 points in its Voronoi will be ascribed a weight of 20 unless you turn of this weighting.  In addition, if you are
running biased simulations and \ref rewweighting these weights will be taken into account when calculating weights in these
analysis algorithm.

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace analysis {

class OutputLandmarks : public AnalysisWithLandmarks {
private:
  std::string ofilename;
  std::string efilename;
  MultiReferenceBase* myembedding;
public:
  static void registerKeywords( Keywords& keys );
  OutputLandmarks( const ActionOptions& );
  ~OutputLandmarks();
  void analyzeLandmarks();
};

PLUMED_REGISTER_ACTION(OutputLandmarks,"OUTPUT_LANDMARKS")

void OutputLandmarks::registerKeywords( Keywords& keys ){
  AnalysisWithLandmarks::registerKeywords( keys );
  keys.add("compulsory","OUTPUT_FILE","dont output","this file will contain the selected landmarks in pdb format");
  keys.add("compulsory","LIST_FILE","dont output","this file contains a list of landmark points.  This can only be used when the distances are calculated from "
                                                  "the values of arguments.  That is to say when the metric does not involve atom positions");
}

OutputLandmarks::OutputLandmarks( const ActionOptions& ao ):
Action(ao),
AnalysisWithLandmarks(ao)
{
  myembedding = new MultiReferenceBase( getMetricName(), false );
  setDataToAnalyze( myembedding );

  parseOutputFile("OUTPUT_FILE",efilename);
  parseOutputFile("LIST_FILE",ofilename);
}

OutputLandmarks::~OutputLandmarks(){
  delete myembedding;
}

void OutputLandmarks::analyzeLandmarks(){

  if( ofilename!="dont output"){
      // Output the embedding as long lists of data
      OFile gfile; gfile.link(*this); 
      gfile.setBackupString("analysis");
      gfile.fmtField(getOutputFormat()+" ");
      gfile.open( ofilename.c_str() );

      // Can't print out all landmark data if we have reference atom positions
      ReferenceAtoms* myat=dynamic_cast<ReferenceAtoms*>( myembedding->getFrame(0) );
      plumed_assert( !myat );
      
      // Print embedding coordinates
      for(unsigned i=0;i<myembedding->getNumberOfReferenceFrames();++i){
          ReferenceArguments* myref=dynamic_cast<ReferenceArguments*>( myembedding->getFrame(i) );
          plumed_assert( myref );
          for(unsigned j=0;j<myref->getReferenceArguments().size();++j){
              gfile.printField( myref->getArgumentNames()[j], myref->getReferenceArgument(j) );
          }
          gfile.printField();
      }  
      gfile.close();
  }

  // Output the embedding in plumed format
  if( efilename!="dont output"){
     OFile afile; afile.link(*this); afile.setBackupString("analysis"); std::size_t psign=getOutputFormat().find("%");
     afile.open( efilename.c_str() ); std::string descr="REMARK WEIGHT=%-" + getOutputFormat().substr(psign+1);
     for(unsigned j=0;j<myembedding->getNumberOfReferenceFrames();++j){
         afile.printf("DESCRIPTION: landmark configuration from %s performed at time %f",getLabel().c_str(),getTime() );
         afile.printf(descr.c_str(),(myembedding->getFrame(j))->getWeight() ); 
         if( plumed.getAtoms().usingNaturalUnits() ) (myembedding->getFrame(j))->print( 1.0, afile, getOutputFormat() );
         else (myembedding->getFrame(j))->print( plumed.getAtoms().getUnits().getLength()/0.1, afile, getOutputFormat() );
     }
     afile.close();
  }
}

}
}
