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
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionShortcut.h"
#include "core/ActionWithValue.h"

//+PLUMEDOC FUNCTION MAHALANOBIS_DISTANCE
/*
Calculate the mahalanobis distance between two points in CV space

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace refdist {

class MahalanobisDistance : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit MahalanobisDistance(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(MahalanobisDistance,"MAHALANOBIS_DISTANCE")

void MahalanobisDistance::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords(keys);
  keys.add("compulsory","ARG1","The point that we are calculating the distance from");
  keys.add("compulsory","ARG2","The point that we are calculating the distance to");
  keys.add("compulsory","METRIC","The inverse covariance matrix that should be used when calculating the distance");
  keys.addFlag("SQUARED",false,"The squared distance should be calculated");
  keys.addFlag("VON_MISSES",false,"Compute the mahalanobis distance in a way that is more sympathetic to the periodic boundary conditions");
  keys.setValueDescription("the Mahalanobis distances between the input vectors");
  keys.needsAction("DISPLACEMENT"); keys.needsAction("CUSTOM"); keys.needsAction("OUTER_PRODUCT");
  keys.needsAction("TRANSPOSE"); keys.needsAction("MATRIX_PRODUCT_DIAGONAL"); keys.needsAction("CONSTANT");
  keys.needsAction("MATRIX_VECTOR_PRODUCT"); keys.needsAction("MATRIX_PRODUCT"); keys.needsAction("COMBINE");
}

MahalanobisDistance::MahalanobisDistance( const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao)
{
  std::string arg1, arg2, metstr; parse("ARG1",arg1); parse("ARG2",arg2); parse("METRIC",metstr);
  // Check on input metric
  ActionWithValue* mav=plumed.getActionSet().selectWithLabel<ActionWithValue*>( metstr );
  if( !mav ) error("could not find action named " + metstr + " to use for metric");
  if( mav->copyOutput(0)->getRank()!=2 ) error("metric has incorrect rank");

  readInputLine( getShortcutLabel() + "_diff: DISPLACEMENT ARG1=" + arg1 + " ARG2=" + arg2 );
  readInputLine( getShortcutLabel() + "_diffT: TRANSPOSE ARG=" + getShortcutLabel() + "_diff");
  bool von_miss, squared; parseFlag("VON_MISSES",von_miss); parseFlag("SQUARED",squared);
  if( von_miss ) {
    unsigned nrows = mav->copyOutput(0)->getShape()[0];
    if( mav->copyOutput(0)->getShape()[1]!=nrows ) error("metric is not symmetric");
    // Create a matrix that can be used to compute the off diagonal elements
    std::string valstr, nrstr;
    Tools::convert( mav->copyOutput(0)->get(0), valstr ); Tools::convert( nrows, nrstr );
    std::string diagmet = getShortcutLabel() + "_diagmet: CONSTANT VALUES=" + valstr;
    std::string offdiagmet = getShortcutLabel() + "_offdiagmet: CONSTANT NROWS=" + nrstr + " NCOLS=" + nrstr + " VALUES=0";
    for(unsigned i=0; i<nrows; ++i) {
      for(unsigned j=0; j<nrows; ++j) {
        Tools::convert( mav->copyOutput(0)->get(i*nrows+j), valstr );
        if( i==j && i>0 ) { offdiagmet += ",0"; diagmet += "," + valstr; }
        else if( i!=j ) { offdiagmet += "," + valstr; }
      }
    }
    readInputLine( diagmet ); readInputLine( offdiagmet );
    // Compute distances scaled by periods
    ActionWithValue* av=plumed.getActionSet().selectWithLabel<ActionWithValue*>( getShortcutLabel() + "_diff" ); plumed_assert( av );
    if( !av->copyOutput(0)->isPeriodic() ) error("VON_MISSES only works with periodic variables");
    std::string min, max; av->copyOutput(0)->getDomain(min,max);
    readInputLine( getShortcutLabel() + "_scaled: CUSTOM ARG=" + getShortcutLabel() + "_diffT FUNC=2*pi*x/(" + max +"-" + min + ") PERIODIC=NO");
    // We start calculating off-diagonal elements by computing the sines of the scaled differences (this is a column vector)
    readInputLine( getShortcutLabel() + "_sinediffT: CUSTOM ARG=" + getShortcutLabel() + "_scaled FUNC=sin(x) PERIODIC=NO");
    // Transpose sines to get a row vector
    readInputLine( getShortcutLabel() + "_sinediff: TRANSPOSE ARG=" + getShortcutLabel() + "_sinediffT");
    // Compute the off diagonal elements
    readInputLine( getShortcutLabel() + "_matvec: MATRIX_PRODUCT ARG=" + getShortcutLabel() + "_offdiagmet," + getShortcutLabel() +"_sinediffT");
    readInputLine( getShortcutLabel() + "_offdiag: MATRIX_PRODUCT_DIAGONAL ARG=" + getShortcutLabel() + "_sinediff," + getShortcutLabel() +"_matvec");
    // Sort out the metric for the diagonal elements
    std::string metstr2 = getShortcutLabel() + "_diagmet";
    // If this is a matrix we need create a matrix to multiply by
    if( av->copyOutput(0)->getShape()[0]>1 ) {
      // Create some ones
      std::string ones=" VALUES=1"; for(unsigned i=1; i<av->copyOutput(0)->getShape()[0]; ++i ) ones += ",1";
      readInputLine( getShortcutLabel() + "_ones: CONSTANT " + ones );
      // Now do some multiplication to create a matrix that can be multiplied by our "inverse variance" vector
      readInputLine( getShortcutLabel() + "_" + metstr + ": OUTER_PRODUCT ARG=" + metstr2 + "," + getShortcutLabel() + "_ones");
      metstr2 = getShortcutLabel() + "_" + metstr;
    }
    // Compute the diagonal elements
    readInputLine( getShortcutLabel() + "_prod: CUSTOM ARG=" + getShortcutLabel() + "_scaled," + metstr2 + " FUNC=2*(1-cos(x))*y PERIODIC=NO");
    std::string ncstr; Tools::convert( nrows, ncstr ); Tools::convert( av->copyOutput(0)->getShape()[0], nrstr );
    std::string ones=" VALUES=1"; for(unsigned i=1; i<av->copyOutput(0)->getNumberOfValues(); ++i) ones += ",1";
    readInputLine( getShortcutLabel() + "_matones: CONSTANT NROWS=" + nrstr + " NCOLS=" + ncstr + ones );
    readInputLine( getShortcutLabel() + "_diag: MATRIX_PRODUCT_DIAGONAL ARG=" + getShortcutLabel() + "_matones," + getShortcutLabel() + "_prod");
    // Sum everything
    if( !squared ) readInputLine( getShortcutLabel() + "_2: COMBINE ARG=" + getShortcutLabel() + "_offdiag," + getShortcutLabel() + "_diag PERIODIC=NO");
    else readInputLine( getShortcutLabel() + ": COMBINE ARG=" + getShortcutLabel() + "_offdiag," + getShortcutLabel() + "_diag PERIODIC=NO");
  } else {
    ActionWithValue* av=plumed.getActionSet().selectWithLabel<ActionWithValue*>( getShortcutLabel() + "_diffT" ); plumed_assert( av && av->getNumberOfComponents()==1 );
    if( (av->copyOutput(0))->getRank()==1 ) readInputLine( getShortcutLabel() + "_matvec: MATRIX_VECTOR_PRODUCT ARG=" + metstr + "," + getShortcutLabel() +"_diffT");
    else readInputLine( getShortcutLabel() + "_matvec: MATRIX_PRODUCT ARG=" + metstr + "," + getShortcutLabel() +"_diffT");
    std::string olab = getShortcutLabel(); if( !squared ) olab += "_2";
    readInputLine( olab + ": MATRIX_PRODUCT_DIAGONAL ARG=" + getShortcutLabel() + "_diff," + getShortcutLabel() +"_matvec");
  }
  if( !squared ) readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + getShortcutLabel() + "_2 FUNC=sqrt(x) PERIODIC=NO");
}

}
}
