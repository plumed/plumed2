/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015,2016 The plumed team
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
#include "AlignedMatrixBase.h"
#include "core/ActionRegister.h"
#include "tools/Matrix.h"

//+PLUMEDOC MATRIX ALIGNED_MATRIX 
/*
Adjacency matrix in which two molecule are adjacent if they are within a certain cutoff and if they have the same orientation.

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace adjmat {

class ContactAlignedMatrix : public AlignedMatrixBase {
private:
   Matrix<SwitchingFunction> sf;
public:
   /// 
   static void registerKeywords( Keywords& keys );
   ///
   explicit ContactAlignedMatrix(const ActionOptions&);
   void readOrientationConnector( const unsigned& i, const unsigned& j, const std::vector<std::string>& desc );
   double computeVectorFunction( const unsigned& iv, const unsigned& jv,
                                 const Vector& conn, const std::vector<double>& vec1, const std::vector<double>& vec2,
                                 Vector& dconn, std::vector<double>& dvec1, std::vector<double>& dvec2 ) const ;
};

PLUMED_REGISTER_ACTION(ContactAlignedMatrix,"ALIGNED_MATRIX")

void ContactAlignedMatrix::registerKeywords( Keywords& keys ){
   AlignedMatrixBase::registerKeywords( keys );
   keys.add("numbered","ORIENTATION_SWITCH","A switching function that transforms the dot product of the input vectors.");
}

ContactAlignedMatrix::ContactAlignedMatrix( const ActionOptions& ao ):
Action(ao),
AlignedMatrixBase(ao)
{
   unsigned nrows, ncols, ig; retrieveTypeDimensions( nrows, ncols, ig ); 
   sf.resize( nrows, ncols );
   parseConnectionDescriptions("ORIENTATION_SWITCH",false,0);
}

void ContactAlignedMatrix::readOrientationConnector( const unsigned& i, const unsigned& j, const std::vector<std::string>& desc ){
   plumed_assert( desc.size()==1 ); std::string errors; sf(j,i).set(desc[0],errors);
   if( j!=i ) sf(i,j).set(desc[0],errors);
   log.printf("  vectors in %u th and %u th groups must have a dot product that is greater than %s \n",i+1,j+1,(sf(i,j).description()).c_str() ); 
}

double ContactAlignedMatrix::computeVectorFunction( const unsigned& iv, const unsigned& jv, 
                                                    const Vector& conn, const std::vector<double>& vec1, const std::vector<double>& vec2,
                                                    Vector& dconn, std::vector<double>& dvec1, std::vector<double>& dvec2 ) const {
   double dot_df, dot=0; dconn.zero();
   for(unsigned k=2;k<vec1.size();++k) dot+=vec1[k]*vec2[k]; 
   double f_dot = sf(iv,jv).calculate( dot, dot_df ); 
   for(unsigned k=2;k<vec1.size();++k){ dvec1[k]=dot_df*vec2[k]; dvec2[k]=dot_df*vec1[k]; }
   return f_dot;
}

}
}



