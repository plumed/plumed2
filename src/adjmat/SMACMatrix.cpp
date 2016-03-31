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
#include "tools/KernelFunctions.h"
#include "tools/Torsion.h"
#include "tools/Matrix.h"

//+PLUMEDOC MATRIX SMAC_MATRIX 
/*
Adjacency matrix in which two molecules are adjacent if they are within a certain cutoff and if the angle between them is within certain ranges.


\bug Contribution to virial is wrong

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace adjmat {

class SMACMatrix : public AlignedMatrixBase {
private:
   Matrix<std::vector<KernelFunctions> > kernels;
public:
   /// 
   static void registerKeywords( Keywords& keys );
   ///
   explicit SMACMatrix(const ActionOptions&);
   void readOrientationConnector( const unsigned& i, const unsigned& j, const std::vector<std::string>& desc );
   double computeVectorFunction( const unsigned& iv, const unsigned& jv, 
                                 const Vector& conn, const std::vector<double>& vec1, const std::vector<double>& vec2,
                                 Vector& dconn, std::vector<double>& dvec1, std::vector<double>& dvec2 ) const ;
};

PLUMED_REGISTER_ACTION(SMACMatrix,"SMAC_MATRIX")

void SMACMatrix::registerKeywords( Keywords& keys ){
   AlignedMatrixBase::registerKeywords( keys );
   keys.add("numbered","KERNEL","The various kernels that are used to determine whether or not the molecules are aligned");
}

SMACMatrix::SMACMatrix( const ActionOptions& ao ):
Action(ao),
AlignedMatrixBase(ao)
{
   unsigned nrows, ncols, ig; retrieveTypeDimensions( nrows, ncols, ig );
   kernels.resize( nrows, ncols ); parseConnectionDescriptions("KERNEL",true,0);
}

void SMACMatrix::readOrientationConnector( const unsigned& iv, const unsigned& jv, const std::vector<std::string>& desc ){
   for(int i=0;i<desc.size();i++){
      KernelFunctions mykernel( desc[i], false );
      kernels(iv,jv).push_back( mykernel );
      if( jv!=iv ) kernels(jv,iv).push_back( mykernel );
   }
   if( kernels(iv,jv).size()==0 ) error("no kernels defined");
}

double SMACMatrix::computeVectorFunction( const unsigned& iv, const unsigned& jv,
                                          const Vector& conn, const std::vector<double>& vec1, const std::vector<double>& vec2,
                                          Vector& dconn, std::vector<double>& dvec1, std::vector<double>& dvec2 ) const {
  double dot=0; Vector v1, v2;
  for(unsigned k=2;k<vec1.size();++k){
      dot+=vec1[k]*vec2[k];
      v1[k-2]=vec1[k]; v2[k-2]=vec2[k];
  }
  Vector dv1, dv2; Torsion t;
  double angle = t.compute( v1, conn, v2, dv1, dconn, dv2 );
  
  std::vector<Value*> pos; pos.push_back( new Value() ); std::vector<double> deriv(1);
  pos[0]->setDomain( "-pi", "pi" ); pos[0]->set( angle );
  double ans=0, df=0; 
  for(unsigned i=0;i<kernels(iv,jv).size();++i){
      ans += kernels(iv,jv)[i].evaluate( pos, deriv );
      df += deriv[0];           
  }
  dconn*=df; for(unsigned k=2;k<vec1.size();++k){ dvec1[k]=df*dv1[k-2]; dvec2[k]=df*dv2[k-2]; }
  delete pos[0]; return ans;
}

}
}



