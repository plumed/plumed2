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
#include "Colvar.h"
#include "ActionRegister.h"
#include "PlumedMain.h"
#include "Matrix.h"

#include <string>
#include <cmath>
#include <cassert>

using namespace std;

namespace PLMD{

//+PLUMEDOC COLVAR RGYR
/*
Calculate the radius of gyration for a chain of atoms.

The radius of gyration is calculated using:

\f[
s_{\rm Gyr}=\Big ( \frac{\sum_i^{n}  
 \vert {r}_i -{r}_{\rm COM} \vert ^2 }{\sum_i^{n} m_i} \Big)^{1/2} 
\f]

with the position of the center of mass \f${r}_{\rm COM}\f$ given by:

\f[
{r}_{\rm COM}=\frac{\sum_i^{n} {r}_i\ m_i }{\sum_i^{n} m_i}
\f]

\bug This was a very quick implementation of RGYR for a project that I am working on. It has very little of the functionality that is available in plumed 1.0. 

\par Examples

The following input tells plumed to print the radius of gyration of the 
chain containing atoms 10 to 20.
\verbatim
GYRATION TYPE=RADIUS ATOMS=10-20 LABEL=rg
PRINT ARG=rg STRIDE=1 FILE=colvar 
\endverbatim
(See also \ref PRINT)

*/
//+ENDPLUMEDOC

class ColvarRGYR : public Colvar {
private:
  std::string Type;
  enum CV_TYPE {RADIUS, TRACE, GTPC_1, GTPC_2, GTPC_3, ASPHERICITY, ACYLINDRICITY, KAPPA2, RGYR_3, RGYR_2, RGYR_1, TOT};
  int  rg_type;
  bool use_masses;
public:
  static void registerKeywords( Keywords& keys );
  ColvarRGYR(const ActionOptions&);
  virtual void calculate();
};

PLUMED_REGISTER_ACTION(ColvarRGYR,"GYRATION")

void ColvarRGYR::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords( keys );
  keys.add("compulsory","TYPE","RADIUS","The type of calculation relative to the Gyration Tensor you want to perform");
  keys.add("atoms","ATOMS","the group of atoms that you are calculating the Gyration Tensor for");
  keys.addFlag("NOT_MASS_WEIGHTED",false,"set the masses of all the atoms equal to one");
}

ColvarRGYR::ColvarRGYR(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
use_masses(true)
{
  std::vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  if(atoms.size()==0) error("no atoms specified");
  bool not_use_masses=!use_masses;
  parseFlag("NOT_MASS_WEIGHTED",not_use_masses);
  use_masses=!not_use_masses;
  parse("TYPE",Type);
  checkRead();

  if(Type=="RADIUS") rg_type=RADIUS;
  else if(Type=="TRACE") rg_type=TRACE;
  else if(Type=="GTPC_1") rg_type=GTPC_1;
  else if(Type=="GTPC_2") rg_type=GTPC_2;
  else if(Type=="GTPC_3") rg_type=GTPC_3;
  else if(Type=="ASPHERICITY") rg_type=ASPHERICITY;
  else if(Type=="ACYLINDRICITY") rg_type=ACYLINDRICITY;
  else if(Type=="KAPPA2") rg_type=KAPPA2;
  else if(Type=="RGYR_3") rg_type=RGYR_3;
  else if(Type=="RGYR_2") rg_type=RGYR_2;
  else if(Type=="RGYR_1") rg_type=RGYR_1;
  else error("Unknown GYRATION type");

  switch(rg_type)
  {
    case RADIUS:   log.printf("  GYRATION RADIUS (Rg);"); break;
    case TRACE:  log.printf("  TRACE OF THE GYRATION TENSOR;"); break;
    case GTPC_1: log.printf("  THE LARGEST PRINCIPAL MOMENT OF THE GYRATION TENSOR (S'_1);"); break;
    case GTPC_2: log.printf("  THE MIDDLE PRINCIPAL MOMENT OF THE GYRATION TENSOR (S'_2);");  break;
    case GTPC_3: log.printf("  THE SMALLEST PRINCIPAL MOMENT OF THE GYRATION TENSOR (S'_3);"); break;
    case ASPHERICITY: log.printf("  THE ASPHERICITY (b');"); break;
    case ACYLINDRICITY: log.printf("  THE ACYLINDRICITY (c');"); break; 
    case KAPPA2: log.printf("  THE RELATIVE SHAPE ANISOTROPY (kappa^2);"); break;
    case RGYR_3: log.printf("  THE SMALLEST PRINCIPAL RADIUS OF GYRATION (r_g3);"); break;
    case RGYR_2: log.printf("  THE MIDDLE PRINCIPAL RADIUS OF GYRATION (r_g2);"); break;
    case RGYR_1: log.printf("  THE LARGEST PRINCIPAL RADIUS OF GYRATION (r_g1);"); break;
  }
  if(rg_type>TRACE) log<<"  Bibliography "<<plumed.cite("Jirí Vymetal and Jirí Vondrasek, J. Phys. Chem. A 115, 11455 (2011)"); log<<"\n";

  log.printf("  atoms involved : ");
  for(unsigned i=0;i<atoms.size();++i) log.printf("%d ",atoms[i].serial());
  log.printf("\n");

  addValueWithDerivatives(); setNotPeriodic();
  requestAtoms(atoms);
}

void ColvarRGYR::calculate(){

  std::vector<Vector> derivatives( getNumberOfAtoms() );
  Tensor virial; virial.zero();
  double totmass = 0.; 
  double d=0., rgyr=0.;
  Vector pos0, com, diff;
  Matrix<double> gyr_tens(3,3);
  for(unsigned i=0;i<3;i++)  for(unsigned j=0;j<3;j++) gyr_tens(i,j)=0.;

  // Find the center of mass
  if( use_masses ) totmass += getMass(0); 
  else totmass += 1.0;
  pos0=getPosition(0); 
  com.zero();
  for(unsigned i=1;i<getNumberOfAtoms();i++){
     diff=delta( pos0, getPosition(i) );
     if( use_masses ){
         totmass += getMass(i);
         com += getMass(i)*diff;
     } else {
         totmass += 1.0;
         com += diff;
     }
  }
  com = com / totmass + pos0;

  // Now compute radius of gyration
  for(unsigned i=0;i<getNumberOfAtoms();i++){
    diff=delta( com, getPosition(i) );
    d=diff.modulo();
    switch(rg_type) {
      case RADIUS:
      case TRACE:
        if( use_masses ){
          rgyr += getMass(i)*d*d;
          derivatives[i]=diff*getMass(i);
        } else {
          rgyr += d*d;
          derivatives[i]=diff;
        }
        break;
      default: //calculate gyration tensor
      {
        if( use_masses ) {
          gyr_tens[0][0]+=getMass(i)*diff[0]*diff[0];
          gyr_tens[1][1]+=getMass(i)*diff[1]*diff[1];
          gyr_tens[2][2]+=getMass(i)*diff[2]*diff[2];
          gyr_tens[0][1]+=getMass(i)*diff[0]*diff[1];
          gyr_tens[0][2]+=getMass(i)*diff[0]*diff[2];
          gyr_tens[1][2]+=getMass(i)*diff[1]*diff[2];
        } else {
          gyr_tens[0][0]+=diff[0]*diff[0];
          gyr_tens[1][1]+=diff[1]*diff[1];
          gyr_tens[2][2]+=diff[2]*diff[2];
          gyr_tens[0][1]+=diff[0]*diff[1];
          gyr_tens[0][2]+=diff[0]*diff[2];
          gyr_tens[1][2]+=diff[1]*diff[2];
        }
        break;
      }
    }
  }

  switch(rg_type) {
    case RADIUS:
    {
      rgyr=sqrt(rgyr/totmass);
      for(unsigned i=0;i<getNumberOfAtoms();i++){
        derivatives[i] /= rgyr*totmass;
        setAtomsDerivatives(i,derivatives[i]);
        virial=virial+(-1.0*Tensor(getPosition(i),derivatives[i]));
      }
      break;
    }
    case TRACE:
    {
      rgyr *= 2.;
      for(unsigned i=0;i<getNumberOfAtoms();i++) {
        derivatives[i] *= 4.;  
        setAtomsDerivatives(i,derivatives[i]);
        virial=virial+(-1.0*Tensor(getPosition(i),derivatives[i]));
      }
      break;
    }
    default:
    {
      // first make the matrix symmetric

      gyr_tens[1][0] = gyr_tens[0][1];
      gyr_tens[2][0] = gyr_tens[0][2];
      gyr_tens[2][1] = gyr_tens[1][2];

      Matrix<double> ttransf(3,3), transf(3,3);
      std::vector<double> princ_comp(3), prefactor(3);
      prefactor[0]=prefactor[1]=prefactor[2]=0.;

      //diagonalize gyration tensor
      diagMat(gyr_tens, princ_comp, ttransf);
      transpose(ttransf, transf);
      //sort eigenvalues and eigenvectors
      if (princ_comp[0]<princ_comp[1]){
        double tmp=princ_comp[0]; princ_comp[0]=princ_comp[1]; princ_comp[1]=tmp;
        for (unsigned i=0; i<3; i++){tmp=transf[i][0]; transf[i][0]=transf[i][1]; transf[i][1]=tmp;}
      }
      if (princ_comp[1]<princ_comp[2]){
        double tmp=princ_comp[1]; princ_comp[1]=princ_comp[2]; princ_comp[2]=tmp;
        for (unsigned i=0; i<3; i++){tmp=transf[i][1]; transf[i][1]=transf[i][2]; transf[i][2]=tmp;}
      }
      if (princ_comp[0]<princ_comp[1]){
        double tmp=princ_comp[0]; princ_comp[0]=princ_comp[1]; princ_comp[1]=tmp;
        for (unsigned i=0; i<3; i++){tmp=transf[i][0]; transf[i][0]=transf[i][1]; transf[i][1]=tmp;}      
      }
      //calculate determinant of transformation matrix
      double det;  
      det=transf[0][0]*transf[1][1]*transf[2][2]+transf[0][1]*transf[1][2]*transf[2][0]+transf[0][2]*transf[1][0]*transf[2][1]- transf[0][2]*transf[1][1]*transf[2][0]-transf[0][1]*transf[1][0]*transf[2][2]-transf[0][0]*transf[1][2]*transf[2][1]; 
      // trasformation matrix for rotation must have positive determinant, otherwise multiply one column by (-1)
      if (det<0) { for (unsigned j=0;j<3;j++) transf[j][2]=-transf[j][2];
                   det=transf[0][0]*transf[1][1]*transf[2][2]+transf[0][1]*transf[1][2]*transf[2][0]+transf[0][2]*transf[1][0]*transf[2][1]- transf[0][2]*transf[1][1]*transf[2][0]-transf[0][1]*transf[1][0]*transf[2][2]-transf[0][0]*transf[1][2]*transf[2][1]; 
      }
      if (fabs(det-1.)>0.0001) error("Plumed Error: Cannot diagonalize gyration tensor\n"); //check again, if det(transf)!=1 something is wrong, die

      switch(rg_type) {
        case GTPC_1:
        case GTPC_2:
        case GTPC_3:
        {
          int pc_index = rg_type-2; //index of principal component
          rgyr=sqrt(princ_comp[pc_index]/totmass);
          if(rgyr*totmass>1e-6) prefactor[pc_index]=1.0/(totmass*rgyr); //some parts of derivate
          break;
        }
	case RGYR_3:        //the smallest principal radius of gyration
        {
          rgyr=sqrt((princ_comp[1]+princ_comp[2])/totmass);
	  if (rgyr*totmass>1e-6){
            prefactor[1]=1.0/(totmass*rgyr);
            prefactor[2]=1.0/(totmass*rgyr);
	  }
          break;
        }
	case RGYR_2:       //the midle principal radius of gyration
        {
          rgyr=sqrt((princ_comp[0]+princ_comp[2])/totmass);
	  if (rgyr*totmass>1e-6){
            prefactor[0]=1.0/(totmass*rgyr);
            prefactor[2]=1.0/(totmass*rgyr);
	  }
          break;
        }
	case RGYR_1:      //the largest principal radius of gyration
        {
          rgyr=sqrt((princ_comp[0]+princ_comp[1])/totmass);
	  if (rgyr*totmass>1e-6){
            prefactor[0]=1.0/(totmass*rgyr);
            prefactor[1]=1.0/(totmass*rgyr);
	  }
          break;             
        }
	case ASPHERICITY:
        {
          rgyr=sqrt((princ_comp[0]-0.5*(princ_comp[1]+princ_comp[2]))/totmass); 
	  if (rgyr*totmass>1e-6){   //avoid division by zero 
            prefactor[0]= 1.0/(totmass*rgyr);
            prefactor[1]=-0.5/(totmass*rgyr);
            prefactor[2]=-0.5/(totmass*rgyr);
	  }
	  break;
        }
	case ACYLINDRICITY:
        {
          rgyr=sqrt((princ_comp[1]-princ_comp[2])/totmass); 
	  if (rgyr*totmass>1e-6){   //avoid division by zero  
            prefactor[1]= 1.0/(totmass*rgyr);
            prefactor[2]=-1.0/(totmass*rgyr);
	  }
          break;
        }
	case KAPPA2: // relative shape anisotropy
        {
          double trace = princ_comp[0]+princ_comp[1]+princ_comp[2];
          double tmp=princ_comp[0]*princ_comp[1]+ princ_comp[1]*princ_comp[2]+ princ_comp[0]*princ_comp[2]; 
          rgyr=1.0-3*(tmp/(trace*trace));
	  if (rgyr>1e-6){
            prefactor[0]= -3*((princ_comp[1]+princ_comp[2])-2*tmp/trace)/(trace*trace) *2;
            prefactor[1]= -3*((princ_comp[0]+princ_comp[2])-2*tmp/trace)/(trace*trace) *2;
            prefactor[2]= -3*((princ_comp[0]+princ_comp[1])-2*tmp/trace)/(trace*trace) *2;
	  }
          break;
        }
      }      
      for(unsigned i=0;i<getNumberOfAtoms();i++){
        Vector tX;
        diff=delta( com,getPosition(i) );
        //project atomic postional vectors to diagonalized frame
        for (unsigned j=0;j<3;j++) tX[j]=transf[0][j]*diff[0]+transf[1][j]*diff[1]+transf[2][j]*diff[2];  
        if( use_masses ) { 
          for (unsigned j=0;j<3;j++) 
            derivatives[i][j]=getMass(i)*(prefactor[0]*transf[j][0]*tX[0]+prefactor[1]*transf[j][1]*tX[1]+prefactor[2]*transf[j][2]*tX[2]);
        } else {                   
          for (unsigned j=0;j<3;j++) 
            derivatives[i][j]=(prefactor[0]*transf[j][0]*tX[0]+prefactor[1]*transf[j][1]*tX[1]+prefactor[2]*transf[j][2]*tX[2]);
        }
        setAtomsDerivatives(i,derivatives[i]);
        virial=virial+(-1.0*Tensor(getPosition(i),derivatives[i]));
      }
      break;
    }
  }
  setValue(rgyr);
  setBoxDerivatives(virial);
}

}
