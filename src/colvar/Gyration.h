/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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
#ifndef __PLUMED_colvar_Gyration_h
#define __PLUMED_colvar_Gyration_h
#include "Colvar.h"
#include "ColvarInput.h"
#include "MultiColvarTemplate.h"
#include "core/PlumedMain.h"

// The number of types of gyration are defined here
// The enum then contains double the number of types defined here as you can do every one of them with
// masses as weights or with unit weights
#define NUMBER_OF_GYRATION_TYPES 11

namespace PLMD {
namespace colvar {

template <typename T=double>
class Gyration : public Colvar {
private:
  enum class CV_TYPE:unsigned {RADIUS=0, TRACE=1, GTPC_1=2, GTPC_2=3, GTPC_3=4, ASPHERICITY=5, ACYLINDRICITY=6, KAPPA2=7, GYRATION_3=8, GYRATION_2=9, GYRATION_1=10,
                         RADIUS_MASS=NUMBER_OF_GYRATION_TYPES+RADIUS,TRACE_MASS=NUMBER_OF_GYRATION_TYPES+TRACE,
                         GTPC_1_MASS=NUMBER_OF_GYRATION_TYPES+GTPC_1, GTPC_2_MASS=NUMBER_OF_GYRATION_TYPES+GTPC_2, GTPC_3_MASS=NUMBER_OF_GYRATION_TYPES+GTPC_3,
                         ASPHERICITY_MASS=NUMBER_OF_GYRATION_TYPES+ASPHERICITY, ACYLINDRICITY_MASS=NUMBER_OF_GYRATION_TYPES+ACYLINDRICITY, KAPPA2_MASS=NUMBER_OF_GYRATION_TYPES+KAPPA2,
                         GYRATION_3_MASS=NUMBER_OF_GYRATION_TYPES+GYRATION_3, GYRATION_2_MASS=NUMBER_OF_GYRATION_TYPES+GYRATION_2, GYRATION_1_MASS=NUMBER_OF_GYRATION_TYPES+GYRATION_1
                        } rg_type{RADIUS};
  bool nopbc;
  std::vector<precision> value;
  std::vector<precision> derivs;
public:
  //move this on top of the class so that the two vectors above will know what precision is
  using precision = T;
  static void registerKeywords(Keywords& keys);
  explicit Gyration(const ActionOptions&);
  static void parseAtomList( const int& num,
                             std::vector<AtomNumber>& t,
                             ActionAtomistic* aa );
  static unsigned getModeAndSetupValues( ActionWithValue* av );
  void calculate() override;
  static void calculateCV( const ColvarInput<T>& cvin, ColvarOutput<T>& cvout );
};

template <typename T>
void Gyration<T>::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords(keys);
  keys.setDisplayName("GYRATION");
  keys.add("atoms","ATOMS","the group of atoms that you are calculating the Gyration Tensor for");
  keys.add("compulsory","TYPE","RADIUS","The type of calculation relative to the Gyration Tensor you want to perform");
  keys.addFlag("MASS_WEIGHTED",false,"use the masses of the atoms as the weights when calculating the gyration tensor");
  keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
  keys.setValueDescription("scalar","the radius of gyration");
}

template <typename T>
Gyration<T>::Gyration(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  nopbc(false),
  value(1) {
  std::vector<AtomNumber> atoms;
  parseAtomList(-1,atoms,this);
  if(atoms.size()==0) {
    error("no atoms specified");
  }

  parseFlag("NOPBC",nopbc);

  rg_type = static_cast<CV_TYPE>( getModeAndSetupValues(this) );
  checkRead();

  log.printf("  atoms involved : ");
  for(unsigned i=0; i<atoms.size(); ++i) {
    log.printf("%d ",atoms[i].serial());
  }
  log.printf("\n");

  if(nopbc) {
    log<<"  PBC will be ignored\n";
  } else {
    log<<"  broken molecules will be rebuilt assuming atoms are in the proper order\n";
  }

  requestAtoms(atoms);
}

template <typename T>
void Gyration<T>::parseAtomList( const int& num,
                                 std::vector<AtomNumber>& t,
                                 ActionAtomistic* aa ) {
  aa->parseAtomList("ATOMS",num,t);
}

template <typename T>
unsigned Gyration<T>::getModeAndSetupValues( ActionWithValue* av ) {
  bool use_masses;
  av->parseFlag("MASS_WEIGHTED",use_masses);
  std::string Type;
  av->parse("TYPE",Type);
  unsigned rg_type;

  if(Type=="RADIUS") {
    rg_type=CV_TYPE::RADIUS;
  } else if(Type=="TRACE") {
    rg_type=CV_TYPE::TRACE;
  } else if(Type=="GTPC_1") {
    rg_type=CV_TYPE::GTPC_1;
  } else if(Type=="GTPC_2") {
    rg_type=CV_TYPE::GTPC_2;
  } else if(Type=="GTPC_3") {
    rg_type=CV_TYPE::GTPC_3;
  } else if(Type=="ASPHERICITY") {
    rg_type=CV_TYPE::ASPHERICITY;
  } else if(Type=="ACYLINDRICITY") {
    rg_type=CV_TYPE::ACYLINDRICITY;
  } else if(Type=="KAPPA2") {
    rg_type=CV_TYPE::KAPPA2;
  } else if(Type=="RGYR_3") {
    rg_type=CV_TYPE::GYRATION_3;
  } else if(Type=="RGYR_2") {
    rg_type=CV_TYPE::GYRATION_2;
  } else if(Type=="RGYR_1") {
    rg_type=CV_TYPE::GYRATION_1;
  } else {
    av->error("Unknown GYRATION type");
  }
  if( use_masses ) {
    rg_type=static_cast<CV_TYPE>( rg_type + NUMBER_OF_GYRATION_TYPES );
    av->log.printf("  computing mass weighted averages \n");
  } else {
    av->log.printf("  atomic positions all have unit weights in averages \n");
  }

  switch(rg_type) {
  case CV_TYPE::RADIUS:
  case CV_TYPE::RADIUS_MASS:
    av->log.printf("  GYRATION RADIUS (Rg);");
    break;
  case CV_TYPE::TRACE:
  case CV_TYPE::TRACE_MASS:
    av->log.printf("  TRACE OF THE GYRATION TENSOR;");
    break;
  case CV_TYPE::GTPC_1:
  case CV_TYPE::GTPC_1_MASS:
    av->log.printf("  THE LARGEST PRINCIPAL MOMENT OF THE GYRATION TENSOR (S'_1);");
    break;
  case CV_TYPE::GTPC_2:
  case CV_TYPE::GTPC_2_MASS:
    av->log.printf("  THE MIDDLE PRINCIPAL MOMENT OF THE GYRATION TENSOR (S'_2);");
    break;
  case CV_TYPE::GTPC_3:
  case CV_TYPE::GTPC_3_MASS:
    av->log.printf("  THE SMALLEST PRINCIPAL MOMENT OF THE GYRATION TENSOR (S'_3);");
    break;
  case CV_TYPE::ASPHERICITY:
  case CV_TYPE::ASPHERICITY_MASS:
    av->log.printf("  THE ASPHERICITY (b');");
    break;
  case CV_TYPE::ACYLINDRICITY:
  case CV_TYPE::ACYLINDRICITY_MASS:
    av->log.printf("  THE ACYLINDRICITY (c');");
    break;
  case CV_TYPE::KAPPA2:
  case CV_TYPE::KAPPA2_MASS:
    av->log.printf("  THE RELATIVE SHAPE ANISOTROPY (kappa^2);");
    break;
  case CV_TYPE::GYRATION_3:
  case CV_TYPE::GYRATION_3_MASS:
    av->log.printf("  THE SMALLEST PRINCIPAL RADIUS OF GYRATION (r_g3);");
    break;
  case CV_TYPE::GYRATION_2:
  case CV_TYPE::GYRATION_2_MASS:
    av->log.printf("  THE MIDDLE PRINCIPAL RADIUS OF GYRATION (r_g2);");
    break;
  case CV_TYPE::GYRATION_1:
  case CV_TYPE::GYRATION_1_MASS:
    av->log.printf("  THE LARGEST PRINCIPAL RADIUS OF GYRATION (r_g1);");
    break;
  }
  av->log<<"\n";
  av->addValueWithDerivatives();
  av->setNotPeriodic();
  return rg_type;
}

template <typename T>
void Gyration<T>::calculate() {

  if(!nopbc) {
    makeWhole();
  }
  using CVInput= ColvarInput<double>;
  using CVOutput=ColvarOutput<double>;
  auto cvout =  CVOutput::createColvarOutput(value, derivs, this);
  Gyration<double>::calculateCV( CVInput::createColvarInput( rg_type, getPositions(), this ), cvout );
  for(unsigned i=0; i<getNumberOfAtoms(); ++i) {
    setAtomsDerivatives(i,cvout.getAtomDerivatives(0,i) );
  }
  setBoxDerivatives(cvout.virial[0]);
  setValue           (value[0]);
}

template <typename T>
void Gyration<T>::calculateCV( const ColvarInput<T>& cvin, ColvarOutput<T>& cvout ) {

  Vector com;
  double totmass = 0.;
  if( cvin.mode>=NUMBER_OF_GYRATION_TYPES ) {
    for(unsigned i=0; i<cvin.pos.size(); i++) {
      totmass+=cvin.mass[i];
      com+=cvin.mass[i]*Vector(cvin.pos[i][0],cvin.pos[i][1],cvin.pos[i][2]);
    }
  } else {
    totmass = static_cast<double>(cvin.pos.size());
    for(unsigned i=0; i<cvin.pos.size(); i++) {
      com+=Vector(cvin.pos[i][0],cvin.pos[i][1],cvin.pos[i][2]);
    }
  }
  com /= totmass;

  double rgyr=0.;

  if(cvin.mode==CV_TYPE::RADIUS||cvin.mode==CV_TYPE::TRACE||cvin.mode==CV_TYPE::RADIUS_MASS||cvin.mode==CV_TYPE::TRACE_MASS) {
    if( cvin.mode>=NUMBER_OF_GYRATION_TYPES ) {
      for(unsigned i=0; i<cvin.pos.size(); i++) {
        auto diff = delta( com, Vector(cvin.pos[i][0],cvin.pos[i][1],cvin.pos[i][2]) );
        rgyr          += cvin.mass[i]*diff.modulo2();
        cvout.derivs[0][i] = cvin.mass[i]*diff;
      }
    } else {
      for(unsigned i=0; i<cvin.pos.size(); i++) {
        const Vector diff = delta( com, Vector(cvin.pos[i][0],cvin.pos[i][1],cvin.pos[i][2]) );
        rgyr          += diff.modulo2();
        cvout.derivs[0][i] = diff;
      }
    }
    double fact;
    if(cvin.mode==CV_TYPE::RADIUS||cvin.mode==CV_TYPE::RADIUS_MASS) {
      rgyr = std::sqrt(rgyr/totmass);
      fact = 1./(rgyr*totmass);
    } else {
      rgyr = 2.*rgyr;
      fact = 4;
    }
    cvout.values[0]=rgyr;
    for(unsigned i=0; i<cvin.pos.size(); i++) {
      for(unsigned j=0; j<3; ++j) {
        cvout.derivs[0][i][j] = fact*cvout.derivs[0][i][j];
      }
    }
    ColvarInput<T>::setBoxDerivativesNoPbc( cvin, cvout );
    return;
  }


  Tensor3d gyr_tens;
  //calculate gyration tensor
  if( cvin.mode>=NUMBER_OF_GYRATION_TYPES ) {
    for(unsigned i=0; i<cvin.pos.size(); i++) {
      const Vector diff=delta( com, Vector(cvin.pos[i][0],cvin.pos[i][1],cvin.pos[i][2]) );
      gyr_tens[0][0]+=cvin.mass[i]*diff[0]*diff[0];
      gyr_tens[1][1]+=cvin.mass[i]*diff[1]*diff[1];
      gyr_tens[2][2]+=cvin.mass[i]*diff[2]*diff[2];
      gyr_tens[0][1]+=cvin.mass[i]*diff[0]*diff[1];
      gyr_tens[0][2]+=cvin.mass[i]*diff[0]*diff[2];
      gyr_tens[1][2]+=cvin.mass[i]*diff[1]*diff[2];
    }
  } else {
    for(unsigned i=0; i<cvin.pos.size(); i++) {
      const Vector diff=delta( com, Vector(cvin.pos[i][0],cvin.pos[i][1],cvin.pos[i][2]) );
      gyr_tens[0][0]+=diff[0]*diff[0];
      gyr_tens[1][1]+=diff[1]*diff[1];
      gyr_tens[2][2]+=diff[2]*diff[2];
      gyr_tens[0][1]+=diff[0]*diff[1];
      gyr_tens[0][2]+=diff[0]*diff[2];
      gyr_tens[1][2]+=diff[1]*diff[2];
    }
  }

  // first make the matrix symmetric
  gyr_tens[1][0] = gyr_tens[0][1];
  gyr_tens[2][0] = gyr_tens[0][2];
  gyr_tens[2][1] = gyr_tens[1][2];
  Tensor3d ttransf,transf;
  Vector princ_comp,prefactor;
  //diagonalize gyration tensor
  diagMatSym(gyr_tens, princ_comp, ttransf);
  transf=transpose(ttransf);
  //sort eigenvalues and eigenvectors
  if (princ_comp[0]<princ_comp[1]) {
    double tmp=princ_comp[0];
    princ_comp[0]=princ_comp[1];
    princ_comp[1]=tmp;
    for (unsigned i=0; i<3; i++) {
      tmp=transf[i][0];
      transf[i][0]=transf[i][1];
      transf[i][1]=tmp;
    }
  }
  if (princ_comp[1]<princ_comp[2]) {
    double tmp=princ_comp[1];
    princ_comp[1]=princ_comp[2];
    princ_comp[2]=tmp;
    for (unsigned i=0; i<3; i++) {
      tmp=transf[i][1];
      transf[i][1]=transf[i][2];
      transf[i][2]=tmp;
    }
  }
  if (princ_comp[0]<princ_comp[1]) {
    double tmp=princ_comp[0];
    princ_comp[0]=princ_comp[1];
    princ_comp[1]=tmp;
    for (unsigned i=0; i<3; i++) {
      tmp=transf[i][0];
      transf[i][0]=transf[i][1];
      transf[i][1]=tmp;
    }
  }
  //calculate determinant of transformation matrix
  double det = determinant(transf);
  // transformation matrix for rotation must have positive determinant, otherwise multiply one column by (-1)
  if(det<0) {
    for(unsigned j=0; j<3; j++) {
      transf[j][2]=-transf[j][2];
    }
    det = -det;
  }
  if(std::abs(det-1.)>0.0001) {
    plumed_merror("Plumed Error: Cannot diagonalize gyration tensor\n");
  }
  switch(cvin.mode) {
  case RADIUS:
  case RADIUS_MASS:
  case TRACE:
  case TRACE_MASS: {
    break;
  }
  case GTPC_1:
  case GTPC_2:
  case GTPC_3:
  case GTPC_1_MASS:
  case GTPC_2_MASS:
  case GTPC_3_MASS: {
    int pc_index = cvin.mode-2; //index of principal component
    if( pc_index>=NUMBER_OF_GYRATION_TYPES ) {
      pc_index = pc_index - NUMBER_OF_GYRATION_TYPES;
    }
    rgyr=std::sqrt(princ_comp[pc_index]/totmass);
    double rm = rgyr*totmass;
    if(rm>1e-6) {
      prefactor[pc_index]=1.0/rm;  //some parts of derivate
    }
    break;
  }
  case GYRATION_3:
  case GYRATION_3_MASS: {      //the smallest principal radius of gyration
    rgyr=std::sqrt((princ_comp[1]+princ_comp[2])/totmass);
    double rm = rgyr*totmass;
    if (rm>1e-6) {
      prefactor[1]=1.0/rm;
      prefactor[2]=1.0/rm;
    }
    break;
  }
  case GYRATION_2:
  case GYRATION_2_MASS: {     //the midle principal radius of gyration
    rgyr=std::sqrt((princ_comp[0]+princ_comp[2])/totmass);
    double rm = rgyr*totmass;
    if (rm>1e-6) {
      prefactor[0]=1.0/rm;
      prefactor[2]=1.0/rm;
    }
    break;
  }
  case GYRATION_1:
  case GYRATION_1_MASS: {    //the largest principal radius of gyration
    rgyr=std::sqrt((princ_comp[0]+princ_comp[1])/totmass);
    double rm = rgyr*totmass;
    if (rm>1e-6) {
      prefactor[0]=1.0/rm;
      prefactor[1]=1.0/rm;
    }
    break;
  }
  case ASPHERICITY:
  case ASPHERICITY_MASS: {
    rgyr=std::sqrt((princ_comp[0]-0.5*(princ_comp[1]+princ_comp[2]))/totmass);
    double rm = rgyr*totmass;
    if (rm>1e-6) {
      prefactor[0]= 1.0/rm;
      prefactor[1]=-0.5/rm;
      prefactor[2]=-0.5/rm;
    }
    break;
  }
  case ACYLINDRICITY:
  case ACYLINDRICITY_MASS: {
    rgyr=std::sqrt((princ_comp[1]-princ_comp[2])/totmass);
    double rm = rgyr*totmass;
    if (rm>1e-6) {  //avoid division by zero
      prefactor[1]= 1.0/rm;
      prefactor[2]=-1.0/rm;
    }
    break;
  }
  case KAPPA2:
  case KAPPA2_MASS: { // relative shape anisotropy
    double trace = princ_comp[0]+princ_comp[1]+princ_comp[2];
    double tmp=princ_comp[0]*princ_comp[1]+ princ_comp[1]*princ_comp[2]+ princ_comp[0]*princ_comp[2];
    rgyr=1.0-3*(tmp/(trace*trace));
    if (rgyr>1e-6) {
      prefactor[0]= -3*((princ_comp[1]+princ_comp[2])-2*tmp/trace)/(trace*trace) *2;
      prefactor[1]= -3*((princ_comp[0]+princ_comp[2])-2*tmp/trace)/(trace*trace) *2;
      prefactor[2]= -3*((princ_comp[0]+princ_comp[1])-2*tmp/trace)/(trace*trace) *2;
    }
    break;
  }
  }

  if( cvin.mode>=NUMBER_OF_GYRATION_TYPES ) {
    for(unsigned i=0; i<cvin.pos.size(); i++) {
      Vector tX;
      const Vector diff=delta( com,Vector(cvin.pos[i][0],cvin.pos[i][1],cvin.pos[i][2]) );
      //project atomic postional vectors to diagonalized frame
      for(unsigned j=0; j<3; j++) {
        tX[j]=transf[0][j]*diff[0]+transf[1][j]*diff[1]+transf[2][j]*diff[2];
      }
      for(unsigned j=0; j<3; j++)
        cvout.derivs[0][i][j]=cvin.mass[i]*(prefactor[0]*transf[j][0]*tX[0]+
                                            prefactor[1]*transf[j][1]*tX[1]+
                                            prefactor[2]*transf[j][2]*tX[2]);
    }
  } else {
    for(unsigned i=0; i<cvin.pos.size(); i++) {
      Vector tX;
      const Vector diff=delta( com, Vector(cvin.pos[i][0],cvin.pos[i][1],cvin.pos[i][2]) );
      //project atomic postional vectors to diagonalized frame
      for(unsigned j=0; j<3; j++) {
        tX[j]=transf[0][j]*diff[0]+transf[1][j]*diff[1]+transf[2][j]*diff[2];
      }
      for(unsigned j=0; j<3; j++)
        cvout.derivs[0][i][j]=prefactor[0]*transf[j][0]*tX[0]+
                              prefactor[1]*transf[j][1]*tX[1]+
                              prefactor[2]*transf[j][2]*tX[2];
    }
  }

  cvout.values[0]=rgyr;
  ColvarInput<T>::setBoxDerivativesNoPbc( cvin, cvout );
}

}
}
#endif
