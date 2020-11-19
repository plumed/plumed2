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
#include "Colvar.h"
#include "core/PlumedMain.h"
#include "ActionRegister.h"
#include "tools/Torsion.h"

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR PUCKERING
/*
 Calculate sugar pseudorotation coordinates.

 This command can be used to calculate ring's pseudorotations in sugars (puckers). It works for both
 5-membered and 6-membered rings. Notice that there are two different implementations depending if
 one passes 5 or 6 atoms in the ATOMS keyword.

 For 5-membered rings the implementation is the one discussed in \cite huang2014improvement .
 This implementation is simple and can be used in RNA to distinguish C2'-endo and C3'-endo conformations.
 Both the polar coordinates (phs and amp) and the Cartesian coordinates (Zx and Zy) are provided.
 C2'-endo conformations have negative Zx, whereas C3'-endo conformations have positive Zy.
 Notation is consistent with \cite huang2014improvement .
 The five atoms should be provided as C4',O4',C1',C2',C3'.
 Notice that this is the same order that can be obtained using the \ref MOLINFO syntax (see example below).

 For 6-membered rings the implementation is the general Cremer-Pople one \cite cremer1975general
 as also discussed in \cite biarnes2007conformational .
 This implementation provides both a triplet with Cartesian components (qx, qy, and qz)
 and a triplet of polar components (amplitude, phi, and theta).
 Applications of this particular implementation are to be published (paper in preparation).

 \note The 6-membered ring implementation distributed with previous versions of PLUMED lead to
 qx and qy values that had an opposite sign with respect to those originally defined in \cite cremer1975general.
 The bug is fixed in version 2.5.

 Components of this action are:

 \par Examples

 This input tells plumed to print the puckering phase angle of the 3rd nucleotide of a RNA molecule on file COLVAR.
 \plumedfile
 #SETTINGS MOLFILE=regtest/basic/rt65/AA.pdb
 MOLINFO STRUCTURE=rna.pdb MOLTYPE=rna
 PUCKERING ATOMS=@sugar-2 LABEL=puck
 PRINT ARG=puck.phs FILE=COLVAR
 \endplumedfile

*/
//+ENDPLUMEDOC

class Puckering : public Colvar {

public:
  explicit Puckering(const ActionOptions&);
  void calculate() override;
  static void registerKeywords(Keywords& keys);
  void calculate5m();
  void calculate6m();
};

PLUMED_REGISTER_ACTION(Puckering,"PUCKERING")

void Puckering::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords( keys );
  keys.remove("NOPBC");
  keys.add("atoms","ATOMS","the five or six atoms of the sugar ring in the proper order");
  keys.addOutputComponent("phs","default","Pseudorotation phase (5 membered rings)");
  keys.addOutputComponent("amp","default","Pseudorotation amplitude (5 membered rings)");
  keys.addOutputComponent("Zx","default","Pseudorotation x Cartesian component (5 membered rings)");
  keys.addOutputComponent("Zy","default","Pseudorotation y Cartesian component (5 membered rings)");
  keys.addOutputComponent("phi","default","Pseudorotation phase (6 membered rings)");
  keys.addOutputComponent("theta","default","Theta angle (6 membered rings)");
  keys.addOutputComponent("amplitude","default","Pseudorotation amplitude (6 membered rings)");
  keys.addOutputComponent("qx","default","Cartesian component x (6 membered rings)");
  keys.addOutputComponent("qy","default","Cartesian component y (6 membered rings)");
  keys.addOutputComponent("qz","default","Cartesian component z (6 membered rings)");
}

Puckering::Puckering(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao)
{
  std::vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  if(atoms.size()!=5 && atoms.size()!=6) error("only for 5 or 6-membered rings");
  checkRead();

  if(atoms.size()==5) {
    log.printf("  between atoms %d %d %d %d %d\n",atoms[0].serial(),atoms[1].serial(),atoms[2].serial(),atoms[3].serial(),atoms[4].serial());
  } else if(atoms.size()==6) {
    log.printf("  between atoms %d %d %d %d %d %d\n",atoms[0].serial(),atoms[1].serial(),atoms[2].serial(),atoms[3].serial(),atoms[4].serial(),atoms[5].serial());
  } else error("ATOMS should specify 5 atoms");

  if(atoms.size()==5) {
    addComponentWithDerivatives("phs"); componentIsPeriodic("phs","-pi","pi");
    addComponentWithDerivatives("amp"); componentIsNotPeriodic("amp");
    addComponentWithDerivatives("Zx"); componentIsNotPeriodic("Zx");
    addComponentWithDerivatives("Zy"); componentIsNotPeriodic("Zy");
  } else if(atoms.size()==6) {
    addComponentWithDerivatives("qx"); componentIsNotPeriodic("qx");
    addComponentWithDerivatives("qy"); componentIsNotPeriodic("qy");
    addComponentWithDerivatives("qz"); componentIsNotPeriodic("qz");
    addComponentWithDerivatives("phi"); componentIsPeriodic("phi","0","2pi");
    addComponentWithDerivatives("theta"); componentIsNotPeriodic("theta");
    addComponentWithDerivatives("amplitude"); componentIsNotPeriodic("amplitude");
  }

  log<<"  Bibliography ";
  if(atoms.size()==5) log<<plumed.cite("Huang, Giese, Lee, York, J. Chem. Theory Comput. 10, 1538 (2014)");
  if(atoms.size()==6) log<<plumed.cite("Cremer and Pople, J. Am. Chem. Soc. 97, 1354 (1975)");
  log<<"\n";

  requestAtoms(atoms);
}

// calculator
void Puckering::calculate() {
  makeWhole();
  if(getNumberOfAtoms()==5) calculate5m();
  else calculate6m();
}

void Puckering::calculate5m() {

  Vector d0,d1,d2,d3,d4,d5;

  d0=delta(getPosition(2),getPosition(1));
  d1=delta(getPosition(3),getPosition(2));
  d2=delta(getPosition(4),getPosition(3));
  d3=delta(getPosition(4),getPosition(3));
  d4=delta(getPosition(0),getPosition(4));
  d5=delta(getPosition(1),getPosition(0));

  Vector dd0,dd1,dd2,dd3,dd4,dd5;

  PLMD::Torsion t;

  double v1=t.compute(d0,d1,d2,dd0,dd1,dd2);
  double v3=t.compute(d3,d4,d5,dd3,dd4,dd5);

  double Zx=(v1+v3)/(2.0*cos(4.0*pi/5.0));
  double Zy=(v1-v3)/(2.0*sin(4.0*pi/5.0));
  double phase=atan2(Zy,Zx);
  double amplitude=std::sqrt(Zx*Zx+Zy*Zy);

  Vector dZx_dR[5];
  Vector dZy_dR[5];

  dZx_dR[0]=(dd5-dd4);
  dZx_dR[1]=(dd0-dd5);
  dZx_dR[2]=(dd1-dd0);
  dZx_dR[3]=(dd2+dd3-dd1);
  dZx_dR[4]=(dd4-dd3-dd2);

  dZy_dR[0]=(dd4-dd5);
  dZy_dR[1]=(dd0+dd5);
  dZy_dR[2]=(dd1-dd0);
  dZy_dR[3]=(dd2-dd3-dd1);
  dZy_dR[4]=(dd3-dd4-dd2);

  for(unsigned j=0; j<5; j++) dZx_dR[j]*=(1.0/(2.0*cos(4.0*pi/5.0)));
  for(unsigned j=0; j<5; j++) dZy_dR[j]*=(1.0/(2.0*sin(4.0*pi/5.0)));

  Vector dphase_dR[5];
  for(unsigned j=0; j<5; j++) dphase_dR[j]=(1.0/(Zx*Zx+Zy*Zy))*(-Zy*dZx_dR[j] + Zx*dZy_dR[j]);

  Vector damplitude_dR[5];
  for(unsigned j=0; j<5; j++) damplitude_dR[j]=(1.0/amplitude)*(Zx*dZx_dR[j] + Zy*dZy_dR[j]);

  Value* vzx=getPntrToComponent("Zx");
  vzx->set(Zx);
  setAtomsDerivatives (vzx,0, dZx_dR[0]);
  setAtomsDerivatives (vzx,1, dZx_dR[1]);
  setAtomsDerivatives (vzx,2, dZx_dR[2]);
  setAtomsDerivatives (vzx,3, dZx_dR[3]);
  setAtomsDerivatives (vzx,4, dZx_dR[4]);
  setBoxDerivativesNoPbc(vzx);

  Value* vzy=getPntrToComponent("Zy");
  vzy->set(Zy);
  setAtomsDerivatives (vzy,0, dZy_dR[0]);
  setAtomsDerivatives (vzy,1, dZy_dR[1]);
  setAtomsDerivatives (vzy,2, dZy_dR[2]);
  setAtomsDerivatives (vzy,3, dZy_dR[3]);
  setAtomsDerivatives (vzy,4, dZy_dR[4]);
  setBoxDerivativesNoPbc(vzy);


  Value* vph=getPntrToComponent("phs");
  vph->set(phase);
  setAtomsDerivatives (vph,0, dphase_dR[0]);
  setAtomsDerivatives (vph,1, dphase_dR[1]);
  setAtomsDerivatives (vph,2, dphase_dR[2]);
  setAtomsDerivatives (vph,3, dphase_dR[3]);
  setAtomsDerivatives (vph,4, dphase_dR[4]);
  setBoxDerivativesNoPbc(vph);

  Value* vam=getPntrToComponent("amp");
  vam->set(amplitude);
  setAtomsDerivatives (vam,0, damplitude_dR[0]);
  setAtomsDerivatives (vam,1, damplitude_dR[1]);
  setAtomsDerivatives (vam,2, damplitude_dR[2]);
  setAtomsDerivatives (vam,3, damplitude_dR[3]);
  setAtomsDerivatives (vam,4, damplitude_dR[4]);
  setBoxDerivativesNoPbc(vam);


}

void Puckering::calculate6m() {

  std::vector<Vector> r(6);
  for(unsigned i=0; i<6; i++) r[i]=getPosition(i);

  std::vector<Vector> R(6);
  Vector center;
  for(unsigned j=0; j<6; j++) center+=r[j]/6.0;
  for(unsigned j=0; j<6; j++) R[j]=(r[j]-center);

  Vector Rp,Rpp;
  for(unsigned j=0; j<6; j++) Rp +=R[j]*sin(2.0/6.0*pi*j);
  for(unsigned j=0; j<6; j++) Rpp+=R[j]*cos(2.0/6.0*pi*j);

  Vector n=crossProduct(Rp,Rpp);
  Vector nhat=n/modulo(n);

  Tensor dn_dRp=dcrossDv1(Rp,Rpp);
  Tensor dn_dRpp=dcrossDv2(Rp,Rpp);

  Tensor dnhat_dn= 1.0/modulo(n)*( Tensor::identity() - extProduct(nhat,nhat));
  Tensor dnhat_dRp=matmul(dnhat_dn,dn_dRp);
  Tensor dnhat_dRpp=matmul(dnhat_dn,dn_dRpp);

  std::vector<double> z(6);
  for(unsigned j=0; j<6; j++) z[j]=dotProduct(R[j],nhat);

  std::vector<std::vector<Vector> > dz_dR(6);
  for(unsigned j=0; j<6; j++) dz_dR[j].resize(6);

  for(unsigned i=0; i<6; i++) for(unsigned j=0; j<6; j++) {
      if(i==j) dz_dR[i][j]+=nhat;
      dz_dR[i][j]+=matmul(R[i],dnhat_dRp)*sin(2.0/6.0*pi*j);
      dz_dR[i][j]+=matmul(R[i],dnhat_dRpp)*cos(2.0/6.0*pi*j);
    }

  double B=0.0;
  for(unsigned j=0; j<6; j++) B+=z[j]*cos(4.0/6.0*pi*j);

  std::vector<Vector> dB_dR(6);
  for(unsigned i=0; i<6; i++) for(unsigned j=0; j<6; j++) {
      dB_dR[i]+=dz_dR[j][i]*cos(4.0/6.0*pi*j);
    }
  Vector Bsum;
  for(unsigned j=0; j<6; j++) Bsum+=dB_dR[j];
  for(unsigned j=0; j<6; j++) dB_dR[j]-=Bsum/6.0;;

  double A=0.0;
  for(unsigned j=0; j<6; j++) A+=z[j]*sin(4.0/6.0*pi*j);

  std::vector<Vector> dA_dR(6);
  for(unsigned i=0; i<6; i++) for(unsigned j=0; j<6; j++) {
      dA_dR[i]+=dz_dR[j][i]*sin(4.0/6.0*pi*j);
    }
  Vector Asum;
  for(unsigned j=0; j<6; j++) Asum+=dA_dR[j];
  for(unsigned j=0; j<6; j++) dA_dR[j]-=Asum/6.0;;

  double C=0.0;
  for(unsigned j=0; j<6; j++) C+=z[j]*Tools::fastpow(-1.0,(j));

  std::vector<Vector> dC_dR(6);
  for(unsigned i=0; i<6; i++) for(unsigned j=0; j<6; j++) {
      dC_dR[i]+=dz_dR[j][i]*Tools::fastpow(-1.0,(j));
    }

  Vector Csum;
  for(unsigned j=0; j<6; j++) Csum+=dC_dR[j];
  for(unsigned j=0; j<6; j++) dC_dR[j]-=Csum/6.0;;


// qx
  double qx = -A/std::sqrt(3);

// qx derivaties
  std::vector<Vector> dqx_dR(6);
  for(unsigned j=0; j<6; j++) {
    dqx_dR[j]=-1/std::sqrt(3) * dA_dR[j];
  }

  Value* vqx=getPntrToComponent("qx");
  vqx->set(qx);
  setAtomsDerivatives (vqx,0, dqx_dR[0] );
  setAtomsDerivatives (vqx,1, dqx_dR[1] );
  setAtomsDerivatives (vqx,2, dqx_dR[2] );
  setAtomsDerivatives (vqx,3, dqx_dR[3] );
  setAtomsDerivatives (vqx,4, dqx_dR[4] );
  setAtomsDerivatives (vqx,5, dqx_dR[5] );
  setBoxDerivativesNoPbc(vqx);

// qy
  double qy = B/std::sqrt(3);

// qy derivatives
  std::vector<Vector> dqy_dR(6);
  for(unsigned j=0; j<6; j++) {
    dqy_dR[j]=1/std::sqrt(3) * dB_dR[j];
  }

  Value* vqy=getPntrToComponent("qy");
  vqy->set(qy);
  setAtomsDerivatives (vqy,0, dqy_dR[0] );
  setAtomsDerivatives (vqy,1, dqy_dR[1] );
  setAtomsDerivatives (vqy,2, dqy_dR[2] );
  setAtomsDerivatives (vqy,3, dqy_dR[3] );
  setAtomsDerivatives (vqy,4, dqy_dR[4] );
  setAtomsDerivatives (vqy,5, dqy_dR[5] );
  setBoxDerivativesNoPbc(vqy);

// qz
  double qz = C/std::sqrt(6);

// qz derivatives
  std::vector<Vector> dqz_dR(6);
  for(unsigned j=0; j<6; j++) {
    dqz_dR[j]=1/std::sqrt(6) * dC_dR[j];
  }

  Value* vqz=getPntrToComponent("qz");
  vqz->set(qz);
  setAtomsDerivatives (vqz,0, dqz_dR[0] );
  setAtomsDerivatives (vqz,1, dqz_dR[1] );
  setAtomsDerivatives (vqz,2, dqz_dR[2] );
  setAtomsDerivatives (vqz,3, dqz_dR[3] );
  setAtomsDerivatives (vqz,4, dqz_dR[4] );
  setAtomsDerivatives (vqz,5, dqz_dR[5] );
  setBoxDerivativesNoPbc(vqz);


// PHASE
  double phi=atan2(-A,B);

// PHASE DERIVATIVES
  std::vector<Vector> dphi_dR(6);
  for(unsigned j=0; j<6; j++) {
    dphi_dR[j]=1.0/(A*A+B*B) * (-B*dA_dR[j] + A*dB_dR[j]);
  }

  Value* vphi=getPntrToComponent("phi");
  vphi->set(phi);
  setAtomsDerivatives (vphi,0, dphi_dR[0] );
  setAtomsDerivatives (vphi,1, dphi_dR[1] );
  setAtomsDerivatives (vphi,2, dphi_dR[2] );
  setAtomsDerivatives (vphi,3, dphi_dR[3] );
  setAtomsDerivatives (vphi,4, dphi_dR[4] );
  setAtomsDerivatives (vphi,5, dphi_dR[5] );
  setBoxDerivativesNoPbc(vphi);

//  AMPLITUDE
  double amplitude=std::sqrt((2*(A*A+B*B)+C*C)/6);

//  AMPLITUDE DERIVATIES
  std::vector<Vector> damplitude_dR(6);
  for (unsigned j=0; j<6; j++) {
    damplitude_dR[j]=0.5*std::sqrt(2.0/6.0)/(std::sqrt(A*A+B*B+0.5*C*C)) * (2*A*dA_dR[j] + 2*B*dB_dR[j] + C*dC_dR[j]);
  }

  Value* vamplitude=getPntrToComponent("amplitude");
  vamplitude->set(amplitude);
  setAtomsDerivatives (vamplitude,0, damplitude_dR[0] );
  setAtomsDerivatives (vamplitude,1, damplitude_dR[1] );
  setAtomsDerivatives (vamplitude,2, damplitude_dR[2] );
  setAtomsDerivatives (vamplitude,3, damplitude_dR[3] );
  setAtomsDerivatives (vamplitude,4, damplitude_dR[4] );
  setAtomsDerivatives (vamplitude,5, damplitude_dR[5] );
  setBoxDerivativesNoPbc(vamplitude);

//  THETA
  double theta=acos( C / std::sqrt(2.*(A*A+B*B) +C*C ) );

//  THETA DERIVATIVES
  std::vector<Vector> dtheta_dR(6);
  for(unsigned j=0; j<6; j++) {
    dtheta_dR[j]=1.0/(3.0*std::sqrt(2)*amplitude*amplitude) * (C/(std::sqrt(A*A+B*B)) * (A*dA_dR[j] + B*dB_dR[j]) - std::sqrt(A*A+B*B)*dC_dR[j]);
  }
  Value* vtheta=getPntrToComponent("theta");
  vtheta->set(theta);
  setAtomsDerivatives (vtheta,0, dtheta_dR[0] );
  setAtomsDerivatives (vtheta,1, dtheta_dR[1] );
  setAtomsDerivatives (vtheta,2, dtheta_dR[2] );
  setAtomsDerivatives (vtheta,3, dtheta_dR[3] );
  setAtomsDerivatives (vtheta,4, dtheta_dR[4] );
  setAtomsDerivatives (vtheta,5, dtheta_dR[5] );
  setBoxDerivativesNoPbc(vtheta);
}


}
}



