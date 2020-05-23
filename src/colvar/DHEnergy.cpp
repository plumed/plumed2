/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
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
#include "CoordinationBase.h"
#include "tools/SwitchingFunction.h"
#include "ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"

#include <iostream>

#include <string>

using namespace std;

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR DHENERGY
/*
Calculate Debye-Huckel interaction energy among GROUPA and GROUPB.

This variable calculates the electrostatic interaction among GROUPA and GROUPB
using a Debye-Huckel approximation defined as
\f[
\frac{1}{4\pi\epsilon_r\epsilon_0}
\sum_{i\in A} \sum_{j \in B} q_i q_j
\frac{e^{-\kappa |{\bf r}_{ij}|}}{|{\bf r}_{ij}|}
\f]

This collective variable can be used to analyze or induce electrostatically driven reactions \cite do13jctc.
Notice that the value of the DHENERGY is returned in plumed units (see \ref UNITS).

If GROUPB is empty, it will sum the N*(N-1)/2 pairs in GROUPA. This avoids computing
twice permuted indexes (e.g. pair (i,j) and (j,i)) thus running at twice the speed.

Notice that if there are common atoms between GROUPA and GROUPB their interaction is discarded.


\par Examples

\plumedfile
# this is printing the electrostatic interaction between two groups of atoms
dh: DHENERGY GROUPA=1-10 GROUPB=11-20 EPSILON=80.0 I=0.1 TEMP=300.0
PRINT ARG=dh
\endplumedfile

*/
//+ENDPLUMEDOC

class DHEnergy : public CoordinationBase {
  double k; // Inverse Debye screening length
  double constant;
  double epsilon;

public:
  explicit DHEnergy(const ActionOptions&);
// active methods:
  static void registerKeywords( Keywords& keys );
  double pairing(double distance,double&dfunc,unsigned i,unsigned j)const override;
};

PLUMED_REGISTER_ACTION(DHEnergy,"DHENERGY")

void DHEnergy::registerKeywords( Keywords& keys ) {
  CoordinationBase::registerKeywords(keys);
  keys.add("compulsory","I","1.0","Ionic strength (M)");
  keys.add("compulsory","TEMP","300.0","Simulation temperature (K)");
  keys.add("compulsory","EPSILON","80.0","Dielectric constant of solvent");
}

/*
Global constants in SI unit used in this calculation:
      N_A = 6.0221412927 * 10^(23) mol^(-1) : Avogadro number
      q = 1.60217656535 * 10^(-19) C : proton charge
      e_0 = 8.854187817620 * 10^(-12) C^2/(N*m^2) : vacuum's dielectric constant
      k_B = 1.380648813 * 10^(-23) N*m/K : Boltzmann constant
In SI unit, Debye Huckel CV is defined as:
      DHen = \sum_i\sum_j (q_i*q_j*q^2*N_A)/(4*pi*eps*e_0) * exp(-k*|f_ij|)/(|f_ij|)
             + \sum_i\sum_j (q_i*q_j*q^2*N_A)/(4*pi*epp*e_0) * (1/|r_ij| - 1/|f_ij|)
           = (q^2*N_A)/(4*pi*e_0) * \sum_i\sum_j q_i*q_j * (exp(-k*|f_ij|)/(eps*|f_ij|) + 1/epp*(1/|r_ij| - 1/|f_ij|))
(in which |f_ij| = \sqrt(|r_ij|^2+\sigma_i*\sigma_j*exp(-|r_ij|^2/4*\sigma_i*\sigma_j)),
 \sigma_i and \sigma_j are the effective Born radius.)
For an efficient calculation, we group constants and variables into groups:
      constant = (q^2*N_A)/(4*pi*e_0)
      tmp = 1/eps*exp(-k*|f_ij|)/(|f_ij|) + 1/epp*(1/|r_ij| - 1/|f_ij|)

To speed up the loop calculation, constant can be modified as followed:
      constant= (q^2*N_A)/(4*pi*e_0*10^(-9))*10^(-3) (kJ/mol)
              = ((1.60217656535*10^(-19))^2*6.0221412927*10^(23)*10^(-3))/(4*3.14159265*8.854187817620*10^(-12)*10^(-9))
              = 138.935458111 (kJ/mol)

*/

DHEnergy::DHEnergy(const ActionOptions&ao):
  Action(ao),
  CoordinationBase(ao),
  k(0.0),
  constant(0.0)
{
  double I,T;
  parse("I",I);
  parse("TEMP",T);
  parse("EPSILON",epsilon);
  checkRead();
  if( plumed.getAtoms().usingNaturalUnits() ) error("DHENERGY cannot be used for calculations performed with natural units");
  constant=138.935458111/atoms.getUnits().getEnergy()/atoms.getUnits().getLength()*atoms.getUnits().getCharge()*atoms.getUnits().getCharge();
  k=sqrt(I/(epsilon*T))*502.903741125*atoms.getUnits().getLength();
  checkRead();
  log<<"  with solvent dielectric constant "<<epsilon<<"\n";
  log<<"  at temperature "<<T<<" K\n";
  log<<"  at ionic strength "<<I<< "M\n";
  log<<"  these parameters correspond to a screening length of "<<(1.0/k)<<"\n";
  log<<"  Bibliography "<<plumed.cite("Do, Carloni, Varani and Bussi, J. Chem. Theory Comput. 9, 1720 (2013)")<<" \n";
}

double DHEnergy::pairing(double distance2,double&dfunc,unsigned i,unsigned j)const {
  double distance=std::sqrt(distance2);
  if(getAbsoluteIndex(i)==getAbsoluteIndex(j)) {
    dfunc=0.0;
    return 0.0;
  }
  double invdistance=1.0/distance;
  double tmp=exp(-k*distance)*invdistance*constant*getCharge(i)*getCharge(j)/epsilon;
  double dtmp=-(k+invdistance)*tmp;
  dfunc=dtmp*invdistance;
  return tmp;
}

}

}
