/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2019 The plumed team
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
#include "ActionRegister.h"
#include "core/PlumedMain.h"

#include <string>
#include <cmath>
#include <cassert>
#include <iostream>
#include <vector>

using namespace std;

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR DIMER
/*
This CV computes the dimer interaction energy for a collection of dimers.

Each dimer represents an atom, as described in the dimer paper \cite dimer-metad.
A system of N atoms is thus represented with N dimers, each
Dimer being composed of two beads and eventually a virtual site representing its center of mass.

A typical configuration for a dimerized system has the following ordering of atoms:

1    TAG1 X Y Z          N atoms representing the first bead of each Dimer

2    TAG2 X Y Z

...

N    TAGN X Y Z          N atoms representing the second bead of each Dimer

N+1  TAG1 X Y Z

N+2  TAG2 X Y Z

...

2N   TAGN X Y Z          Optional: N atoms representing the center of mass of each Dimer

2N+1 TAG1 X Y Z

2N+2 TAG2 X Y Z

...

3N   TAGN X Y Z          The configuration might go on with un-dimerized atoms (like a solvent)

3N+1

3N+2

...


The Dimer interaction energy is defined between atoms x and N+x, for x=1,...,N and is
characterized by two parameters Q and DSIGMA. These are passed as mandatory arguments along with
the temperature of the system.

\par Examples

This line tells Plumed to compute the Dimer interaction energy for every dimer in the system.

\plumedfile
dim: DIMER TEMP=300 Q=0.5 ALLATOMS DSIGMA=0.002
\endplumedfile

If the simulation doesn't use virtual sites for the dimers centers of mass,
Plumed has to know in order to determine correctly the total number of dimers from
the total number of atoms:
\plumedfile
dim: DIMER TEMP=300 Q=0.5 ALLATOMS DSIGMA=0.002 NOVSITES
\endplumedfile

The NOVSITES flag is not required if one provides the atom serials of each Dimer. These are
defined through two lists of atoms provided __instead__ of the ALLATOMS keyword.
For example, the Dimer interaction energy of dimers specified by beads (1;23),(5;27),(7;29) is:
\plumedfile
dim: DIMER TEMP=300 Q=0.5 ATOMS1=1,5,7 ATOMS2=23,27,29 DSIGMA=0.002
\endplumedfile

Note that the ATOMS1,ATOMS2 keywords can support atom groups and
interval notation as defined in \ref GROUP.


In a Replica Exchange simulation the keyword DSIGMA can be used in two ways:
if a plumed.n.dat file is provided for each replica, then DSIGMA is passed as a single value,
like in the previous examples, and each replica will read its own DSIGMA value. If
a unique plumed.dat is given, DSIGMA has to be a list containing a value for each replica.
For 4 replicas:
\plumedfile
#SETTINGS NREPLICAS=4
dim: DIMER TEMP=300 Q=0.5 ATOMS1=1,5,7 ATOMS2=23,27,29 DSIGMA=0.002,0.002,0.004,0.01
\endplumedfile


\par Usage of the CV

The dimer interaction is not coded in the driver program and has to be inserted
in the Hamiltonian of the system as a linear RESTRAINT (see \ref RESTRAINT):
\plumedfile
dim: DIMER TEMP=300 Q=0.5 ALLATOMS DSIGMA=0.002
RESTRAINT ARG=dim AT=0 KAPPA=0 SLOPE=1 LABEL=dimforces
\endplumedfile

In a replica exchange, Metadynamics (see \ref METAD) can be used on the Dimer CV to reduce
the number of replicas. Just keep in mind that METAD SIGMA values should be tuned
in the standard way for each replica according to the value of DSIGMA.
*/
//+ENDPLUMEDOC

class Dimer : public Colvar {
public:
  static void registerKeywords( Keywords& keys);
  explicit Dimer(const ActionOptions&);
  void calculate() override;
protected:
  bool trimer,useall;
  int myrank, nranks;
  double qexp,temperature,beta,dsigma;
  vector<double> dsigmas;
private:
  void consistencyCheck();
  vector<AtomNumber> usedatoms1;
  vector<AtomNumber> usedatoms2;

};

PLUMED_REGISTER_ACTION(Dimer, "DIMER")



void Dimer::registerKeywords( Keywords& keys) {
  Colvar::registerKeywords(keys);

  keys.add("compulsory","DSIGMA","The interaction strength of the dimer bond.");
  keys.add("compulsory", "Q", "The exponent of the dimer potential.");
  keys.add("compulsory", "TEMP", "The temperature (in Kelvin) of the simulation.");
  keys.add("atoms", "ATOMS1", "The list of atoms representing the first bead of each Dimer being considered by this CV. Used if ALLATOMS flag is missing");
  keys.add("atoms", "ATOMS2", "The list of atoms representing the second bead of each Dimer being considered by this CV. Used if ALLATOMS flag is missing");
  keys.addFlag("ALLATOMS", false, "Use EVERY atom of the system. Overrides ATOMS keyword.");
  keys.addFlag("NOVSITES", false, "If present the configuration is without virtual sites at the centroid positions.");

}



Dimer::Dimer(const ActionOptions& ao):
  PLUMED_COLVAR_INIT(ao)
{

  log<<" Bibliography "<<plumed.cite("M Nava, F. Palazzesi, C. Perego and M. Parrinello, J. Chem. Theory Comput. 13, 425(2017)")<<"\n";
  parseVector("DSIGMA",dsigmas);
  parse("Q",qexp);
  parse("TEMP",temperature);


  vector<AtomNumber> atoms;
  parseFlag("ALLATOMS",useall);
  trimer=true;
  bool notrim;
  parseFlag("NOVSITES",notrim);
  trimer=!notrim;

  nranks=multi_sim_comm.Get_size();
  myrank=multi_sim_comm.Get_rank();
  if(dsigmas.size()==1)
    dsigma=dsigmas[0];
  else
    dsigma=dsigmas[myrank];




  if(useall)
  {
    // go with every atom in the system but not the virtuals...
    int natoms;
    if(trimer)
      natoms= 2*getTotAtoms()/3;
    else
      natoms=getTotAtoms()/2;

    for(unsigned int i=0; i<((unsigned int)natoms); i++)
    {
      AtomNumber ati;
      ati.setIndex(i);
      atoms.push_back(ati);
    }
  }
  else  // serials for the first beads of each dimer are given
  {
    parseAtomList("ATOMS1",usedatoms1);
    parseAtomList("ATOMS2",usedatoms2);

    int isz1 = usedatoms1.size();

    for(unsigned int i=0; i<isz1; i++)
    {
      AtomNumber ati;
      ati.setIndex(usedatoms1[i].index());
      atoms.push_back(ati);
    }

    int isz2 = usedatoms2.size();
    for(unsigned int i=0; i<isz2; i++)
    {
      AtomNumber atip2;
      atip2.setIndex(usedatoms2[i].index());
      atoms.push_back(atip2);
    }

  }
  consistencyCheck();
  checkRead();
  beta = 1./(kBoltzmann*temperature);

  addValueWithDerivatives();  // allocate
  requestAtoms(atoms);
  setNotPeriodic();
}

void Dimer::calculate()
{
  double cv_val=0;
  Tensor virial;
  vector<Vector> derivatives;
  vector<Vector> my_pos=getPositions();
  int atms = my_pos.size();
  vector<Vector> der_b2;
  for(int i=0; i<atms/2; i++)
  {
    Vector dist;
    dist = pbcDistance(my_pos[i],my_pos[i+atms/2]);
    double distquad=0;
    for(int j=0; j<3; j++)
      distquad += dist[j]*dist[j];

    double dsigquad = dsigma*dsigma;
    double fac1 = 1.0 + distquad/(2*qexp*dsigquad);
    double fac1qm1 = pow(fac1,qexp-1);


    cv_val += (fac1*fac1qm1-1.0)/beta;
    Vector der_val;
    Vector mder_val;
    for(int j=0; j<3; j++)
    {
      der_val[j] = -fac1qm1*dist[j]/(dsigquad*beta);
      mder_val[j]=-der_val[j];
    }
    derivatives.push_back(der_val);
    der_b2.push_back(mder_val);

    // virial part: each dimer contributes -x_{ij}*ds/dx_{ij}  (s is the CV)
    double dfunc = fac1qm1/(beta*dsigquad);
    Vector dd(dfunc*dist);
    Tensor vv(dd,dist);
    virial -= vv;

  }

  derivatives.insert(derivatives.end(), der_b2.begin(), der_b2.end());

  for(unsigned int i=0; i<derivatives.size(); i++)
    setAtomsDerivatives(i,derivatives[i]);

  setValue(cv_val);
  setBoxDerivatives(virial);

}



/*****************
There are some conditions that a valid input should satisfy.
These are checked here and PLUMED error handlers are (eventually) called.
******************/
void Dimer::consistencyCheck()
{
  if(!useall && usedatoms1.size()!=usedatoms2.size())
    error("The provided atom lists are of different sizes.");

  if(qexp<0.5 || qexp>1)
    warning("Dimer CV is meant to be used with q-exponents between 0.5 and 1. We are not responsible for any black hole. :-)");

  if(dsigma<0)
    error("Please use positive sigma values for the Dimer strength constant");

  if(temperature<0)
    error("Please, use a positive value for the temperature...");

  // if dsigmas has only one element means that either
  // you are using different plumed.x.dat files or a plumed.dat with a single replica
  if(dsigmas.size()!=nranks && dsigmas.size()!=1)
    error("Mismatch between provided sigmas and number of replicas");

}


}
}

