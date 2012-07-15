#include "ActionWithVirtualAtom.h"
#include "ActionRegister.h"
#include "Vector.h"
#include "PlumedException.h"

using namespace std;

namespace PLMD{

//+PLUMEDOC VATOM GHOST 
/*
Calculate the absolute position of a ghost atom with fixed coordinates
in the local reference frame formed by three atoms. 
The computed ghost atom is stored as a virtual atom that can be accessed in
an atom list through the the label for the GHOST action that creates it.

\par Examples
The following input instructs plumed to print the distance between the
ghost atom and the center of mass for atoms 15,20:
\verbatim
GHOST ATOMS=1,5,10 COORDINATES=10.0,10.0,10.0 LABEL=c1
COM ATOMS=15,20       LABEL=c2
DISTANCE ATOMS=c1,c2  LABEL=d1
PRINT ARG=d1
\endverbatim
(See also \ref DISTANCE and \ref PRINT).

*/
//+ENDPLUMEDOC


class GenericGhostAtom:
  public ActionWithVirtualAtom
{
  vector<double> coord;
public:
  GenericGhostAtom(const ActionOptions&ao);
  void calculate();
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(GenericGhostAtom,"GHOST")

void GenericGhostAtom::registerKeywords(Keywords& keys){
  ActionWithVirtualAtom::registerKeywords(keys);
  keys.add("atoms","COORDINATES","coordinates of the ghost atom in the local reference frame");
}

GenericGhostAtom::GenericGhostAtom(const ActionOptions&ao):
  Action(ao),
  ActionWithVirtualAtom(ao)
{
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  plumed_massert(atoms.size()==3,"ATOMS should contain a list of three atoms");

  parseVector("COORDINATES",coord);
  plumed_massert(coord.size()==3,"COORDINATES should be a list of three real numbers");

  checkRead();
  log.printf("  of atoms");
  for(unsigned i=0;i<atoms.size();++i) log.printf(" %d",atoms[i].serial());
  log.printf("\n");
  requestAtoms(atoms);
}

void GenericGhostAtom::calculate(){
  Vector pos;
  vector<Tensor> deriv(getNumberOfAtoms());
  vector<Vector> n;

// first versor 
  Vector n01 = delta(getPosition(0), getPosition(1));
  n.push_back(n01/n01.modulo());

// auxiliary vector
  Vector n02 = delta(getPosition(0), getPosition(2));
  double n02_norm = n02.modulo();
  n02 /= n02_norm;

// second versor 
  n.push_back(crossProduct(n[0],n02));

// third versor 
  n.push_back(crossProduct(n[0],n[1]));

// origin of the reference system
  pos = getPosition(0);

  for(unsigned i=0;i<3;++i){
    pos += coord[i] * n[i];
  }

  setPosition(pos);
  setMass(1.0);
  setCharge(0.0);

// some useful tensors for derivatives 
  Tensor dn0d0  = (-Tensor::identity()+Tensor(n[0],n[0]))/n01.modulo();
  Tensor dn0d1  = (+Tensor::identity()-Tensor(n[0],n[0]))/n01.modulo();
  Tensor dn02d0 = (-Tensor::identity()+Tensor(n02,n02))/n02_norm;
  Tensor dn02d2 = (+Tensor::identity()-Tensor(n02,n02))/n02_norm;

// derivative of n1 = n0 x n02
  Tensor dn1d0, dn1d1, dn1d2;

  for(unsigned j=0;j<3;++j){
// derivative of n1 = n0 x n02 with respect to point 0
   Vector tmp00  = Vector( dn0d0(0,j),  dn0d0(1,j),  dn0d0(2,j));
   Vector tmp020 = Vector(dn02d0(0,j), dn02d0(1,j), dn02d0(2,j));
   Vector tmp0   = crossProduct(tmp00,n02) + crossProduct(n[0],tmp020);
// derivative of n1 = n0 x n02 with respect to point 1
   Vector tmp01  = Vector( dn0d1(0,j),  dn0d1(1,j),  dn0d1(2,j));
   Vector tmp1   = crossProduct(tmp01,n02);
// derivative of n1 = n0 x n02 with respect to point 2
   Vector tmp022 = Vector(dn02d2(0,j), dn02d2(1,j), dn02d2(2,j));
   Vector tmp2   = crossProduct(n[0],tmp022);
   for(unsigned i=0;i<3;++i) {
    dn1d0(i,j) = tmp0[i];
    dn1d1(i,j) = tmp1[i];
    dn1d2(i,j) = tmp2[i];
   }
  }

// Derivative of the last versor n2 = n0 x ( n0 x n02 ) = n0( n0 n02 ) - n02 
// Scalar product and derivatives
  double n0_n02 = dotProduct(n[0],n02); 
  Vector dn0_n02d0, dn0_n02d1, dn0_n02d2;

  for(unsigned i=0;i<3;++i){
   for(unsigned j=0;j<3;++j){
    dn0_n02d0[i] += dn0d0(j,i)*n02[j] + n[0][j]*dn02d0(j,i);
    dn0_n02d1[i] += dn0d1(j,i)*n02[j];
    dn0_n02d2[i] +=                     n[0][j]*dn02d2(j,i);
   }
  }

  Tensor dn2d0, dn2d1, dn2d2;

  for(unsigned i=0;i<3;++i){
   for(unsigned j=0;j<3;++j){
    dn2d0(i,j) = dn0d0(i,j) * n0_n02 + n[0][i] * dn0_n02d0[j] - dn02d0(i,j);
    dn2d1(i,j) = dn0d1(i,j) * n0_n02 + n[0][i] * dn0_n02d1[j];
    dn2d2(i,j) =                       n[0][i] * dn0_n02d2[j] - dn02d2(i,j);
   }
  }

// Finally, the derivative tensor
  deriv[0] = Tensor::identity() + coord[0]*dn0d0 + coord[1]*dn1d0 + coord[2]*dn2d0;
  deriv[1] =                      coord[0]*dn0d1 + coord[1]*dn1d1 + coord[2]*dn2d1;
  deriv[2] =                                       coord[1]*dn1d2 + coord[2]*dn2d2;

  setAtomsDerivatives(deriv);
}

}
