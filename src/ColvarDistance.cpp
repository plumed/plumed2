#include "ColvarWithModifiers.h"
#include "ActionRegister.h"

#include <string>
#include <cmath>
#include <cassert>

using namespace std;

namespace PLMD{

//+PLUMEDOC COLVAR DISTANCE
/**
Calculate distances between atoms.  To calculate minimum distances, the number of distances less or more than a given value etc you should use this
keyword in conjuction with the colvar modifiers described in \ref Colvar. 

\par Example

The following calculates the distance between atoms 3 and 5 and stores the value on d1.value0. 
\verbatim
DISTANCE ATOMS=3,5 LABEL=d1
\endverbatim

The following calculates two distances.  That between atoms 3 and 5 and that between atoms 2 and 4.  These two values
are stored on d1.value0 and d1.value1
\verbatim
DISTANCE ATOMS1=3,5 ATOMS2=2,4 LABEL=d1
\endverbatim
 
The following calculates all the distance between the atoms in the group - i.e. the distances between atoms 3 and 4, between
3 and 5 and between 4 and 5.

\verbatim
DISTANCE GROUP=3,4,5 LABEL=d1
\endverbatim

Lastly, this calculates all the distances between the atoms in the two groups - i.e. the distances between 3 and 4 and 3 and 5.

\verbatim
DISTNACE GROUP1=3 GROUP2=4,5 LABEL=d1
\endverbatim

*/
//+ENDPLUMEDOC
   
class ColvarDistance : public ColvarWithModifiers {
  int component;
public:
  ColvarDistance(const ActionOptions&);
// active methods:
  double compute( const std::vector<unsigned>& indexes, std::vector<Vector>& derivatives, Tensor& virial );
  void interpretGroupsKeyword( const std::vector<unsigned>& boundaries, unsigned& maxatoms );
};

PLUMED_REGISTER_ACTION(ColvarDistance,"DISTANCE")

ColvarDistance::ColvarDistance(const ActionOptions&ao):
ColvarWithModifiers(ao),
component(-1)
{
  setNeighbourListStyle("skipAll");
  registerKeyword(2,"ATOMS","calculate the distance between this pair of atoms.  To calculate multiple distances use ATOMS1, ATOMS2, ...");
  registerKeyword(2,"GROUP","calculate the distances between every pair of atoms in the group.  Alternatively use GROUP1, GROUP2 to calculate all the distances between the atoms in group 1 and the atoms in group 2."); 
  registerKeyword(0,"COMPONENT","use this if you only want the X,Y or Z component of the distance");

  readActionAtomistic();
  // Read in the atoms
  int maxatoms=2; 
  readAtomsKeyword(maxatoms);
  readGroupsKeyword(maxatoms);

  // Setup the neighbour list
  std::vector< std::pair<unsigned, unsigned> > tmp_pair(1); 
  tmp_pair[0].first=0; tmp_pair[0].second=1;
  setupNeighbourList( tmp_pair ); 
  
  std::string lab="none"; parse("COMPONENT",lab);
  if( lab=="X" || lab=="x" ){
     component=0;
  } else if( lab=="Y" || lab=="y" ){
     component=1;
  } else if( lab=="Z" || lab=="z" ){
     component=2;
  } else if( lab!="none") {
      error( lab + " is not a valid argument for the COMPONENT keyword use X, Y or Z");
  } 
  finishColvarSetup( 0, 0 );
}

void ColvarDistance::interpretGroupsKeyword( const std::vector<unsigned>& boundaries, unsigned& maxatoms ){
   AtomicNeighbourList nlist( dynamic_cast<ActionAtomistic*>(this) );
   if( boundaries.size()==1 ){
       assert(boundaries[0]==getNumberOfAtoms());
       for(unsigned i=1;i<getNumberOfAtoms();++i){
           for(unsigned j=0;j<i;++j){ 
               nlist.clear(); nlist.addAtom(i); nlist.addAtom(j); addColvar( nlist );
           }
       }
   } else if ( boundaries.size()==2 ){
       assert(boundaries[1]==getNumberOfAtoms());
       for(unsigned i=0;i<boundaries[0];++i){
           for(unsigned j=boundaries[0];j<boundaries[1];++j){ 
               nlist.clear(); nlist.addAtom(i); nlist.addAtom(j); addColvar( nlist );
           }
       }
   } else {
      assert(false);
   }
   maxatoms=2;
}

double ColvarDistance::compute( const std::vector<unsigned>& indexes, std::vector<Vector>& derivatives, Tensor& virial ){
  assert( indexes.size()==2 && derivatives.size()==2 );
  Vector distance=getSeparation( indexes[0], indexes[1] ); 

  if ( component<0 ){
    const double value=distance.modulo();
    const double invvalue=1.0/value;
    derivatives[0]=-invvalue*distance;
    derivatives[1]=invvalue*distance;
    virial=-invvalue*Tensor(distance,distance);
    return value;
  } else {
    const double value=distance[component];
    const double invvalue=1.0/value;
    derivatives[0][component]=-1.0; 
    derivatives[1][component]=1.0; 
    virial=Tensor( distance,derivatives[0] );
    return value;
  }
}

}



