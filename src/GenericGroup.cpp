#include "ActionRegister.h"
#include "Group.h"

using namespace std;

namespace PLMD {

//+PLUMEDOC GROUP STATIC_GROUP
/**
Define a group of atoms

\par Example
The following contains a static group containing atoms 1-20.  Wherever the label
of the group appears after the GROUP keyword the specified list of atom will be used
to calculate the colvar.  
\verbatim
GROUP LABEL=label ATOMS=1-20
\endverbatim

*/
//+ENDPLUMEDOC

class GenericGroup : public Group {
public:
  GenericGroup(const ActionOptions&ao);
  virtual double compute( const std::vector<Vector>& positions, std::vector<double>& contributions, std::vector<Vector>& derivatives, Tensor& virial );
};

PLUMED_REGISTER_ACTION(GenericGroup,"STATIC_GROUP")

GenericGroup::GenericGroup(const ActionOptions&ao):
Group(ao)
{
  readGroup();
  checkRead();
}

double GenericGroup::compute( const std::vector<Vector>& positions, std::vector<double>& contributions, std::vector<Vector>& derivatives, Tensor& virial ){
  assert( positions.size()==contributions.size() && derivatives.size()==positions.size() );
  for(unsigned i=0;i<positions.size();++i){ contributions[i]=1.0; derivatives[i][0]=derivatives[i][1]=derivatives[i][2]=0.0; }
  virial.clear();
}

} 
