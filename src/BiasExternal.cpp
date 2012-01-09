#include "Bias.h"
#include "ActionRegister.h"
#include "Grid.h"
#include "Tools.h"
#include "PlumedMain.h"
#include <cassert>


using namespace std;


namespace PLMD{

//+PLUMEDOC BIAS EXTERNAL 
/**
Restraint from external potential defined on a grid

\par Syntax
\verbatim
EXTERNAL ...
  ARG=x1,x2,... 
  FILE=filename
  [NOSPLINE]
  [SPARSE]
... EXTERNAL 
\endverbatim
FILE the name of the file where the external potential is stored, 
NOSPLINE to disable the use of splines, SPARSE to use a sparse grid.
\par Example
The following input is for a calculation with an external potential
defined in the file bias.dat and acting on the distance between atoms 3 and 5.
\verbatim
DISTANCE ATOMS=3,5 LABEL=d1
EXTERNAL ARG=d1 FILENAME=bias.dat LABEL=external 
\endverbatim
(See also \ref DISTANCE \ref PRINT).

*/
//+ENDPLUMEDOC

class BiasExternal : public Bias{

private:
  Grid* BiasGrid_;
  
public:
  BiasExternal(const ActionOptions&);
  ~BiasExternal();
  void calculate();

};

PLUMED_REGISTER_ACTION(BiasExternal,"EXTERNAL")

BiasExternal::~BiasExternal(){
  delete BiasGrid_;
}

BiasExternal::BiasExternal(const ActionOptions& ao):
PLUMED_BIAS_INIT(ao),
BiasGrid_(NULL)
{
  string filename;
  parse("FILE",filename);
  assert(filename.length()>0);
  bool sparsegrid=false;
  parseFlag("SPARSE",sparsegrid);
  bool nospline=false;
  parseFlag("NOSPLINE",nospline);
  bool spline=!nospline;

  checkRead();

  log.printf("  External potential from file %s\n",filename.c_str());
  if(spline){log.printf("  External potential uses spline interpolation\n");}
  if(sparsegrid){log.printf("  External potential uses sparse grid\n");}
  
  addValue("bias");

// read grid
  FILE* gridfile=fopen(filename.c_str(),"r");  
  BiasGrid_=Grid::create(gridfile,sparsegrid,spline,true);
  fclose(gridfile);
  assert(BiasGrid_->getDimension()==getNumberOfArguments());
  for(unsigned i=0;i<getNumberOfArguments();++i){
   assert(getArguments()[i]->isPeriodic()==BiasGrid_->getIsPeriodic()[i]);
  } 
}

void BiasExternal::calculate()
{
  unsigned ncv=getNumberOfArguments();
  vector<double> cv(ncv), der(ncv);

  for(unsigned i=0;i<ncv;++i){cv[i]=getArgument(i);}

  double ene=BiasGrid_->getValueAndDerivatives(cv,der);

  Value* value=getValue("bias");
  setValue(value,ene);

// set Forces 
  for(unsigned i=0;i<ncv;++i){
   const double f=-der[i];
   setOutputForces(i,f);
  }
}

}
