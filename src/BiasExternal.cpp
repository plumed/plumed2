#include "Bias.h"
#include "ActionRegister.h"
#include "Grid.h"
#include "PlumedMain.h"
#include "PlumedException.h"


using namespace std;


namespace PLMD{

//+PLUMEDOC BIAS EXTERNAL 
/**
Calculate a restraint that is defined on a grid that is read in at start up

\par Examples
The following is an input for a calculation with an external potential that is
defined in the file bias.dat and that acts on the distance between atoms 3 and 5.
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
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(BiasExternal,"EXTERNAL")

void BiasExternal::registerKeywords(Keywords& keys){
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","FILE","the name of the file containing the external potential");
  keys.addFlag("NOSPLINE",false,"specifies that no spline interpolation is to be used when calculating the energy and forces due to the external potential");
  keys.addFlag("SPARSE",false,"specifies that the external potential uses a sparse grid");
}

BiasExternal::~BiasExternal(){
  delete BiasGrid_;
}

BiasExternal::BiasExternal(const ActionOptions& ao):
PLUMED_BIAS_INIT(ao),
BiasGrid_(NULL)
{
  string filename;
  parse("FILE",filename);
  plumed_assert(filename.length()>0);
  bool sparsegrid=false;
  parseFlag("SPARSE",sparsegrid);
  bool nospline=false;
  parseFlag("NOSPLINE",nospline);
  bool spline=!nospline;

  checkRead();

  log.printf("  External potential from file %s\n",filename.c_str());
  if(spline){log.printf("  External potential uses spline interpolation\n");}
  if(sparsegrid){log.printf("  External potential uses sparse grid\n");}
  
  addComponent("bias");

// read grid
  FILE* gridfile=fopen(filename.c_str(),"r");  
  BiasGrid_=Grid::create(gridfile,sparsegrid,spline,true);
  fclose(gridfile);
  if(BiasGrid_->getDimension()!=getNumberOfArguments()) error("mismatch between dimensionality of input grid and number of arguments");
  for(unsigned i=0;i<getNumberOfArguments();++i){
    if( getPntrToArgument(i)->isPeriodic()!=BiasGrid_->getIsPeriodic()[i] ) error("periodicity mismatch between arguments and input bias"); 
  } 
}

void BiasExternal::calculate()
{
  unsigned ncv=getNumberOfArguments();
  vector<double> cv(ncv), der(ncv);

  for(unsigned i=0;i<ncv;++i){cv[i]=getArgument(i);}

  double ene=BiasGrid_->getValueAndDerivatives(cv,der);

  getPntrToComponent("bias")->set(ene);

// set Forces 
  for(unsigned i=0;i<ncv;++i){
   const double f=-der[i];
   setOutputForces(i,f);
  }
}

}
