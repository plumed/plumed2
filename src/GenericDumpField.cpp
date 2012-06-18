#include "ActionPilot.h"
#include "ActionRegister.h"
#include "ActionWithDistribution.h"
#include "PlumedMain.h"
#include "ActionSet.h"
#include "Field.h"
#include "Grid.h"

namespace PLMD {

//+PLUMEDOC ANALYSIS DUMPFIELD
/**
Print out the instantaneous of a field that would be used to represent the instantaneous
state of the system in a method such as \ref FIELD_METAD 

\par Examples
The following input writes out the histogram of distances between the atoms in the 
specified group every 10 steps to file called field.1, field.2...
\verbatim
DISTANCES GROUP=1-7 FIELD=(MIN=0.5 MAX=3.0 NSPLINE=20 SIGMA=0.1) LABEL=f1
DUMPFIELD FIELD=f1 NGRID=200 STRIDE=50 FILE_BASE=field
\endverbatim
(See also \ref DISTANCES)

*/
//+ENDPLUMEDOC

class GenericDumpField :
public ActionPilot
{
private:
  int fnum;
  Field* myfield;
  std::vector<unsigned> ngrid;
  std::string fbase;
public:
  GenericDumpField(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
  void calculate(){};
  void apply(){};
  void update();
  bool retrieveForces( std::vector<double>& force ){ return false; }
};

PLUMED_REGISTER_ACTION(GenericDumpField,"DUMPFIELD")

void GenericDumpField::registerKeywords(Keywords& keys){
  Action::registerKeywords(keys);
  ActionPilot::registerKeywords(keys);
  keys.add("compulsory","FIELD","the field we are writing out");
  keys.add("compulsory","STRIDE","the frequency with which to output the field");
  keys.add("compulsory","FILE_BASE","field","the base name to use to output the field");
  keys.add("compulsory","NGRID","number of grid points to use in each direction");
}

GenericDumpField::GenericDumpField(const ActionOptions& ao):
Action(ao),
ActionPilot(ao)
{

  // Find the field we are using
  std::string ll; parse("FIELD",ll);
  ActionWithDistribution* field=plumed.getActionSet().selectWithLabel<ActionWithDistribution*>(ll);
  if(!field) error("cannot find action named " + ll);
  myfield=field->getField();
  if(!myfield) error("action " + ll + " calculates a colvar and not a field");
  addDependency(field);

  ngrid.resize( myfield->get_Ndx() ); 
  parseVector("NGRID",ngrid);

  parse("FILE_BASE",fbase);
  log.printf("  printing field on files named %s\n",fbase.c_str() );
  log.printf("  printing %d dimensional field %s on a grid with dimensions ",ngrid.size(), ll.c_str() );
  for(unsigned i=0;i<ngrid.size();++i) log.printf("%d ",ngrid[i]);
  log.printf("\n");
  fnum=1;
  checkRead();
}

void GenericDumpField::update(){
  std::vector<double> min, max; myfield->retrieveBoundaries( min, max );
  plumed_assert( min.size()==ngrid.size() && max.size()==ngrid.size() );
  std::vector<bool> pbc(min.size(), false ); 
  Grid gg( min, max, ngrid, pbc, false, false );

  // Check field is properly normalized
  double norm=0; std::vector<double> pp( myfield->get_Ndx() );
  for(unsigned i=0;i<gg.getSize();++i){
      gg.getPoint( i, pp ); gg.setValue( i, myfield->calculateField( pp ) );
      norm+=gg.getValue(i);
  } 
  norm*=gg.getBinVolume(); 

  double norm2=0.0;
  for(unsigned i=0;i<gg.getSize();++i){
      gg.getPoint( i, pp ); gg.setValue( i, gg.getValue(i)/norm );
      norm2+=gg.getValue(i);
  } 
  norm2*=gg.getBinVolume();

  if(comm.Get_rank()==0){
     std::string num, fname; 
     Tools::convert( fnum, num );
     fname = fbase + "." + num;
     FILE* gfile=fopen( fname.c_str(), "w" );
     fprintf(gfile,"#! DEFINITE INTEGRAL %f\n",norm2);
     gg.writeToFile( gfile );
     fclose(gfile);
  }
  fnum++;
}

}
