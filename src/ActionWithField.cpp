#include "ActionWithField.h"
#include "PlumedMain.h"
#include "ActionSet.h"

using namespace std;
using namespace PLMD;

void ActionWithField::registerKeywords(Keywords& keys){
  Action::registerKeywords(keys);
  ActionWithValue::registerKeywords(keys);
  ActionPilot::registerKeywords(keys);
  keys.add("compulsory","STRIDE","1","the frequency with which to output the field");
  keys.add("compulsory","FIELD","the input for this action is the field calculated during one of the other actions.");
  keys.add("compulsory","NGRID","number of grid points to use in each direction");
  keys.addFlag("SERIAL", false, "do the calculation in serial");
}

ActionWithField::ActionWithField(const ActionOptions&ao):
  Action(ao),
  ActionWithValue(ao),
  ActionPilot(ao),
  serial(false)
{
  if( checkNumericalDerivatives() ) warning("cannot parallelize field with numerical derivatives over actions");
  parseFlag("SERIAL",serial);
  if(serial) log.printf("  doing calculation in serial\n");

  // Find the field we are using
  std::string ll; parse("FIELD",ll);
  ActionWithDistribution* field=plumed.getActionSet().selectWithLabel<ActionWithDistribution*>(ll);
  addDependency(field);
  if(!field) error("cannot find action named " + ll);
  myfield=field->getField();
  if(!myfield) error("action " + ll + " calculates a colvar and not a field");
  apply_action=dynamic_cast<Action*>( field );

  // Read how we descritize the integrals
  std::vector<unsigned> ngrid( myfield->get_Ndx() ); parseVector("NGRID",ngrid);
  log.printf("  using field %s and descritizing %d dimensional integrals over ",ll.c_str(), ngrid.size() );
  for(unsigned i=0;i<ngrid.size();++i) log.printf("%d ",ngrid[i]);
  log.printf("points\n");

  // Create the grid where we store the bias
  std::vector<double> min, max; 
  myfield->retrieveBoundaries( min, max );
  plumed_assert( min.size()==ngrid.size() && max.size()==ngrid.size() );
  std::vector<bool> pbc(min.size(), false );
  bias=new Grid( min, max, ngrid, pbc, false, false );
  buffer.resize( bias->getSize() + 2 );

  // Setup the blocks for parallelizing the grid calculations
  unsigned stride=comm.Get_size();
  if(serial) stride=1;

  blocks.resize( stride+1 );
  unsigned nn=std::floor( bias->getSize() / stride );
  unsigned nrem=bias->getSize() - nn*stride;

  blocks[0]=0;
  for(unsigned i=1;i<blocks.size();++i){
      for(unsigned j=0;j<i;++j) blocks[i]+=blocks[j];
      if( i<=nrem ) blocks[i]+=nn + 1;
      else blocks[i]+=nn;
  }
  plumed_assert( blocks[blocks.size()-1]==bias->getSize() );

  // Add something for the bias and set up the forces
  addComponentWithDerivatives("bias"); 
  getPntrToComponent("bias")->resizeDerivatives( 1 );
}

void ActionWithField::clearBias(){
  std::vector<unsigned> ngrid; ngrid=bias->getNbin();
  delete bias;
  std::vector<double> min, max; myfield->retrieveBoundaries( min, max );
  plumed_assert( min.size()==ngrid.size() && max.size()==ngrid.size() );
  std::vector<bool> pbc(min.size(), false );
  bias=new Grid( min, max, ngrid, pbc, false, false );
}

void ActionWithField::calculate(){
  unsigned rank=comm.Get_rank();
  if(serial) rank=0;

  unsigned nder=myfield->get_NdX();
  if( derivatives.size()!=nder ){
      derivatives.resize( nder );
      getPntrToComponent("bias")->resizeDerivatives( nder );
  }

  // This calculates the current field if we are doing numerical derivatives
  if( checkNumericalDerivatives() ) apply_action->calculate(); 

  // The loop for the bias
  buffer.assign( buffer.size(), 0.0 );
  std::vector<double> pp( myfield->get_Ndx() );
  for(unsigned i=blocks[rank];i<blocks[rank+1];++i){
      bias->getPoint( i, pp );
      buffer[i+2]=myfield->calculateField( pp );
      buffer[0]+=buffer[i+2];
      buffer[1]+=buffer[i+2]*bias->getValue( i );
  }
  buffer[0]*=bias->getBinVolume();   
  buffer[1]*=bias->getBinVolume() / buffer[0]; 
  if(!serial) comm.Sum( &buffer[0], buffer.size() );
  getPntrToComponent("bias")->set(buffer[1]);

  if( checkNumericalDerivatives() ) return ;
 
  // The loop for the derivatives 
  derivatives.assign( derivatives.size(), 0.0 );
  std::vector<double> tmpforce( derivatives.size() );
  for(unsigned i=blocks[rank];i<blocks[rank+1];++i){
      bias->getPoint( i, pp );
      myfield->calculateFieldDerivatives( pp, tmpforce ); 
      for(unsigned j=0;j<derivatives.size();++j){
          derivatives[j] += tmpforce[j] * ( bias->getValue( i ) - buffer[1] ); 
      }
  }
  for(unsigned j=0;j<derivatives.size();++j) derivatives[j]*=bias->getBinVolume() / buffer[0];
  if(!serial) comm.Sum( &derivatives[0], derivatives.size() );
  Value* bb=getPntrToComponent("bias");
  for(unsigned j=0;j<derivatives.size();++j) bb->addDerivative( j, derivatives[j] );
}

void ActionWithField::calculateNumericalDerivatives( ActionWithValue* a ){
  apply_action->calculateNumericalDerivatives( this );
}

void ActionWithField::addFieldToBias( const double& hh ){
  for(unsigned i=0;i<bias->getSize();++i){
      bias->addValue( i, hh*buffer[i+2]/buffer[0] );
  }
}

void ActionWithField::apply(){
  if(onStep()){
     for(unsigned j=0;j<derivatives.size();++j) derivatives[j]*=-1.0*getStride();
     myfield->addForces( derivatives ); 
  }
}

ActionWithField::~ActionWithField(){
  delete bias; 
}


