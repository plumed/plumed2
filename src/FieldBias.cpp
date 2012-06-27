#include "FieldBias.h"
#include "PlumedMain.h"
#include "ActionSet.h"

using namespace std;
using namespace PLMD;

void FieldBias::registerKeywords(Keywords& keys){
  Action::registerKeywords(keys);
  ActionWithValue::registerKeywords(keys);
  ActionPilot::registerKeywords(keys);
  keys.add("compulsory","STRIDE","1","the frequency with which to output the field");
  keys.add("compulsory","FIELD","the input for this action is the field calculated during one of the other actions.");
  keys.add("compulsory","NORM","the normalization to use");
  keys.add("optional","NGRID","number of grid points to use in each direction - if you are using function interpolation");
  keys.add("optional","START_BIAS","the bias at the start of the simulation");
  keys.addFlag("SERIAL", false, "do the calculation in serial");
  keys.addFlag("DEBUG_DERIVATIVES",false,"used to debug the derivatives of the bias");
}

FieldBias::FieldBias(const ActionOptions&ao):
  Action(ao),
  ActionWithValue(ao),
  ActionPilot(ao),
  myfield(NULL),
  serial(false),
  debug(false)
{
  if( checkNumericalDerivatives() ) warning("cannot parallelize field with numerical derivatives over actions");
  parseFlag("SERIAL",serial);
  if(serial) log.printf("  doing calculation in serial\n");
  parseFlag("DEBUG_DERIVATIVES",debug);
  if(debug) log.printf("  storing derivatives in a value so they can be output with dump derivatives");

  // Find the field we are using
  std::string ll; parse("FIELD",ll);
  ActionWithDistribution* field=plumed.getActionSet().selectWithLabel<ActionWithDistribution*>(ll);
  if(field){
     myfield=field->getField();
     if(!myfield) error("input action " + ll + " does not calcualte a field");
     apply_action=dynamic_cast<ActionWithValue*>( field );
  } else {
     apply_action=plumed.getActionSet().selectWithLabel<ActionWithValue*>(ll);
     apply_action->checkFieldsAllowed();
     for(unsigned i=0;i<apply_action->getNumberOfComponents();++i) f_arg.push_back(apply_action->copyOutput(i));
  }
  addDependency(apply_action);

  // Read how we descritize the integrals
  std::vector<unsigned> ngrid; 
  if( f_arg.size()==0 ) ngrid.resize( myfield->get_Ndx() ); 
  parseVector("NGRID",ngrid);
  if( !myfield ){ ngrid.resize( 1 ); ngrid[0]=f_arg.size()-1; }   // myfield->get_nspline( ngrid ); }
  else if( ngrid.size()!=myfield->get_Ndx() ) error("grid size is wrong");
  log.printf("  using field %s and descritizing %d dimensional integrals over ",ll.c_str(), ngrid.size() );
  for(unsigned i=0;i<ngrid.size();++i) log.printf("%d ",ngrid[i]);
  // Work out the norm we are using
  unsigned nn; parse("NORM",nn); norm=static_cast<double>(nn);
  log.printf("points. Normalizing with %d norm\n",nn);

  // Create the grid where we store the bias
  std::vector<double> min, max;
  if( myfield ) myfield->retrieveBoundaries( min, max );
  else{ min.resize(1); max.resize(1); min[0]=0; max[0]=static_cast<double>( f_arg.size()-1 ); }
  plumed_assert( min.size()==ngrid.size() && max.size()==ngrid.size() );

  std::string sbias; parse("START_BIAS",sbias);
  if( sbias.length()==0 ){
    std::vector<bool> pbc(min.size(), false );
    bias=new Grid( min, max, ngrid, pbc, false, false );
  } else {
    log.printf("  reading initial bias from file named %s\n",sbias.c_str());
    FILE* bfile=fopen( sbias.c_str(),"r");
    bias=Grid::create( bfile, false, false, false );
    fclose( bfile );
    for(unsigned i=0;i<ngrid.size();++i){
       if( bias->getMin()[i]!=min[i] ) error("minimum in input grid does not match minimum of field");
       if( ( bias->getMax()[i]- (max[i] + bias->getDx()[i]) )>0.00001 ) error("maximum in input grid does not match maximum of field");
       if( bias->getNbin()[i]!=(ngrid[i]+1) ) error("number of bins in input grid does not match plumed input");
       if( bias->getIsPeriodic()[i]==true ) error("periodic input grid is not compatible with field overlap");
    }
  }
  // Prepare the buffers for the calculation
  buffer.resize( bias->getSize() + 2 );

  // Setup the blocks for parallelizing the grid calculations
  unsigned stride=comm.Get_size();
  if(serial) stride=1;

  blocks.resize( stride+1 );
  nn=std::floor( bias->getSize() / stride );
  unsigned nrem=bias->getSize() - nn*stride;

  blocks[0]=0;
  for(unsigned i=1;i<blocks.size();++i){
      blocks[i]=blocks[i-1];
      if( blocks[i]<bias->getSize() ){
         if( i<=nrem ) blocks[i]+=nn + 1;
         else blocks[i]+=nn;
      } else {
         warning("number of integration points is less than number of nodes so some nodes are idle : consider running in serial");
      }
  }
  if(!serial){
     for(unsigned i=0;i<stride;++i){
         log.printf("  node %d is doing integration of %d points\n",i,blocks[i+1]-blocks[i]);
     }
  }
  plumed_assert( blocks[blocks.size()-1]==bias->getSize() );

  // Add something for the bias and set up the forces
  addComponentWithDerivatives("bias"); 
  if( debug ) getPntrToComponent("bias")->resizeDerivatives( 1 );
}

void FieldBias::clearBias(){
  std::vector<unsigned> ngrid; ngrid=bias->getNbin();
  delete bias;
  std::vector<double> min, max; 
  if( myfield ) myfield->retrieveBoundaries( min, max );
  else{ min.resize(1); max.resize(1); min[0]=0; max[0]=static_cast<double>( f_arg.size()-1 ); }
  plumed_assert( min.size()==ngrid.size() && max.size()==ngrid.size() );
  std::vector<bool> pbc(min.size(), false );
  bias=new Grid( min, max, ngrid, pbc, false, false );
}

void FieldBias::calculate(){
  unsigned rank=comm.Get_rank();
  if(serial) rank=0;

  unsigned nder;
  if( myfield ){ nder=myfield->get_NdX(); }
  else{ nder=f_arg.size(); }
  if( derivatives.size()!=nder ){
      derivatives.resize( nder ); 
      if( debug ) getPntrToComponent("bias")->resizeDerivatives( apply_action->getNumberOfDerivatives() ); 
  }
  // This calculates the current field if we are doing numerical derivatives
  if( checkNumericalDerivatives() ) apply_action->calculate(); 

  // The loop for the bias
  buffer.assign( buffer.size(), 0.0 );
  unsigned nldim;
  if( myfield ) nldim=myfield->get_Ndx(); else nldim=1;
  std::vector<double> pp( nldim );
  for(unsigned i=blocks[rank];i<blocks[rank+1];++i){
      bias->getPoint( i, pp );
      if( myfield ) buffer[i+2]=myfield->calculateField( pp );
      else buffer[i+2]=f_arg[i]->get();
      buffer[0]+=pow(buffer[i+2], norm);
      buffer[1]+=buffer[i+2]*bias->getValue( i );
  }
  buffer[0]*=bias->getBinVolume();   
  if(!serial) comm.Sum( &buffer[0], buffer.size() );
  double normali=pow( buffer[0], 1./static_cast<double>(norm) );
  buffer[1]*=bias->getBinVolume() / normali;
  getPntrToComponent("bias")->set( buffer[1]  );

  if( checkNumericalDerivatives() ) return ;

  // The loop for the derivatives 
  derivatives.assign( derivatives.size(), 0.0 );
  std::vector<double> tmpforce( derivatives.size() );
  for(unsigned i=blocks[rank];i<blocks[rank+1];++i){
      bias->getPoint( i, pp );
      if( myfield ) myfield->calculateFieldDerivatives( pp, tmpforce ); 
      else { 
        for(unsigned j=0;j<tmpforce.size();++j){ tmpforce[j]=0.0; }
        tmpforce[i]=1.0;
      }
      for(unsigned j=0;j<derivatives.size();++j){
          derivatives[j] += tmpforce[j] * ( bias->getValue( i ) / normali - (buffer[1]/buffer[0])*pow( buffer[i+2], norm-1) ); 
      }
  }
  for(unsigned j=0;j<derivatives.size();++j) derivatives[j]*=bias->getBinVolume();
  if(!serial) comm.Sum( &derivatives[0], derivatives.size() );
  if( debug ) apply_action->mergeFieldDerivatives( derivatives, getPntrToComponent("bias") );
}

void FieldBias::calculateNumericalDerivatives( ActionWithValue* a ){
  apply_action->calculateNumericalDerivatives( this );
}

void FieldBias::apply(){
  if(onStep()){
     for(unsigned j=0;j<derivatives.size();++j) derivatives[j]*=-1.0*getStride();
     if( myfield ){
        myfield->addForces( derivatives ); 
     } else {
        for(unsigned i=0;i<f_arg.size();++i) f_arg[i]->addForce( derivatives[i] );
     }
  }
}

FieldBias::~FieldBias(){
  delete bias; 
}


