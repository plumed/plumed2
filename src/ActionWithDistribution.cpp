#include "ActionWithDistribution.h"

using namespace std;
using namespace PLMD;

void ActionWithDistribution::registerKeywords(Keywords& keys){
  keys.add("optional","NL_STRIDE","the frequency with which the neighbor list should be updated.");
  keys.add("optional","NL_TOL","when accumulating sums quantities that contribute less than this will be ignored.");
}

void ActionWithDistribution::autoParallelize(Keywords& keys){
  keys.addFlag("SERIAL",false,"do the calculation in serial.  Do not parallelize over collective variables");
}

ActionWithDistribution::ActionWithDistribution(const ActionOptions&ao):
  Action(ao),
  read(false),
  all_values(true),
  use_field(false),
  serial(true),
  updateFreq(0),
  lastUpdate(0),
  reduceAtNextStep(false),
  tolerance(0)
{
  if( keywords.exists("SERIAL") ) parseFlag("SERIAL",serial);
  if(serial)log.printf("  doing calculation in serial\n");
  if( keywords.exists("NL_STRIDE") ) parse("NL_STRIDE",updateFreq);
  if( keywords.exists("NL_TOL") ){
      if(updateFreq>0 && keywords.exists("NL_TOL") ){ 
         tolerance=epsilon; parse("NL_TOL",tolerance); 
         log.printf("  updating calculation every %d steps.  Ignoring contributions less than %lf\n",updateFreq,tolerance);
      } else if( updateFreq>0 ){ 
         log.printf("  updating calculation every %d steps.\n");  
      }
  }
}

ActionWithDistribution::~ActionWithDistribution(){
  for(unsigned i=0;i<functions.size();++i) delete functions[i];
}

void ActionWithDistribution::addDistributionFunction( std::string name, DistributionFunction* fun ){
  if(all_values) all_values=false;  // Possibly will add functionality to delete all values here

  // Check function is good 
  if( !fun->check() ){
     log.printf("ERROR for keyword %s in action %s with label %s : %s \n \n",name.c_str(), getName().c_str(), getLabel().c_str(), ( fun->errorMessage() ).c_str() );
     fun->printKeywords( log ); 
     plumed_merror("ERROR for keyword " + name + " in action "  + getName() + " with label " + getLabel() + " : " + fun->errorMessage() );
     exit(1);
  }

  // Add a value
  ActionWithValue*a=dynamic_cast<ActionWithValue*>(this);
  a->addComponentWithDerivatives( fun->getLabel() );
  unsigned fno=a->getNumberOfComponents()-1;
  final_values.push_back( a->copyOutput( fno ) );

  // Add the function   
  plumed_massert( fno==functions.size(), "Number of functions does not match number of values" );
  functions.push_back( fun );
  log.printf("  value %s contains %s\n",( (a->copyOutput( fno ))->getName() ).c_str(),( functions[fno]->message() ).c_str() );
}

void ActionWithDistribution::requestDistribution(){
  read=true; 
  ActionWithValue*a=dynamic_cast<ActionWithValue*>(this);
  plumed_massert(a,"can only do distribution on ActionsWithValue");

  if(all_values){
      error("No function has been specified");
  } else if( !use_field ){
      plumed_massert( functions.size()==final_values.size(), "number of functions does not match number of values" );
  }
  // This sets up the dynamic list that holds what we are calculating
  for(unsigned i=0;i<getNumberOfFunctionsInDistribution();++i){ members.addIndexToList(i); }
  members.activateAll(); members.updateActiveMembers();
  // We prepare the first step as if we are doing a neighbor list update
  prepareForNeighborListUpdate(); reduceAtNextStep=true;
}

void ActionWithDistribution::setupField( unsigned ldim ){
  plumed_massert( keywords.exists("FIELD"), "you have chosen to use fields but have not written documentation for the FIELD keyword"); 
  plumed_massert( ldim<=2 , "fields don't work with more than two dimensions" );
  
  std::string fieldin;
  parse("FIELD",fieldin);
  if( fieldin.size()==0 ) return ;
  all_values=false; use_field=true;
  std::vector<std::string> data=Tools::getWords(fieldin);

  std::vector<unsigned> nspline; std::vector<double> min,max; double sigma;
  bool found_s=Tools::parseVector( data,"NSPLINE", nspline ); 
  if(!found_s) field_error("did not find NPLINE keyword");
  if(nspline.size()!=ldim){
     std::string ll,ww; 
     Tools::convert(ldim,ll);
     Tools::convert(nspline.size(),ww);
     field_error("found " + ww + " values for NSPLINE when expecting only " + ll);
  }
  bool found_min=Tools::parseVector( data, "MIN", min );
  if(!found_min) field_error("did not find MIN keyword");
  if(min.size()!=ldim){ 
     std::string ll,ww; 
     Tools::convert(ldim,ll);
     Tools::convert(min.size(),ww);
     field_error("found " + ww + " values for MIN when expecting only " + ll);
  }
  bool found_max=Tools::parseVector( data, "MAX", max );
  if(!found_max) field_error("did not find MAX keyword");
  if(max.size()!=ldim){
     std::string ll,ww;
     Tools::convert(ldim,ll);
     Tools::convert(max.size(),ww);
     field_error("found " + ww + " values for MAX when expecting only " + ll);
  }
  bool found_sigma=Tools::parse( data, "SIGMA", sigma );
  if(!found_sigma) field_error("did not find SIGMA keyword");

  if( data.size()!=0 ){
      std::string err="found the following rogue keywords : ";
      for(unsigned i=0;i<data.size();++i) err=err + data[i] + ", ";
      field_error(err); 
  }
  // Setup the field stuff in derived class
  derivedFieldSetup( sigma );

  // Setup the field
  unsigned np=1; for(unsigned i=0;i<nspline.size();++i) np*=nspline[i]; 
  ActionWithValue*a=dynamic_cast<ActionWithValue*>(this);
  plumed_massert(a,"can only do fields on ActionsWithValue");
  myfield.setup( getNumberOfFunctionsInDistribution(), np, a->getNumberOfComponents(), ldim );
}

void ActionWithDistribution::prepare(){
 if(reduceAtNextStep){
    completeNeighborListUpdate();
    // Setup the functions by declaring enough space to hold the derivatives
    for(unsigned i=0;i<functions.size();++i) functions[i]->setNumberOfDerivatives( final_values[0]->getNumberOfDerivatives() );
    // Setup the buffers for mpi gather
    if(!serial){
       if( use_field ){
         std::vector<unsigned> cv_sizes( getNumberOfFunctionsInDistribution() );
         for(unsigned i=0;i<getNumberOfFunctionsInDistribution();++i){ cv_sizes[i]=getThisFunctionsNumberOfDerivatives(i); }
         myfield.resizeBaseQuantityBuffers( cv_sizes ); 
       } else {
         unsigned bufsize=0;
         for(unsigned i=0;i<functions.size();++i) bufsize+=functions[i]->requiredBufferSpace();
         buffer.resize( bufsize ); 
       }
    }
    reduceAtNextStep=false;
 }
 if( updateFreq>0 && (getStep()-lastUpdate)>=updateFreq ){
    members.activateAll();
    members.updateActiveMembers();
    prepareForNeighborListUpdate();
    // Setup the functions by declaring enough space to hold the derivatives
    for(unsigned i=0;i<functions.size();++i) functions[i]->setNumberOfDerivatives( final_values[0]->getNumberOfDerivatives() );
    // Setup the buffers for mpi gather
    if(!serial){
       if( use_field ){
         std::vector<unsigned> cv_sizes( getNumberOfFunctionsInDistribution() );
         for(unsigned i=0;i<getNumberOfFunctionsInDistribution();++i){ cv_sizes[i]=getThisFunctionsNumberOfDerivatives(i); }
         myfield.resizeBaseQuantityBuffers( cv_sizes ); 
       } else {
         unsigned bufsize=0;
         for(unsigned i=0;i<functions.size();++i) bufsize+=functions[i]->requiredBufferSpace();
         buffer.resize( bufsize ); 
       }
    }
    reduceAtNextStep=true;
    lastUpdate=getStep();
 }
}

void ActionWithDistribution::calculate(){
  plumed_massert( read, "you must have a call to requestDistribution somewhere" );

  if( use_field ) calculateField();
  else calculateFunctions();
}

void ActionWithDistribution::calculateField(){
  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  if(serial){ stride=1; rank=0; }

  // Set everything in the field to zero
  myfield.clear();
  unsigned kk; std::vector<Value> aux; Value* tmpvalue=new Value();
  // A loop over the functions in the distribution
  for(unsigned i=rank;i<members.getNumberActive();i+=stride){
      // Retrieve the function we are calculating from the dynamic list
      kk=members[i];
      // Make sure we have enough derivatives in this value
      unsigned nder=getThisFunctionsNumberOfDerivatives(kk);
      if( tmpvalue->getNumberOfDerivatives()!=nder ){ tmpvalue->resizeDerivatives(nder); }
      
      // Calculate the value of this particular function 
      calculateThisFunction( kk, tmpvalue, aux );
      // Transfer the value to the field buffers
      myfield.setBaseQuantity( kk, tmpvalue ); 

      // Skip if we are not calculating this particular value
      if( reduceAtNextStep && !tmpvalue->valueHasBeenSet() ){ members.deactivate(kk); deactivate(kk); }

      // Reset everything ready for next step
      tmpvalue->clearDerivatives();
  }
  if(!serial){ myfield.gatherBaseQuantities( comm ); }
  // Set the output values for this quantity (we use these to chain rule)
  for(unsigned i=0;i<myfield.get_NdX();++i){
     unsigned nder=getThisFunctionsNumberOfDerivatives(i);
     if( tmpvalue->getNumberOfDerivatives()!=nder ){ tmpvalue->resizeDerivatives(nder); }
     myfield.extractBaseQuantity( i, tmpvalue );
     setFieldOutputValue( i, tmpvalue );
  }
  delete tmpvalue;
}

void ActionWithDistribution::calculateFunctions(){
  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  if(serial){ stride=1; rank=0; }

  std::vector<Value> aux;
  // Reset everything
  for(unsigned j=0;j<functions.size();++j) functions[j]->reset();
  // Create a value to store stuff in 
  Value* tmpvalue=new Value();

  unsigned kk; double contribution; bool keep;
  for(unsigned i=rank;i<members.getNumberActive();i+=stride){
      // Retrieve the function we are calculating from the dynamic list
      kk=members[i];
      // Make sure we have enough derivatives in this value
      unsigned nder=getThisFunctionsNumberOfDerivatives(kk);
      if( tmpvalue->getNumberOfDerivatives()!=nder ) tmpvalue->resizeDerivatives( nder );
      // Retrieve the periodicity of this value
      if( isPeriodic(kk) ){ 
         double min, max; retrieveDomain( kk, min, max );
         tmpvalue->setDomain( min, max ); 
      } else { tmpvalue->setNotPeriodic(); }
 
      // Calculate the value of this particular function 
      calculateThisFunction( kk, tmpvalue, aux );

      // Skip if we are not calculating this particular value
      if( reduceAtNextStep && !tmpvalue->valueHasBeenSet() ){ members.deactivate(kk); deactivate(kk); continue; }

      // Now incorporate the derivative of the function into the derivatives for the min etc
      if( updateFreq>0 ){ keep=false; } else { keep=true; }
      for(unsigned j=0;j<functions.size();++j){
         functions[j]->calculate( tmpvalue, aux );
         if( functions[j]->sizableContribution( tolerance ) ){ 
             keep=true; functions[j]->mergeDerivatives( kk, *this );
         }
      }
      tmpvalue->clearDerivatives();
      for(unsigned i=0;i<aux.size();++i) aux[i].clearDerivatives();
      // If the contribution of this quantity is very small at neighbour list time ignore it
      // untill next neighbour list time
      if( reduceAtNextStep && !keep ){ members.deactivate(kk); deactivate(kk); } 
  }
  // Update the dynamic list 
  if(reduceAtNextStep){ members.mpi_gatherActiveMembers( comm ); }
  // MPI Gather everything
  if(!serial){ 
     unsigned bufsize=0;
     for(unsigned i=0;i<functions.size();++i) functions[i]->copyDataToBuffers( bufsize, buffer );
     plumed_assert( bufsize==buffer.size() ); 
     comm.Sum( &buffer[0],buffer.size() ); 
     bufsize=0;
     for(unsigned i=0;i<functions.size();++i) functions[i]->retrieveDataFromBuffers( bufsize, buffer );
     plumed_assert( bufsize==buffer.size() );
  }

  // Delete the tmpvalues
  delete tmpvalue; 
  // Set the final value of the function
  for(unsigned j=0;j<final_values.size();++j) functions[j]->finish( final_values[j] ); 
}

void ActionWithDistribution::retrieveDomain( const unsigned nn, double& min, double& max ){
  plumed_massert(0, "If your function is periodic you need to add a retrieveDomain function so that ActionWithDistribution can retrieve the domain");
}

void ActionWithDistribution::field_error( const std::string msg ){
  Keywords field_keys;
  field_keys.add("compulsory","NPSPLINE","the number of points in each direction at which to calculate the value of the field");
  field_keys.add("compulsory","MIN","the minimum value at which the value of the field should be calculated in each direction");
  field_keys.add("compulsory","MAX","the maximum value at which the value of the field should be calculated in each direction");
  field_keys.add("compulsory","SIGMA","the value of the sigma parameter in the field");

  log.printf("ERROR for keyword FIELD in action %s with label %s : %s \n \n",getName().c_str(), getLabel().c_str(), msg.c_str() );
  field_keys.print( log );
  plumed_merror("ERROR for keyword FIELD in action "  + getName() + " with label " + getLabel() + " : " + msg );
  exit(1); 
}

void ActionWithDistribution::FieldClass::setup( const unsigned nfunc, const unsigned np, const unsigned D, const unsigned d ){
// Set everything for base quantities 
  baseq_nder.resize(nfunc); baseq_starts.resize(nfunc); 
// Set everything for grid
  npoints=np; ndX=D; ndx=d; nper=(ndX+1)*(ndx+1); grid_buffer.resize( npoints*nper );
}

void ActionWithDistribution::FieldClass::clear(){
  baseq_buffer.assign( baseq_buffer.size(), 0.0 );
  grid_buffer.assign( grid_buffer.size(), 0.0 );
}

void ActionWithDistribution::FieldClass::resizeBaseQuantityBuffers( const std::vector<unsigned>& cv_sizes ){
  plumed_assert( cv_sizes.size()==baseq_nder.size() && cv_sizes.size()==baseq_starts.size() );
  unsigned nn=0;
  for(unsigned i=0;i<cv_sizes.size();++i){
      baseq_starts[i]=nn; baseq_nder[i]=cv_sizes[i]; nn+=cv_sizes[i]+1;
  }
  baseq_buffer.resize(nn); // And resize the buffers
}

void ActionWithDistribution::FieldClass::setBaseQuantity( const unsigned nn, Value* val ){
  plumed_assert( nn<baseq_nder.size() ); 
  plumed_assert( val->getNumberOfDerivatives()==baseq_nder[nn] );

  unsigned kk=baseq_starts[nn];
  baseq_buffer[kk]=val->get(); kk++;
  for(unsigned i=0;i<val->getNumberOfDerivatives();++i){ baseq_buffer[kk]=val->getDerivative(i); kk++; }
  if( (nn+1)==baseq_starts.size() ){ plumed_assert( kk==baseq_buffer.size() ); }
  else{ plumed_assert( kk==baseq_starts[nn+1] ); }
}

void ActionWithDistribution::FieldClass::gatherBaseQuantities( PlumedCommunicator& comm ){
  comm.Sum( &baseq_buffer[0],baseq_buffer.size() );
}

void ActionWithDistribution::FieldClass::extractBaseQuantity( const unsigned nn, Value* val ){
  plumed_assert( nn<baseq_nder.size() ); 
  if( baseq_nder[nn]!=val->getNumberOfDerivatives() ) val->resizeDerivatives(baseq_nder[nn]); 

  val->clearDerivatives();
  unsigned kk=baseq_starts[nn];
  val->set( baseq_buffer[kk] ); kk++;
  for(unsigned i=0;i<val->getNumberOfDerivatives();++i){ val->addDerivative( i, baseq_buffer[kk] ); kk++; }
  if( (nn+1)==baseq_starts.size() ){ plumed_assert( kk==baseq_buffer.size() ); }
  else{ plumed_assert( kk==baseq_starts[nn+1] ); }
}
