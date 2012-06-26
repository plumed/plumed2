#include "TargetDist.h"

namespace PLMD {

void TargetDist::read( const PDB& pdb, std::vector<Value*> ar ){
  // Clear values in target actions
  for(unsigned i=0;i<ar.size();++i){
     (ar[i]->getPntrToAction())->clearInputForces();
     (ar[i]->getPntrToAction())->clearDerivatives();
  }

  // Caclulate target actions from input in PDB file
  std::vector<double> targ( ar.size() );
  for(unsigned i=0;i<ar.size();++i){
      if( ar[i]->valueHasBeenSet() ){ 
         targ[i]=ar[i]->get();
      } else {
         ActionWithValue* vv=ar[i]->getPntrToAction();
         (ar[i]->getPntrToAction())->calculateFromPDB( pdb );
         targ[i]=ar[i]->get();
      }
  }
  read( targ, ar );
}

void TargetDist::read( const std::vector<double>& targ, std::vector<Value*> ar ){
  plumed_assert( targ.size()==ar.size() );

  target.resize( ar.size() ); args.resize( ar.size() );
  log.printf("  distance from this point in cv space : ");
  for(unsigned i=0;i<target.size();++i){ log.printf("%f ", targ[i]); target[i]=targ[i]; args[i]=ar[i]; }
  log.printf("\n");
}

double TargetDist::calculate( std::vector<double>& derivs ){
  plumed_assert( derivs.size()==args.size() );
  double dist=0, tmp;
  for(unsigned i=0;i<args.size();++i){
      tmp=args[i]->difference( target[i], args[i]->get() );
      derivs[i]=tmp; dist+=tmp*tmp; 
  }
  dist=sqrt(dist);
  for(unsigned i=0;i<args.size();++i) derivs[i]/=dist;
  return dist;
}

}
