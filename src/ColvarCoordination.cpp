#include "ColvarCoordination.h"

namespace PLMD{ 

ColvarCoordination::ColvarCoordination(const ActionOptions&ao):
Colvar(ao)
{
  allowKeyword("GROUP"); isCSphereF=true;
}

void ColvarCoordination::interpretGroupsKeyword( const unsigned& natoms, const std::string& atomGroupName, const std::vector<std::vector<unsigned> >& groups ){
  assert( natoms==0 );

  if( atomGroupName!=getLabel() ){
      log.printf("  using atoms specified in group %s\n", atomGroupName.c_str() );
  } else {
      for(unsigned i=0;i<groups.size();++i){
          log.printf("  atoms in group %d : ", i+1 );
          for(unsigned j=0;j<groups[i].size();++j) log.printf("%s ", plumed.getAtoms().interpretIndex( groups[i][j] ).c_str() );
          log.printf("\n"); 
      }
  }
 
  if( groups.size()==1 ){
     unsigned k;
     central.resize( groups[0].size() ); sphere.resize( groups[0].size() );
     for(unsigned i=0;i<groups[0].size();++i){
        k=0; central[i]=i; sphere[i].resize( groups[0].size()-1 );
        for(unsigned j=0;j<groups[0].size();++j){
            if(i!=j){ sphere[i].index[k]=k; k++; } 
        }
     }
     neighbours.resize( groups[0].size() -1 ); derivatives.resize( groups[0].size() );
  } else if ( groups.size()==2 ){
     central.resize( groups[0].size() ); sphere.resize( groups[0].size() );
     for(unsigned i=0;i<groups[0].size();++i){
        central[i]=i; sphere[i].resize( groups[1].size() );
        for(unsigned j=0;j<groups[1].size();++j){
            if(i!=j){ sphere[i].index[j]=groups[0].size()+j; }
        }
     }
     neighbours.resize( groups[1].size() ); derivatives.resize( 1+groups[1].size() );
  } else {
     assert(false);
  }
  for(unsigned i=0;i<sphere.size();++i){
     for(unsigned j=0;j<sphere[i].skipto.size();++j)  sphere[i].skipto[j]=j+1; 
  }
} 

void ColvarCoordination::interpretAtomsKeyword( const std::vector<std::vector<unsigned> >& flist ){ assert(false); }

void ColvarCoordination::updateNeighbourList( const double& cutoff, std::vector<bool>& skips ){
  unsigned n=0;
  
  std::vector<unsigned> tmpskip( central.size() );
  for(unsigned i=0;i<central.size();++i){
     if( skips[ central[i] ] ){ skipAllColvarFrom(n,i); tmpskip[n]=i; n=i; }
  }

  std::vector<bool> required_atoms(skips.size(),false); unsigned kk,maxvecs=0;
  for(unsigned i=0;i<central.size();i=tmpskip[i]){
     required_atoms[ central[i] ] = required_atoms[ sphere[i].index[0] ] = true; kk=0;
     for(unsigned n=0;n<sphere[i].skipto.size();++n){
        if( getSeparation( central[i], sphere[i].index[n] ).modulo()<=cutoff ){ 
            required_atoms[ sphere[i].index[n] ] = true; sphere[i].skipto[kk]=n; kk=n; 
        } 
     }
  } 
  for(unsigned i=0;i<skips.size();++i){ if( !required_atoms[i] ) skips[i]=true; }
}

}
