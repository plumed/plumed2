#include "ColvarDistinguishable.h"

using namespace PLMD;

ColvarDistinguishable::ColvarDistinguishable(const ActionOptions&ao) :
Colvar(ao)
{
  allowKeyword("ATOMS");
}

void ColvarDistinguishable::addIndexes( const unsigned& astart, const std::vector<unsigned>& new_indexes ){
  if( function_indexes.size()==0 ) derivatives.resize( new_indexes.size() );
  else if( new_indexes.size()!=function_indexes[0].size() ) error("mismatch for number of atoms in colvar");

  unsigned accum=astart; std::vector<unsigned> tmplist; tmplist.push_back( accum ); 
  log.printf("  a colvar will be calculated from the positions of the following set of atoms : ( %s",plumed.getAtoms().interpretIndex( new_indexes[0] ).c_str() );
  for(unsigned j=1;j<new_indexes.size();++j){
      tmplist.push_back( accum + j );
      log.printf( ", %s",plumed.getAtoms().interpretIndex( new_indexes[j] ).c_str() );
  }
  log.printf(" ) \n");
  function_indexes.push_back( tmplist );
}

void ColvarDistinguishable::interpretGroupsKeyword( const unsigned& natoms, const std::string& atomGroupName, const std::vector<std::vector<unsigned> >& groups ){
  std::vector<unsigned> tmplist; 
  if( groups.size()>2 ) error("you can only use groups for indistinguishable colvars if the number of atoms in each colvar is equal to 2");

  derivatives.resize( natoms );
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
      for(unsigned i=1;i<groups[0].size();++i){
          for(unsigned j=0;j<i;++j){ 
              tmplist.resize(0); tmplist.push_back(j); 
              tmplist.push_back(i); function_indexes.push_back( tmplist ); 
          } 
      }
  } else if( groups.size()==2 ){
      for(unsigned i=0;i<groups[0].size();++i){
          for(unsigned j=0;j<groups[1].size();++j){
              tmplist.resize(0); tmplist.push_back(i); 
              tmplist.push_back(groups[0].size()+j); function_indexes.push_back( tmplist ); 
          }
      }  
  }
}

void ColvarDistinguishable::interpretAtomsKeyword( const std::vector<std::vector<unsigned> >& flist ){
  unsigned accum=0; std::vector<unsigned> tmplist;
  for(unsigned i=0;i<flist.size();++i){ 
    addIndexes( accum, flist[i] ); accum+=flist[i].size(); 
  }
}

void ColvarDistinguishable::updateDynamicContent( const double& cutoff, std::vector<bool>& skips ){
  bool calcfunc; unsigned n=0;

  std::vector<unsigned> tmpskip( function_indexes.size() ), tmpblocks( comm.Get_size() + 1);
  for(unsigned i=0;i<function_indexes.size();++i){ tmpskip[i]=i+1; }
  if( isParallel() ) {
     comm.splitList( tmpskip, tmpblocks );
  } else {
     tmpblocks[0]=0;
     for(unsigned i=1;i<tmpblocks.size();++i){ tmpblocks[i]=tmpskip.size(); }
  }
  for(unsigned i=0;i<function_indexes.size();++i){ tmpskip[i]=0; }
  
  unsigned rank=comm.Get_rank();
  //std::vector<unsigned> tmpskip( function_indexes.size() );
  for(unsigned i=tmpblocks[rank];i<tmpblocks[rank+1];++i){
      calcfunc=true;
      for(unsigned j=1;j<function_indexes[i].size();++j){
         for(unsigned k=0;k<j;++k){
            if( skips[ function_indexes[i][j] ] || skips[ function_indexes[i][k] ] ){
                calcfunc=false;
            } else if( getSeparation( function_indexes[i][j], function_indexes[i][k] ).modulo()>cutoff ){
                calcfunc=false;
            }
         }
      } 
      if( calcfunc ) { tmpskip[n]=i; n=i; }  
  } 

  // PARALLEL Do an allgather on tmpskip

  setSkips( tmpskip );
    
  std::vector<bool> required_atoms(skips.size(),false);
  for(unsigned i=0;i<function_indexes.size();i=tmpskip[i]){ 
     for(unsigned n=0;n<function_indexes[i].size();++n) required_atoms[ function_indexes[i][n] ] = true;
  }
  for(unsigned i=0;i<skips.size();++i){ if( !required_atoms[i] ) skips[i]=true; }
  updateParallelLoops(); // And finally update the parallelized loops
}


