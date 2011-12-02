#include "ColvarWithModifiers.h"
#include "ColvarModifier.h"
#include "ColvarModifierFunctions.h"

namespace PLMD {

ColvarWithModifiers::ColvarWithModifiers(const ActionOptions& ao) :
Colvar(ao),
readInput(false)
{
  forbidKeyword("STRIDE");
  registerKeyword(0, "MIN", "modifier");
  registerKeyword(0, "MAX", "modifier");
  registerKeyword(0, "SUM", "modifier");
  registerKeyword(0, "AVERAGE", "modifier");
  registerKeyword(0, "LESS_THAN", "modifier");
  registerKeyword(0, "MORE_THAN", "modifier");
  registerKeyword(0, "HISTOGRAM", "modifier");
}

void ColvarWithModifiers::readAtomsKeyword( int& natoms ){
  if( testForKey("ATOMS") ){
      if(readInput) error("cannot mix input keywords");
      maxatoms=natoms; 

      if( testForNumberedKeys("ATOMS") ){
         std::string nn, report; AtomicNeighbourList nlist( dynamic_cast<ActionAtomistic*>(this) );
         for(int i=1;;++i ){
            Tools::convert(i, nn); report= nn + "th cv calculated from atoms";
            if( !parseAtomList( "ATOMS", i ) ) break;
            if( i==1 && maxatoms<0 ){
                maxatoms=getNumberOfAtoms();
            } else if ( (getNumberOfAtoms() - (i-1)*maxatoms)!=maxatoms ){
                std::string nn, ss, tt; Tools::convert( i, nn ); 
                Tools::convert( static_cast<int>( getNumberOfAtoms() - (i-1)*maxatoms ), ss ); Tools::convert( maxatoms, tt );
                error("in ATOMS" + nn + " keyword found " + ss + " atoms when there should only be " + tt + " atoms");
            }
            nlist.clear(); unsigned nat=getNumberOfAtoms() - maxatoms;
            for(unsigned j=0;j<maxatoms;++j) nlist.addAtom( nat+j ); 
            addColvar( nlist );
        }
      } else {
         AtomicNeighbourList nlist( dynamic_cast<ActionAtomistic*>(this) );
         parseAtomList("ATOMS", 0 );
         for(unsigned j=1;j<maxatoms;++j){
             for(unsigned k=0;k<j;++k){ nlist.addPair( j, k ); }
         }   
         addColvar( nlist );
         if( maxatoms<0 ){
            maxatoms=getNumberOfAtoms();
         } else if ( getNumberOfAtoms()!=maxatoms ) {
            std::string ss, tt; 
            Tools::convert( getNumberOfAtoms(), ss ); Tools::convert( maxatoms, tt );
            error("in ATOMS keyword found " + ss + " atoms when there should only be " + tt + " atoms");
         }
      }
      readInput=true; derivatives.resize( maxatoms );
  }
}

void ColvarWithModifiers::readGroupsKeyword( int& maxgroups ){
  if( testForKey("GROUP") ){ 
     if(readInput) error("cannot mix input keywords");

     std::vector<unsigned> boundaries;
     if( testForNumberedKeys("GROUP") ){
         std::string nn, report;
         for(int i=1;i<=maxgroups;++i ){
            Tools::convert(i, nn); report= nn + "th group contains atoms";
            if( !parseAtomList( "GROUP", i ) ) break;
            boundaries.push_back( getNumberOfAtoms() );
         }
     } else {
         parseAtomList("GROUP", 0 );  
         boundaries.push_back( getNumberOfAtoms() );
     }
     unsigned maxatoms=0;
     interpretGroupsKeyword( boundaries, maxatoms );
     if( maxatoms==0 ) error("Your interpretGroupsKeyword didn't specify the number of atoms");
     derivatives.resize( maxatoms ); 
     readInput=true;
  }
}

// Note you can overwrite this routine by writing a new one inside your class
void ColvarWithModifiers::interpretGroupsKeyword( const std::vector<unsigned>& boundaries, unsigned& maxatoms ){
  if( !checkNeighbourListType(2) ) error("the neighbour list is of the wrong type to use interpretGroupsKeyword function");

  unsigned nsphere;
  AtomicNeighbourList nlist( dynamic_cast<ActionAtomistic*>(this) );
  if( boundaries.size()==1 ){
      assert( boundaries[0]==getNumberOfAtoms() );
      for(unsigned i=0;i<getNumberOfAtoms();++i){
          nlist.clear(); nlist.addAtom(i);
          for(unsigned j=0;j<getNumberOfAtoms();++j){
              if(i!=j) nlist.addAtom(j);
          }
          addColvar( nlist );
      }
      nsphere=boundaries[0]-1;
      maxatoms=getNumberOfAtoms();
  } else if ( boundaries.size()==2 ){
      assert( boundaries[1]==getNumberOfAtoms() );
      for(unsigned i=0;i<boundaries[0];++i){
         nlist.clear(); nlist.addAtom(i);
         for(unsigned j=boundaries[0];j<getNumberOfAtoms();++j) nlist.addAtom(j);
         addColvar( nlist );
      }
      nsphere=boundaries[1]-boundaries[0];
      maxatoms=nsphere + 1;
  } else {
      assert(false);
  }

  // Now create the neighbour list
  std::vector< std::pair<unsigned,unsigned> > pair_list;
  for(unsigned i=1;i<=nsphere;++i){
      pair_list.push_back( std::pair<unsigned,unsigned>( 0, i ) );
  }
  setupNeighbourList( pair_list );
}

void ColvarWithModifiers::finishColvarSetup( const double& min, const double& max ){
  assert(readInput);  // Check input has been read in
  
  readActionColvar(); 
  std::vector<double> domain(2); domain[0]=min; domain[1]=max;
  readActionWithExternalArguments( 3*getNumberOfAtoms()+9, domain ); 

  // Make it so that we can get number of colvars for modifiers
  int ncol=getNumberOfColvars(); 
  std::string sncol; Tools::convert(ncol,sncol);
  mod_params.push_back( std::pair<std::string,std::string>("NCOLVAR",sncol) );

  bool doall=true, dothis; 
  // First all the colvar modifiers that must be done alone
  dothis=false; parseFlag("MIN",dothis);
  if (dothis){
     doall=false;
     ColvarModifier* mod=new ColvarModifierMin(this);
     modifiers.push_back( mod );
  }
//  dothis=false; parseFlag("MAX",dothis);
//  if (dothis){
//     doall=false;
//     log.printf("  using maximum value \n");
//     abort();    /// Must implement max value at some stage GAT
//  }
  dothis=false; parseFlag("SUM",dothis);
  if( dothis ){
     doall=false;
     ColvarModifier* mod=new ColvarModifierSum(this);
     modifiers.push_back( mod );
  }
  dothis=false; parseFlag("AVERAGE",dothis);
  if( dothis ){
     doall=false;
     ColvarModifier* mod=new ColvarModifierMean(this);
     modifiers.push_back( mod );
  }
  std::string r_0="none"; parse("LESS_THAN",r_0);
  if( r_0!="none" ){
     doall=false; 
     mod_params.push_back( std::pair<std::string,std::string>("LESS_THAN",r_0) );
     ColvarModifier* mod=new ColvarModifierLess(this);
     modifiers.push_back( mod );
  }
  r_0="none"; parse("MORE_THAN",r_0);
  if( r_0!="none" ){
     doall=false;
     mod_params.push_back( std::pair<std::string,std::string>("MORE_THAN",r_0) );
     ColvarModifier* mod=new ColvarModifierMore(this);
     modifiers.push_back( mod );
  }
  std::string hist_input="none"; 
  dothis=false; parseFlag("HISTOGRAM",dothis); 
  if (dothis){
      doall=false;
      ColvarModifier* mod=new ColvarModifierHistogram(this);
      modifiers.push_back( mod );
  }

  // And lastly the default - do everything
  if( doall ){
     // This checks we are not doing things with neighbour lists that are forbidden
     checkForBadNeighbourLists();
     // This adds all the values
     std::string n;
     for(unsigned i=0;i<getNumberOfColvars();++i){
        Tools::convert(i,n); addValue("value" + n, false, true );
     }
  }
  mod_params.resize(0);  // Get rid of this stuff as it is only there for input
  checkRead();
}

void ColvarWithModifiers::calculate() {

   Tensor virial; unsigned nat=getNumberOfAtoms(); 
   double value; std::vector<unsigned> indexes;
   for(unsigned i=0;i<getNumberOfColvars();++i){
       if( !getColvarAtoms( i , indexes ) ) error("something has gone wrong with neighbour lists"); 
       assert( indexes.size()==derivatives.size() );    
       value=compute( indexes, derivatives, virial );
       if( modifiers.size()==0 ){
           for(unsigned j=0;j<indexes.size();++j){
               addDerivative( i, 3*indexes[j] + 0, derivatives[j][0] );
               addDerivative( i, 3*indexes[j] + 1, derivatives[j][1] );
               addDerivative( i, 3*indexes[j] + 2, derivatives[j][2] );
           } 
           addDerivative( i, 3*nat + 0, virial(0,0) );
           addDerivative( i, 3*nat + 1, virial(0,1) );
           addDerivative( i, 3*nat + 2, virial(0,2) );
           addDerivative( i, 3*nat + 3, virial(1,0) );
           addDerivative( i, 3*nat + 4, virial(1,1) );
           addDerivative( i, 3*nat + 5, virial(1,2) );
           addDerivative( i, 3*nat + 6, virial(2,0) );
           addDerivative( i, 3*nat + 7, virial(2,1) );
           addDerivative( i, 3*nat + 8, virial(2,2) );
           setValue( i, value, 1.0 );
       } else {
           for(unsigned j=0;j<modifiers.size();++j) modifiers[j]->mergeDerivatives( indexes, value, derivatives, virial );
       }
   }
   for(unsigned j=0;j<modifiers.size();++j) modifiers[j]->finishCalculation();    
}

}
