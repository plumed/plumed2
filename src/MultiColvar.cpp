#include "MultiColvar.h"
#include "PlumedMain.h"
#include <vector>
#include <string>

using namespace std;
using namespace PLMD;

void MultiColvar::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  keys.addFlag("PBC",true,"use the periodic boundary conditions when calculating distances");
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  keys.reserve("atoms","ATOMS","the atoms involved in each of the collective variables you wish to calculate.  To compute a single CV use ATOMS.  If you use ATOMS1, ATOMS2, ATOMS3... multiple CVs will be calculated - one for each ATOM keyword you specify (all ATOM keywords should define the same number of atoms).  The eventual number of quantities calculated by this action will depend on what functions of the distribution you choose to calculate."); 
  ActionWithDistribution::registerKeywords( keys );
} 

MultiColvar::MultiColvar(const ActionOptions&ao):
Action(ao),
ActionAtomistic(ao),
ActionWithValue(ao),
ActionWithDistribution(ao),
readatoms(false),
usepbc(true)
{
  if( keywords.style("NOPBC", "flag") ){ 
    bool nopbc=!usepbc; parseFlag("NOPBC",nopbc);
    usepbc=!nopbc; parseFlag("PBC",usepbc);
  }
}

void MultiColvar::readAtoms( int& natoms ){
  if( keywords.exists("ATOMS") ) readAtomsKeyword( natoms );

  if( !readatoms ) error("No atoms have been read in");

  // Now set up all the lists for the first time
  for(unsigned i=0;i<colvar_atoms.size();++i){
     colvar_atoms[i].activateAll();
     colvar_atoms[i].updateActiveMembers();
     activateLinks( colvar_atoms[i], all_atoms );
  }
  all_atoms.updateActiveMembers();

  readDistributionKeywords();  // And read the ActionWithDistributionKeywords
  requestAtoms();              // Request the atoms in ActionAtomistic and set up the value sizes
}

void MultiColvar::readAtomsKeyword( int& natoms ){ 
  if( readatoms) return; 

  std::vector<AtomNumber> t;
  parseAtomList("ATOMS",t); 
  if( t.size()!=0 ){
     readatoms=true;
     if( natoms>0 && t.size()!=natoms ){
        std::string nat; Tools::convert(natoms, nat );
        error("ATOMS keyword did not specify " + nat  + " atoms.");
     } else {
        natoms=t.size();
     }
     DynamicList newlist;
     for(unsigned i=0;i<natoms;++i){ 
        newlist.addIndexToList(i);
        all_atoms.addIndexToList( t[i].index() ); 
     }
     colvar_atoms.push_back( newlist );
     log.printf("  Colvar 1 is calculated from atoms : ");
     for(unsigned i=0;i<t.size();++i) log.printf("%d ",t[i].serial() );
     log.printf("\n");
  } else {
     bool readone=false; DynamicList newlist;
     for(int i=1;;++i ){
        parseAtomList("ATOMS", i, t );
        if( t.size()==0 ) break;

        log.printf("  Colvar %d is calculated from atoms : ", i);
        for(unsigned i=0;i<t.size();++i) log.printf("%d ",t[i].serial() );
        log.printf("\n"); 

        if( i==1 && natoms<0 ) natoms=t.size();
        if( t.size()!=natoms ){
            std::string ss; Tools::convert(i,ss); 
            error("ATOMS" + ss + " keyword has the wrong number of atoms"); 
        }
        for(unsigned j=0;j<natoms;++j){ 
           newlist.addIndexToList( natoms*(i-1)+j ); 
           all_atoms.addIndexToList( t[j].index() );
        }
        t.resize(0); colvar_atoms.push_back( newlist );
        newlist.clear(); readatoms=true;
     }
  }
}

void MultiColvar::requestAtoms(){
   unsigned natoms=all_atoms.getNumberActive();
   std::vector<AtomNumber> a(natoms); 
   for(unsigned i=0;i<natoms;++i) a[i].setIndex( all_atoms(i) );

//   printf("Requesting atoms : ");
//   for(unsigned i=0;i<natoms;++i) printf(" %d",a[i].serial() );
//   printf("\n");

   ActionAtomistic::requestAtoms(a);
   if( usingDistributionFunctions() ){
       for(unsigned i=0;i<getNumberOfComponents();++i) getPntrToComponent(i)->resizeDerivatives(3*a.size()+9);
   } else {
       for(unsigned i=0;i<getNumberOfComponents();++i) getPntrToComponent(i)->resizeDerivatives( getThisFunctionsNumberOfDerivatives(i) );
   }
}

// This is the preparation for neighbour list update
//void MultiColvar::prepare(){
//  if( updateFreq>0 && (getStep()-lastUpdate)>=updateFreq ){
//     for(unsigned i=0;i<index_translator.size();++i) index_translator[i]=i;
//     // update the atoms we require
//     requestAtoms();
//  }
//}
//
//void MultiColvar::startUpdateStep(){
//  if( updateFreq>0 && (getStep()-lastUpdate)>=updateFreq ){
//      for(unsigned i=0;i<ColvarAtoms.size();++i) ColvarAtoms[i].resetCurrentAtoms(); 
//  }
//}
//
//void MultiColvar::completeUpdateStep(){
//  if( updateFreq>0 && (getStep()-lastUpdate)>=updateFreq ){
//      // Currently we get no atoms
//      for(unsigned i=0;i<index_translator.size();++i) index_translator[i]=-1; 
//
//      // Update the neighbour list
//      std::vector<Vector> pos(natoms); 
//      unsigned a1, a2; double d; bool skip;
//      for(unsigned i=0;i<ColvarAtoms.size();++i){
//          ColvarAtoms[i].getAtoms( pos ); skip=false;
//          for(unsigned j=0;j<ColvarAtoms[i].getNpairs();++j){
//              a1=pairs[j].first; a2=pairs[j].second;
//              d=getSeparation( pos[a1], pos[a2] );
//              if( sum && d<=nl_cut ){ 
//                  index_translator[ ColvarAtoms[i].getAtomNumber(a1) ]=1;
//                  index_translator[ ColvarAtoms[i].getAtomNumber(a2) ]=1;
//              } else if( !sum && d>nl_cut ){
//                  skip=true; break;
//              }
//          }
//          // Collect atoms if we are not skipping
//          if( !skip ){
//              unsigned natoms=ColvarAtoms[i].getNumberOfAtoms();
//              for(unsigned j=0;j<natoms;++j){
//                  a1=ColvarAtoms[i].getAtomNumber(j);
//                  index_translator[a1]=1;
//              }
//          } 
//      }
//
//      // And store the last neighour list update time
//      lastUpdate=getStep();
//      // update the atoms array and so on
//      requestAtoms();
//  }
//}

void MultiColvar::retrieveAtoms( const unsigned& j, std::vector<Vector>& pos ){
  unsigned natoms=colvar_atoms[j].getNumberActive();
  plumed_assert( pos.size()==natoms );
  for(unsigned i=0;i<natoms;++i){
      pos[i]=getPosition( colvar_atoms[j](i) );
  }
}

void MultiColvar::calculateThisFunction( const unsigned& j, Value* value_in ){
  Tensor vir; unsigned natoms=colvar_atoms[j].getNumberActive();
  std::vector<Vector> pos(natoms), der(natoms); 
  // Retrieve the atoms
  retrieveAtoms( j, pos );  
  // Compute the derivatives
  double value=compute( pos, der, vir );
  // Put all this in the value we are passing back
  for(unsigned i=0;i<natoms;++i){
      value_in->addDerivative( 3*i+0,der[i][0] );
      value_in->addDerivative( 3*i+1,der[i][1] );
      value_in->addDerivative( 3*i+2,der[i][2] );
  }
  value_in->addDerivative( 3*natoms+0, vir(0,0) );
  value_in->addDerivative( 3*natoms+1, vir(0,1) );
  value_in->addDerivative( 3*natoms+2, vir(0,2) );
  value_in->addDerivative( 3*natoms+3, vir(1,0) );
  value_in->addDerivative( 3*natoms+4, vir(1,1) );
  value_in->addDerivative( 3*natoms+5, vir(1,2) );
  value_in->addDerivative( 3*natoms+6, vir(2,0) );
  value_in->addDerivative( 3*natoms+7, vir(2,1) );
  value_in->addDerivative( 3*natoms+8, vir(2,2) );

  // And store the value
  value_in->set(value);
}

void MultiColvar::mergeDerivatives( const unsigned j, Value* value_in, Value* value_out ){    

  unsigned thisatom; unsigned innat=colvar_atoms[j].getNumberActive();
  for(unsigned i=0;i<innat;++i){
     thisatom=linkIndex( i, colvar_atoms[j], all_atoms );
     plumed_assert( thisatom>=0 ); 
     value_out->addDerivative( 3*thisatom+0, value_in->getDerivative(3*i+0) );
     value_out->addDerivative( 3*thisatom+1, value_in->getDerivative(3*i+1) );
     value_out->addDerivative( 3*thisatom+2, value_in->getDerivative(3*i+2) ); 
  }

  // Easy to merge the virial
  unsigned outnat=getNumberOfAtoms(); 
  value_out->addDerivative( 3*outnat+0, value_in->getDerivative(3*innat+0) );
  value_out->addDerivative( 3*outnat+1, value_in->getDerivative(3*innat+1) );
  value_out->addDerivative( 3*outnat+2, value_in->getDerivative(3*innat+2) );
  value_out->addDerivative( 3*outnat+3, value_in->getDerivative(3*innat+3) );
  value_out->addDerivative( 3*outnat+4, value_in->getDerivative(3*innat+4) );
  value_out->addDerivative( 3*outnat+5, value_in->getDerivative(3*innat+5) );
  value_out->addDerivative( 3*outnat+6, value_in->getDerivative(3*innat+6) );
  value_out->addDerivative( 3*outnat+7, value_in->getDerivative(3*innat+7) );
  value_out->addDerivative( 3*outnat+8, value_in->getDerivative(3*innat+8) );
}

Vector MultiColvar::getSeparation( const Vector& vec1, const Vector& vec2 ) const {
  if(usepbc){ return pbcDistance( vec1, vec2 ); }
  else{ return delta( vec1, vec2 ); }
}

void MultiColvar::apply(){
  vector<Vector>&   f(modifyForces());
  Tensor&           v(modifyVirial());

  for(unsigned i=0;i<f.size();i++){
    f[i][0]=0.0;
    f[i][1]=0.0;
    f[i][2]=0.0;
  }
  v.clear();

  unsigned nat=getNumberOfAtoms(); std::vector<double> forces; unsigned nder;
  if( usingDistributionFunctions() ) forces.resize(3*getNumberOfAtoms()+9);

  for(int i=0;i<getNumberOfComponents();++i){
     nder=getThisFunctionsNumberOfDerivatives(i);
     if( !usingDistributionFunctions() && forces.size()!=nder ) forces.resize(nder);
 
     if( getPntrToComponent(i)->applyForce( forces ) ){
        for(unsigned j=0;j<nat;++j){
           f[j][0]+=forces[3*j+0];
           f[j][1]+=forces[3*j+1];
           f[j][2]+=forces[3*j+2];
        }
        v(0,0)+=forces[3*nat+0];
        v(0,1)+=forces[3*nat+1];
        v(0,2)+=forces[3*nat+2];
        v(1,0)+=forces[3*nat+3];
        v(1,1)+=forces[3*nat+4];
        v(1,2)+=forces[3*nat+5];
        v(2,0)+=forces[3*nat+6];
        v(2,1)+=forces[3*nat+7];
        v(2,2)+=forces[3*nat+8];
     }
  }
}





