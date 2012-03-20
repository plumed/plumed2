#include "MultiColvar.h"
#include "PlumedMain.h"
#include "DistributionFunctions.h"
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
  keys.reserve("atoms","ATOMS","the atoms involved in each of the collective variables you wish to calculate. "
                               "Keywords like ATOMS1, ATOMS2, ATOMS3,... should be listed and one CV will be "
                               "calculated for each ATOM keyword you specify (all ATOM keywords should "
                               "define the same number of atoms).  The eventual number of quantities calculated by this "
                               "action will depend on what functions of the distribution you choose to calculate."); 
  keys.reserve("atoms","GROUP","this keyword is used for colvars that are calculated from a pair of atoms. "
                               "One colvar is calculated for each distinct pair of atoms in the group.");
  keys.reserve("atoms","GROUPA","this keyword is used for colvars that are calculated from a pair of atoms and must appaer with the keyword GROUPB. "
                                "Every pair of atoms which involves one atom from GROUPA and one atom from GROUPB defines one colvar");
  keys.reserve("atoms","GROUPB","this keyword is used for colvars that are calculate from a pair of atoms and must appaer with the keyword GROUPA. "
                                "Every pair of atoms which involves one atom from GROUPA and one atom from GROUPB defines one colvar");
  keys.reserve("atoms","SPECIES","this keyword is used for colvars such as coordination number. In that context it specifies that plumed should calculate "
                                 "one coordination number for each of the atoms specified.  Each of these coordination numbers specifies how many of the "
                                 "other specified atoms are within a certain cutoff of the central atom.");
  keys.reserve("atoms","SPECIESA","this keyword is used for colvars such as the coordination number.  In that context it species that plumed should calculate "
                                  "one coordination number for each of the atoms specified in SPECIESA.  Each of these cooordination numbers specifies how many "
                                  "of the atoms specifies using SPECIESB is within the specified cutoff");
  keys.reserve("atoms","SPECIESB","this keyword is used for colvars such as the coordination number.  It must appear with SPECIESA.  For a full explanation see " 
                                  "the documentation for that keyword");
  keys.add("optional","FIELD","create a field cv from these collective variables.");
  keys.add("optional","MIN","calculate the minimum value and store it in a value called min. " + min::documentation() );
//  keys.add("optional","MAX", "take the maximum value from these variables");
  keys.addFlag("AVERAGE",false,"take the average value of these variables and store it in value called average.");
  keys.add("optional","LESS_THAN", "take the number of variables less than the specified target and store it in a value called lt<target>. " + less_than::documentation() );
  keys.add("optional","MORE_THAN", "take the number of variables more than the specified target and store it in a value called gt<target>. " + more_than::documentation() ); 
  keys.add("optional","HISTOGRAM", "create a discretized histogram of the distribution of collective variables.  " + HistogramBead::histodocs(false) );
  keys.add("numbered", "WITHIN", "calculate the number variabels that are within a certain range and store it in a value called between<lowerbound>&<upperbound>. " + within::documentation() );
  keys.reserve("numbered", "SUBCELL", "calculate the average value of the CV within a portion of the box and store it in a value called subcell. " + cvdens::documentation() );
  ActionWithDistribution::registerKeywords( keys );
} 

MultiColvar::MultiColvar(const ActionOptions&ao):
Action(ao),
ActionAtomistic(ao),
ActionWithValue(ao),
ActionWithDistribution(ao),
readatoms(false),
setperiods(false),
needsCentralAtomPosition(false),
usepbc(true)
{
  if( keywords.style("NOPBC", "flag") ){ 
    bool nopbc=!usepbc; parseFlag("NOPBC",nopbc);
    usepbc=!nopbc; parseFlag("PBC",usepbc);
  }
}

void MultiColvar::readAtoms( int& natoms ){
  if( keywords.exists("ATOMS") ) readAtomsKeyword( natoms );
  if( keywords.exists("GROUP") ) readGroupsKeyword( natoms );
  if( keywords.exists("SPECIES") ) readSpeciesKeyword( natoms );

  if( !readatoms ) error("No atoms have been read in");

  // -- Now read in distribution keywords -- //
  bool dothis; std::string params;

  // Read AVERAGE keyword
  if( keywords.exists("AVERAGE") ) parseFlag("AVERAGE",dothis);
  else dothis=false;
  if( dothis ){
     Tools::convert(getNumberOfFunctionsInDistribution(),params);
     addDistributionFunction("AVERAGE", new mean(params) );
     params.clear();
  }
  // Read MIN keyword
  if( keywords.exists("MIN") ) parse("MIN",params);
  if( params.size()!=0 ){
     addDistributionFunction("MIN", new min(params) );
     params.clear();
  }
  // Read Less_THAN keyword
  if( keywords.exists("LESS_THAN") ) parse("LESS_THAN",params);
  if( params.size()!=0 ){
     addDistributionFunction("LESS_THAN", new less_than(params) );
     params.clear();
  }
  // Read MORE_THAN keyword
  if( keywords.exists("MORE_THAN") ) parse("MORE_THAN",params);
  if( params.size()!=0 ){
     addDistributionFunction("MORE_THAN", new more_than(params) );
     params.clear();
  }
  // Read HISTOGRAM keyword
  if( keywords.exists("HISTOGRAM") ) parse("HISTOGRAM",params);
  if( params.size()!=0 ){
       std::vector<std::string> bins;
       HistogramBead::generateBins( params, "", bins );
       for(unsigned i=0;i<bins.size();++i){
           addDistributionFunction( "HISTOGRAM", new within( bins[i] ) );
       }
       params.clear();
  }
  // Read within keywords
  if( keywords.exists("WITHIN") ){
     parse("WITHIN",params);
     if( params.size()!=0 ){
         addDistributionFunction( "WITHIN", new within(params) );
         params.clear();
     } else if( keywords.exists("WITHIN") ){
         for(unsigned i=1;;++i){
            if( !parseNumbered("WITHIN",i,params) ) break;
            std::string ss; Tools::convert(i,ss);
            addDistributionFunction( "WITHIN" + ss, new within(params) );
            params.clear();
         }
     }
  }

  // Establish whether or not this is a density ( requires special method )
  std::string dens; 
  if( keywords.exists("SPECIES") && !keywords.exists("SPECIESA") && !keywords.exists("SPECIESB") ){ dens="density"; }

  // Read CV_DENSITY keywords
  if( keywords.exists("SUBCELL") ){
      parse("SUBCELL",params); params=params + dens;
      if( params.size()!=0 ){
          needsCentralAtomPosition=true; 
          addDistributionFunction( "subcell", new cvdens(params) );
          params.clear();
      } else {
         for(unsigned i=1;;++i){
            if( !parseNumbered("SUBCELL",i,params) ) break;
            params=params+dens;
            needsCentralAtomPosition=true;
            std::string num; Tools::convert(i, num);
            addDistributionFunction( "subcell"+num, new cvdens(params) );
            params.clear();
         }
      }
  } 

  if( !usingDistributionFunctions() ) setupField(1);  // Setup the field if we are not using any of the above 
  requestDistribution();  // And setup the ActionWithDistribution
  requestAtoms();         // Request the atoms in ActionAtomistic and set up the value sizes
}

void MultiColvar::readAtomsKeyword( int& natoms ){ 
  if( readatoms) return; 

  std::vector<AtomNumber> t; DynamicList<unsigned> newlist;
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
        all_atoms.addIndexToList( t[j] );
     }
     t.resize(0); colvar_atoms.push_back( newlist );
     newlist.clear(); readatoms=true;
  }
}

void MultiColvar::readGroupsKeyword( int& natoms ){
  if( readatoms ) return;

  if( natoms==2 ){
      if( !keywords.exists("GROUPA") ) error("use GROUPA and GROUPB keywords as well as GROUP");
      if( !keywords.exists("GROUPB") ) error("use GROUPA and GROUPB keywords as well as GROUP");
  } else {
      error("Cannot use groups keyword unless the number of atoms equals 2");
  }
  
  std::vector<AtomNumber> t;
  parseAtomList("GROUP",t);
  if( t.size()!=0 ){
      readatoms=true;
      for(unsigned i=0;i<t.size();++i) all_atoms.addIndexToList( t[i] );
      DynamicList<unsigned> newlist; 
      for(unsigned i=1;i<t.size();++i){ 
          for(unsigned j=0;j<i;++j){ 
             newlist.addIndexToList(i); newlist.addIndexToList(j);
             colvar_atoms.push_back( newlist ); newlist.clear();
             log.printf("  Colvar %d is calculated from atoms : %d %d \n", colvar_atoms.size(), t[i].serial(), t[j].serial() ); 
          }
      }
  } else {
      std::vector<AtomNumber> t1,t2; 
      parseAtomList("GROUPA",t1);
      if( t1.size()!=0 ){
         readatoms=true;
         parseAtomList("GROUPB",t2);
         if ( t2.size()==0 ) error("GROUPB keyword defines no atoms or is missing. Use either GROUPA and GROUPB or just GROUP"); 
         for(unsigned i=0;i<t1.size();++i) all_atoms.addIndexToList( t1[i] ); 
         for(unsigned i=0;i<t2.size();++i) all_atoms.addIndexToList( t2[i] ); 
         DynamicList<unsigned> newlist;
         for(unsigned i=0;i<t1.size();++i){
             for(unsigned j=0;j<t2.size();++j){
                 newlist.addIndexToList(i); newlist.addIndexToList( t1.size() + j );
                 colvar_atoms.push_back( newlist ); newlist.clear();
                 log.printf("  Colvar %d is calculated from atoms : %d %d \n", colvar_atoms.size(), t1[i].serial(), t2[j].serial() );
             }
         }
      }
  }
}

void MultiColvar::readSpeciesKeyword( int& natoms ){
  if( readatoms ) return ;

//  if( !keywords.exists("SPECIESA") ) error("use SPECIESA and SPECIESB keywords as well as SPECIES");
//  if( !keywords.exists("SPECIESB") ) error("use SPECIESA and SPECIESB keywords as well as SPECIES");

  std::vector<AtomNumber> t;
  parseAtomList("SPECIES",t);
  if( t.size()!=0 ){
      readatoms=true; natoms=t.size();
      for(unsigned i=0;i<t.size();++i) all_atoms.addIndexToList( t[i] );
      DynamicList<unsigned> newlist;
      if( keywords.exists("SPECIESA") && keywords.exists("SPECIESB") ){
          for(unsigned i=0;i<t.size();++i){
              newlist.addIndexToList(i);
              log.printf("  Colvar %d involves central atom %d and atoms : ", colvar_atoms.size()+1,t[i].serial() );
              for(unsigned j=0;j<t.size();++j){
                  if(i!=j){ newlist.addIndexToList(j); log.printf("%d ",t[j].serial() ); }
              }
              log.printf("\n");
              colvar_atoms.push_back( newlist ); newlist.clear();
          }
      } else if( !( keywords.exists("SPECIESA") && keywords.exists("SPECIESB") ) ){
          DynamicList<unsigned> newlist;
          log.printf("  involving atoms : ");
          for(unsigned i=0;i<t.size();++i){ 
             newlist.addIndexToList(i); log.printf(" %d",t[i].serial() ); 
             colvar_atoms.push_back( newlist ); newlist.clear();
          }
          log.printf("\n");  
      } else {
          plumed_massert(0,"SPECIES keyword is not for density or coordination like CV");
      }
  } else if( keywords.exists("SPECIESA") && keywords.exists("SPECIESB") ) {
      std::vector<AtomNumber> t1,t2;
      parseAtomList("SPECIESA",t1);
      if( t1.size()!=0 ){
         readatoms=true; 
         parseAtomList("SPECIESB",t2);
         if ( t2.size()==0 ) error("SPECIESB keyword defines no atoms or is missing. Use either SPECIESA and SPECIESB or just SPECIES");
         natoms=1+t2.size();
         for(unsigned i=0;i<t1.size();++i) all_atoms.addIndexToList( t1[i] );
         for(unsigned i=0;i<t2.size();++i) all_atoms.addIndexToList( t2[i] );
         DynamicList<unsigned> newlist;
         for(unsigned i=0;i<t1.size();++i){
            newlist.addIndexToList(i);
            log.printf("  Colvar %d involves central atom %d and atoms : ", colvar_atoms.size()+1,t1[i].serial() );
            for(unsigned j=0;j<t2.size();++j){
                newlist.addIndexToList( t1.size() + j ); log.printf("%d ",t2[j].serial() ); 
            }
            log.printf("\n");
            colvar_atoms.push_back( newlist ); newlist.clear(); 
         }
      }
  } 
}

void MultiColvar::checkRead(){
  plumed_massert(setperiods, "You must set the periodicity of the various component functions");
  Action::checkRead();
}

void MultiColvar::setNotPeriodic(){
  setperiods=true;
  if( !usingDistributionFunctions() ){
      std::string num;
      for(unsigned i=0;i<getNumberOfComponents();++i){
          Tools::convert(i+1,num);
          componentIsNotPeriodic( "fval"+num );
      }
  }
}

void MultiColvar::setPeriodicDomain( const double& min, const double max ){
  setperiods=true;
  if( !usingDistributionFunctions() ){
     std::string num;
     for(unsigned i=0;i<getNumberOfComponents();++i){ 
         Tools::convert(i+1,num);
         componentIsPeriodic( "fval"+num , min, max );
     }
  }     
}


void MultiColvar::prepareForNeighborListUpdate(){
   for(unsigned i=0;i<colvar_atoms.size();++i){
      colvar_atoms[i].activateAll(); colvar_atoms[i].updateActiveMembers();
   }
   all_atoms.activateAll(); 
   all_atoms.updateActiveMembers();
   requestAtoms(); 
}

void MultiColvar::completeNeighborListUpdate(){
   for(unsigned i=0;i<colvar_atoms.size();++i){
      colvar_atoms[i].mpi_gatherActiveMembers( comm );
      activateLinks( colvar_atoms[i], all_atoms );
   }
   all_atoms.updateActiveMembers(); requestAtoms();
}

void MultiColvar::requestAtoms(){
   ActionAtomistic::requestAtoms( all_atoms.retrieveActiveList() );
   if( usingDistributionFunctions() ){
       for(unsigned i=0;i<getNumberOfComponents();++i) getPntrToComponent(i)->resizeDerivatives(3*getNumberOfAtoms()+9);
   } else {
       for(unsigned i=0;i<getNumberOfComponents();++i){ getPntrToComponent(i)->resizeDerivatives( getThisFunctionsNumberOfDerivatives(i) ); }
   }
}

void MultiColvar::getCentralAtom( const std::vector<Vector>& pos, std::vector<Value>& cpos){
   plumed_massert(0,"gradient and related cv distribution functions are not available in this colvar");
}

void MultiColvar::calculateThisFunction( const unsigned& j, Value* value_in, std::vector<Value>& aux ){
  unsigned natoms=colvar_atoms[j].getNumberActive();

  if ( natoms==0 ) return;   // Do nothing if there are no active atoms in the colvar
  // Retrieve the atoms
  Tensor vir; std::vector<Vector> pos(natoms), der(natoms); 
  for(unsigned i=0;i<natoms;++i){ pos[i]=getPosition( colvar_atoms[j][i] ); der.clear(); }
  
  // Compute the derivatives
  stopcondition=false; current=j;
  double value=compute( pos, der, vir );

  // This is the end of neighbor list update
  if(stopcondition){
     plumed_massert(isTimeForNeighborListUpdate(), "found stop but not during neighbor list step");
     return;
  }

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

  if(needsCentralAtomPosition){
     if( aux.size()!=3 ){ 
        aux.resize(3); 
        for(unsigned i=0;i<3;++i) aux[i].resizeDerivatives( value_in->getNumberOfDerivatives() );
     }
     getCentralAtom( pos, aux ); 
  }
}

void MultiColvar::mergeDerivatives( const unsigned j, Value* value_in, Value* value_out ){    

  int thisatom; unsigned innat=colvar_atoms[j].getNumberActive();
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

void MultiColvar::derivedFieldSetup( const double sigma ){
  const double pi=3.141592653589793238462643383279502884197169399375105820974944592307;
  log.printf("  generating field cv from histogram in which each component CV is represented by a Gaussian of width %f\n",sigma);
  fsigma2=2*sigma*sigma;
  fsigma4=fsigma2*fsigma2;
  fnorm = 1.0 / ( sqrt(2*pi)*sigma );
  for(unsigned i=0;i<getNumberOfFunctionsInDistribution();++i){
       std::string ss; Tools::convert(i+1,ss);
       addComponentWithDerivatives("fval"+ss); 
  }
}

void MultiColvar::setFieldOutputValue( const unsigned& j, Value* value_in ){
  plumed_assert( value_in->getNumberOfDerivatives()==getPntrToComponent(j)->getNumberOfDerivatives() );
  Value* value_out=getPntrToComponent(j);
  for(unsigned i=0;i<value_in->getNumberOfDerivatives();++i){ value_out->addDerivative( i, value_in->getDerivative(i) ); }
  // And store the value
  value_out->set( value_in->get() );
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

  unsigned vstart=3*nat-9; unsigned thisatom;
  for(int i=0;i<getNumberOfComponents();++i){

     if( !usingDistributionFunctions() ){
       nder=getThisFunctionsNumberOfDerivatives(i); 
       if (forces.size()!=nder ) forces.resize(nder);
       if( getPntrToComponent(i)->applyForce( forces ) ){
           for(unsigned j=0;j<forces.size();++j){
               thisatom=linkIndex( i, colvar_atoms[j], all_atoms );
               plumed_assert( thisatom>=0 );
               f[thisatom][0]+=forces[3*j+0];
               f[thisatom][1]+=forces[3*j+1];
               f[thisatom][2]+=forces[3*j+2];
           }
           vstart=nder-9;
       }
     } else if( getPntrToComponent(i)->applyForce( forces ) ){
        for(unsigned j=0;j<nat;++j){
           f[j][0]+=forces[3*j+0];
           f[j][1]+=forces[3*j+1];
           f[j][2]+=forces[3*j+2];
        }
     }
     v(0,0)+=forces[vstart+0];
     v(0,1)+=forces[vstart+1];
     v(0,2)+=forces[vstart+2];
     v(1,0)+=forces[vstart+3];
     v(1,1)+=forces[vstart+4];
     v(1,2)+=forces[vstart+5];
     v(2,0)+=forces[vstart+6];
     v(2,1)+=forces[vstart+7];
     v(2,2)+=forces[vstart+8];
  }
}





