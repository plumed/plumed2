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
  keys.add("nohtml","MIN","take the minimum value from these variables.  The continuous version of the minimum described above is calculated and beta must be specified in input");
//  keys.add("optional","MAX", "take the maximum value from these variables");
  keys.reserve("nohtml","SUM", "take the sum of these variables");
  keys.add("nohtml","AVERAGE", "take the average value of these variables");
  keys.add("nohtml","LESS_THAN", "take the number of variables less than the specified target and store it in a value called lt<target>. " + SwitchingFunction::documentation() );
  keys.add("nohtml","MORE_THAN", "take the number of variables more than the specified target and store it in a value called gt<target>. " + SwitchingFunction::documentation() ); 
  keys.add("nohtml","HISTOGRAM", "create a discretized histogram of the distribution.  This keyword's input should consist of one or two numbers");
  keys.add("nohtml","RANGE", "the range in which to calculate the histogram");
  keys.add("nohtml", "WITHIN", "The number of values within a certain range.  This keyword's input should consist of one or two numbers.");
  keys.reserve("nohtml","CV_DENSITY_X","Takes 2/3 params specifying the upper and lower boundaries of a bin in x in fractional coordinates and the ammount of smearing to use");
  keys.reserve("nohtml","CV_DENSITY_Y","Takes 2/3 params specifying the upper and lower boundaries of a bin in y in fractional coordinates and the ammount of smearing to use");
  keys.reserve("nohtml","CV_DENSITY_Z","Takes 2/3 params specifying the upper and lower boundaries of a bin in z in fractional coordinates and the ammount of smearing to use");
  ActionWithDistribution::registerKeywords( keys );
} 

MultiColvar::MultiColvar(const ActionOptions&ao):
Action(ao),
ActionAtomistic(ao),
ActionWithValue(ao),
ActionWithDistribution(ao),
readatoms(false),
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
  bool dothis; std::vector<std::string> params;

  // Read SUM keyword
  if( keywords.exists("SUM") ) parseFlag("SUM",dothis);
  else dothis=false;
  if( dothis ){
     params.resize(1); Tools::convert(getNumberOfFunctionsInDistribution(),params[0]);
     addDistributionFunction( "sum", new sum(params) );
  }
  // Read AVERAGE keyword
  if( keywords.exists("AVERAGE") ) parseFlag("AVERAGE",dothis);
  else dothis=false;
  if( dothis ){
     params.resize(1); Tools::convert(getNumberOfFunctionsInDistribution(),params[0]);
     addDistributionFunction( "average", new mean(params) );
  }
  // Read MIN keyword
  if( keywords.exists("MIN") ) parseVector("MIN",params);
  else params.resize(0);
  if( params.size()!=0 ){
     addDistributionFunction( "min", new min(params) );
  }
  // Read Less_THAN keyword
  std::string tmpparam="none";
  if( keywords.exists("LESS_THAN") ) parse("LESS_THAN",tmpparam);
  if( tmpparam!="none" ){
     std::vector<std::string> data=Tools::getWords(tmpparam);
     std::string r0; Tools::parse(data,"R_0",r0); 
     params.resize(0); params.push_back(tmpparam);
     addDistributionFunction( "lt" + r0, new less_than(params) );
     params.resize(0);
  }
  // Read MORE_THAN keyword
  tmpparam="none";
  if( keywords.exists("MORE_THAN") ) parse("MORE_THAN",tmpparam);
  if( tmpparam!="none" ){
     std::vector<std::string> data=Tools::getWords(tmpparam);     
     std::string r0; Tools::parse(data,"R_0",r0);
     params.resize(0); params.push_back(tmpparam);
     addDistributionFunction( "gt" + r0, new more_than(params) );
     params.resize(0);
  }
  // Read HISTOGRAM keyword
  if( keywords.exists("HISTOGRAM") ) parseVector("HISTOGRAM",params);
  else params.resize(0);
  if( params.size()!=0 ){
      std::vector<double> range(2); parseVector("RANGE",range);
      if(range[1]<range[0]) error("range is meaningless");
      int nbins; Tools::convert(params[0],nbins);
      std::vector<std::string> hparams(2);
      if(params.size()==2){
          hparams.resize(3);
          hparams[2]=params[1];
      } else if(params.size()!=1){
          error("Histogram keyword should either specify just the number"
                " of bins or the number of bins and the ammount of smearing");
      }
      double lb,ub,delr=(range[1]-range[0])/static_cast<double>(nbins);
      for(int i=0;i<nbins;++i){
          lb=range[0]+i*delr; Tools::convert( lb, hparams[0] );
          ub=range[0]+(i+1)*delr; Tools::convert( ub, hparams[1] );
          addDistributionFunction( "between" + hparams[0] + "&" +hparams[1], new within(hparams) );
      }
  }
  // Read within keywords
  if( keywords.exists("WITHIN") ) parseVector("WITHIN",params);
  else params.resize(0);
  if( params.size()!=0 ){
      addDistributionFunction( "between" + params[0] + "&" +params[1], new within(params) );
  } else if( keywords.exists("WITHIN") ){
      for(unsigned i=1;;++i){
         if( !parseNumberedVector("WITHIN",i,params) ) break;
         addDistributionFunction( "between" + params[0] + "&" +params[1], new within(params) );
      }
  }

  // Establish whether or not this is a density ( requires special method )
  std::string dens="notdensity"; 
  if( keywords.exists("SPECIES") && !keywords.exists("SPECIESA") && !keywords.exists("SPECIESB") ){ dens="isdensity"; }

  // Read CV_DENSITY keywords
  if( keywords.exists("CV_DENSITY_X") ){
      plumed_massert( keywords.exists("CV_DENSITY_Y") && keywords.exists("CV_DENSITY_Z"), "please use CV_DENSITY_Y and CV_DENSITY_Z");
      std::vector<double> xbin,ybin,zbin;
      parseVector("CV_DENSITY_X",xbin); parseVector("CV_DENSITY_Y",ybin); parseVector("CV_DENSITY_Z",zbin);
      if( xbin.size()!=0 || ybin.size()!=0 || zbin.size()!=0 ){
          needsCentralAtomPosition=true;
          // Sort out xbin
          if( xbin.size()==0 ){ 
              xbin.push_back(0.0); xbin.push_back(1.0); xbin.push_back(20.0);
          } else if( xbin.size()==2){
              xbin.push_back(0.5);
          } else if( xbin.size()!=3){
              error("CV_DENSITY_X should have two or three arguments");
          }
          // Sort out ybin
          if( ybin.size()==0 ){ 
              ybin.push_back(0.0); ybin.push_back(1.0); ybin.push_back(20.0);
          } else if( ybin.size()==2){
              ybin.push_back(0.5);
          } else if( ybin.size()!=3){
              error("CV_DENSITY_Y should have two or three arguments");
          }
          // Sort out zbin
          if( zbin.size()==0 ){ 
              zbin.push_back(0.0); zbin.push_back(1.0); zbin.push_back(20.0);
          } else if( zbin.size()==2){
              zbin.push_back(0.5);
          } else if( zbin.size()!=3){
              error("CV_DENSITY_Z should have two or three arguments");
          }
          params.resize(10);
          Tools::convert( xbin[0], params[0] );
          Tools::convert( xbin[1], params[1] );
          Tools::convert( xbin[2], params[2] );
          Tools::convert( ybin[0], params[3] ); 
          Tools::convert( ybin[1], params[4] );
          Tools::convert( ybin[2], params[5] );
          Tools::convert( zbin[0], params[6] ); 
          Tools::convert( zbin[1], params[7] );
          Tools::convert( zbin[2], params[8] );
          params[9]=dens;
          addDistributionFunction( "meanCVinBox" , new cvdens(params) ); 
      } else {
          std::vector<double> xbin,ybin,zbin; params.resize(10);
          for(unsigned i=1;;++i){
              parseNumberedVector("CV_DENSITY_X",i,xbin); 
              parseNumberedVector("CV_DENSITY_Y",i,ybin);
              parseNumberedVector("CV_DENSITY_Z",i,zbin);
              if( xbin.size()==0 && ybin.size()==0 && zbin.size()==0 ) break;
              needsCentralAtomPosition=true;
              // Sort out xbin
              if( xbin.size()==0 ){ 
                  xbin.push_back(0.0); xbin.push_back(1.0); xbin.push_back(20.0);
              } else if( xbin.size()==2){
                  xbin.push_back(0.5);
              } else if( xbin.size()!=3){
                  error("CV_DENSITY_X should have two or three arguments");
              }
              // Sort out ybin
              if( ybin.size()==0 ){
                  ybin.push_back(0.0); ybin.push_back(1.0); ybin.push_back(20.0);
              } else if( ybin.size()==2){
                  ybin.push_back(0.5);
              } else if( ybin.size()!=3){
                  error("CV_DENSITY_Y should have two or three arguments");
              }
              // Sort out zbin
              if( zbin.size()==0 ){
                  zbin.push_back(0.0); zbin.push_back(1.0); zbin.push_back(20.0);
              } else if( zbin.size()==2){
                  zbin.push_back(0.5);
              } else if( zbin.size()!=3){
                  error("CV_DENSITY_Z should have two or three arguments");
              }
              Tools::convert( xbin[0], params[0] );
              Tools::convert( xbin[1], params[1] );
              Tools::convert( xbin[2], params[2] );
              Tools::convert( ybin[0], params[3] );
              Tools::convert( ybin[1], params[4] );
              Tools::convert( ybin[2], params[5] );
              Tools::convert( zbin[0], params[6] );
              Tools::convert( zbin[1], params[7] );
              Tools::convert( zbin[2], params[8] );
              params[9]=dens;
              std::string lab; Tools::convert(i,lab);
              addDistributionFunction( "meanCVforXin" + lab , new cvdens(params) ); 
              xbin.resize(0); ybin.resize(0); zbin.resize(0);
          }
      }
  }

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
       for(unsigned i=0;i<getNumberOfComponents();++i) getPntrToComponent(i)->resizeDerivatives( getThisFunctionsNumberOfDerivatives(i) );
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

     if( !usingDistributionFunctions() ){
       nder=getThisFunctionsNumberOfDerivatives(i); 
       if (forces.size()!=nder ) forces.resize(nder);
     } 

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





