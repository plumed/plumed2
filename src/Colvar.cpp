#include "Colvar.h"
#include "PlumedMain.h"
#include <vector>
#include <string>
#include <cassert>

using namespace std;
using namespace PLMD;

Colvar::Colvar(const ActionOptions&ao) :
ActionAtomistic(ao),
doall(true),
domin(false),
domax(false),
dototal(false),
domean(false),
dolt(false),
domt(false)
{
  forbidKeyword("STRIDE");
  registerKeyword(0, "UPDATE", "frequency for updates of neighbour lists and dynamic groups");
  registerKeyword(0, "NL_CUT", "distance cutoff for neighbour lists");
  registerKeyword(0, "MIN", "calculate the minimum for the defined colvars");
  registerKeyword(0, "MAX", "calculate the maximum for the defined colvars");
  registerKeyword(0, "BETA", "(default=50) value used to create smooth derivatives in calculations of MIN/MAX"); 
  registerKeyword(0, "SUM", "calculate the sum of all the defined colvars");
  registerKeyword(0, "AVERAGE", "compute the average of the defined colvars");
  registerKeyword(0, "LESS_THAN", "compute the number of colvars that that are less than than a particular value using a smooth switching function");
  registerKeyword(0, "MORE_THAN", "compute the number of colvars that that are more than than a particular value using a smooth switching function");
  registerKeyword(0, "LOGIC_NN", "(default=6) value of NN in switching functions for MORE_THAN/LESS_THAN");
  registerKeyword(0, "LOGIC_MM", "(default=12) value of MM in switching functions for MORE_THAN/LESS_THAN");
}

void Colvar::readActionColvar( int natoms, const std::vector<double>& domain ){
  unsigned ngrp=static_cast<unsigned>(natoms);
  readActionAtomistic( natoms, ngrp );
  readActionWithExternalArguments( 3*getNumberOfAtoms()+9, domain );

  // Setup everything for calculation of individual colvars
  skipto.resize( function_indexes.size() ); derivatives.resize( natoms ); 
  for(unsigned i=0;i<function_indexes.size();++i) skipto[i]=i+1;

  // Resize stuff for applying forces
  f.resize( getNumberOfAtoms() ); forces.resize( 3*getNumberOfAtoms()+9 );

  // Read in everything that tells us what sort of calculation we are doing
  parseFlag("MIN",domin);
  if (domin){
     doall=false;
     beta=50.0; parse("BETA",beta);
     addValue("min", true, true);
     log.printf("  using minimum value \n");
  }
  parseFlag("MAX",domax);
  if (domax){
     doall=false;
     error("Dont know how to implement this yet");
     addValue("max", true, true);
     log.printf("  using maximum value \n");
  }
  parseFlag("SUM",dototal);
  if( dototal ){
     doall=false;
     addValue("sum", true, true);
     log.printf("  using sum of all values \n");
  }
  parseFlag("AVERAGE",domean);
  if( domean ){
     doall=false;
     addValue("average", true, true);
     log.printf("  using average value \n");
  }

  /*
  if ( testForKey("LESS_THAN") ) {
     doall=false; dolt=true;
     std::vector<double> r_0; parseVector("LESS_THAN",r_0);
     if ( r_0.size()!=1 ) error("Input for LESS_THAN makes no sense");
     int nn=6; parse("LOGIC_NN",nn);
     int mm=12; parse("LOGIC_MM",mm);
     ltswitch.set(nn, mm, r_0[0], 0.0);
     log.printf("  number of values less than %f.  Switching function paramers are %d %d \n", r_0[0], nn, mm );
     addValue("less_than", true, true);
  }

  if ( testForKey("MORE_THAN") ) {
     doall=false; domt=true;
     std::vector<double> r_0; parseVector("MORE_THAN",r_0);
     if ( r_0.size()!=1 ) error("Input for MORE_THAN makes no sense");
     int nn=6; parse("LOGIC_NN",nn);
     int mm=12; parse("LOGIC_MM",mm);
     mtswitch.set(nn, mm, r_0[0], 0.0);
     log.printf("  number of values greater than %f.  Switching function paramers are %d %d \n", r_0[0], nn, mm );
     addValue("more_than", true, true);
  }
  */

  if( doall ){
     std::string n;
     for(unsigned i=0;i<function_indexes.size();++i){
        Tools::convert(i,n);
        addValue("value" + n, false, true );
     }
  }
}

void Colvar::interpretGroupsKeyword( const unsigned& natoms, const std::string& atomGroupName, const std::vector<std::vector<unsigned> >& groups ){
  std::vector<unsigned> tmplist;

  if( atomGroupName!=getLabel() ){
      log.printf("  using atoms specified in group %s\n", atomGroupName.c_str() );
  } else {
      for(unsigned i=0;i<groups.size();++i){
          log.printf("  atoms in group %d : ", i+1 );
          for(unsigned j=0;j<groups[i].size();++j) log.printf("%s ", plumed.getAtoms().interpretIndex( groups[i][j] ).c_str() );
          log.printf("\n");
      }
  }

  if( natoms==0 ){
      unsigned accum=0;
      for(unsigned i=0;i<groups.size();++i){
         tmplist.resize(0);
         for(unsigned j=0;j<groups[i].size();++j) tmplist.push_back(accum+j);
         function_indexes.push_back( tmplist ); accum+=groups[i].size();
      }
  } else {
      if( natoms!=2 ){
          error("you can only use groups for indistinguishable colvars if the number of atoms in each colvar is equal to 2");
          assert( groups.size()<=2 );
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
}

void Colvar::interpretAtomsKeyword( const std::vector<std::vector<unsigned> >& flist ){
  unsigned accum=0; std::vector<unsigned> tmplist;
  log.printf("  one colvar will be calculated for each of the following sets of atoms ");
  for(unsigned i=0;i<flist.size();++i){ 
    log.printf("( %s",plumed.getAtoms().interpretIndex( flist[i][0] ).c_str() );
    tmplist.resize(0); tmplist.push_back( accum );
    for(unsigned j=1;j<flist[i].size();++j){
        tmplist.push_back( accum + j );
        log.printf( ", %s",plumed.getAtoms().interpretIndex( flist[i][j] ).c_str() );
    }
    log.printf(" ) : "); accum+=flist[i].size();
    function_indexes.push_back( tmplist );
  }
  log.printf("\n");
}

void Colvar::calculate(){
  double value;
  for (unsigned i=0; i<skipto.size(); i=skipto[i] ) {
     value=calcFunction( function_indexes[i], derivatives, virial );
     if (doall) {
        mergeFunctions( i, 1.0 );
        setValue( i, value, 1.0 );
     } else {
        assert(false);
     }
  }

}

//void ActionAtomistic::calculate(){

  // This will eventually update dynamic content (neighbour lists and so on)   -- Gareth we must sort prepare
  //if( nl_st!=0 && getStep()-nl_last>=nl_st ){
  //   for(unsigned i=0;i<skips.size();++i) skips[i]=false;  

  //   if( atomGroupName!=getLabel() ){ 
  //      // Find the action
  //       GenericDynamicGroup* action=plumed.getActionSet().selectWithLabel<GenericDynamicGroup*>(atomGroupName);
  //       assert(action);
  //       // Run a routine that selects the atoms in the group
  //       action.updateSelection( skips );
  //   }
  //   updateNeighbourLists( skips );
  //   plumed.getAtoms().updateAtomsInSelection( atomGroupName, skips );
  //}

  // MPI this eventually
/*  double mintotal, ttotal, atotal, maxtotal,lttotal,mttotal;
  mintotal=maxtotal=ttotal=atotal=lttotal=mttotal=0.0; 
  for (unsigned i=0; i<skipto.size(); i=skipto[i] ) {
     value=calcFunction( i, derivatives, virial );
     if (doall) {
        mergeFunctions(i, 1.0, derivatives, virial);
        setOutputValue( i, value, 1.0 );
     } else {
       if ( domin ) {
           // Minimum
           tmp=exp( beta/value ); mintotal+=tmp; df=tmp/(value*value);
           mergeFunctions("min", df, derivatives, virial );
        } else if ( domax ) {
           // Maximum
        } else if ( dototal ) {
           ttotal+=value;
           mergeFunctions("sum", 1.0, derivatives, virial );
        } else if ( domean ) {
           // The mean
           atotal+=value;
           mergeFunctions("average", 1.0, derivatives, virial );
        } else if ( dolt ) {
           // Less than
           double tmp=ltswitch.calculate(value, df);
           lttotal+=tmp;
           mergeFunctions("less_than", df*value, derivatives, virial );
        } else if ( domt ) {
           // More than
           mttotal+=1.0 - mtswitch.calculate(value, df);
           mergeFunctions("more_than", -df*value, derivatives, virial );
        }
     }
  }
  if ( domin ){ double dist=beta/std::log(mintotal); setOutputValue("min", dist, dist*dist/mintotal ); }
  if ( domax ) { }
  if ( dototal ) { setOutputValue("sum", ttotal, 1.0 ); }
  if ( domean ) { setOutputValue("average", atotal/skipto.size(), 1.0/skipto.size() ); }
  if ( dolt ) { setOutputValue("less_than", lttotal, 1.0 ); }
  if ( domt ) { setOutputValue("more_than", mttotal, 1.0 ); }
}
*/

void Colvar::apply(){

  const unsigned nat=f.size();
  for(unsigned i=0;i<nat;i++){
    f[i][0]=0.0; f[i][1]=0.0; f[i][2]=0.0;
  }

  Tensor v; v.clear();
  for(int i=0;i<getNumberOfValues();++i){
    if( getForces( i, forces ) ){
       for(unsigned j=0;j<nat;++j){
          f[j][0]+=forces[3*j+0];
          f[j][0]+=forces[3*j+1];
          f[j][0]+=forces[3*j+2];
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
  applyForces( f, v );
}





