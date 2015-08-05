/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013,2014 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "core/ActionRegister.h"
#include "tools/KernelFunctions.h"
#include "tools/IFile.h"
#include "MultiColvarFunction.h"

//+PLUMEDOC MCOLVARF PAMM
/*
Probabilistic analysis of molecular mofifs.

This is an implementation of the methods discussed in \cite{pamm}

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace multicolvar {

class PAMM : public MultiColvarFunction {
private:
  double regulariser;
  std::vector<KernelFunctions> kernels;
  std::vector<Value*> pos;
public:
  static void registerKeywords( Keywords& keys );
  PAMM(const ActionOptions&);
  ~PAMM();
/// We have to overwrite this here
  unsigned getNumberOfQuantities();
/// Calculate the weight of this object ( average of input weights )
  void calculateWeight( AtomValuePack& myatoms );
/// Actually do the calculation
  double compute( const unsigned& tindex, AtomValuePack& myatoms ) const ;
/// This returns the position of the central atom
  Vector getCentralAtom();
/// Is the variable periodic
  bool isPeriodic(){ return false; }
};

PLUMED_REGISTER_ACTION(PAMM,"PAMM")

void PAMM::registerKeywords( Keywords& keys ){
  MultiColvarFunction::registerKeywords( keys ); 
  keys.add("compulsory","CLUSTERS","the name of the file that contains the definitions of all the clusters");
  keys.add("compulsory","REGULARISE","0.001","don't allow the denominator to be smaller then this value");
  keys.use("MEAN"); keys.use("MORE_THAN"); keys.use("SUM"); keys.use("LESS_THAN"); 
}

PAMM::PAMM(const ActionOptions& ao):
Action(ao),
MultiColvarFunction(ao)
{
   // Check for reasonable input
   for(unsigned i=0;i<getNumberOfBaseMultiColvars();++i){
      if( getBaseMultiColvar(i)->getNumberOfQuantities()!=2 ) error("cannot use PAMM with " + getBaseMultiColvar(i)->getName() );
   }
   // This builds the lists
   buildSets(); 

   std::string filename; parse("CLUSTERS",filename);
   IFile ifile; 
   if( !ifile.FileExist(filename) ) error("could not find file named " + filename);
   ifile.open(filename); ifile.allowIgnoredFields(); double weight, wij;
   unsigned ncolv = getNumberOfBaseMultiColvars(); 
   std::vector<double> center( ncolv );
   std::vector<double> sig( 0.5*ncolv*(ncolv-1) + ncolv );   
   while( ifile.scanField("weight",weight) ) {
       for(unsigned i=0;i<center.size();++i){
          ifile.scanField("center-"+getBaseMultiColvar(i)->getLabel(),center[i]);
       }   
       unsigned k=0;
       for(unsigned i=0;i<center.size();++i){
           for(unsigned j=0;j<center.size();++j){
              ifile.scanField("width-"+getBaseMultiColvar(i)->getLabel()+"-" + getBaseMultiColvar(j)->getLabel(),wij );
              if( i<=j) sig[k++]=wij;
           }
       }  
       ifile.scanField();
       // std::string ktype; ifile.scanField("kerneltype",ktype);
       kernels.push_back( KernelFunctions( center, sig, "GAUSSIAN", true, weight, true ) );
   } 

 //  std::string kernelinpt;
 //  parse("TARGET-KERNEL",kernelinpt);
 //  KernelFunctions target( kernelinpt, true ); 
 //  if( target.ndim()!=getNumberOfBaseMultiColvars() ) error("mismatch for dimensionality of target kernel");
 //  kernels.push_back( target );
 //  for(int i=1;;i++){
 //     if( !parseNumbered("MIX-KERNEL",i,kernelinpt) ) break;
 //     KernelFunctions mykernel( kernelinpt, true );
 //     if( mykernel.ndim()!=getNumberOfBaseMultiColvars() ) error("mismatch for dimensionality of mix kernel");
 //     kernels.push_back( mykernel );
 //  }
 //  if( kernels.size()==1 ) error("no mix kernels defined");

   // I think we can put this back for anything with not usespecies
   // Now need to pass cutoff information to underlying distance multicolvars
   // for(unsigned i=0;i<getNumberOfBaseMultiColvars();++i){
   //     Distances* mydist = dynamic_cast<Distances*>( getBaseMultiColvar(i) );
   //     for(unsigned k=0;k<kernels.size();++k){
   //         mydist->setLinkCellCutoff( kernels[k].getCenter()[i]+kernels[k].getContinuousSupport()[i] );
   //     }
   // }

   for(unsigned i=0;i<getNumberOfBaseMultiColvars();++i){
       pos.push_back( new Value() );
       if( !getBaseMultiColvar(i)->isPeriodic() ){
           pos[i]->setNotPeriodic();
       } else {
           std::string min, max;
           getBaseMultiColvar(i)->retrieveDomain( min, max );
           pos[i]->setDomain( min, max );
       }
   }
}

PAMM::~PAMM(){
   for(unsigned i=0;i<pos.size();++i) delete pos[i];
}

unsigned PAMM::getNumberOfQuantities(){
   return 1 + kernels.size();    // getBaseMultiColvar(0)->getNumberOfQuantities();
}

void PAMM::calculateWeight( AtomValuePack& myatoms ){
   unsigned nvars = getNumberOfBaseMultiColvars();
   // Weight of point is average of weights of input colvars?
   std::vector<double> tval(2); double ww=0;
   for(unsigned i=0;i<nvars;++i){
       getVectorForTask( myatoms.getIndex(i), false, tval ); ww+=tval[i];
   }
   myatoms.setValue( 0, ww / static_cast<double>( nvars ) );

   if(!doNotCalculateDerivatives() ){
      double pref = 1.0 / static_cast<double>( nvars );
      MultiValue myder( getBaseMultiColvar(0)->getNumberOfQuantities(), getBaseMultiColvar(0)->getNumberOfDerivatives() );
      for(unsigned ivar=0;ivar<nvars;++ivar){
          // Resize the holder that gets the derivatives from the base class
          if( ivar>0 ) myder.resize( getBaseMultiColvar(ivar)->getNumberOfQuantities(), getBaseMultiColvar(ivar)->getNumberOfDerivatives() );
          // Get the values of derivatives
          getVectorDerivatives( myatoms.getIndex(ivar), false, myder );
          for(unsigned j=0;j<myder.getNumberActive();++j){
              unsigned jder=myder.getActiveIndex(j);
              myatoms.addDerivative( 0, jder, pref*myder.getDerivative( 1, jder ) );
          }
      }
   }
}

double PAMM::compute( const unsigned& tindex, AtomValuePack& myatoms ) const {
   unsigned nvars = getNumberOfBaseMultiColvars();    
   std::vector<std::vector<double> > tderiv( kernels.size() );
   for(unsigned i=0;i<kernels.size();++i) tderiv[i].resize( nvars );
   std::vector<double> tval(2), vals( kernels.size() ), dderiv( kernels.size(), 0 );

   for(unsigned i=0;i<nvars;++i){
       getVectorForTask( myatoms.getIndex(i), false, tval ); pos[i]->set( tval[1] );
   }

   // Evaluate the set of kernels 
   double denom=regulariser;
   for(unsigned i=0;i<kernels.size();++i){ 
       vals[i]=kernels[i].evaluate( pos, tderiv[i] );
       denom+=vals[i]; 
       for(unsigned j=0;j<nvars;++j) dderiv[j] += tderiv[i][j];
   }
   // Now set all values other than the first one 
   // This is because of some peverse choices in multicolvar
   for(unsigned i=1;i<kernels.size();++i) myatoms.setValue( 1+i, vals[i]/denom );

   if( !doNotCalculateDerivatives() ){
       double denom2=denom*denom; std::vector<double> mypref( 1 + kernels.size() );
       MultiValue myder( 2, getBaseMultiColvar(0)->getNumberOfDerivatives() );
       for(unsigned ivar=0;ivar<nvars;++ivar){
           // Resize the holder that gets the derivatives from the base class
           if( ivar>0 ){
              if( getBaseMultiColvar(ivar)->getNumberOfDerivatives()!=getBaseMultiColvar(ivar-1)->getNumberOfDerivatives() ){
                 myder.resize( 2, getBaseMultiColvar(ivar)->getNumberOfDerivatives() );
              }
           }
           // Get the values of the derivatives
           getVectorDerivatives( myatoms.getIndex(ivar), false, myder );
           // And calculate the derivatives
           for(unsigned i=0;i<kernels.size();++i){
               mypref[1+i] = tderiv[i][ivar]/denom - vals[i]*dderiv[ivar]/denom2;
               // This is basically doing the chain rule to get the final derivatives
               mergeVectorDerivatives( 1+i, 1+i, 2+i, myatoms.getIndex(ivar), mypref, myder, myatoms );
               myder.clearAll();
           }
       }
   }

   return vals[0] / denom;
}

Vector PAMM::getCentralAtom(){
   // Who knows how this should work
   plumed_error();
   return Vector(1.0,0.0,0.0);
}

}
}
