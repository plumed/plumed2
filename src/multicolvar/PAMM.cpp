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

Probabilistic analysis of molecular motifs (PAMM) was introduced in this paper \cite{pamm}.
The essence of this approach involves calculating some large set of collective variables 
for a set of atoms in a short trajectory and fitting this data using a Gaussian Mixture Model.
The idea is that modes in these distributions can be used to identify features such as hydrogen bonds or 
secondary structure types.

The assumption within this implementation is that the fitting of the Gaussian mixture model has been 
done elsewhere by a separate code.  You thus provide an input file to this action which contains the 
means, covariances and weights for a set of Gaussian kernels, \f$\{ \phi \}\f$.  The values and
derivatives for the following set of quantities is then computed:

\f[
s_k = \frac{ \phi_k}{ \sum_i \phi_i }
\f]

Each of the \f$\phi_k\f$ is a Gaussian function that acts on a set of quantities calculated within 
a \ref mcolv.  These might be \ref TORSIONS, \ref DISTANCES, \ref ANGLES or any one of the many
symmetry functions that are available within \ref mcolv actions.  These quantities are then inserted into 
the set of \f$n\f$ kernels that are in the the input file.   This will be done for multiple sets of values
for the input quantities and a final quantity will be calculated by summing the above \f$s_k\f$ values or
some transformation of the above.  This sounds less complicated than it is and is best understood by 
looking through the example given below.

\warning Mixing periodic and aperiodic \ref mcolv actions has not been tested

\par Examples

In this example I will explain in detail what the following input is computing:

\verbatim
MOLINFO MOLTYPE=protein STRUCTURE=M1d.pdb
psi: TORSIONS ATOMS1=@psi-2 ATOMS2=@psi-3 ATOMS3=@psi-4 
phi: TORSIONS ATOMS1=@phi-2 ATOMS2=@phi-3 ATOMS3=@phi-4
p: PAMM DATA=phi,psi CLUSTERS=clusters.dat MEAN1={COMPONENT=1} MEAN2={COMPONENT=2} 
PRINT ARG=p.mean-1,mean-2 FILE=colvar
\endverbatim

The best place to start our explanation is to look at the contents of the clusters.dat file

\verbatim
#! FIELDS weight center-phi center-psi width-phi-phi width-phi-psi width-psi-phi width-psi-psi      
      0.4     -1.0      -1.0      0.2     -0.1    -0.1    0.2   
      0.6      1.0      +1.0      0.1     -0.03   -0.03   0.1   
\endverbatim

This files contains the parameters of two two-dimensional Gaussian functions.  Each of these Gaussians has a weight, \f$w_k\f$,
a vector that specifies the position of its centre, \f$\mathbf{c}_\f$, and a covariance matrix, \f$\Sigma_k\f$.  The \f$\phi_k\f$ functions that 
we use to calculate our PAMM components are thus:

\f[
\phi_k = \frac{w_k}{N_k} \exp\left( -(\mathbf{s} - \mathbf{c}_k)^T \Sigma^{-1}_k (\mathbf{s} - \mathbf{c}_k) \right)
\f]

In the above \f$N_k\f$ is a normalisation factor that is calculated based on \f$\Sigma\f$.  The vector \f$\mathbf{s}\f$ is a vector of quantities
that are calculated by the \ref TORSIONS actions.  This vector must be two dimensional and in this case each component is the value of a 
torsion angle.  If we look at the two \ref TORSIONS actions in the above we are calculating the \f$\phi\f$ and \f$\psi\f$ backbone torsional
angles in a protein (Note the use of \ref MOLFINTO to make specification of atoms straightforward).  We thus calculate the values of our 
2 \f$ \{ \phi \} \f$  kernels 3 times.  The first time we use the \f$\phi\f$ and \f$\psi\f$ angles in the 2nd resiude of the protein,
the second time it is the \f$\phi\f$ and \f$\psi\f$ angles of the 3rd residue of the protein and the third time it is the \f$\phi\f$ and \f$\psi\f$ angles
of the 4th residue in the protein.  The final two quantities that are output by the print command, p.mean-1 and p.mean-2, are the averages 
over these three residues for the quantities:
\f[
s_1 = \frac{ \phi_1}{ \phi_1 + \phi_2 }
\f]
and 
\f[
s_2 = \frac{ \phi_2}{ \phi_1 + \phi_2 }
\f]
There is a great deal of flexibility in this input.  We can work with, and examine, any number of components, we can use any set of collective variables
and compute these PAMM variables and we can transform the PAMM variables themselves in a large number of different ways when computing these sums.
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
  keys.setComponentsIntroduction("When the label of this action is used as the input for a second you are not referring to a scalar quantity as you are in "
                                 "regular collective variables.  The label is used to reference the full set of quantities calculated by "
                                 "the action.  This is usual when using \\ref multicolvarfunction. Generally when doing this the set of PAMM variables "
                                 "will be referenced using the DATA keyword rather than ARG.\n\n"
                                 "This Action can be used to calculate the following scalar quantities directly from the underlying set of PAMM variables. "  
                                 "These quantities are calculated by employing the keywords listed below and they can be referenced elsewhere in the input "
                                 "file by using this Action's label followed by a dot and the name of the quantity.  The particular PAMM variable that should "
                                 "be averaged in a MEAN command or transformed by a swiching function in a LESS_THAN command is specified using the COMPONENT "
                                 "keyword. COMPONENT=1 refers to the PAMM variable in which the first kernel in your input file is on the numerator, COMPONENT=2 refers to "
                                 "PAMM variable in which the second kernel in the input file is on the numerator and so on.  The same quantity can be calculated "
                                 "multiple times for different PAMM components by a single PAMM action.  In this case the relevant keyword must appear multiple "
                                 "times on the input line followed by a numerical identifier i.e. MEAN1, MEAN2, ...  The quantities calculated when multiple "
                                 "MEAN commands appear on the input line can be referenece elsewhere in the input file by using the name of the quantity followed " 
                                 "followed by a numerical identifier e.g. <em>label</em>.lessthan-1, <em>label</em>.lessthan-2 etc.  Alternatively, you can "
                                 "customize the labels of the quantities by using the LABEL keyword in the description of the keyword.");
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

   bool mixed=getBaseMultiColvar(0)->isPeriodic();
   for(unsigned i=0;i<getNumberOfBaseMultiColvars();++i){
       pos.push_back( new Value() );
       if( getBaseMultiColvar(i)->isPeriodic()!=mixed ) warning("mixing of periodic and aperiodic base variables in pamm is untested");
       if( !getBaseMultiColvar(i)->isPeriodic() ){
           pos[i]->setNotPeriodic();
       } else {
           std::string min, max;
           getBaseMultiColvar(i)->retrieveDomain( min, max );
           pos[i]->setDomain( min, max );
       }
   }

   std::string filename; parse("CLUSTERS",filename);
   unsigned ncolv = getNumberOfBaseMultiColvars();
   IFile ifile; Matrix<double> covar( ncolv, ncolv ), icovar( ncolv, ncolv );
   if( !ifile.FileExist(filename) ) error("could not find file named " + filename);
   ifile.open(filename); ifile.allowIgnoredFields(); double weight, wij;
   std::vector<double> center( ncolv );
   std::vector<double> sig( 0.5*ncolv*(ncolv-1) + ncolv );   
   while( ifile.scanField("weight",weight) ) {
       for(unsigned i=0;i<center.size();++i){
          ifile.scanField("center-"+getBaseMultiColvar(i)->getLabel(),center[i]);
       }   
       for(unsigned i=0;i<center.size();++i){
           for(unsigned j=0;j<center.size();++j){
              ifile.scanField("width-"+getBaseMultiColvar(i)->getLabel()+"-" + getBaseMultiColvar(j)->getLabel(),covar(i,j) );
           }
       }  
       Invert( covar, icovar );
       unsigned k=0; 
       for(unsigned i=0;i<ncolv;i++){
           for(unsigned j=i;j<ncolv;j++){ sig[k]=icovar(i,j); k++; }
       }
       ifile.scanField(); 
       // std::string ktype; ifile.scanField("kerneltype",ktype);
       kernels.push_back( KernelFunctions( center, sig, "GAUSSIAN", "VON-MISES", weight ) );
       kernels[kernels.size()-1].normalize( pos );
   } 

   // I think we can put this back for anything with not usespecies
   // Now need to pass cutoff information to underlying distance multicolvars
   // for(unsigned i=0;i<getNumberOfBaseMultiColvars();++i){
   //     Distances* mydist = dynamic_cast<Distances*>( getBaseMultiColvar(i) );
   //     for(unsigned k=0;k<kernels.size();++k){
   //         mydist->setLinkCellCutoff( kernels[k].getCenter()[i]+kernels[k].getContinuousSupport()[i] );
   //     }
   // }
}

PAMM::~PAMM(){
   for(unsigned i=0;i<pos.size();++i) delete pos[i];
}

unsigned PAMM::getNumberOfQuantities(){
   return 1 + kernels.size();    
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
           for(unsigned i=0;i<kernels.size();++i) mypref[1+i] = tderiv[i][ivar]/denom - vals[i]*dderiv[ivar]/denom2;
           // This is basically doing the chain rule to get the final derivatives
           superChainRule( 1, 1, 1+kernels.size(), myatoms.getIndex(ivar), mypref, myder, myatoms );
           // And clear the derivatives
           myder.clearAll();
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
