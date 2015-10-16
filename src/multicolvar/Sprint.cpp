/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2015 The plumed team
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
#include "AdjacencyMatrixAction.h"
#include "core/ActionRegister.h"

//+PLUMEDOC MCOLVARF SPRINT
/*
Calculate SPRINT topological variables.

The SPRINT topological variables are calculated from the largest eigenvalue, \f$\lambda\f$ of
an \f$n\times n\f$ adjacency matrix and its corresponding eigenvector, \f$\mathbf{V}\f$, using:

\f[
s_i = \sqrt{n} \lambda v_i
\f]

You can use different quantities to measure whether or not two given atoms/molecules are
adjacent or not in the adjacency matrix.  The simplest measure of adjacency is is whether
two atoms/molecules are within some cutoff of each other.  Further complexity can be added by 
insisting that two molecules are adjacent if they are within a certain distance of each 
other and if they have similar orientations.

\par Examples

This example input calculates the 7 SPRINT coordinates for a 7 atom cluster of Lennard-Jones
atoms and prints their values to a file.  In this input the SPRINT coordinates are calculated
in the manner described in ?? so two atoms are adjacent if they are within a cutoff:

\verbatim
DENSITY SPECIES=1-7 LABEL=d1
SPRINT ARG=d1 SWITCH={RATIONAL R_0=0.1} LABEL=ss
PRINT ARG=ss.* FILE=colvar 
\endverbatim

This example input calculates the 14 SPRINT coordinates foa a molecule composed of 7 hydrogen and
7 carbon atoms.  Once again two atoms are adjacent if they are within a cutoff:

\verbatim
DENSITY SPECIES=1-7 LABEL=c
DENSITY SPECIES=8-14 LABEL=h

SPRINT ...
 ARG=c,h
 SWITCH11={RATIONAL R_0=2.6 NN=6 MM=12}
 SWITCH12={RATIONAL R_0=2.2 NN=6 MM=12}
 SWITCH22={RATIONAL R_0=2.2 NN=6 MM=12}
 LABEL=ss
... SPRINT

PRINT ARG=ss.* FILE=colvar
\endverbatim

*/
//+ENDPLUMEDOC


namespace PLMD {
namespace multicolvar {

class Sprint : public AdjacencyMatrixAction {
private:
/// Square root of number of atoms
  double sqrtn;
/// Vector that stores eigenvalues
  std::vector<double> eigvals;
/// Vector that stores max eigenvector
  std::vector< std::pair<double,int> > maxeig;
/// Adjacency matrix
  Matrix<double> mymatrix;
/// Matrix that stores eigenvectors
  Matrix<double> eigenvecs;
public:
/// Create manual
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit Sprint(const ActionOptions&);
/// Do the matrix calculation
  void completeCalculation();
/// Sprint needs its only apply routine as it creates values
  void apply();
};

PLUMED_REGISTER_ACTION(Sprint,"SPRINT")

void Sprint::registerKeywords( Keywords& keys ){
  AdjacencyMatrixAction::registerKeywords( keys );
  componentsAreNotOptional(keys);
  keys.addOutputComponent("coord","default","all \\f$n\\f$ sprint coordinates are calculated and then stored in increasing order. "
                                            "the smallest sprint coordinate will be labelled <em>label</em>.coord-1, "
                                            "the second smallest will be labelleled <em>label</em>.coord-1 and so on");
}

Sprint::Sprint(const ActionOptions&ao):
Action(ao),
AdjacencyMatrixAction(ao),
eigvals( getFullNumberOfBaseTasks() ),
maxeig( getFullNumberOfBaseTasks() ),
mymatrix( getFullNumberOfBaseTasks(), getFullNumberOfBaseTasks() ),
eigenvecs( getFullNumberOfBaseTasks(), getFullNumberOfBaseTasks() )
{
   // Check on setup
   if( getNumberOfVessels()!=1 ) error("there should be no vessel keywords");
   // Check for bad colvar input
   for(unsigned i=0;i<getNumberOfBaseMultiColvars();++i){
      if( !getBaseMultiColvar(i)->hasDifferentiableOrientation() ) error("cannot use multicolvar of type " + getBaseMultiColvar(i)->getName() );
   }

   // Create all the values
   sqrtn = sqrt( static_cast<double>( getFullNumberOfBaseTasks() ) );
   for(unsigned i=0;i<getFullNumberOfBaseTasks();++i){
      std::string num; Tools::convert(i,num);
      addComponentWithDerivatives("coord-"+num);
      componentIsNotPeriodic("coord-"+num);
      getPntrToComponent(i)->resizeDerivatives( getNumberOfDerivatives() );
   }
}

void Sprint::completeCalculation(){
   // Get the adjacency matrix
   retrieveMatrix( mymatrix ); 
   // Diagonalize it
   diagMat( mymatrix, eigvals, eigenvecs );
   // Get the maximum eigevalue
   double lambda = eigvals[ getFullNumberOfBaseTasks()-1 ];
   // Get the corresponding eigenvector
   for(unsigned j=0;j<maxeig.size();++j){
       maxeig[j].first = fabs( eigenvecs( getFullNumberOfBaseTasks()-1, j ) );
       maxeig[j].second = j;
       // Must make all components of principle eigenvector +ve
       eigenvecs( getFullNumberOfBaseTasks()-1, j ) = maxeig[j].first;
   }

   // Reorder each block of eigevectors
   unsigned startnum=0;
   for(unsigned j=0;j<getNumberOfBaseMultiColvars();++j){
       unsigned nthis = getBaseMultiColvar(j)->getFullNumberOfTasks();
       // Sort into ascending order
       std::sort( maxeig.begin() + startnum, maxeig.begin() + startnum + nthis );
       // Used so we can do sorting in blocks 
       startnum += nthis;
   }
   // Set the sprint coordinates
   for(unsigned icomp=0;icomp<getNumberOfComponents();++icomp){
      getPntrToComponent(icomp)->set( sqrtn*lambda*maxeig[icomp].first );
   }

   // Parallelism
   unsigned rank, stride;
   if( serialCalculation() ){ stride=1; rank=0; }
   else { rank=comm.Get_rank(); stride=comm.Get_size(); } 

   // Derivatives
   MultiValue myvals( 2, getNumberOfDerivatives() );
   Matrix<double> mymat_ders( getNumberOfComponents(), getNumberOfDerivatives() );  
   std::vector<unsigned> catoms(2); unsigned nval = getFullNumberOfBaseTasks(); mymat_ders=0; 
   for(unsigned i=rank;i<getNumberOfActiveMatrixElements();i+=stride){
      decodeIndexToAtoms( getTaskCode(getActiveMatrixElement(i)), catoms ); unsigned j=catoms[0], k=catoms[1];
      double tmp1 = 2 * eigenvecs(nval-1,j)*eigenvecs(nval-1,k);
      for(unsigned icomp=0;icomp<getNumberOfComponents();++icomp){
          double tmp2 = 0.; 
          for(unsigned n=0;n<nval-1;++n){  // Need care on following line
              tmp2 += eigenvecs(n,maxeig[icomp].second) * ( eigenvecs(n,j)*eigenvecs(nval-1,k) + eigenvecs(n,k)*eigenvecs(nval-1,j) ) / ( lambda - eigvals[n] );
          }
          double prefactor=sqrtn*( tmp1*maxeig[icomp].first + tmp2*lambda );
          getAdjacencyVessel()->retrieveDerivatives( getActiveMatrixElement(i), false, myvals );
          for(unsigned jd=0;jd<myvals.getNumberActive();++jd){
              unsigned ider=myvals.getActiveIndex(jd);
              mymat_ders( icomp, ider ) += prefactor*myvals.getDerivative( 1, ider );
          }
      }
   }
   if( !serialCalculation() ) comm.Sum( mymat_ders );

   for(unsigned j=0;j<getNumberOfComponents();++j){
       Value* val=getPntrToComponent(j); 
       for(unsigned i=0;i<getNumberOfDerivatives();++i) val->addDerivative( i, mymat_ders(j,i) );
   }
}

void Sprint::apply(){
  std::vector<Vector>&   f(modifyForces());
  Tensor&           v(modifyVirial());
  unsigned          nat=getNumberOfAtoms();

  std::vector<double> forces( 3*getNumberOfAtoms() + 9 );
  for(unsigned i=0;i<getNumberOfComponents();++i){
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

}
}
