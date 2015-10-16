/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014,2015 The plumed team
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
#include "DFSClustering.h"
#include "core/ActionRegister.h"

//+PLUMEDOC MCOLVARF DFSCLUSTERING 
/*
Cluster atoms based on their proximities and find the average properties of those atoms
in a cluster.

This collective variable was developed for looking at nucleation phenomena, where you are 
interested in using studying the behavior of atoms in small aggregates or nuclei.  In these sorts of 
problems you might be interested in the degree the atoms in a nucleus have adopted their crystalline 
structure or (in the case of heterogenous nucleation of a solute from a solvent) you might be 
interested in how many atoms are present in the largest cluster.  

This collective variable is a function of a multicolvar and as such it must take a multicolvar as input. 
At first this may seem rather peverse but after a while one hopefully will come to appreciate just how 
many different properties of the cluster simply by using DFSCLUSTERING in tandem with the different multicolvars.  
As examples of some things that could be done you could use DFSCLUSTERING in tandem with \ref COORDINATIONNUMBERS 
to look at the coordination numbers of the atoms in the largest cluster or you could use DFSCLUSTERING in tandem 
with \ref Q6 to look at the crystallinity in the largest cluster.  Regardless of what you do, however, this works
because a multicolvar is essentially calculating a large number of (we will call them) symmetry functions and 
each of these symmetry functions can be ascribed to a particular location in space.  As an example the symmetry
functions calculated by the command \ref COORDINATIONNUMBERS SPECIES=1-10 are the coordination numbers of atoms
1 through 10.  In other words, the first of these symmetry functions measures the number of atoms that are within
a certain cutoff of atom 1, the second measures the number of atoms that are within a cutoff of atom 2 and so on.
In terms of location it seems logical to suppose that the first of these symmetry functions is located at atom 1's 
position, the second is at atom 2's position and so on.  My point is that from the point of view of DFSCLUSTERING
it does not particularly matter what is calculated by the underlying multicolvar.  This action operates on the assumption
that what we have is a set of positions in space that have a scalar (or vector) quantity associated with them. 

Lets summarise thus far: we use multicolvar to convert an atomic configuration into a set of coordinates that have
symmetry function values associated with them.  Usually this means that we calculate something akin to a coordination
number for each of the atoms in our system and that the positions of these ``coordination numbers" are essentially the
positions of the atoms.  The DFS clustering begins here.  We start by calculating an adjacency matrix based on the 
positions of the various atoms.  In other words, we construct an \f$n\f$-symmetry functions by \f$n\f$-symmetry functions matrix
in which element \f$i,j\f$ is either 0 or 1 and measures whether or not symmetry functions \f$i\f$ and 
\f$j\f$ are adjacent.  As with \ref NLINKS or \ref SPRINT this measure of adjacency can measure whether
or not the positions of the two symmetry functions are within a certain cutoff.  Alternatively, if the symmetry functions 
are vector quanties such as \ref MOLECULES or \ref Q4 it can measure whether the two symmetry functions are within a 
certain cutoff of each other and whether or not the vector symmetry functions have the same orientation.  An important difference
between this Action and \ref NLINKS or \ref SPRINT is that in DFCLUSTERING a hard cutoff is used - although the user inputs a 
soft switching function the code sets any number that is greater than the tolerance set using the keyword TOL (which is by 
default machine epsilon) equal to one while any number less than this tolerance is set equal to zero.   As a consequence of this 
we can use the adjacency matrix to determine whether or not a particular pair of atoms are connected or not.  We
can thus use a DFS clustering algorithm to determine the sizes of the various connected components in the underlying graph.  These
connected components should then correspond to the various clusters in our system.

Once we have determined which atoms are in each of the connected components of our system we can then determine what the properties 
of the atoms are in each of the clusters for our system.  In other just as we would in any other multicolvar we can determine the 
average properties of the atoms in the largest cluster.  Alternatively, we can determine how many of the atoms in the largest cluster
in our system have a coordination number that is greater than a particular threshold and so on.  To be clear the mathematics that 
we use in determining these quantities is identical to that used in multicolvar.  The only difference is that now, rather than 
summing these quantities over all the symmetry functions, we sum only those symmetry functions for the atoms in a particular connected
cluster.  

\par Examples

In this example the \ref FCCCUBIC multicolvar is used to determine how similar the environment around each atom is
to an FCC unit cell.  DFS clustering is used with the adjacency matrix constructed using the switching function 
specified by SWITCH={CUBIC D_0=0.4   D_MAX=0.5} - in pratice this means the clustering algorithm views two atoms
as connected as long as they are within 0.5 nm of each other (as the tolerance is machine epsilon).  The final 
quantity calculated measures the number of atoms in the largest cluster that have an FCCCUBIC parameter greater than 0.035.

\verbatim
cubic1: FCCUBIC SPECIES=1-1000 SWITCH={CUBIC D_0=0.4  D_MAX=0.5} 
clust: DFSCLUSTERING DATA=cubic1 CLUSTER=1 SWITCH={CUBIC D_0=0.4   D_MAX=0.5} MORE_THAN={CUBIC D_0=0.035  D_MAX=0.045}
\endverbatim

This second example adds a layer of complexity to the previous one.  Once again the \ref FCCCUBIC multicolvar is used 
to determine how similar the environment around each atom is to an FCC unit cell.  However, in this case the intervening
\ref MFILTER_MORE means that when the clustering is performed by DFSCLUSTERING only those atoms from the \ref FCCUBIC action
that have an fcccubic switching function greater than 0.045 are explicitly clustered using the DFS algorithm.  As described on 
the manual page for \ref MFILTER_MORE this command takes each of the various symmetry functions input by the input multicolvar (in this 
case the FCCCUBIC multicolvar labelled cubic1) and ascribes each of them a weight, \f$w_i\f$ using:
\f[
w_i = 1.0 - s( v_i ) 
\f]
where \f$v_i\f$ is the value of the \f$i\f$th symmetry function and \f$s()\f$ is the switching function specified through the 
SWITCH keyword.  Within the DFSCLUSTERING these continuous switching functions have to be converted into discrete (discontinuous)
switching functions, which is done through the WTOL keyword.  The value given to this keyword ensures that any symmetry function
calculated by cubic1 whose weight (calculated by MFILTER_MORE) is less than 0.03 is ignored during the DFS clustering.  

Once again the final quantity to be calculated here is the number of atoms in the largest cluster that have an FCCCUBIC 
parameter greater than 0.035.  There is a very important and subtle distinction between what is being calculated here and
what was calculated in the first example, however.  The difference is best understood by considering an example.  Imagine we had
a very large extended but largely disordered cluster of atoms.  By virtue of sheer blind luck a few widely-spaced atoms in this disordered
cluster would have a configuration around them that might be very similar to the structure inside an fcc crystal.  This despite the
fact that the overall structure of the cluster does not resemble an fcc crystal.  Suppose that this were clustered using the first 
command above.  The final result you would get would be equal to the number of fcc-like atoms in the disordered cluster as essentially
all the atoms in the system are connected - the fact that there is not an ordered arrangement is not particularly important.  Now 
suppose instead that this structure were analysed using the command below.  The final value of this quantity would most likely be
equal to one as it is unlikely that two of the fcc-like atoms are close together in the cluster.  The point is that by inserting the 
MFILTER_MORE command we have restricted the DFS clustering to only those atoms that are reasonably close in structure to fcc.  These
atoms are not necessarily close together in are disordered agregate and are thus not necessarily connected. 

\verbatim
cubic1: FCCUBIC SPECIES=1-1000 SWITCH={CUBIC D_0=0.4  D_MAX=0.5} TOL=0.03 #MEAN
cf: MFILTER_MORE DATA=cubic1 SWITCH={CUBIC D_0=0.035 D_MAX=0.045}
clust: DFSCLUSTERING DATA=cf WTOL=0.03 CLUSTER=1 SWITCH={CUBIC D_0=0.4   D_MAX=0.5} MORE_THAN={CUBIC D_0=0.035  D_MAX=0.045}
\endverbatim


*/
//+ENDPLUMEDOC

namespace PLMD {
namespace crystallization {

class DFSBasic : public DFSClustering {
private:
/// The cluster we are looking for
  unsigned clustr;
/// The buffer (we keep a copy here to avoid resizing)
  std::vector<double> buffer;
public:
/// Create manual
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit DFSBasic(const ActionOptions&);
///
  void doCalculationOnCluster();
};

PLUMED_REGISTER_ACTION(DFSBasic,"DFSCLUSTERING")

void DFSBasic::registerKeywords( Keywords& keys ){
  DFSClustering::registerKeywords( keys );
  keys.add("compulsory","CLUSTER","1","which cluster would you like to look at 1 is the largest cluster, 2 is the second largest, 3 is the the third largest and so on.");
  keys.use("WTOL"); keys.use("USE_ORIENTATION");
  keys.use("MEAN"); keys.use("MORE_THAN"); keys.use("LESS_THAN");
  if( keys.reserved("VMEAN") ) keys.use("VMEAN");
  if( keys.reserved("VSUM") ) keys.use("VSUM");
  keys.use("BETWEEN"); keys.use("HISTOGRAM"); keys.use("MOMENTS"); keys.use("ALT_MIN"); 
  keys.use("MIN"); keys.use("MAX"); keys.use("SUM"); keys.remove("LOWMEM"); keys.use("HIGHMEM");
  keys.use("LOWEST"); keys.use("HIGHEST");
}

DFSBasic::DFSBasic(const ActionOptions&ao):
Action(ao),
DFSClustering(ao)
{
   // Find out which cluster we want
   parse("CLUSTER",clustr);

   if( clustr<1 ) error("cannot look for a cluster larger than the largest cluster");
   if( clustr>getFullNumberOfBaseTasks() ) error("cluster selected is invalid - too few atoms in system");

   // Setup the various things this will calculate
   readVesselKeywords();
}

void DFSBasic::doCalculationOnCluster(){
   std::vector<unsigned> myatoms; retrieveAtomsInCluster( clustr, myatoms );
   unsigned size=comm.Get_size(), rank=comm.Get_rank();

   // Now calculate properties of the largest cluster 
   ActionWithVessel::doJobsRequiredBeforeTaskList();  // Note we loose adjacency data by doing this
   // Get rid of bogus derivatives
   clearDerivatives(); 

   // Get size for buffer
   getAdjacencyVessel()->setFinishedTrue();    // This ensures buffer size is smaller
   unsigned bsize=0, bufsize=getSizeOfBuffer(bsize); 
   if( buffer.size()!=bufsize ) buffer.resize( bufsize );
   buffer.assign( bufsize, 0.0 );   

   std::vector<double> vals( getNumberOfQuantities() ); std::vector<unsigned> der_index;
   MultiValue myvals( getNumberOfQuantities(), getNumberOfDerivatives() );
   MultiValue bvals( getNumberOfQuantities(), getNumberOfDerivatives() );
   for(unsigned j=rank;j<myatoms.size();j+=size){
       // Note loop above over array containing atoms so this is load balanced
       unsigned i=myatoms[j];
       // Need to copy values from base action
       getVectorForTask( i, false, vals );
       if( !doNotCalculateDerivatives() ) getVectorDerivatives( i, false, myvals );
       for(unsigned k=0;k<vals.size();++k) myvals.setValue( k, vals[k] );
       // Run calculate all vessels
       calculateAllVessels( i, myvals, bvals, buffer, der_index );
       myvals.clearAll();
   }
   // MPI Gather everything
   if( buffer.size()>0 ) comm.Sum( buffer );
   finishComputations( buffer );
}

}
}
