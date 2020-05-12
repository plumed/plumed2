/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

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
#include "MultiColvarBase.h"
#include "AtomValuePack.h"
#include "core/ActionRegister.h"
#include "tools/SwitchingFunction.h"

//+PLUMEDOC MCOLVARF LOCAL_AVERAGE
/*
Calculate averages over spherical regions centered on atoms

As is explained in <a href="http://www.youtube.com/watch?v=iDvZmbWE5ps"> this video </a> certain multicolvars
calculate one scalar quantity or one vector for each of the atoms in the system.  For example
\ref COORDINATIONNUMBER measures the coordination number of each of the atoms in the system and \ref Q4 measures
the fourth order Steinhardt parameter for each of the atoms in the system.  These quantities provide tell us something about
the disposition of the atoms in the first coordination sphere of each of the atoms of interest.  Lechner and Dellago \cite dellago-q6
have suggested that one can probe local order in a system by taking the average value of such symmetry functions over
the atoms within a spherical cutoff of each of these atoms in the systems.  When this is done with Steinhardt parameters
they claim this gives a coordinate that is better able to distinguish solid and liquid configurations of Lennard-Jones atoms.

You can calculate such locally averaged quantities within plumed by using the LOCAL_AVERAGE command.  This command calculates
the following atom-centered quantities:

\f[
s_i = \frac{ c_i + \sum_j \sigma(r_{ij})c_j }{ 1 + \sum_j \sigma(r_{ij}) }
\f]

where the \f$c_i\f$ and \f$c_j\f$ values can be for any one of the symmetry functions that can be calculated using plumed
multicolvars.  The function \f$\sigma( r_{ij} )\f$ is a \ref switchingfunction that acts on the distance between
atoms \f$i\f$ and \f$j\f$.  Lechner and Dellago suggest that the parameters of this function should be set so that it the function is equal to one
when atom \f$j\f$ is in the first coordination sphere of atom \f$i\f$ and is zero otherwise.

The \f$s_i\f$ quantities calculated using the above command can be again thought of as atom-centered symmetry functions.  They
thus operate much like multicolvars.  You can thus calculate properties of the distribution of \f$s_i\f$ values using MEAN, LESS_THAN, HISTOGRAM
and so on.  You can also probe the value of these averaged variables in regions of the box by using the command in tandem with the
\ref AROUND command.

\par Examples

This example input calculates the coordination numbers for all the atoms in the system.  These coordination numbers are then averaged over
spherical regions.  The number of averaged coordination numbers that are greater than 4 is then output to a file.

\plumedfile
COORDINATIONNUMBER SPECIES=1-64 D_0=1.3 R_0=0.2 LABEL=d1
LOCAL_AVERAGE SPECIES=d1 SWITCH={RATIONAL D_0=1.3 R_0=0.2} MORE_THAN={RATIONAL R_0=4} LABEL=la
PRINT ARG=la.* FILE=colvar
\endplumedfile

This example input calculates the \f$q_4\f$ (see \ref Q4) vectors for each of the atoms in the system.  These vectors are then averaged
component by component over a spherical region.  The average value for this quantity is then output to a file.  This calculates the
quantities that were used in the paper by Lechner and Dellago \cite dellago-q6

\plumedfile
Q4 SPECIES=1-64 SWITCH={RATIONAL D_0=1.3 R_0=0.2} LABEL=q4
LOCAL_AVERAGE SPECIES=q4 SWITCH={RATIONAL D_0=1.3 R_0=0.2} MEAN LABEL=la
PRINT ARG=la.* FILE=colvar
\endplumedfile

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace multicolvar {

class LocalAverage : public MultiColvarBase {
private:
/// Cutoff
  double rcut2;
/// The switching function that tells us if atoms are close enough together
  SwitchingFunction switchingFunction;
public:
  static void registerKeywords( Keywords& keys );
  explicit LocalAverage(const ActionOptions&);
/// We have to overwrite this here
  unsigned getNumberOfQuantities() const override;
/// Actually do the calculation
  double compute( const unsigned& tindex, AtomValuePack& myatoms ) const override;
/// We overwrite this in order to have dumpmulticolvar working for local average
  void normalizeVector( std::vector<double>& vals ) const override {}
/// Is the variable periodic
  bool isPeriodic() override { return false; }
};

PLUMED_REGISTER_ACTION(LocalAverage,"LOCAL_AVERAGE")

void LocalAverage::registerKeywords( Keywords& keys ) {
  MultiColvarBase::registerKeywords( keys );
  keys.add("compulsory","NN","6","The n parameter of the switching function ");
  keys.add("compulsory","MM","0","The m parameter of the switching function; 0 implies 2*NN");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous switching function defined above. "
           "The following provides information on the \\ref switchingfunction that are available. "
           "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
  // Use actionWithDistributionKeywords
  keys.use("SPECIES"); keys.use("SPECIESA"); keys.use("SPECIESB");
  keys.remove("LOWMEM"); keys.use("MEAN"); keys.use("MORE_THAN"); keys.use("LESS_THAN");
  keys.use("BETWEEN"); keys.use("HISTOGRAM"); keys.use("MOMENTS");
  keys.addFlag("LOWMEM",false,"lower the memory requirements");
  if( keys.reserved("VMEAN") ) keys.use("VMEAN");
  if( keys.reserved("VSUM") ) keys.use("VSUM");
}

LocalAverage::LocalAverage(const ActionOptions& ao):
  Action(ao),
  MultiColvarBase(ao)
{
  if( getNumberOfBaseMultiColvars()>1 ) error("local average with more than one base colvar makes no sense");
  // Read in the switching function
  std::string sw, errors; parse("SWITCH",sw);
  if(sw.length()>0) {
    switchingFunction.set(sw,errors);
  } else {
    double r_0=-1.0, d_0; int nn, mm;
    parse("NN",nn); parse("MM",mm);
    parse("R_0",r_0); parse("D_0",d_0);
    if( r_0<0.0 ) error("you must set a value for R_0");
    switchingFunction.set(nn,mm,r_0,d_0);
  }
  log.printf("  averaging over central molecule and those within %s\n",( switchingFunction.description() ).c_str() );
  rcut2 = switchingFunction.get_dmax()*switchingFunction.get_dmax();
  setLinkCellCutoff( switchingFunction.get_dmax() );
  std::vector<AtomNumber> all_atoms; setupMultiColvarBase( all_atoms );
}

unsigned LocalAverage::getNumberOfQuantities() const {
  return getBaseMultiColvar(0)->getNumberOfQuantities();
}

double LocalAverage::compute( const unsigned& tindex, AtomValuePack& myatoms ) const {
  double sw, dfunc; MultiValue& myvals = myatoms.getUnderlyingMultiValue();
  std::vector<double> values( getBaseMultiColvar(0)->getNumberOfQuantities() );

  getInputData( 0, false, myatoms, values );
  myvals.addTemporyValue( values[0] );
  if( values.size()>2 ) {
    for(unsigned j=2; j<values.size(); ++j) myatoms.addValue( j, values[0]*values[j] );
  } else {
    myatoms.addValue( 1, values[0]*values[1] );
  }

  if( !doNotCalculateDerivatives() ) {
    MultiValue& myder=getInputDerivatives( 0, false, myatoms );

    // Convert input atom to local index
    unsigned katom = myatoms.getIndex( 0 ); plumed_dbg_assert( katom<atom_lab.size() ); plumed_dbg_assert( atom_lab[katom].first>0 );
    // Find base colvar
    unsigned mmc=atom_lab[katom].first - 1; plumed_dbg_assert( mybasemulticolvars[mmc]->taskIsCurrentlyActive( atom_lab[katom].second ) );
    // Get start of indices for this atom
    unsigned basen=0; for(unsigned i=0; i<mmc; ++i) basen+=mybasemulticolvars[i]->getNumberOfDerivatives() - 9;
    plumed_dbg_assert( basen%3==0 ); // Check the number of atoms is consistent with input derivatives

    unsigned virbas = myvals.getNumberOfDerivatives()-9;
    if( values.size()>2 ) {
      for(unsigned j=0; j<myder.getNumberActive(); ++j) {
        unsigned jder=myder.getActiveIndex(j);
        if( jder<mybasemulticolvars[mmc]->getNumberOfDerivatives()-9 ) {
          unsigned kder=basen+jder;
          for(unsigned k=2; k<values.size(); ++k) {
            myatoms.addDerivative( k, kder, values[0]*myder.getDerivative(k,jder) );
            myatoms.addDerivative( k, kder, values[k]*myder.getDerivative(0,jder) );
          }
        } else {
          unsigned kder=virbas + (jder - mybasemulticolvars[mmc]->getNumberOfDerivatives() + 9);
          for(unsigned k=2; k<values.size(); ++k) {
            myatoms.addDerivative( k, kder, values[0]*myder.getDerivative(k,jder) );
            myatoms.addDerivative( k, kder, values[k]*myder.getDerivative(0,jder) );
          }
        }
      }
    } else {
      for(unsigned j=0; j<myder.getNumberActive(); ++j) {
        unsigned jder=myder.getActiveIndex(j);
        if( jder<mybasemulticolvars[mmc]->getNumberOfDerivatives()-9 ) {
          unsigned kder=basen+jder;
          myatoms.addDerivative( 1, kder, values[0]*myder.getDerivative(1,jder) );
          myatoms.addDerivative( 1, kder, values[1]*myder.getDerivative(0,jder) );
        } else {
          unsigned kder=virbas + (jder - mybasemulticolvars[mmc]->getNumberOfDerivatives() + 9);
          myatoms.addDerivative( 1, kder, values[0]*myder.getDerivative(1,jder) );
          myatoms.addDerivative( 1, kder, values[1]*myder.getDerivative(0,jder) );
        }
      }
    }
    for(unsigned j=0; j<myder.getNumberActive(); ++j) {
      unsigned jder=myder.getActiveIndex(j);
      if( jder<mybasemulticolvars[mmc]->getNumberOfDerivatives()-9 ) {
        unsigned kder=basen+jder;
        myvals.addTemporyDerivative( kder, myder.getDerivative(0, jder) );
      } else {
        unsigned kder=virbas + (jder - mybasemulticolvars[mmc]->getNumberOfDerivatives() + 9);
        myvals.addTemporyDerivative( kder, myder.getDerivative(0, jder) );
      }
    }
    myder.clearAll();
  }

  for(unsigned i=1; i<myatoms.getNumberOfAtoms(); ++i) {
    Vector& distance=myatoms.getPosition(i);  // getSeparation( myatoms.getPosition(0), myatoms.getPosition(i) );
    double d2;
    if ( (d2=distance[0]*distance[0])<rcut2 &&
         (d2+=distance[1]*distance[1])<rcut2 &&
         (d2+=distance[2]*distance[2])<rcut2 &&
         d2>epsilon) {

      sw = switchingFunction.calculateSqr( d2, dfunc );

      getInputData( i, false, myatoms, values );
      if( values.size()>2 ) {
        for(unsigned j=2; j<values.size(); ++j) myatoms.addValue( j, sw*values[0]*values[j] );
      } else {
        myatoms.addValue( 1, sw*values[0]*values[1] );
      }
      myvals.addTemporyValue(sw);

      if( !doNotCalculateDerivatives() ) {
        Tensor vir(distance,distance);
        MultiValue& myder=getInputDerivatives( i, false, myatoms );

        // Convert input atom to local index
        unsigned katom = myatoms.getIndex( i ); plumed_dbg_assert( katom<atom_lab.size() ); plumed_dbg_assert( atom_lab[katom].first>0 );
        // Find base colvar
        unsigned mmc=atom_lab[katom].first - 1; plumed_dbg_assert( mybasemulticolvars[mmc]->taskIsCurrentlyActive( atom_lab[katom].second ) );
        // Get start of indices for this atom
        unsigned basen=0; for(unsigned j=0; j<mmc; ++j) basen+=mybasemulticolvars[j]->getNumberOfDerivatives() - 9;
        plumed_dbg_assert( basen%3==0 ); // Check the number of atoms is consistent with input derivatives

        unsigned virbas = myvals.getNumberOfDerivatives()-9;
        if( values.size()>2 ) {
          for(unsigned j=0; j<myder.getNumberActive(); ++j) {
            unsigned jder=myder.getActiveIndex(j);
            if( jder<mybasemulticolvars[mmc]->getNumberOfDerivatives()-9 ) {
              unsigned kder=basen+jder;
              for(unsigned k=2; k<values.size(); ++k) {
                myatoms.addDerivative( k, kder, sw*values[0]*myder.getDerivative(k,jder) );
                myatoms.addDerivative( k, kder, sw*values[k]*myder.getDerivative(0,jder) );
              }
            } else {
              unsigned kder=virbas + (jder - mybasemulticolvars[mmc]->getNumberOfDerivatives() + 9);
              for(unsigned k=2; k<values.size(); ++k) {
                myatoms.addDerivative( k, kder, sw*values[0]*myder.getDerivative(k,jder) );
                myatoms.addDerivative( k, kder, sw*values[k]*myder.getDerivative(0,jder) );
              }
            }
          }
          for(unsigned k=2; k<values.size(); ++k) {
            addAtomDerivatives( k, 0, (-dfunc)*values[0]*values[k]*distance, myatoms );
            addAtomDerivatives( k, i, (+dfunc)*values[0]*values[k]*distance, myatoms );
            myatoms.addBoxDerivatives( k, (-dfunc)*values[0]*values[k]*vir );
          }
        } else {
          for(unsigned j=0; j<myder.getNumberActive(); ++j) {
            unsigned jder=myder.getActiveIndex(j);
            if( jder<mybasemulticolvars[mmc]->getNumberOfDerivatives()-9 ) {
              unsigned kder=basen+jder;
              myatoms.addDerivative( 1, kder, sw*values[0]*myder.getDerivative(1,jder) );
              myatoms.addDerivative( 1, kder, sw*values[1]*myder.getDerivative(0,jder) );
            } else {
              unsigned kder=virbas + (jder - mybasemulticolvars[mmc]->getNumberOfDerivatives() + 9);
              myatoms.addDerivative( 1, kder, sw*values[0]*myder.getDerivative(1,jder) );
              myatoms.addDerivative( 1, kder, sw*values[1]*myder.getDerivative(0,jder) );
            }
          }
          addAtomDerivatives( 1, 0, (-dfunc)*values[0]*values[1]*distance, myatoms );
          addAtomDerivatives( 1, i, (+dfunc)*values[0]*values[1]*distance, myatoms );
          myatoms.addBoxDerivatives( 1, (-dfunc)*values[0]*values[1]*vir );
        }
        // And the bit we use to average the vector
        addAtomDerivatives( -1, 0, (-dfunc)*values[0]*distance, myatoms );
        addAtomDerivatives( -1, i, (+dfunc)*values[0]*distance, myatoms );
        for(unsigned j=0; j<myder.getNumberActive(); ++j) {
          unsigned jder=myder.getActiveIndex(j);
          if( jder<mybasemulticolvars[mmc]->getNumberOfDerivatives()-9 ) {
            unsigned kder=basen+jder;
            myvals.addTemporyDerivative( kder, sw*myder.getDerivative(0, jder) );
          } else {
            unsigned kder=virbas + (jder - mybasemulticolvars[mmc]->getNumberOfDerivatives() + 9);
            myvals.addTemporyDerivative( kder, sw*myder.getDerivative(0, jder) );
          }
        }
        myatoms.addTemporyBoxDerivatives( (-dfunc)*values[0]*vir );
        myder.clearAll();
      }
    }
  }

  // Set the tempory weight
  updateActiveAtoms( myatoms );
  if( values.size()>2) {
    double norm=0;
    for(unsigned i=2; i<values.size(); ++i) {
      myvals.quotientRule( i, i );
      // Calculate length of vector
      norm+=myvals.get(i)*myvals.get(i);
    }
    norm=sqrt(norm); myatoms.setValue(1, norm); double inorm = 1.0 / norm;
    for(unsigned j=0; j<myvals.getNumberActive(); ++j) {
      unsigned jder=myvals.getActiveIndex(j);
      for(unsigned i=2; i<values.size(); ++i) {
        myvals.addDerivative( 1, jder, myvals.get(i)*inorm*myvals.getDerivative(i,jder) );
      }
    }
  } else {
    myvals.quotientRule( 1, 1 );
  }

  return myatoms.getValue(1);
}

}
}
