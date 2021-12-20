#ifdef __PLUMED_HAS_RASCAL
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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
#include "core/ActionAtomistic.h"
#include "core/ActionWithValue.h"
#include "core/ActionRegister.h"

#include "rascal/representations/calculator_sorted_coulomb.hh"
#include "rascal/representations/calculator_spherical_expansion.hh"
#include "rascal/representations/calculator_spherical_invariants.hh"
#include "rascal/structure_managers/adaptor_center_contribution.hh"
#include "rascal/structure_managers/adaptor_neighbour_list.hh"
#include "rascal/structure_managers/adaptor_strict.hh"
#include "rascal/structure_managers/atomic_structure.hh"
#include "rascal/structure_managers/make_structure_manager.hh"
#include "rascal/structure_managers/structure_manager_centers.hh"
#include "rascal/utils/basic_types.hh"
#include "rascal/utils/utils.hh"

#include <chrono>
#include <cmath>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <list>
#include <string>
#include <cmath>

using namespace rascal;

namespace PLMD {
namespace rascal {

//+PLUMEDOC COLVAR SPHERICAL_INVARIANTS
/*
Interface to librascal code for computing structural descritors such as SOAP.

To use PLUMED and librascal together you must configure as follows:

\verbatim
./configure --enable-cxx=14 --enable-modules=all --enable-asmjit --enable-rascal LDFLAGS="-L<librascal>/build/external/wigxjpf/lib -lwigxjpf -L/<librascal>/build/src -lrascal" CPPFLAGS="-march=native -I/<librascal>/src/ -I/<librascal>/build/external/Eigen3/ -I/<librascal>/build/external/wigxjpf/inc"
\endverbatim

\par Examples

Here is an example input file:

\plumedfile
r: SPHERICAL_INVARIANTS ...
  SPECIES=1-3
  HYPERPARAMS={ 
                "max_radial": 3,
                "max_angular": 2,
                "compute_gradients": true,
                "soap_type": "PowerSpectrum",
                "normalize": true,
                "cutoff_function": {"type": "ShiftedCosine", "cutoff": {"value": 4.0, "unit": "AA"}, "smooth_width": {"value": 0.5, "unit": "AA"}},
                "gaussian_density": {"type": "Constant", "gaussian_sigma": {"value": 0.4, "unit": "AA"}},
                "radial_contribution": {"type": "GTO"}
              }
...

PRINT ARG=rr FILE=colvar
\endplumedfile

*/
//+ENDPLUMEDOC

template <class T>
class RascalSpherical : 
public ActionAtomistic,
public ActionWithValue 
{
private:
/// These are typedefs for RASCAL stuff
  typedef AdaptorStrict<AdaptorCenterContribution<AdaptorNeighbourList<StructureManagerCenters>>> Manager_t;
  typedef CalculatorSphericalInvariants::Property_t<Manager_t> Prop_t;
  typedef CalculatorSphericalInvariants::PropertyGradient_t<Manager_t> PropGrad_t ;
/// These are the adaptors to use within rascal
  json adaptors;
/// This is the structure object that is passed between PLUMED and librascal
  json structure;
/// Tempory vector used to store the number of neighbours in apply
  std::vector<unsigned> neigh;
/// This is a vector containing the forces that act on the system
  std::vector<double> forcesToApply;
/// This is the representation that is being computed using librascal
  std::unique_ptr<T> representation;
/// This generates the json input for RASCAL
  void structureToJson();
public:
  static void registerKeywords(Keywords& keys);
  explicit RascalSpherical(const ActionOptions&);
// active methods:
  unsigned getNumberOfDerivatives() const override;
  void calculate() override;
  void apply() override;
};

typedef RascalSpherical<CalculatorSphericalInvariants> SphericalInvariants;
PLUMED_REGISTER_ACTION(SphericalInvariants,"SPHERICAL_INVARIANTS")
typedef RascalSpherical<CalculatorSphericalExpansion> SphericalExpansion;
PLUMED_REGISTER_ACTION(SphericalExpansion,"SPHERICAL_EXPANSION")

template <class T>
void RascalSpherical<T>::registerKeywords(Keywords& keys) {
  Action::registerKeywords(keys); ActionAtomistic::registerKeywords(keys); ActionWithValue::registerKeywords(keys);
  keys.add("numbered","SPECIES","the atoms in each species type"); keys.reset_style("SPECIES","atoms");
  keys.add("compulsory","HYPERPARAMS","the json input for the librascal hyperparameters");
}

template <class T>
RascalSpherical<T>::RascalSpherical(const ActionOptions&ao):
  Action(ao),
  ActionAtomistic(ao),
  ActionWithValue(ao)
{
  // Read in the hyper parameters for the representation
  std::string adapt, hypers; parse("HYPERPARAMS",hypers);
  json hyper_params=json::parse( "{" + hypers + "}" );
  log<<"   hyper parameters : \n"<<std::setw(4)<<hyper_params<<"\n";
  double cutoff = hyper_params["cutoff_function"]["cutoff"]["value"];
  // Check that we have gradients being computed
  if( !hyper_params["compute_gradients"] ) { 
      warning("resetting compute_gradients to true as PLUMED cannot operate without gradients");
      hyper_params["compute_gradients"]=true;
  }
  if( hyper_params["cutoff_function"]["cutoff"].find("unit")!=hyper_params["cutoff_function"]["cutoff"].end() || 
      hyper_params["cutoff_function"]["smooth_width"].find("unit")!=hyper_params["cutoff_function"]["smooth_width"].end() ||
      hyper_params["gaussian_density"]["gaussian_sigma"].find("unit")!=hyper_params["gaussian_density"]["gaussian_sigma"].end() ) {
      error("remove units keywords from json input to HYPERPARAMS.  The units of length for PLUMED are nm or whatever length unit was specified in the UNIT action and cannot be changed in input to librascal");
  }
  // Create the representation using the hyper parameters 
  representation=std::unique_ptr<T>( new T{hyper_params} );
  // Now read in the adaptors for the representation
  json ad1{{"name", "AdaptorNeighbourList"},
           {"initialization_arguments", {{"cutoff", cutoff}}}};
  json ad1b{{"name", "AdaptorCenterContribution"},
            {"initialization_arguments", {}}};
  json ad2{{"name", "AdaptorStrict"},
           {"initialization_arguments", {{"cutoff", cutoff}}}};
  adaptors.emplace_back(ad1);
  adaptors.emplace_back(ad1b);
  adaptors.emplace_back(ad2);
  // Create structure object.  This will be used to pass information to rascal
  structure["cell"]={{0,0,0},{0,0,0},{0,0,0}}; structure["pbc"]={true,true,true};
  // Now read in atoms that we are using
  std::vector<AtomNumber> catoms, all_atoms; parseAtomList("SPECIES",all_atoms); unsigned nspecies=0;
  if( all_atoms.size()==0 ) {
      std::vector<AtomNumber> t; 
      for(int i=1;; ++i ) {
          parseAtomList("SPECIES", i, t );
          if( t.empty() ) break;

          log.printf("  Species %d includes atoms : ", i);
          for(unsigned j=0; j<t.size(); ++j) { 
              log.printf("%d ",t[j].serial() ); structure["numbers"].push_back(i); 
              structure["positions"].push_back( std::vector<double>(3) );
              all_atoms.push_back( t[j] ); 
          }
          log.printf("\n"); t.resize(0); nspecies++;
      }     
  } else {
      nspecies = 1;
      for(unsigned j=0; j<all_atoms.size(); ++j) {
          structure["numbers"].push_back( 1 );
          structure["positions"].push_back( std::vector<double>(3) );
      }
  }
  // Request the atoms and check we have read in everything
  requestAtoms(all_atoms); forcesToApply.resize( 3*all_atoms.size() + 9 ); neigh.resize( all_atoms.size() );
  // Setup a matrix to hold the soap vectors
  std::vector<unsigned> shape(2); shape[0]=all_atoms.size(); shape[1]=representation->get_num_coefficients( nspecies ); 
  addValue( shape ); setNotPeriodic(); getPntrToOutput(0)->alwaysStoreValues();
  checkRead();
}

template <class T>
unsigned RascalSpherical<T>::getNumberOfDerivatives() const {
  return 3*getNumberOfAtoms()+9;
}

template <class T>
void RascalSpherical<T>::structureToJson() {
  // Set the atoms from the atomic positions
  for(unsigned i=0;i<getNumberOfAtoms();++i) {
      for(unsigned k=0;k<3;++k) structure["positions"][i][k]=getPosition(i)[k];
  }
  // Set the box and the pbc from the cell 
  for(unsigned i=0; i<3; ++i) for(unsigned j=0; j<3; ++j) { structure["cell"][i][j]=getBox()(i,j); }
}

// calculator
template <class T>
void RascalSpherical<T>::calculate() {
  // tranfer strcutre to json input
  structureToJson();
  // Now lets setup the stuff for librascal
  auto manager = make_structure_manager_stack<StructureManagerCenters, AdaptorNeighbourList, AdaptorCenterContribution, AdaptorStrict>( structure, adaptors );
  // And compute everything with librascal
  representation->compute(manager); 

  // Retrieve the properties that are required
  auto && property{ *manager->template get_property<Prop_t>(representation->get_name())}; auto features = property.get_features();
  Value* valout=getPntrToOutput(0); std::vector<unsigned> shape( valout->getShape() );
  if( shape[0]!=features.rows() || shape[1]!=features.cols() ) {
     shape[0]=features.rows(); shape[1]=features.cols(); valout->setShape( shape );
  }

  for(unsigned i=0; i<shape[0]; ++i ) {
      for(unsigned j=0; j<shape[1]; ++j) valout->set( i*shape[1]+j, features(i,j) );
  }
}

template <class T>
void RascalSpherical<T>::apply() {
  // Do nothing if no forces were added
  if( !getPntrToOutput(0)->forcesWereAdded() ) return ;
  // Clear the forces from last time
  std::fill(forcesToApply.begin(),forcesToApply.end(),0);

  // Recalculate SOAPS if forces were added
  structureToJson();
  // Now lets setup the stuff for librascal
  auto manager = make_structure_manager_stack<StructureManagerCenters, AdaptorNeighbourList, AdaptorCenterContribution, AdaptorStrict>( structure, adaptors );
  // And compute everything with librascal
  representation->compute(manager); 

  // Now lets get the gradients
  auto && soap_vector_gradients{*manager->template get_property<PropGrad_t>(representation->get_gradient_name())};
  math::Matrix_t gradients = soap_vector_gradients.get_features_gradient(); 
  // This gets information on the neighbour lists from librascal
  Eigen::Matrix<int, Eigen::Dynamic, 5> ninfo = manager->get_gradients_info();
  // Determine the number of neighbours for each of the environments
  std::fill(neigh.begin(),neigh.end(),0); for(unsigned i=0;i<ninfo.rows();++i) neigh[ ninfo(i,1) ]++;

  // Apply the forces on the atoms
  Value* outval=getPntrToOutput(0); std::vector<unsigned> shape( outval->getShape() ); unsigned base=0;
  // Loop over environments (atoms)
  for(unsigned i=0; i<shape[0]; ++i) {
      // Loop over features
      for(unsigned j=0; j<shape[1];++j) {
          // Get the force on the jth feature in the ith environment.
          double ff = outval->getForce( i*shape[1] + j );
          // Loop over neighbours in ith environment
          for(unsigned k=0; k<neigh[i]; ++k) {
              // Loop over x, y, z
              for(unsigned n=0;n<3;++n) {
                  // This is the force to add to the central atom  
                  forcesToApply[ 3*i+n ] -= ff*gradients(3*(base+k)+n,j);
                  // This is the force to add to the neighbour atom.
                  forcesToApply[ 3*ninfo(base+k,2) + n ] += ff*gradients(3*(base+k)+n,j);
              }
          }
      }
      // And this is some book keeping to make sure we are navigating gradients and info correctly
      base = base + neigh[i];
  }
  unsigned start=0; setForcesOnAtoms( forcesToApply, start );
}

}
}
#endif


