/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016 The plumed team
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
#include "core/ActionRegister.h"
#include "ActionWithIntegral.h"
#include "tools/Grid.h"
#include "tools/File.h"
#include <math.h>
#include <string>
#include <cstring>
#include <iostream>

namespace PLMD {
namespace gridtools {

class KLDiv : public ActionWithIntegral {
private:
  double (get_refgrid_value)(const unsigned &current);
  // pointer to the distance function used
  double (KLDiv::*diffuncPtr)(const unsigned& current, double& der) const = NULL;
  double var_diff(const unsigned& current, double& der) const ;
  double kull_leib(const unsigned& current, double& der) const ;
  // pointer to function returning correct reference value
  double (KLDiv::*getrefvaluePtr)(const unsigned& current) const = NULL;
  double getuniformvalue(const unsigned& current) const ;
  double getrefgridvalue(const unsigned& current) const ;
  Grid* RefGrid=NULL;
  bool userefgrid=false;
  bool usekldiv;
public:
  static void registerKeywords( Keywords& keys );
  explicit KLDiv(const ActionOptions&ao);
  //unsigned getNumberOfDerivatives();
  void compute( const unsigned& current, MultiValue& myvals ) const ;
  //bool isPeriodic(){ return false; }
};

PLUMED_REGISTER_ACTION(KLDiv,"KLDIV")

void KLDiv::registerKeywords( Keywords& keys ) {
  ActionWithIntegral::registerKeywords( keys );
  keys.add("optional","REFERENCE","The file in which the reference grid is stored."
           "This keyword is used if you want to compute the distance from a distribution different than a uniform one.");
  keys.addFlag("USEKLDIV",false,"Use Kullback-Leibler divergence to compute distance between distributions instead of variational difference.");
}

KLDiv::KLDiv(const ActionOptions&ao):
  Action(ao),
  ActionWithIntegral(ao),
  RefGrid(NULL)
{
  // Decide which distance metric to use
  parseFlag("USEKLDIV",usekldiv);
  if(!usekldiv) {
    diffuncPtr=&KLDiv::var_diff;
    log.printf("  Using variational difference as distance metric.\n");
  } else {
    log.printf("  Using Kullback-Leibler divergence as distance metric.\n");
    diffuncPtr=&KLDiv::kull_leib;
  }

  // Read in the reference grid
  std::string refgridfname; parse("REFERENCE",refgridfname);
  if(refgridfname.length()>0) {
    userefgrid=true;
    log.printf("  Computing distance of distribution from input reference one.\n");
    log.printf("  Reading reference grid from file %s\n", refgridfname.c_str());
    // read in refgrid
    IFile refgridfile;
    refgridfile.link(*this);
    if(refgridfile.FileExist(refgridfname)) {
      refgridfile.open(refgridfname);
    } else {
      error("The GRID file you want to read: " + refgridfname + ", cannot be found!");
    }
    // retrieve input parameters and create vector of pointers to value to initialize refgrid
    std::vector<Value*> ingridargs;
    for(unsigned i=0; i<ingrid->getDimension(); i++) {
      Value* p = new Value(NULL, ingrid->getComponentName(i), false);
      ingridargs.push_back(p);
      if(ingrid->isPeriodic(i)) {
        ingridargs[i]->setDomain( ingrid->getMin()[i], ingrid->getMax()[i] );
      } else {
        ingridargs[i]->setNotPeriodic();
      }
    }
    RefGrid=Grid::create(ingrid->getComponentName(ingrid->getDimension()), ingridargs, refgridfile, ingrid->getMin(), ingrid->getMax(), ingrid->getNbin(), false, false, false);
    refgridfile.close();
    // delete Value objects and the vector of pointers
    //for (std::vector< Value >::iterator it = ingridargs.begin() ; it != ingridargs.end(); ++it){
    for (unsigned i =0; i< ingridargs.size(); i++) {
      delete (ingridargs[i]);
    }
    ingridargs.clear();
    //log.printf("  KLdiv: trying to get ingridargs Value name after deletion: %s\n", ingridargs[0]->getName().c_str());

    // Initialize pointer to function that returns value of RefGrid on a given point
    getrefvaluePtr=&KLDiv::getrefgridvalue;
    // If no reference grid is provided reference distribution is the uniform one
  } else {
    log.printf("  Computing distance of distribution from the uniform one.\n");
    //Initialize pointer to function that returns 1./(volume*ingrid->getNumberOfPoints())
    getrefvaluePtr=&KLDiv::getuniformvalue;
  }

}


// This function will be evaluated at each grid point.  In this case we are just taking the
// value of the function at the grid point, getFunctionValue( current ), and multiplying it
// by the volume of a grid cell (we could perhaps use some smarter method for calculating the
// integral here but it is going to make the differentiation harder.  In addition, this way of
// calculating integrals is easy to generalise to arbitrary dimensions).

// This function is called a parallel loop.  The parallel loop (in ActionWithVessel) looks after
// summing all the values (and derivatives) calculated by this function.  If you look here a MultiValue
// Object is passed to the compute method.  This method contains a N dimensional vector of values and a
// N by M matrix of derivatives (M here is the number of derivatives - in this case the number of grid points).
// What is summed in ActionWithVessel is the first component of this vector multiplied by the second component.
// (originally I conceived this bit of code for calculated weighted averages so the first number was a weight
// while the second was the quantity to be averaged).  This is why I do myvals.setValue( 0, 1.0 ).  What you need to
// do is modify the code below so that when the quanities are summed I get the KL divergence or whatever it is you need.
void KLDiv::compute( const unsigned& current, MultiValue& myvals ) const {
  myvals.setValue( 0, 1.0 );
  double value,deriv;
  value = (this->*diffuncPtr)(current, deriv);
  myvals.setValue( 1, value );
  if( !doNotCalculateDerivatives() ) myvals.addDerivative( 1, current, deriv );
}

// variational difference distribution distance
double KLDiv::var_diff( const unsigned& current, double& der) const {
  // Implementation of L1 norm
  //double val = getVolume()*( getFunctionValue( current ) - (this->*getrefvaluePtr)( current ) );
  //int sign =val<0.0?-1:val>0.0;
  //der = getVolume()*sign;
  //return std::abs( val );
  // Implementation of L2 norm
  double val = getFunctionValue( current ) - (this->*getrefvaluePtr)( current ) ;
  der = getVolume() * 2 * val;
  return getVolume() * val * val;
}

// Kullback-Leibler divergence distribution distance
double KLDiv::kull_leib( const unsigned& current, double& der) const {
  double refval = (this->*getrefvaluePtr)( current );
  double pointval = getFunctionValue( current );
  // Correct kl-div definition would be log(getFunctionValue( current )/r) the average used below is
  // a common procedure to regularize the behaviour where getFunctionValue( current )==0
  // Similarly p*log(p) -> 0 when p->0 so
  // sum should be performed only when r>0 (probably also a tolerance on r should be added?)
  // Alternatively one could make kl-div only if the reference is the uniform distribution
  der = - getVolume()/(pointval + refval);
  return getVolume()*( std::log( refval/(0.5*( pointval + refval ))));
}

// returns the current value of the reference
double KLDiv::getrefgridvalue(const unsigned& current) const {
  return RefGrid->getValue(current);
}

//Return the value of a uniform distribution on the grid.
//The unsigned input var is dummy. It is necessary because getrefvaluePtr can point
//either to this or to getrefgridvalue which needs the grid point as input.
//This could be probably done more elegantly using functors
double KLDiv::getuniformvalue(const unsigned& current) const {
  //return the value of the uniform distribution on the grid
  return 1./( ingrid->getNumberOfPoints() * getVolume());
}


}
}

