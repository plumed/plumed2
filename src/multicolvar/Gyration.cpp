/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2018 The plumed team
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
#include "core/ActionWithValue.h"
#include "core/ActionAtomistic.h"
#include "core/ActionShortcut.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace multicolvar {

//+PLUMEDOC COLVAR GYRATION
/*
Calculate the radius of gyration, or other properties related to it.

The different properties can be calculated and selected by the TYPE keyword:
the Radius of Gyration (RADIUS); the Trace of the Gyration Tensor (TRACE);
the Largest Principal Moment of the Gyration Tensor (GTPC_1); the middle Principal Moment of the Gyration Tensor (GTPC_2);
the Smallest Principal Moment of the Gyration Tensor (GTPC_3); the Asphericiry (ASPHERICITY); the Acylindricity (ACYLINDRICITY);
the Relative Shape Anisotropy (KAPPA2); the Smallest Principal Radius Of Gyration (GYRATION_3);
the Middle Principal Radius of Gyration (GYRATION_2); the Largest Principal Radius of Gyration (GYRATION_1).
A derivation of all these different variants can be found in \cite Vymetal:2011gv

The radius of gyration is calculated using:

\f[
s_{\rm Gyr}=\Big ( \frac{\sum_i^{n}
 m_i \vert {r}_i -{r}_{\rm COM} \vert ^2 }{\sum_i^{n} m_i} \Big)^{1/2}
\f]

with the position of the center of mass \f${r}_{\rm COM}\f$ given by:

\f[
{r}_{\rm COM}=\frac{\sum_i^{n} {r}_i\ m_i }{\sum_i^{n} m_i}
\f]

The radius of gyration usually makes sense when atoms used for the calculation
are all part of the same molecule.
When running with periodic boundary conditions, the atoms should be
in the proper periodic image. This is done automatically since PLUMED 2.2,
by considering the ordered list of atoms and rebuilding PBCs with a procedure
that is equivalent to that done in \ref WHOLEMOLECULES . Notice that
rebuilding is local to this action. This is different from \ref WHOLEMOLECULES
which actually modifies the coordinates stored in PLUMED.

In case you want to recover the old behavior you should use the NOPBC flag.
In that case you need to take care that atoms are in the correct
periodic image.


\par Examples

The following input tells plumed to print the radius of gyration of the
chain containing atoms 10 to 20.
\plumedfile
GYRATION TYPE=RADIUS ATOMS=10-20 LABEL=rg
PRINT ARG=rg STRIDE=1 FILE=colvar
\endplumedfile

*/
//+ENDPLUMEDOC

class GyrationTensor : 
public ActionAtomistic,
public ActionWithValue
{
private:
  bool nopbc;
  Value* val_weights;
  bool weight_mass, weight_charge;
  std::vector<double> weights;
  std::vector<double> forcesToApply;
  double getWeight( const unsigned& task_index ) const ;
public:
  static void registerKeywords(Keywords& keys);
  explicit GyrationTensor(const ActionOptions&);
  unsigned getNumberOfDerivatives() const { return 0; }
  virtual void calculate();
  void apply();
};

PLUMED_REGISTER_ACTION(GyrationTensor,"GYRATION_TENSOR")

void GyrationTensor::registerKeywords(Keywords& keys) {
  Action::registerKeywords( keys ); ActionAtomistic::registerKeywords( keys ); ActionWithValue::registerKeywords( keys );
  keys.remove("NUMERICAL_DERIVATIVES");
  keys.add("atoms","ATOMS","the group of atoms that you are calculating the Gyration Tensor for");
  keys.add("atoms","CENTER","the position to use for the center of the gyration tensor");
  keys.add("optional","WEIGHTS","what weights should be used when calculating the center.  If this keyword is not present the geometric center is computed. "
           "If WEIGHTS=@masses is used the center of mass is computed.  If WEIGHTS=@charges the center of charge is computed.  If "
           "the label of an action is provided PLUMED assumes that that action calculates a list of symmetry functions that can be used "
           "as weights. Lastly, an explicit list of numbers to use as weights can be provided");
  keys.addFlag("MASS",false,"calculate the center of mass");
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
}

GyrationTensor::GyrationTensor(const ActionOptions&ao):
  Action(ao),
  ActionAtomistic(ao),
  ActionWithValue(ao),
  nopbc(false),
  forcesToApply(9)
{
  std::vector<AtomNumber> atoms,catom;
  parseAtomList("ATOMS",atoms);
  if(atoms.size()==0) error("no atoms specified");
  parseAtomList("CENTER",catom);
  if(catom.size()!=1) error("should be one central atom only");
  parseFlag("NOPBC",nopbc);
  //copied this bit from ../vatom/Center.cpp  :: Begin 
  bool usemass = false; parseFlag("MASS",usemass);
  if(atoms.size()==0) error("at least one atom should be specified");
  std::vector<std::string> str_weights; parseVector("WEIGHTS",str_weights);
  if( usemass ) {
      if( str_weights.size()>0 ) error("USEMASS is incompatible with WEIGHTS");
      str_weights.resize(1); str_weights[0]="@masses";
  }
  if( str_weights.size()==0 ) {
    log<<"  computing the geometric center of atoms:\n";
    weights.resize( atoms.size() );
    for(unsigned i=0; i<atoms.size(); i++) weights[i] = 1.;
  } else if( str_weights.size()==1 ) {
    if( str_weights[0]=="@masses" ) {
      weight_mass=true;
      log<<"  computing the center of mass of atoms:\n";
    } else if( str_weights[0]=="@charges" ) {
      weight_charge=true;
      log<<"  computing the center of charge of atoms:\n";
    } else {
      std::size_t dot=str_weights[0].find_first_of("."); unsigned nargs=0; std::vector<Value*> args;
      if( dot!=std::string::npos ) {
        ActionWithValue* action=plumed.getActionSet().selectWithLabel<ActionWithValue*>( str_weights[0].substr(0,dot) );
        if( !action ) {
          std::string str=" (hint! the actions in this ActionSet are: ";
          str+=plumed.getActionSet().getLabelList()+")";
          error("cannot find action named " + str_weights[0] +str);
        }
        action->interpretDataLabel( str_weights[0], this, nargs, args );
      } else {
        ActionWithValue* action=plumed.getActionSet().selectWithLabel<ActionWithValue*>( str_weights[0] );
        if( !action ) {
          std::string str=" (hint! the actions in this ActionSet are: ";
          str+=plumed.getActionSet().getLabelList()+")";
          error("cannot find action named " + str_weights[0] +str);
        }
        action->interpretDataLabel( str_weights[0], this, nargs, args );
      }
    }   
  } else {
    log<<" with weights:";
    if( str_weights.size()!=atoms.size() ) error("number of elements in weight vector does not match the number of atoms");
    weights.resize( atoms.size() );
    for(unsigned i=0; i<weights.size(); ++i) {
      if(i%25==0) log<<"\n";
      Tools::convert( str_weights[i], weights[i] ); log.printf(" %f",weights[i]);
    }
    log.printf("\n");
  }
  //copied this bit from ../vatom/Center.cpp  :: End

  checkRead();

  log<<"  Bibliography "<<plumed.cite("Jirí Vymetal and Jirí Vondrasek, J. Phys. Chem. A 115, 11455 (2011)");
  log<<"\n";

  log.printf("  atoms involved : ");
  for(unsigned i=0; i<atoms.size(); ++i) {
    if(i%25==0) log<<"\n";
    log.printf("%d ",atoms[i].serial());
  }
  log.printf("\n");
  log.printf("  central atom : %d \n", catom[0].serial() );

  if(nopbc) {
    log<<"  PBC will be ignored\n";
  } else {
    log<<"  broken molecules will be rebuilt assuming atoms are in the proper order\n";
  }

  std::vector<unsigned> shape(2); shape[0]=shape[1]=3;
  addValue(shape); setNotPeriodic(); 
  atoms.push_back(catom[0]); requestAtoms(atoms); 
}

double GyrationTensor::getWeight( const unsigned& task_index ) const {
  if( weight_mass ) return getMass(task_index);
  if( weight_charge ) {
      if( !plumed.getAtoms().chargesWereSet() ) plumed_merror("cannot calculate center of charge if chrages are unset");
      return getCharge(task_index);
  }
  if( val_weights ) return val_weights->get(task_index);
  plumed_dbg_assert( task_index<weights.size() );
  return weights[task_index];
}

void GyrationTensor::calculate() {

  if(!nopbc) makeWhole();

  Tensor3d gyr_tens; double totmass = 0;
  //calculate gyration tensor
  for(unsigned i=0; i<getNumberOfAtoms()-1; i++){ 
    double w = getWeight(i);
    const Vector diff=delta( getPosition(getNumberOfAtoms()-1), getPosition(i) );
    gyr_tens[0][0]+= w*diff[0]*diff[0];
    gyr_tens[0][1]+= w*diff[0]*diff[1];
    gyr_tens[0][2]+= w*diff[0]*diff[2];
    gyr_tens[1][1]+= w*diff[1]*diff[1];
    gyr_tens[1][2]+= w*diff[1]*diff[2];
    gyr_tens[2][2]+= w*diff[2]*diff[2];
    totmass += w ;
  }

  Value* myval=getPntrToOutput(0); 
  myval->set(0, gyr_tens[0][0] / totmass );
  myval->set(1, gyr_tens[0][1] / totmass );
  myval->set(2, gyr_tens[0][2] / totmass );
  myval->set(3, gyr_tens[0][1] / totmass );
  myval->set(4, gyr_tens[1][1] / totmass );
  myval->set(5, gyr_tens[1][2] / totmass );
  myval->set(6, gyr_tens[0][2] / totmass );
  myval->set(7, gyr_tens[1][2] / totmass );
  myval->set(8, gyr_tens[2][2] / totmass );
}

void GyrationTensor::apply(){
  if( doNotCalculateDerivatives() ) return ;
  
  // Retrieve the forces from the values
  double totmass=0; for(unsigned i=0; i<getNumberOfAtoms()-1; i++) totmass += getWeight(i);
  for(unsigned i=0;i<9;++i) forcesToApply[i] = getPntrToOutput(0)->getForce( i ) / totmass;

  unsigned n=getNumberOfAtoms()-1; Vector ff;
  std::vector<Vector>& f(modifyForces());
  Tensor&              v(modifyVirial());
  for(unsigned i=0; i<getNumberOfAtoms()-1; i++) {
      double w = getWeight(i);
      const Vector diff=delta( getPosition(getNumberOfAtoms()-1), getPosition(i) );
      ff[0] = w*(2*forcesToApply[0]*diff[0] +   forcesToApply[1]*diff[1] +   forcesToApply[2]*diff[2]); 
      ff[1] = w*(  forcesToApply[3]*diff[0] + 2*forcesToApply[4]*diff[1] +   forcesToApply[5]*diff[2]); 
      ff[2] = w*(  forcesToApply[6]*diff[0] +   forcesToApply[7]*diff[1] + 2*forcesToApply[8]*diff[2]);
      f[i] += ff; f[n] -= ff;
      v -= Tensor(getPosition(i),ff);
  }
}

class Gyration : public ActionShortcut {
public:
    static void registerKeywords( Keywords& keys );
    explicit Gyration(const ActionOptions&);    
};

PLUMED_REGISTER_ACTION(Gyration,"GYRATION")

void Gyration::registerKeywords( Keywords& keys ) {
   ActionShortcut::registerKeywords( keys );
   keys.add("atoms","ATOMS","the group of atoms that you are calculating the Gyration Tensor for");
   keys.add("compulsory","TYPE","RADIUS","The type of calculation relative to the Gyration Tensor you want to perform");
   keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances"); 
}

Gyration::Gyration(const ActionOptions& ao):
Action(ao),
ActionShortcut(ao)
{
    std::string atoms; parse("ATOMS",atoms); bool nopbc; parseFlag("NOPBC",nopbc); 
    std::string pbcstr; if(nopbc) pbcstr = " NOPBC";
    // Create the geometric center of the molecule
    readInputLine( getShortcutLabel() + "_cent: CENTER ATOMS=" + atoms + pbcstr );
    // Now compute the gyration tensor
    readInputLine( getShortcutLabel() + "_tensor: GYRATION_TENSOR ATOMS=" + atoms + pbcstr + " CENTER=" + getShortcutLabel() + "_cent");
    std::string gtype; parse("TYPE",gtype);
    if( gtype=="RADIUS") {
        // And now we need the average trace for the gyration radius
        readInputLine( getShortcutLabel() + "_trace: COMBINE ARG=" + getShortcutLabel() + "_tensor.1.1," + 
            	   getShortcutLabel() + "_tensor.2.2," + getShortcutLabel() + "_tensor.3.3 PERIODIC=NO"); 
        // Square root the radius
        readInputLine( getShortcutLabel() + ": MATHEVAL ARG1=" + getShortcutLabel() + "_trace FUNC=sqrt(x) PERIODIC=NO");
    } else if( gtype=="TRACE" ) {
	// Compte the trace of the gyration tensor
	readInputLine( getShortcutLabel() + ": COMBINE ARG=" + getShortcutLabel() + "_tensor.1.1," +
                   getShortcutLabel() + "_tensor.2.2," + getShortcutLabel() + "_tensor.3.3 PERIODIC=NO");
    } else {
	// Diagonalize the gyration tensor
	readInputLine( getShortcutLabel() + "_diag: DIAGONALIZE ARG=" + getShortcutLabel() + "_tensor VECTORS=all" );    
        if( gtype=="ASPHERICITY" ) {
            readInputLine( getShortcutLabel() + ": MATHEVAL ARG1=" + getShortcutLabel() + "_diag.vals-1" +
			   " ARG2=" + getShortcutLabel() + "_diag.vals-2 ARG3=" + getShortcutLabel() + "_diag.vals-3 FUNC=x-0.5*(y+z) PERIODIC=NO");
	} else if( gtype=="GTPC_1" ) {
            readInputLine( getShortcutLabel() + ": MATHEVAL ARG1=" + getShortcutLabel() + "_diag.vals-1 FUNC=sqrt(x) PERIODIC=NO");
	} else if( gtype=="GTPC_2" ) {
            readInputLine( getShortcutLabel() + ": MATHEVAL ARG1=" + getShortcutLabel() + "_diag.vals-2 FUNC=sqrt(x) PERIODIC=NO");
	} else if( gtype=="GTPC_3" ) {
            readInputLine( getShortcutLabel() + ": MATHEVAL ARG1=" + getShortcutLabel() + "_diag.vals-3 FUNC=sqrt(x) PERIODIC=NO");
	} else error( gtype + " is not a valid type for gyration radius");
    }	
}

}
}
