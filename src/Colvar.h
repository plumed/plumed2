#ifndef __PLUMED_Colvar_h
#define __PLUMED_Colvar_h

#include "ActionAtomistic.h"
#include "ActionWithValue.h"
#include <vector>

#define PLUMED_COLVAR_INIT(ao) Action(ao),Colvar(ao)

namespace PLMD {

/**
\ingroup INHERIT
This is the abstract base class to use for implementing new collective variables, within it there is 
information as to how to go about implementing a new CV. 

To implement a CV one you need to create a single cpp file called ColvarNAME.cpp. If you use the following template for this file then the manual and the calls to the CV will be looked after automatically.

\verbatim
#include "Colvar.h"
#include "ActionRegister.h"

#include <string>
#include <cmath>
#include <cassert>

using namespace std;

namespace PLMD{

/+PLUMEDOC COLVAR NAME
/**
\endverbatim

At this point you provide the description of your CV that will appear in the manual along with an description of the input file syntax and an example.  Merging new features of the code into the plumed main branch without proper documentation is punishable by death!  Some instructions as to how to format this information is provided here: \ref usingDoxygen

\verbatim
*/
//+ENDPLUMEDOC

/**** We begin by declaring a class for your colvar.  This class inherits everything from the Colvar class.
      This ensures it has a label, a place to store its value, places to the store the values of the derivatives
      and that it can access the various atoms it will employ.

class ColvarNAME : public Colvar {
\endverbatim

Insert declarations for your colvar's parameters here using plumed's \ref parsing.

\verbatim
public:
 /---- This is the constructor for your colvar.  It is this routine that will do all the reading.
       Hence it takes as input a line from the input file.
  ColvarNAME(const ActionOptions&);
 /---- This is the routine that will be used to calculate the value of the colvar, whenever its calculation is required.
       This routine and the constructor above must be present - if either of them are not the code will not compile.
  virtual void calculate();
};

 /------ The following command inserts your new colvar into plumed by inserting calls to your new
        routines into the parts of plumed where they are required.  This macro takes two arguments:
        The first is the name of your ColvarClass and the second is the keyword for your CV
        (the first word in the input line for your CV).
PLUMED_REGISTER_ACTION(ColvarNAME,"KEYWORD")

 /---- We now write the actual readin (constructor) and calculations routines for the colvar

ColvarNAME::ColvarNAME(const ActionOptions&ao):
 /------ This line sets up various things in the plumed core which colvars rely on.
PLUMED_COLVAR_INIT(ao)
{
 vector<int> atoms;  /----- You almost always have atoms -----/
\endverbatim

Insert code here to read the arguments of the CV here using plumed's parsing functionality.  N.B.  The label is read in already elsewhere.

\verbatim
  checkRead();     /--- This command checks that everything on the input line has been read properly

/--- The following two lines inform the plumed core that we require space to store the value
     of the CV and that the CV will act on a particular list of atoms.
  addValueWithDerivatives("");
  requestAtoms(atoms);

/ --- For a number of the free energy methods in plumed it is necessary to calculate the
      distance between two points in CV space.  Obviously, for periodic CVs one must take
      periodicities into account when calculating distances and use the minimum image
      convention in distance calculations.  Hence, we set the periodicity of the cv using
      the following two lines.
   getValue("")->setPeridodicity(true);  // Set this true if the CV is periodic otherwise set if false.
   getValue("")->setDomain(min,max);     // This routine is only required if the function is periodic.  It sets the minimum and maximum values of the colvar.
}

void ColvarNAME::calculate(){
/--- These are the things you must calculate for any cv ---/
  double cv_val;              /--- The value of the cv ----/
  Tensor boxDerivatives;      /--- The derivative of the cv with respect to the box vectors ----/
  vector<double> derivatives; /--- The derivative of the cv with respect to the atom positions ---/
\endverbatim

Insert the code to calculate your cv, its derivatives and its contribution to the virial here.  N.B. Please use, where possible, the plumed's set of CV calculation tools.

\verbatim
/---- Having calculated the cv, its derivative and the contribution to the virial you now
      transfer this information to the plumed core using the following three commands. 
  for(int i=0;i<derivatives.size();i++){ setAtomsDerivatives(i,derivatives[i]); }
  setBoxDerivatives(boxDerivatives);
  setValue(cv_val);
}
\endverbatim

\section multicvs Mult-component CVs

To avoid code duplication, and in some cases computational expense, plumed has functionality so that a single line in input can calculate be used to calculate multiple components for a CV.  For example, PATH computes the distance along the path,\f$s\f$, and the distance from the path, \f$z\f$.  Alternatively, a distance can give one the \f$x\f$, \f$y\f$ and \f$z\f$ components of the vector connecting the two atoms.  You can make use of this functionality in your own CVs as follows:

- In the constructor we create an additional value for the CV by adding the call PLMD::addValueWithDerivative("new") as well as PLMD::addValueWithDerivatives(””).  In addtion set any periodicity for our component using getValue("new")->setPeridicity() and getValue("new")->setDomain(min,max).  If this CV is called plum in our input file we can now use both plum and plum.new in any of the functions/methods in plumed.
- Obviously in calculate we now must provide functionality to calculate the values, boxDerivatives and the atom derivatives for both our original plum and its component plum.new. Furthermore, all of this data must be transferred to the plumed core.  This is done by the following code:

Here we transfer the value, box derivatives and atomic derivatives for plum.
\verbatim
for(int i=0;i<derivatives.size();i++){ setAtomsDerivatives(i,derivatives[i]); }
setBoxDerivatives(boxDerivatives);
setValue(cv_val);
\endverbatim
Here we transfer the value, box derivatives and atomic derivatives for plum.new.
\verbatim
Value* nvalue=getValue("new");
for(int i=0;i<nderivatives.size();i++){ setAtomsDerivatives(nvalue i,nderivatives[i]); }
setBoxDerivatives(nvalue,nboxDerivatives);
setValue(nvalue,ncv_val);
\endverbatim

Please only use this functionality for CVs that are VERY similar.

*/

class Colvar :
  public ActionAtomistic,
  public ActionWithValue
  {
private:
/// This is used by apply to retrive the forces on the atoms
  std::vector<double> forces;
protected:
  bool isEnergy;
  void requestAtoms(const std::vector<AtomNumber> & a);
// Set the derivatives for a particular atom equal to the input Vector
// This routine is called setAtomsDerivatives because not because you
// are setting the derivative of many atoms but because you are setting
// the derivatives of a particular atom.  The s is an apostrophe s 
// but you can't put apostrophes in function names
  void           setAtomsDerivatives(int,const Vector&);
  void           setAtomsDerivatives(Value*,int,const Vector&);
  void           setBoxDerivatives(const Tensor&);
  void           setBoxDerivatives(Value*,const Tensor&);
  const Tensor & getBoxDerivatives()const;
  const double & getForce()const;
  void apply();
public:
  bool checkIsEnergy(){return isEnergy;};
  Colvar(const ActionOptions&);
  ~Colvar(){};
  static void registerKeywords( Keywords& keys );
};

inline
void Colvar::setAtomsDerivatives(Value*v,int i,const Vector&d){
  v->addDerivative(3*i+0,d[0]);
  v->addDerivative(3*i+1,d[1]);
  v->addDerivative(3*i+2,d[2]);
}


inline
void Colvar::setBoxDerivatives(Value* v,const Tensor&d){
  unsigned nat=getNumberOfAtoms();
  v->addDerivative(3*nat+0,d(0,0));
  v->addDerivative(3*nat+1,d(0,1));
  v->addDerivative(3*nat+2,d(0,2));
  v->addDerivative(3*nat+3,d(1,0));
  v->addDerivative(3*nat+4,d(1,1));
  v->addDerivative(3*nat+5,d(1,2));
  v->addDerivative(3*nat+6,d(2,0));
  v->addDerivative(3*nat+7,d(2,1));
  v->addDerivative(3*nat+8,d(2,2));
}

inline
void Colvar::setAtomsDerivatives(int i,const Vector&d){
  setAtomsDerivatives(getPntrToValue(),i,d);
}

inline
void Colvar::setBoxDerivatives(const Tensor&d){
  setBoxDerivatives(getPntrToValue(),d);
}

}

#endif
