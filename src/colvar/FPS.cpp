#include "Colvar.h"
#include "ActionRegister.h"
#include <string>
#include <cmath>
#include <cassert>
#include <vector>
#include "tools/Tools.h"
#include "tools/PDB.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include <iostream>
#include "tools/RMSD.h"

using namespace std;

namespace PLMD{
namespace colvar{

//+PLUMEDOC COLVAR FPS
/*
This colvar evaluates the position of an atom with respect to a given line.

By default FPS is computed taking into account periodic
boundary conditions. This behavior can be changed with the NOPBC flag.
FPS is divided in 2 components (linepos and linedist) which evaluate
the distance from a given point and the altitude from the line, in order to
inform the Funnel if to apply or not the potential. The point has to be the
COM of the molecule you want to move inside the funnel and it has to be ONE!

Inside the colvar it will be created a rotation matrix using a reference
given by the user. This matrix is applied for every step of metadynamics
and it will align the frame to my reference. Afterwards, linepos and linedist
will be computed and re-rotated again, using the transpose of the previous
matrix, to fit with the correct derivatives.

\par Examples

The following input tells plumed to print the FPS components for the COM of a ligand.
As input we give a reference structure, which is the structure used for alignment,
the atom of the molecule we want to track and 2 points in the cartesian space (x,y,z)
to draw an imaginary line where the funnel will be constructed.
\verbatim
ligand: COM ATOMS=2446,2447,2448,2449,2451
fps: FPS REFERENCE=protein.pdb LIGAND=ligand POINTS= 5.3478,-0.7278,2.4746,7.3785,6.7364,-9.3624
PRINT ARG=fps.lp,fps.ld
\endverbatim
*/
//+ENDPLUMEDOC


class FPS : public Colvar {

  // Those are arrays that I created for the colvar
  std::vector<AtomNumber> ligand_com;
  std::vector<AtomNumber> anchor;
  std::vector<AtomNumber> numbers;
//  bool components;
  bool pbc;
  PLMD::RMSD* alignment;
  PLMD::PDB* pdb;
  bool squared;
private:
  // Placing it in private to not risk to call something with the same name somewhere else
  vector<double> points;
public:
  explicit FPS(const ActionOptions&);
// active methods:
  virtual void calculate();
  static void registerKeywords(Keywords& keys);
  std::vector<double> align;
  std::vector<double> displace;
  ~FPS();
};

using namespace std;

PLUMED_REGISTER_ACTION(FPS,"FPS")

void FPS::registerKeywords(Keywords& keys){
        Colvar::registerKeywords( keys );
        keys.add("compulsory","REFERENCE","a file in pdb format containing the structure you would like to align.");
        keys.add("atoms","LIGAND","This MUST be a single atom, normally the COM of the ligand");
        keys.add("atoms","ANCHOR","Atom picked to maintain the correct ligand during all the simulation");
        keys.add("compulsory","POINTS","6 values that define 2 points where we construct the line.");
        keys.addFlag("SQUARED-ROOT",false,"Maintained to use the function already implemented, but to no use");
        keys.addOutputComponent("lp","default","the position on the funnel line");
        keys.addOutputComponent("ld","default","the distance from the funnel line");
}

FPS::FPS(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
pbc(true),squared(true)
{
  string reference;
  parse("REFERENCE",reference);
  string type;
  type.assign("OPTIMAL");
  parseAtomList("LIGAND",ligand_com);
  if(ligand_com.size()!=1)
    error("Number of specified atoms should be one, normally the COM of the ligand");
  parseVector("POINTS",points);
  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;
  bool sq;  parseFlag("SQUARED-ROOT",sq);
  if (sq){ squared=false; }
  parseAtomList("ANCHOR",anchor);
  checkRead();

  pdb = new PDB();

  // read everything in ang and transform to nm if we are not in natural units
  if( !pdb->read(reference,plumed.getAtoms().usingNaturalUnits(),0.1/ActionAtomistic::atoms.getUnits().getLength()) )
      error("missing input file " + reference );

  alignment = new RMSD();

  bool remove_com=true;
  bool normalize_weights=true;
  // here displace is a simple vector of ones
  align=pdb->getOccupancy();for(unsigned i=0;i<align.size();i++){align[i]=1.;} ;
  displace=pdb->getBeta();for(unsigned i=0;i<displace.size();i++){displace[i]=1.;} ;
  // reset again to reimpose uniform weights (safe to disable this)
  alignment->set(align,displace,pdb->getPositions(),type,remove_com,normalize_weights);



  // Array with inside both the structure to align and the atom to be aligned
  numbers=pdb->getAtomNumbers();
  numbers.push_back(anchor[0]);
  numbers.push_back(ligand_com[0]);

  log.printf("  average from file %s\n",reference.c_str());
  log.printf("  method for alignment : %s \n",type.c_str() );

  if(pbc) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");

  addComponentWithDerivatives("lp"); componentIsNotPeriodic("lp");
  addComponentWithDerivatives("ld"); componentIsNotPeriodic("ld");

  requestAtoms( numbers );

}

FPS::~FPS(){
  delete alignment;
  delete pdb;
}

// calculator
void FPS::calculate(){

    if(pbc) makeWhole();

	Tensor Rotation;
    Matrix<std::vector<Vector> > drotdpos(3,3);
    std::vector<Vector> alignedpos;
    std::vector<Vector> centeredpos;
    std::vector<Vector> centeredref;
    std::vector<Vector> ddistdpos;

    std::vector<Vector> buffer;

    Vector centerreference;
    Vector centerpositions;

    // Created only to give the correct object to calc_Funnelelements
    std::vector<Vector> sourceAllPositions;
    std::vector<Vector> sourcePositions;

    // SourcePositions contains only the coordinates of the COM of the ligand, we need only this point
    sourceAllPositions=getPositions();
    sourcePositions=sourceAllPositions;
    sourcePositions.resize(sourcePositions.size()-2);// just the protein

    // The two points that define the axis : this can be moved in the initialization
    Vector p1 = VectorGeneric<3>(points[0],points[1],points[2]);
    Vector p2 = VectorGeneric<3>(points[3],points[4],points[5]);
    Vector s = p2 - p1;

    // Function slightly different from the real one which is calc_PCAelements
    // note that the rotation matrix provided is positions to reference
    // remember that the only part that needs to be provided is 
    double rmsd=alignment->calc_FitElements( sourcePositions, Rotation , drotdpos , buffer, centerpositions, squared);

    // I cancelled the additional lines in the library of RMSD.h, thus I am missing the center of the reference
    // Creating variable kito to extract only the center of the reference, since no method is calling
    // function getReferenceCenter()
    PLMD::RMSDCoreData* kito; kito = new RMSDCoreData(align,displace,sourcePositions,pdb->getPositions());
    centerreference = kito->getReferenceCenter();
    delete kito;

    // DEBUG
/*    log.printf(" RMSD: %13.6lf\n",rmsd );
    log.printf(" cpos: %13.6lf %13.6lf %13.6lf\n",centerpositions[0],centerpositions[1],centerpositions[2] );
    log.printf(" Rotation: %13.6lf %13.6lf %13.6lf\n",Rotation[0][0],Rotation[0][1],Rotation[0][2] );
    log.printf("           %13.6lf %13.6lf %13.6lf\n",Rotation[1][0],Rotation[1][1],Rotation[1][2] );
    log.printf("           %13.6lf %13.6lf %13.6lf\n",Rotation[2][0],Rotation[2][1],Rotation[2][2] );
*/

    // the position is   Rot(ligand-com_prot_md)+ com_prot_ref
    Vector ligand_centered =getPositions().back()-centerpositions;
    Vector ligand_aligned =matmul(Rotation,ligand_centered);
    Vector v = ligand_aligned +centerreference -p1;

    // DEBUG
//    log.printf(" ligando: %13.6lf %13.6lf %13.6lf\n",getPositions().back()[0],getPositions().back()[1],getPositions().back()[2] );

    //Projection vector v onto s

    Vector prj = (dotProduct(s,v)/dotProduct(s,s))*s; 
    const double prj_length = prj.modulo() ;  
    const double inv_prj_length = 1.0/prj_length;

    Vector height = v - prj;
    const double prj_height = height.modulo() ;
    const double inv_prj_height = 1.0/prj_height;

    // derivative of the prj: only on the com of the ligand
    Vector der_prj;
    der_prj=s/s.modulo();

    // derivative of the height: only on the com of the ligand
    Vector der_height(inv_prj_height*(height[0]-(s[0]/s.modulo2())*dotProduct(height,s)),
	     	          inv_prj_height*(height[1]-(s[1]/s.modulo2())*dotProduct(height,s)),
	  	              inv_prj_height*(height[2]-(s[2]/s.modulo2())*dotProduct(height,s)));

    Value* valuelp=getPntrToComponent("lp");
    Value* valueld=getPntrToComponent("ld");
    valuelp->set(dotProduct(s,v)/s.modulo()); // this includes the sign 
    valueld->set(prj_height);

    // DEBUG
//    log.printf(" Dopo: %13.6lf  %13.6lf\n",dotProduct(s,v)/s.modulo(),prj_height );

    setAtomsDerivatives(valuelp,getNumberOfAtoms()-1,matmul(transpose(Rotation),der_prj));
    setAtomsDerivatives(valueld,getNumberOfAtoms()-1,matmul(transpose(Rotation),der_height));
    
    Vector der_h;
    Vector der_l;
    double weight=1./float(getNumberOfAtoms()-2.);

    for(unsigned iat=0;iat<getNumberOfAtoms()-2;iat++){
		der_h.zero();
		der_l.zero();
		for(unsigned a=0;a<3;a++){
    			for(unsigned b=0;b<3;b++){
  			    	for(unsigned g=0;g<3;g++){
					der_h[a]+=der_height[b]*drotdpos[b][g][iat][a]*ligand_centered[g];
					der_l[a]+=der_prj[b]*drotdpos[b][g][iat][a]*ligand_centered[g];
				}	
				der_h[a]-=der_height[b]*Rotation[b][a]*weight;
				der_l[a]-=der_prj[b]*Rotation[b][a]*weight;
	        	}		
	        }		
		setAtomsDerivatives(valuelp,iat,der_l);	
		setAtomsDerivatives(valueld,iat,der_h);	
    }

    setBoxDerivativesNoPbc(valuelp);
    setBoxDerivativesNoPbc(valueld);

}

}
}



