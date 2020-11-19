/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2019-2020

   This file is part of funnel code module.

   The FM code respects the CC BY-NC license.
   Users are free to download, adapt and use the code as long as it is not for commercial purposes.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "colvar/Colvar.h"
#include "colvar/ActionRegister.h"
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

namespace PLMD {
namespace funnel {

//+PLUMEDOC FUNNELMOD_COLVAR FUNNEL_PS
/*
FUNNEL_PS implements the Funnel-Metadynamics (FM) technique in PLUMED 2.

Please read the FM \cite limongelli2013funnel \cite raniolo2020ligand papers to better understand the notions hereby reported.

This colvar evaluates the position of a ligand of interest with respect to a given line, built from the two points
A and B, and should be used together with the \ref FUNNEL bias.
The constructed line represents the axis of the funnel-shape restraint potential, which should be placed so
as to include the portion of a macromolecule (i.e., protein, DNA, etc.) that should be explored.
Since it is important that the position of the line is updated based on the motion of the macromolecule during
the simulation, this colvar incorporates an alignment method. The latter uses the TYPE=OPTIMAL option to remove
motions due to rotation and translation of the macromolecule with respect to a reference structure, which is
provided by the user. In order to accomplish the task, an optimal alignment matrix is calculated using the
Kearsley \cite kearsley algorithm.
The reference structure should be given as a pdb file, containing only the atoms of the macromolecule or a
selection of them (e.g., the protein CA atoms of secondary structures). In contrast to the methods reported in
the \ref dists, the values reported in the occupancy and beta-factor columns of the pdb file are not important
since they will be overwritten and replaced by the value 1.00 during the procedure. It is important to understand
that all atoms in the file will be used for the alignment, even if they display 0.00 in the occupancy column.

The ligand can be represented by one single atom or the center of mass (COM) of a group of atoms that should be
provided by the user.

By default FUNNEL_PS is computed taking into account periodic boundary conditions. Since PLUMED 2.5, molecules are
rebuilt using a procedure that is equivalent to that done in \ref WHOLEMOLECULES. We note that this action is local
to this colvar, thus it does not modify the coordinates stored in PLUMED. Moreover, FUNNEL_PS requires an ANCHOR atom
to be specified in order to facilitate the reconstruction of periodic boundary conditions. This atom should be the
closest macromolecule's atom to the ligand and it should reduce the risk of ligand "warping" in the simulation box.
Nevertheless, we highly recommend to add to the PLUMED input file a custom line of \ref WHOLEMOLECULES, in order to
be sure of reconstructing the ligand together with the macromolecule (please look the examples). In this case, the user
can use the NOPBC flag to turn off the internal periodic boundary condition reconstruction.

FUNNEL_PS is divided in two components (fps.lp and fps.ld) which evaluate the projection of the ligand along the funnel line
and the distance from it, respectively. The values attributed to these two components are then used together with the
potential file created by the \ref FUNNEL bias to define if the ligand is within or not in the funnel-shape restraint
potential. In the latter case, the potential will force the ligand to enter within the funnel boundaries.

\par Examples

The following input tells plumed to print the FUNNEL_PS components for the COM of a ligand. The inputs are a reference structure,
which is the structure used for the alignment, the COM of the molecule we want to track, and 2 points in the Cartesian space
(i.e., x, y, and z) to draw the line representing the funnel axis.
\plumedfile
lig: COM ATOMS=2446,2447,2448,2449,2451
fps: FUNNEL_PS REFERENCE=protein.pdb LIGAND=lig POINTS=5.3478,-0.7278,2.4746,7.3785,6.7364,-9.3624
PRINT ARG=fps.lp,fps.ld
\endplumedfile

It is recommended to add a line to force the reconstruction of the periodic boundary conditions. In the following example,
\ref WHOLEMOLECULES was added to make sure that the ligand was reconstructed together with the protein. The list contains
all the atoms reported in the start.pdb file followed by the ANCHOR atom and the ligand. All atoms should be contained in the
same entity and the correct order.
\plumedfile
WHOLEMOLECULES ENTITY0=54,75,212,228,239,258,311,328,348,372,383,402,421,463,487,503,519,657,674,690,714,897,914,934,953,964,
974,985,1007,1018,1037,1247,1264,1283,1302,1324,1689,1708,1727,1738,1962,1984,1994,2312,2329,2349,2467,2490,2500,2517,2524,2536,
2547,2554,2569,2575,2591,2607,2635,2657,2676,2693,2700,2719,2735,2746,2770,2777,2788,2795,2805,2815,2832,2854,2868,2898,2904,
2911,2927,2948,2962,2472,3221,3224,3225,3228,3229,3231,3233,3235,3237
lig: COM ATOMS=3221,3224,3225,3228,3229,3231,3233,3235,3237
fps: FUNNEL_PS LIGAND=lig REFERENCE=start.pdb ANCHOR=2472 POINTS=4.724,5.369,4.069,4.597,5.721,4.343
PRINT ARG=fps.lp,fps.ld
\endplumedfile

*/
//+ENDPLUMEDOC


class FUNNEL_PS : public Colvar {

  // Those are arrays that I created for the colvar
  std::vector<AtomNumber> ligand_com;
  std::vector<AtomNumber> anchor;
  std::vector<AtomNumber> numbers;
  bool pbc;
  PLMD::RMSD* alignment;
  PLMD::PDB* pdb;
  bool squared;
private:
  vector<double> points;
public:
  explicit FUNNEL_PS(const ActionOptions&);
// active methods:
  virtual void calculate();
  static void registerKeywords(Keywords& keys);
// I need a method in RMSDCoreCalc and these were requested
  std::vector<double> align;
  std::vector<double> displace;
// It is written no more desctructors, but an expert said it's necessary for imported variables (pdb and alignment) or else memory leak
  ~FUNNEL_PS();
};

using namespace std;

PLUMED_REGISTER_ACTION(FUNNEL_PS,"FUNNEL_PS")

void FUNNEL_PS::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords( keys );
  keys.add("compulsory","REFERENCE","a file in pdb format containing the structure you would like to align.");
  keys.add("atoms","LIGAND","This MUST be a single atom, normally the COM of the ligand");
  keys.add("atoms","ANCHOR","Closest protein atom to the ligand, picked to avoid pbc problems during the simulation");
  keys.add("compulsory","POINTS","6 values defining x, y, and z of the 2 points used to construct the line. The order should be A_x,A_y,A_z,B_x,B_y,B_z.");
  keys.addFlag("SQUARED-ROOT",false,"Used to initialize the creation of the alignment variable");
  keys.addOutputComponent("lp","default","the position along the funnel line");
  keys.addOutputComponent("ld","default","the distance from the funnel line");
}

FUNNEL_PS::FUNNEL_PS(const ActionOptions&ao):
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
  bool sq;
  parseFlag("SQUARED-ROOT",sq);
  if (sq) {
    squared=false;
  }
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
  align=pdb->getOccupancy();
  for(unsigned i=0; i<align.size(); i++) {
    align[i]=1.;
  } ;
  displace=pdb->getBeta();
  for(unsigned i=0; i<displace.size(); i++) {
    displace[i]=1.;
  } ;
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

  addComponentWithDerivatives("lp");
  componentIsNotPeriodic("lp");
  addComponentWithDerivatives("ld");
  componentIsNotPeriodic("ld");

  requestAtoms( numbers );

}

FUNNEL_PS::~FUNNEL_PS() {
  delete alignment;
  delete pdb;
}

// calculator
void FUNNEL_PS::calculate() {

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

  // Created only to give the correct object to calc_FitElements
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

  // I call the method calc_FitElements that initializes all feature that I need
  // except for centerreference that I need to calculate from scratch
  // Buffer has no meaning but I had to fulfill the requirements of calc_FitElements
  double rmsd = alignment->calc_FitElements( sourcePositions, Rotation, drotdpos, buffer, centerpositions, squared);

  // To Plumed developers: it would be interesting to make the functions to calculate centers of mass public or protected
  centerreference.zero();
  for(unsigned i=0; i<pdb->size(); i++) {
    centerreference+=pdb->getPositions()[i]*align[i]/align.size();
  }

  /*
  // I cancelled the additional lines in the library of RMSD.h, thus I am missing the center of the reference
  // Creating variable kito to extract only the center of the reference, since no method is calling
  // function getReferenceCenter()
  PLMD::RMSDCoreData* kito; kito = new RMSDCoreData(align,displace,sourcePositions,pdb->getPositions());
  centerreference = kito->getReferenceCenter();
  delete kito;
  */

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

  for(unsigned iat=0; iat<getNumberOfAtoms()-2; iat++) {
    der_h.zero();
    der_l.zero();
    for(unsigned a=0; a<3; a++) {
      for(unsigned b=0; b<3; b++) {
        for(unsigned g=0; g<3; g++) {
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



