/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2019 The plumed team
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
#include "Colvar.h"
#include "core/Atoms.h"
#include "core/PlumedMain.h"
#include "ActionRegister.h"
#include "tools/PDB.h"
#include "tools/RMSD.h"
#include "tools/Tools.h"
#include <memory>

using namespace std;

namespace PLMD {
namespace colvar {

class PCARMSD : public Colvar {

  std::unique_ptr<PLMD::RMSD> rmsd;
  bool squared;
  bool nopbc;
  std::vector< std::vector<Vector> > eigenvectors;
  std::vector<PDB> pdbv;
  std::vector<string> pca_names;
public:
  explicit PCARMSD(const ActionOptions&);
  void calculate() override;
  static void registerKeywords(Keywords& keys);
};


using namespace std;

//+PLUMEDOC DCOLVAR PCARMSD
/*
Calculate the PCA components ( see \cite Sutto:2010 and \cite spiwok )  for a number of provided eigenvectors and an average structure. Performs optimal alignment at every step and reports the rmsd so you know if you are far or close from the average structure.
It takes the average structure and eigenvectors in form of a pdb.
Note that beta and occupancy values in the pdb are neglected and all the weights are placed to 1 (differently from the RMSD colvar for example)

\par Examples

\plumedfile
PCARMSD AVERAGE=file.pdb EIGENVECTORS=eigenvectors.pdb
\endplumedfile

The input is taken so to be compatible with the output you get from g_covar utility of gromacs (suitably adapted to have a pdb input format).
The reference configuration (file.pdb) will thus be in a file that looks something like this:

\auxfile{file.pdb}
TITLE     Average structure
MODEL        1
ATOM      1  CL  ALA     1       1.042  -3.070   0.946  1.00  0.00
ATOM      5  CLP ALA     1       0.416  -2.033   0.132  1.00  0.00
ATOM      6  OL  ALA     1       0.415  -2.082  -0.976  1.00  0.00
ATOM      7  NL  ALA     1      -0.134  -1.045   0.677  1.00  0.00
ATOM      9  CA  ALA     1      -0.774   0.053   0.003  1.00  0.00
TER
ENDMDL
\endauxfile

while the eigenvectors will be in a pdb file (eigenvectors.pdb) that looks something like this:

\auxfile{eigenvectors.pdb}
TITLE     frame t= -1.000
MODEL        1
ATOM      1  CL  ALA     1       1.194  -2.988   0.724  1.00  0.00
ATOM      5  CLP ALA     1      -0.996   0.042   0.144  1.00  0.00
ATOM      6  OL  ALA     1      -1.246  -0.178  -0.886  1.00  0.00
ATOM      7  NL  ALA     1      -2.296   0.272   0.934  1.00  0.00
ATOM      9  CA  ALA     1      -0.436   2.292   0.814  1.00  0.00
TER
ENDMDL
TITLE     frame t= 0.000
MODEL        1
ATOM      1  CL  ALA     1       1.042  -3.070   0.946  1.00  0.00
ATOM      5  CLP ALA     1      -0.774   0.053   0.003  1.00  0.00
ATOM      6  OL  ALA     1      -0.849  -0.166  -1.034  1.00  0.00
ATOM      7  NL  ALA     1      -2.176   0.260   0.563  1.00  0.00
ATOM      9  CA  ALA     1       0.314   1.825   0.962  1.00  0.00
TER
ENDMDL

\endauxfile

*/
//+ENDPLUMEDOC

PLUMED_REGISTER_ACTION(PCARMSD,"PCARMSD")

void PCARMSD::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords(keys);
  keys.add("compulsory","AVERAGE","a file in pdb format containing the reference structure and the atoms involved in the CV.");
  keys.add("compulsory","EIGENVECTORS","a file in pdb format containing the reference structure and the atoms involved in the CV.");
  keys.addOutputComponent("eig","default","the projections on each eigenvalue are stored on values labeled eig-1, eig-2, ...");
  keys.addOutputComponent("residual","default","the distance of the present configuration from the configuration supplied as AVERAGE in terms of mean squared displacement after optimal alignment ");
  keys.addFlag("SQUARED_ROOT",false," This should be set if you want RMSD instead of mean squared displacement ");
}

PCARMSD::PCARMSD(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  squared(true),
  nopbc(false)
{
  string f_average;
  parse("AVERAGE",f_average);
  string type;
  type.assign("OPTIMAL");
  string f_eigenvectors;
  parse("EIGENVECTORS",f_eigenvectors);
  bool sq;  parseFlag("SQUARED_ROOT",sq);
  if (sq) { squared=false; }
  parseFlag("NOPBC",nopbc);
  checkRead();

  PDB pdb;

  // read everything in ang and transform to nm if we are not in natural units
  if( !pdb.read(f_average,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength()) )
    error("missing input file " + f_average );

  rmsd.reset( new RMSD() );
  bool remove_com=true;
  bool normalize_weights=true;
  // here align and displace are a simple vector of ones
  std::vector<double> align; align=pdb.getOccupancy(); for(unsigned i=0; i<align.size(); i++) {align[i]=1.;} ;
  std::vector<double> displace;  displace=pdb.getBeta(); for(unsigned i=0; i<displace.size(); i++) {displace[i]=1.;} ;
  // reset again to reimpose unifrom weights (safe to disable this)
  rmsd->set(align,displace,pdb.getPositions(),type,remove_com,normalize_weights);
  requestAtoms( pdb.getAtomNumbers() );

  addComponentWithDerivatives("residual"); componentIsNotPeriodic("residual");

  log.printf("  average from file %s\n",f_average.c_str());
  log.printf("  which contains %d atoms\n",getNumberOfAtoms());
  log.printf("  with indices : ");
  for(unsigned i=0; i<pdb.getAtomNumbers().size(); ++i) {
    if(i%25==0) log<<"\n";
    log.printf("%d ",pdb.getAtomNumbers()[i].serial());
  }
  log.printf("\n");
  log.printf("  method for alignment : %s \n",type.c_str() );
  if(nopbc) log.printf("  without periodic boundary conditions\n");
  else      log.printf("  using periodic boundary conditions\n");

  log<<"  Bibliography "<<plumed.cite("Spiwok, Lipovova and Kralova, JPCB, 111, 3073 (2007)  ");
  log<<" "<<plumed.cite( "Sutto, D'Abramo, Gervasio, JCTC, 6, 3640 (2010)");

  // now get the eigenvectors
  // open the file
  FILE* fp=fopen(f_eigenvectors.c_str(),"r");
  std::vector<AtomNumber> aaa;
  unsigned neigenvects;
  neigenvects=0;
  if (fp!=NULL)
  {
    log<<"  Opening the eigenvectors file "<<f_eigenvectors.c_str()<<"\n";
    bool do_read=true;
    while (do_read) {
      PDB mypdb;
      // check the units for reading this file: how can they make sense?
      do_read=mypdb.readFromFilepointer(fp,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength());
      if(do_read) {
        neigenvects++;
        if(mypdb.getAtomNumbers().size()==0) error("number of atoms in a frame should be more than zero");
        unsigned nat=mypdb.getAtomNumbers().size();
        if(nat!=mypdb.getAtomNumbers().size()) error("frames should have the same number of atoms");
        if(aaa.empty()) aaa=mypdb.getAtomNumbers();
        if(aaa!=mypdb.getAtomNumbers()) error("frames should contain same atoms in same order");
        log<<"  Found eigenvector: "<<neigenvects<<" containing  "<<mypdb.getAtomNumbers().size()<<" atoms\n";
        pdbv.push_back(mypdb);
        eigenvectors.push_back(mypdb.getPositions());
      } else {break ;}
    }
    fclose (fp);
    log<<"  Found total "<<neigenvects<< " eigenvectors in the file "<<f_eigenvectors.c_str()<<" \n";
    if(neigenvects==0) error("at least one eigenvector is expected");
  }
  // the components
  for(unsigned i=0; i<neigenvects; i++) {
    std::string num; Tools::convert( i, num );
    string name; name=string("eig-")+num;
    pca_names.push_back(name);
    addComponentWithDerivatives(name); componentIsNotPeriodic(name);
  }
  turnOnDerivatives();

}

// calculator
void PCARMSD::calculate() {
  if(!nopbc) makeWhole();
  Tensor rotation,invrotation;
  Matrix<std::vector<Vector> > drotdpos(3,3);
  std::vector<Vector> alignedpos;
  std::vector<Vector> centeredpos;
  std::vector<Vector> centeredref;
  std::vector<Vector> ddistdpos;
  double r=rmsd->calc_PCAelements( getPositions(), ddistdpos, rotation,  drotdpos, alignedpos,centeredpos, centeredref,squared);
  invrotation=rotation.transpose();

  Value* verr=getPntrToComponent("residual");
  verr->set(r);
  for(unsigned iat=0; iat<getNumberOfAtoms(); iat++) {
    setAtomsDerivatives (verr,iat,ddistdpos[iat]);
  }

  std::vector< Vector > der;
  der.resize(getNumberOfAtoms());


  for(unsigned i=0; i<eigenvectors.size(); i++) {
    Value* value=getPntrToComponent(pca_names[i].c_str());
    double val; val=0.;
    for(unsigned iat=0; iat<getNumberOfAtoms(); iat++) {
      val+=dotProduct(alignedpos[iat]-centeredref[iat],eigenvectors[i][iat]);	der[iat].zero();
    }
    value->set(val);
    // here the loop is reversed to better suit the structure of the derivative of the rotation matrix
    double tmp1;
    for(unsigned a=0; a<3; a++) {
      for(unsigned b=0; b<3; b++) {
        tmp1=0.;
        for(unsigned n=0; n<getNumberOfAtoms(); n++) {
          tmp1+=centeredpos[n][b]*eigenvectors[i][n][a];
        }
        for(unsigned iat=0; iat<getNumberOfAtoms(); iat++) {
          der[iat]+=drotdpos[a][b][iat]*tmp1;
        }
      }
    }
    Vector v1;
    for(unsigned n=0; n<getNumberOfAtoms(); n++) {
      v1+=(1./getNumberOfAtoms())*matmul(invrotation,eigenvectors[i][n]);
    }
    for(unsigned iat=0; iat<getNumberOfAtoms(); iat++) {
      der[iat]+=matmul(invrotation,eigenvectors[i][iat])-v1;
      setAtomsDerivatives (value,iat,der[iat]);
    }
  }

  for(int i=0; i<getNumberOfComponents(); ++i) setBoxDerivativesNoPbc( getPntrToComponent(i) );

}

}
}



