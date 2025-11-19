/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2023 The plumed team
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
#include "core/PlumedMain.h"
#include "core/ActionRegister.h"
#include "tools/PDB.h"
#include "tools/RMSD.h"
#include "tools/Tools.h"

namespace PLMD {
namespace colvar {

class PCARMSD : public Colvar {

  std::unique_ptr<PLMD::RMSD> rmsd;
  bool squared;
  bool nopbc;
  std::vector< std::vector<Vector> > eigenvectors;
  std::vector<PDB> pdbv;
  std::vector<std::string> pca_names;
public:
  explicit PCARMSD(const ActionOptions&);
  void calculate() override;
  static void registerKeywords(Keywords& keys);
};

//+PLUMEDOC DCOLVAR PCARMSD
/*
Calculate the PCA components for a number of provided eigenvectors and an average structure.

Information about this method can be found in the reference papers in the bibliography below.  An example input is provided below:

```plumed
#SETTINGS INPUTFILES=regtest/trajectories/pca/average.pdb,regtest/trajectories/pca/eigenvec.pdb
PCARMSD ...
  AVERAGE=regtest/trajectories/pca/average.pdb
  EIGENVECTORS=regtest/trajectories/pca/eigenvec.pdb
...
```

This input performs optimal alignment at every step and reports the rmsd so you know if you are far or close from the average structure.
It takes the average structure and eigenvectors in form of a pdb.
Note that beta and occupancy values in the pdb are neglected and all the weights are placed to 1 (differently from the RMSD colvar for example)

*/
//+ENDPLUMEDOC

PLUMED_REGISTER_ACTION(PCARMSD,"PCARMSD")

void PCARMSD::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords(keys);
  keys.add("compulsory","AVERAGE","a file in pdb format containing the reference structure and the atoms involved in the CV.");
  keys.add("compulsory","EIGENVECTORS","a file in pdb format containing the reference structure and the atoms involved in the CV.");
  keys.addOutputComponent("eig","default","scalar","the projections on each eigenvalue are stored on values labeled eig-1, eig-2, ...");
  keys.addOutputComponent("residual","default","scalar","the distance of the present configuration from the configuration supplied as AVERAGE in terms of mean squared displacement after optimal alignment ");
  keys.addFlag("SQUARED_ROOT",false," This should be set if you want RMSD instead of mean squared displacement ");
  keys.addDOI("10.1021/ct100413b");
  keys.addDOI("10.1021/jp068587c");
}

PCARMSD::PCARMSD(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  squared(true),
  nopbc(false) {
  std::string f_average;
  parse("AVERAGE",f_average);
  std::string type;
  type.assign("OPTIMAL");
  std::string f_eigenvectors;
  parse("EIGENVECTORS",f_eigenvectors);
  bool sq;
  parseFlag("SQUARED_ROOT",sq);
  if (sq) {
    squared=false;
  }
  parseFlag("NOPBC",nopbc);
  checkRead();

  PDB pdb;

  // read everything in ang and transform to nm if we are not in natural units
  if( !pdb.read(f_average,usingNaturalUnits(),0.1/getUnits().getLength()) ) {
    error("missing input file " + f_average );
  }

  rmsd=Tools::make_unique<RMSD>();
  bool remove_com=true;
  bool normalize_weights=true;
  // here align and displace are a simple vector of ones
  std::vector<double> align;
  align=pdb.getOccupancy();
  for(unsigned i=0; i<align.size(); i++) {
    align[i]=1.;
  } ;
  std::vector<double> displace;
  displace=pdb.getBeta();
  for(unsigned i=0; i<displace.size(); i++) {
    displace[i]=1.;
  } ;
  // reset again to reimpose unifrom weights (safe to disable this)
  rmsd->set(align,displace,pdb.getPositions(),type,remove_com,normalize_weights);
  requestAtoms( pdb.getAtomNumbers() );

  addComponentWithDerivatives("residual");
  componentIsNotPeriodic("residual");

  log.printf("  average from file %s\n",f_average.c_str());
  log.printf("  which contains %d atoms\n",getNumberOfAtoms());
  log.printf("  with indices : ");
  for(unsigned i=0; i<pdb.getAtomNumbers().size(); ++i) {
    if(i%25==0) {
      log<<"\n";
    }
    log.printf("%d ",pdb.getAtomNumbers()[i].serial());
  }
  log.printf("\n");
  log.printf("  method for alignment : %s \n",type.c_str() );
  if(nopbc) {
    log.printf("  without periodic boundary conditions\n");
  } else {
    log.printf("  using periodic boundary conditions\n");
  }

  log<<"  Bibliography "<<plumed.cite("Spiwok, Lipovova and Kralova, JPCB, 111, 3073 (2007)  ");
  log<<" "<<plumed.cite( "Sutto, D'Abramo, Gervasio, JCTC, 6, 3640 (2010)");

  unsigned neigenvects;
  neigenvects=0;
  // now get the eigenvectors
  // open the file
  if (FILE* fp=this->fopen(f_eigenvectors.c_str(),"r")) {
// call fclose when exiting this block
    auto deleter=[this](FILE* f) {
      this->fclose(f);
    };
    std::unique_ptr<FILE,decltype(deleter)> fp_deleter(fp,deleter);

    std::vector<AtomNumber> aaa;
    log<<"  Opening the eigenvectors file "<<f_eigenvectors.c_str()<<"\n";
    bool do_read=true;
    unsigned nat=0;
    while (do_read) {
      PDB mypdb;
      // check the units for reading this file: how can they make sense?
      do_read=mypdb.readFromFilepointer(fp,usingNaturalUnits(),0.1/getUnits().getLength());
      if(do_read) {
        neigenvects++;
        if(mypdb.getAtomNumbers().size()==0) {
          error("number of atoms in a frame should be more than zero");
        }
        if(nat==0) {
          nat=mypdb.getAtomNumbers().size();
        }
        if(nat!=mypdb.getAtomNumbers().size()) {
          error("frames should have the same number of atoms");
        }
        if(aaa.empty()) {
          aaa=mypdb.getAtomNumbers();
        }
        if(aaa!=mypdb.getAtomNumbers()) {
          error("frames should contain same atoms in same order");
        }
        log<<"  Found eigenvector: "<<neigenvects<<" containing  "<<mypdb.getAtomNumbers().size()<<" atoms\n";
        pdbv.push_back(mypdb);
        eigenvectors.push_back(mypdb.getPositions());
      } else {
        break ;
      }
    }
    log<<"  Found total "<<neigenvects<< " eigenvectors in the file "<<f_eigenvectors.c_str()<<" \n";
    if(neigenvects==0) {
      error("at least one eigenvector is expected");
    }
  }
  // the components
  for(unsigned i=0; i<neigenvects; i++) {
    std::string num;
    Tools::convert( i, num );
    std::string compName=std::string("eig-")+num;
    pca_names.push_back(compName);
    addComponentWithDerivatives(compName);
    componentIsNotPeriodic(compName);
  }
  turnOnDerivatives();

}

// calculator
void PCARMSD::calculate() {
  if(!nopbc) {
    makeWhole();
  }
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
    double val;
    val=0.;
    for(unsigned iat=0; iat<getNumberOfAtoms(); iat++) {
      val+=dotProduct(alignedpos[iat]-centeredref[iat],eigenvectors[i][iat]);
      der[iat].zero();
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

  for(unsigned i=0; i<getNumberOfComponents(); ++i) {
    setBoxDerivativesNoPbc( getPntrToComponent(i) );
  }

}

}
}



