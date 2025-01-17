/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2021, Andrea Arsiccio

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not
   claim that you wrote the original software. If you use this software
   in a product, an acknowledgment in the product documentation would be
   appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#include "Sasa.h"
#include "ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/GenericMolInfo.h"
#include "core/ActionSet.h"
#include <cstdio>
#include <iostream>
#include <string>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <fstream>


using namespace std;

namespace PLMD {
namespace sasa {

//+PLUMEDOC SASAMOD_COLVAR SASA_LCPO
/*
Calculates the solvent accessible surface area (SASA) of a protein molecule, or other properties related to it.

The atoms for which the SASA is desired should be indicated with the keyword ATOMS, and a pdb file of the protein must be provided in input with the MOLINFO keyword. The LCPO algorithm is used for the calculation (please, read and cite \cite Weiser1999). The radius of the solvent is assumed to be 0.14 nm, which is the radius of water molecules. Using the keyword NL_STRIDE it is also possible to specify the frequency with which the neighbor list for the calculation of the SASA is updated (the default is every 10 steps).

Different properties can be calculated and selected using the TYPE keyword:

1) the total SASA (TOTAL);

2) the free energy of transfer for the protein according to the transfer model (TRANSFER). This keyword can be used, for instance, to compute the transfer of a protein to different temperatures, as detailed in \cite Arsiccio-T-SASA-2021, or to different pressures, as detailed in \cite Arsiccio-P-SASA-2021, or to different osmolyte solutions, as detailed in \cite Arsiccio-C-SASA-2022.


When the TRANSFER keyword is used, a file with the free energy of transfer values for the sidechains and backbone atoms should be provided (using the keyword DELTAGFILE). Such file should have the following format:

\verbatim
----------------Sample DeltaG.dat file---------------------
ALA	0.711019999999962
ARG	-2.24832799999996
ASN	-2.74838799999999
ASP	-2.5626376
CYS	3.89864000000006
GLN	-1.76192
GLU	-2.38664400000002
GLY	0
HIS	-3.58152799999999
ILE	2.42634399999986
LEU	1.77233599999988
LYS	-1.92576400000002
MET	-0.262827999999956
PHE	1.62028800000007
PRO	-2.15598800000001
SER	-1.60934800000004
THR	-0.591559999999987
TRP	1.22936000000027
TYR	0.775547999999958
VAL	2.12779200000011
BACKBONE	1.00066920000002
-----------------------------------------------------------
\endverbatim

where the second column is the free energy of transfer for each sidechain/backbone, in kJ/mol.

A Python script for the computation of free energy of transfer values to describe the effect of osmolyte concentration, temperature and pressure (according to \cite Arsiccio-C-SASA-2022, \cite Arsiccio-T-SASA-2021 and \cite Arsiccio-P-SASA-2021) is freely available at https://github.com/andrea-arsiccio/DeltaG-calculation. The script automatically outputs a DeltaG.dat file compatible with this SASA module.

If the DELTAGFILE is not provided, the program computes the free energy of transfer values as if they had to take into account the effect of temperature according to approaches 2 or 3 in the paper \cite Arsiccio-T-SASA-2021. Please read and cite this paper if using the transfer model for computing the effect of temperature in implicit solvent simulations. For this purpose, the keyword APPROACH should be added, and set to either 2 or 3.

The SASA usually makes sense when atoms used for the calculation are all part of the same molecule. When running with periodic boundary conditions, the atoms should be in the proper periodic image. This is done automatically since PLUMED 2.2, by considering the ordered list of atoms and rebuilding the broken entities using a procedure that is equivalent to that done in \ref WHOLEMOLECULES. Notice that rebuilding is local to this action. This is different from \ref WHOLEMOLECULES which actually modifies the coordinates stored in PLUMED.

In case you want to recover the old behavior you should use the NOPBC flag.
In that case you need to take care that atoms are in the correct periodic image.

The SASA may also be computed using the SASA_HASEL collective variable, which makes use of the algorithm described in \cite Hasel1988. SASA_HASEL is less accurate then SASA_LCPO, but the computation is faster.



\par Examples

The following input tells plumed to print the total SASA for atoms 10 to 20 in a protein chain.
\plumedfile
SASA_LCPO TYPE=TOTAL ATOMS=10-20 NL_STRIDE=10 LABEL=sasa
PRINT ARG=sasa STRIDE=1 FILE=colvar
\endplumedfile


The following input tells plumed to compute the transfer free energy for the protein chain containing atoms 10 to 20. Such transfer free energy is then used as a bias in the simulation (e.g., implicit solvent simulations). The free energy of transfer values are read from a file called DeltaG.dat.

\plumedfile
SASA_LCPO TYPE=TRANSFER ATOMS=10-20 NL_STRIDE=10 DELTAGFILE=DeltaG.dat LABEL=sasa

bias: BIASVALUE ARG=sasa

PRINT ARG=sasa,bias.* STRIDE=1 FILE=colvar
\endplumedfile

The following input tells plumed to compute the transfer free energy for the protein chain containing atoms 10 to 20. Such transfer free energy is then used as a bias in the simulation (e.g., implicit solvent simulations). The free energy of transfer values are computed according to \cite Arsiccio-T-SASA-2021, and take into account the effect of temperature using approach 2 as described in the paper.

\plumedfile
SASA_LCPO TYPE=TRANSFER ATOMS=10-20 NL_STRIDE=10 APPROACH=2 LABEL=sasa

bias: BIASVALUE ARG=sasa

PRINT ARG=sasa,bias.* STRIDE=1 FILE=colvar
\endplumedfile

*/
//+ENDPLUMEDOC

class SASA_LCPO : public Colvar {
private:
  enum CV_TYPE {TOTAL, TRANSFER};
  int sasa_type;
  bool nopbc;
  double rs;
  string DeltaGValues;
  int approach;
  unsigned stride;
  unsigned nl_update;
  int firstStepFlag;
  double Ti;
  // cppcheck-suppress duplInheritedMember
  std::vector<AtomNumber> atoms;
  vector < vector < std::string > > AtomResidueName;
  vector < vector < double > > LCPOparam;
  unsigned natoms;
  vector < vector < double > > MaxSurf;
  vector < vector < double > > DeltaG;
  vector < vector < int > > Nlist;
public:
  static void registerKeywords(Keywords& keys);
  explicit SASA_LCPO(const ActionOptions&);
  void readPDB();
  map<string, vector<double> > setupLCPOparam();
  void readLCPOparam();
  void calcNlist();
  map<string, vector<double> > setupMaxSurfMap();
  void readMaxSurf();
  void readDeltaG();
  void computeDeltaG();
  virtual void calculate();
};

PLUMED_REGISTER_ACTION(SASA_LCPO,"SASA_LCPO")

void SASA_LCPO::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords(keys);
  keys.add("atoms","ATOMS","the group of atoms that you are calculating the SASA for");
  keys.add("compulsory","TYPE","TOTAL","The type of calculation you want to perform. Can be TOTAL or TRANSFER");
  keys.add("compulsory", "NL_STRIDE", "The frequency with which the neighbor list is updated.");
  keys.add("optional","DELTAGFILE","a file containing the free energy values for backbone and sidechains. Necessary only if TYPE = TRANSFER. A Python script for the computation of free energy of transfer values to describe the effect of osmolyte concentration, temperature and pressure is freely available at https://github.com/andrea-arsiccio/DeltaG-calculation. The script automatically outputs a DeltaG.dat file compatible with this SASA module. If TYPE = TRANSFER and no DELTAGFILE is provided, the free energy values are those describing the effect of temperature, and are computed using the temperature value passed by the MD engine");
  keys.add("optional","APPROACH","either approach 2 or 3. Necessary only if TYPE = TRANSFER and no DELTAGFILE is provided. If TYPE = TRANSFER and no DELTAGFILE is provided, the free energy values are those describing the effect of temperature, and the program must know if approach 2 or 3 (as described in Arsiccio and Shea, Protein Cold Denaturation in Implicit Solvent Simulations: A Transfer Free Energy Approach, J. Phys. Chem. B, 2021) needs to be used to compute them");
}


SASA_LCPO::SASA_LCPO(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  nopbc(false),
  DeltaGValues("absent"),
  Ti(0),
  stride(10),
  nl_update(0),
  firstStepFlag(0)
{
  rs = 0.14;
  parse("DELTAGFILE",DeltaGValues);
  parse("APPROACH", approach);
  parseAtomList("ATOMS",atoms);
  if(atoms.size()==0) error("no atoms specified");
  std::string Type;
  parse("TYPE",Type);
  parse("NL_STRIDE", stride);
  parseFlag("NOPBC",nopbc);
  checkRead();

  if(Type=="TOTAL") sasa_type=TOTAL;
  else if(Type=="TRANSFER") sasa_type=TRANSFER;
  else error("Unknown SASA type");

  switch(sasa_type)
  {
  case TOTAL:   log.printf("  TOTAL SASA;"); break;
  case TRANSFER: log.printf("  TRANSFER MODEL;"); break;
  }

  log.printf("  atoms involved : ");
  for(unsigned i=0; i<atoms.size(); ++i) {
    if(i%25==0) log<<"\n";
    log.printf("%d ",atoms[i].serial());
  }
  log.printf("\n");

  if(nopbc) {
    log<<"  PBC will be ignored\n";
  } else {
    log<<"  broken molecules will be rebuilt assuming atoms are in the proper order\n";
  }


  addValueWithDerivatives(); setNotPeriodic();
  requestAtoms(atoms);

  natoms = getNumberOfAtoms();
  AtomResidueName.resize(2);
  LCPOparam.resize(natoms);
  MaxSurf.resize(natoms);
  DeltaG.resize(natoms+1);
  Nlist.resize(natoms);


}


//splits strings into tokens. Used to read into LCPO parameters file and into reference pdb file
template <class Container>
void split(const std::string& str, Container& cont)
{
  std::istringstream iss(str);
  std::copy(std::istream_iterator<std::string>(iss),
            std::istream_iterator<std::string>(),
            std::back_inserter(cont));
}


//reads input PDB file provided with MOLINFO, and assigns atom and residue names to each atom involved in the CV

void SASA_LCPO::readPDB() {
  auto* moldat = plumed.getActionSet().selectLatest<GenericMolInfo*>(this);
  if( ! moldat ) error("Unable to find MOLINFO in input");
  AtomResidueName[0].clear();
  AtomResidueName[1].clear();

  for(unsigned i=0; i<natoms; i++) {
    string Aname = moldat->getAtomName(atoms[i]);
    string Rname = moldat->getResidueName(atoms[i]);
    AtomResidueName[0].push_back (Aname);
    AtomResidueName[1].push_back (Rname);
  }

}

//LCPO parameters database
map<string, vector<double> > SASA_LCPO::setupLCPOparam() {
  map<string, vector<double> > lcpomap;
  lcpomap["ALA_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["ALA_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["ALA_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["ALA_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["ALA_CB"] = { 1.7,  0.77887,  -0.28063,  -0.0012968,  0.00039328};
  lcpomap["ASP_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["ASP_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["ASP_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["ASP_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["ASP_CB"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["ASP_CG"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["ASP_OD1"] = { 1.6,  0.77914,  -0.25262,  -0.0016056,  0.00035071};
  lcpomap["ASP_OD2"] = { 1.6,  0.77914,  -0.25262,  -0.0016056,  0.00035071};
  lcpomap["ASN_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["ASN_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["ASN_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["ASN_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["ASN_CB"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["ASN_CG"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["ASN_OD1"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["ASN_ND2"] = { 1.65,  0.73511,  -0.22116,  -0.00089148,  0.0002523};
  lcpomap["ARG_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["ARG_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["ARG_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["ARG_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["ARG_CB"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["ARG_CG"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["ARG_CD"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["ARG_NE"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["ARG_NH1"] = { 1.65,  0.73511,  -0.22116,  -0.00089148,  0.0002523};
  lcpomap["ARG_NH2"] = { 1.65,  0.73511,  -0.22116,  -0.00089148,  0.0002523};
  lcpomap["ARG_CZ"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["CYS_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["CYS_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["CYS_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["CYS_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["CYS_CB"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["CYS_SG"] = { 1.9,  0.54581,  -0.19477,  -0.0012873,  0.00029247};
  lcpomap["GLU_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["GLU_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["GLU_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["GLU_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["GLU_CB"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["GLU_CG"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["GLU_CD"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["GLU_OE1"] = { 1.6,  0.77914,  -0.25262,  -0.0016056,  0.00035071};
  lcpomap["GLU_OE2"] = { 1.6,  0.77914,  -0.25262,  -0.0016056,  0.00035071};
  lcpomap["GLN_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["GLN_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["GLN_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["GLN_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["GLN_CB"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["GLN_CG"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["GLN_CD"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["GLN_OE1"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["GLN_NE2"] = { 1.65,  0.73511,  -0.22116,  -0.00089148,  0.0002523};
  lcpomap["GLY_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["GLY_CA"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["GLY_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["GLY_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["HIS_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["HIS_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["HIS_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["HIS_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["HIS_CB"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["HIS_CG"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["HIS_ND1"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["HIS_CE1"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["HIS_NE2"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["HIS_CD2"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["ILE_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["ILE_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["ILE_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["ILE_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["ILE_CB"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["ILE_CG2"] = { 1.7,  0.77887,  -0.28063,  -0.0012968,  0.00039328};
  lcpomap["ILE_CG1"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["ILE_CD1"] = { 1.7,  0.77887,  -0.28063,  -0.0012968,  0.00039328};
  lcpomap["LEU_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["LEU_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["LEU_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["LEU_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["LEU_CB"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["LEU_CG"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["LEU_CD1"] = { 1.7,  0.77887,  -0.28063,  -0.0012968,  0.00039328};
  lcpomap["LEU_CD2"] = { 1.7,  0.77887,  -0.28063,  -0.0012968,  0.00039328};
  lcpomap["LYS_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["LYS_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["LYS_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["LYS_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["LYS_CB"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["LYS_CG"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["LYS_CD"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["LYS_CE"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["LYS_NZ"] = { 1.65,  0.73511,  -0.22116,  -0.00089148,  0.0002523};
  lcpomap["MET_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["MET_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["MET_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["MET_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["MET_CB"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["MET_CG"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["MET_SD"] = { 1.9,  0.54581,  -0.19477,  -0.0012873,  0.00029247};
  lcpomap["MET_CE"] = { 1.7,  0.77887,  -0.28063,  -0.0012968,  0.00039328};
  lcpomap["PHE_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["PHE_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["PHE_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["PHE_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["PHE_CB"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["PHE_CG"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["PHE_CD1"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["PHE_CE1"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["PHE_CZ"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["PHE_CE2"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["PHE_CD2"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["PRO_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["PRO_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["PRO_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["PRO_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["PRO_CB"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["PRO_CG"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["PRO_CD"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["SER_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["SER_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["SER_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["SER_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["SER_CB"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["SER_OG"] = { 1.6,  0.77914,  -0.25262,  -0.0016056,  0.00035071};
  lcpomap["THR_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["THR_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["THR_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["THR_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["THR_CB"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["THR_CG2"] = { 1.7,  0.77887,  -0.28063,  -0.0012968,  0.00039328};
  lcpomap["THR_OG1"] = { 1.6,  0.77914,  -0.25262,  -0.0016056,  0.00035071};
  lcpomap["TRP_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["TRP_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["TRP_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["TRP_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["TRP_CB"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["TRP_CG"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["TRP_CD1"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["TRP_NE1"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["TRP_CE2"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["TRP_CZ2"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["TRP_CH2"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["TRP_CZ3"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["TRP_CE3"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["TRP_CD2"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["TYR_N"] = { 1.65,  0.062577,  -0.017874,  -8.312e-05,  1.9849e-05};
  lcpomap["TYR_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["TYR_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["TYR_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["TYR_CB"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["TYR_CG"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["TYR_CD1"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["TYR_CE1"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["TYR_CZ"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["TYR_OH"] = { 1.6,  0.77914,  -0.25262,  -0.0016056,  0.00035071};
  lcpomap["TYR_CE2"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["TYR_CD2"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["VAL_N"] = { 1.65,  0.062577,  -0.017874,  -8.312e-05,  1.9849e-05};
  lcpomap["VAL_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["VAL_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["VAL_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["VAL_CB"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["VAL_CG1"] = { 1.7,  0.77887,  -0.28063,  -0.0012968,  0.00039328};
  lcpomap["VAL_CG2"] = { 1.7,  0.77887,  -0.28063,  -0.0012968,  0.00039328};
  return lcpomap;
}

//assigns LCPO parameters to each atom reading from database
void SASA_LCPO::readLCPOparam() {

  for(unsigned i=0; i<natoms; i++) {
    LCPOparam[i].clear();
  }

  map<string, vector<double> > lcpomap;
  lcpomap = setupLCPOparam();
  vector<double> LCPOparamVector;
  string identifier;
  for(unsigned i=0; i<natoms; i++) {
    if ((AtomResidueName[0][i][0]=='O') || (AtomResidueName[0][i][0]=='N') || (AtomResidueName[0][i][0]=='C') || (AtomResidueName[0][i][0]=='S')) {
      identifier = AtomResidueName[1][i]+"_"+AtomResidueName[0][i];
      LCPOparamVector = lcpomap.at(identifier);
      LCPOparam[i].push_back (LCPOparamVector[0]+rs*10);
      LCPOparam[i].push_back (LCPOparamVector[1]);
      LCPOparam[i].push_back (LCPOparamVector[2]);
      LCPOparam[i].push_back (LCPOparamVector[3]);
      LCPOparam[i].push_back (LCPOparamVector[4]);
    }
  }


  for(unsigned i=0; i<natoms; i++) {
    if (LCPOparam[i].size()==0 ) {
      if ((AtomResidueName[0][i][0]=='O') || (AtomResidueName[0][i][0]=='N') || (AtomResidueName[0][i][0]=='C') || (AtomResidueName[0][i][0]=='S')) {
        cout << "Could not find LCPO paramaters for atom " << AtomResidueName[0][i] << " of residue " << AtomResidueName[1][i] << endl;
        error ("missing LCPO parameters\n");
      }
    }
  }

  if (AtomResidueName[0][0] == "N") {
    LCPOparam[0][1] = 7.3511e-01;
    LCPOparam[0][2] = -2.2116e-01;
    LCPOparam[0][3] = -8.9148e-04;
    LCPOparam[0][4] = 2.5230e-04;
  }

  if (AtomResidueName[0][natoms-1] == "O") {
    LCPOparam[natoms-1][1] = 8.8857e-01;
    LCPOparam[natoms-1][2] = -3.3421e-01;
    LCPOparam[natoms-1][3] = -1.8683e-03;
    LCPOparam[natoms-1][4] = 4.9372e-04;
  }
}


//Max Surf values, used only if TYPE=TRANSFER
map<string, vector<double> > SASA_LCPO::setupMaxSurfMap() {
  // Max Surface Area for backbone and sidechain, in nm2
  map<string, vector<double> > valuemap;
  valuemap ["ALA"]= {0.671446,0.420263,};
  valuemap ["ARG"]= {0.578582,1.95454,};
  valuemap ["ASN"]= {0.559411,1.07102,};
  valuemap ["ASP"]= {0.558363,1.03971,};
  valuemap ["CYS"]= {0.609271,0.657612,};
  valuemap ["GLN"]= {0.565651,1.433031,};
  valuemap ["GLU"]= {0.572399,1.41848,};
  valuemap ["GLY"]= {0.861439,0,};
  valuemap ["HIS"]= {0.559929,1.143827,};
  valuemap ["ILE"]= {0.553491,1.25334,};
  valuemap ["LEU"]= {0.570103,1.260459,};
  valuemap ["LYS"]= {0.580304,1.687487,};
  valuemap ["MET"]= {0.5856,1.498954,};
  valuemap ["PHE"]= {0.54155,1.394861,};
  valuemap ["PRO"]= {0.456048,0.849461,};
  valuemap ["SER"]= {0.59074,0.714538,};
  valuemap ["THR"]= {0.559116,0.951644,};
  valuemap ["TRP"]= {0.55536,1.59324,};
  valuemap ["TYR"]= {0.451171,1.566918,};
  valuemap ["VAL"]= {0.454809,0.928685,};
  return valuemap;
}



//reads maximum surface values per residue type and assigns values to each atom, only used if sasa_type = TRANSFER

void SASA_LCPO::readMaxSurf() {
  map<string, vector<double> > valuemap;
  valuemap = setupMaxSurfMap();
  vector<double> MaxSurfVector;

  for(unsigned i=0; i<natoms; i++) {
    MaxSurf[i].clear();
    MaxSurfVector = valuemap.at(AtomResidueName[1][i]);
    MaxSurf[i].push_back (MaxSurfVector[0]*100);
    MaxSurf[i].push_back (MaxSurfVector[1]*100);
  }
}

//reads file with free energy values for sidechains and for the backbone, and assigns values to each atom. Only used if sasa_type = TRANSFER

void SASA_LCPO::readDeltaG() {

  for(unsigned i=0; i<natoms; i++) {
    DeltaG[i].clear();
  }

  string DeltaGline;
  fstream DeltaGFile;
  DeltaGFile.open(DeltaGValues);
  if (DeltaGFile) {
    int backboneflag = 0;
    while(getline(DeltaGFile, DeltaGline)) {
      if(!DeltaGline.empty()) {
        std::vector<std::string> DeltaGtoken;
        split(DeltaGline, DeltaGtoken);
        for(unsigned i=0; i<natoms; i++) {
          if (DeltaGtoken[0].compare(AtomResidueName[1][i])==0 ) {
            DeltaG[i].push_back (std::atof(DeltaGtoken[1].c_str()));
          }
        }
        if (DeltaGtoken[0].compare("BACKBONE")==0 ) {
          backboneflag = 1;
          DeltaG[natoms].push_back (std::atof(DeltaGtoken[1].c_str()));
        }
      }
    }
    if ( backboneflag == 0) error("Cannot find backbone value in Delta G parameters file\n");
  }
  else error("Cannot open DeltaG file");

  for(unsigned i=0; i<natoms; i++) {
    if (DeltaG[i].size()==0 ) {
      cout << "Delta G value for residue " << AtomResidueName[1][i] << " not found " << endl;
      error ("missing Delta G parameter\n");
    }
  }

}

//computes free energy values for the sidechains and for the backbone, and assigns values to each atom. Only used if sasa_type = TRANSFER, and if no DELTAGFILE is provided. In this case, the free energy values are those describing the effect of temperature, and the program must know if approach 2 or 3 (as described in Arsiccio and Shea, Protein Cold Denaturation in Implicit Solvent Simulations: A Transfer Free Energy Approach, JPCB, 2021) needs to be used to compute them.

void SASA_LCPO::computeDeltaG() {

  for(unsigned i=0; i<natoms; i++) {
    DeltaG[i].clear();
  }

  double T;
  T = plumed.getAtoms().getKbT()/plumed.getAtoms().getKBoltzmann();

  if (T != Ti) {
    for(unsigned i=0; i<natoms; i++) {
      if (approach==2) {
        if (AtomResidueName[1][i].compare("ALA")==0) {
          DeltaG[i].push_back (-2.995/1000*std::pow(T,2)+1.808*T-272.895);
        }
        if (AtomResidueName[1][i].compare("ARG")==0) {
          DeltaG[i].push_back (-3.182/1000*std::pow(T,2)+1.894*T-282.032);
        }
        if (AtomResidueName[1][i].compare("ASN")==0) {
          DeltaG[i].push_back (-1.047/1000*std::pow(T,2)+0.6068*T-87.846);
        }
        if (AtomResidueName[1][i].compare("ASP")==0) {
          DeltaG[i].push_back (-0.1794/1000*std::pow(T,2)+0.1091*T-16.526);
        }
        if (AtomResidueName[1][i].compare("CYS")==0) {
          DeltaG[i].push_back (-3.09/1000*std::pow(T,2)+1.835*T-272.26);
        }
        if (AtomResidueName[1][i].compare("GLN")==0) {
          DeltaG[i].push_back (-2.23/1000*std::pow(T,2)+1.335*T-199.707);
        }
        if (AtomResidueName[1][i].compare("GLU")==0) {
          DeltaG[i].push_back (-1.511/1000*std::pow(T,2)+0.8904*T-131.168);
        }
        if (AtomResidueName[1][i].compare("GLY")==0) {
          DeltaG[i].push_back (0);
        }
        if (AtomResidueName[1][i].compare("HIS")==0) {
          DeltaG[i].push_back (-3.482/1000*std::pow(T,2)+2.084*T-311.694);
        }
        if (AtomResidueName[1][i].compare("ILE")==0) {
          DeltaG[i].push_back (-6.364/1000*std::pow(T,2)+3.8*T-567.444);
        }
        if (AtomResidueName[1][i].compare("LEU")==0) {
          DeltaG[i].push_back (-7.466/1000*std::pow(T,2)+4.417*T-653.394);
        }
        if (AtomResidueName[1][i].compare("LYS")==0) {
          DeltaG[i].push_back (-2.091/1000*std::pow(T,2)+1.2458*T-185.549);
        }
        if (AtomResidueName[1][i].compare("MET")==0) {
          DeltaG[i].push_back (-3.807/1000*std::pow(T,2)+2.272*T-339.007);
        }
        if (AtomResidueName[1][i].compare("PHE")==0) {
          DeltaG[i].push_back (-7.828/1000*std::pow(T,2)+4.644*T-688.874);
        }
        if (AtomResidueName[1][i].compare("PRO")==0) {
          DeltaG[i].push_back (-2.347/1000*std::pow(T,2)+1.411*T-212.059);
        }
        if (AtomResidueName[1][i].compare("SER")==0) {
          DeltaG[i].push_back (1.813/1000*std::pow(T,2)-1.05*T+151.957);
        }
        if (AtomResidueName[1][i].compare("THR")==0) {
          DeltaG[i].push_back (-2.64/1000*std::pow(T,2)+1.591*T-239.516);
        }
        if (AtomResidueName[1][i].compare("TRP")==0) {
          DeltaG[i].push_back (-11.66/1000*std::pow(T,2)+6.916*T-1025.293);
        }
        if (AtomResidueName[1][i].compare("TYR")==0) {
          DeltaG[i].push_back (-7.513/1000*std::pow(T,2)+4.478*T-667.261);
        }
        if (AtomResidueName[1][i].compare("VAL")==0) {
          DeltaG[i].push_back (-4.902/1000*std::pow(T,2)+2.921*T-435.309);
        }
        DeltaG[natoms].push_back (-0.6962/1000*std::pow(T,2)+0.4426*T-70.015);
      }
      if (approach==3) {
        if (AtomResidueName[1][i].compare("ALA")==0) {
          DeltaG[i].push_back (-2.995/1000*std::pow(T,2)+1.808*T-272.105);
        }
        if (AtomResidueName[1][i].compare("ARG")==0) {
          DeltaG[i].push_back (-3.182/1000*std::pow(T,2)+1.894*T-284.086);
        }
        if (AtomResidueName[1][i].compare("ASN")==0) {
          DeltaG[i].push_back (-1.047/1000*std::pow(T,2)+0.6068*T-90.597);
        }
        if (AtomResidueName[1][i].compare("ASP")==0) {
          DeltaG[i].push_back (-0.1794/1000*std::pow(T,2)+0.1091*T-19.143);
        }
        if (AtomResidueName[1][i].compare("CYS")==0) {
          DeltaG[i].push_back (-3.09/1000*std::pow(T,2)+1.835*T-268.527);
        }
        if (AtomResidueName[1][i].compare("GLN")==0) {
          DeltaG[i].push_back (-2.23/1000*std::pow(T,2)+1.335*T-201.559);
        }
        if (AtomResidueName[1][i].compare("GLU")==0) {
          DeltaG[i].push_back (-1.511/1000*std::pow(T,2)+0.8904*T-133.543);
        }
        if (AtomResidueName[1][i].compare("GLY")==0) {
          DeltaG[i].push_back (0);
        }
        if (AtomResidueName[1][i].compare("HIS")==0) {
          DeltaG[i].push_back (-3.482/1000*std::pow(T,2)+2.084*T-315.398);
        }
        if (AtomResidueName[1][i].compare("ILE")==0) {
          DeltaG[i].push_back (-6.364/1000*std::pow(T,2)+3.8*T-564.825);
        }
        if (AtomResidueName[1][i].compare("LEU")==0) {
          DeltaG[i].push_back (-7.466/1000*std::pow(T,2)+4.417*T-651.483);
        }
        if (AtomResidueName[1][i].compare("LYS")==0) {
          DeltaG[i].push_back (-2.091/1000*std::pow(T,2)+1.2458*T-187.485);
        }
        if (AtomResidueName[1][i].compare("MET")==0) {
          DeltaG[i].push_back (-3.807/1000*std::pow(T,2)+2.272*T-339.242);
        }
        if (AtomResidueName[1][i].compare("PHE")==0) {
          DeltaG[i].push_back (-7.828/1000*std::pow(T,2)+4.644*T-687.134);
        }
        if (AtomResidueName[1][i].compare("PRO")==0) {
          DeltaG[i].push_back (-2.347/1000*std::pow(T,2)+1.411*T-214.211);
        }
        if (AtomResidueName[1][i].compare("SER")==0) {
          DeltaG[i].push_back (1.813/1000*std::pow(T,2)-1.05*T+150.289);
        }
        if (AtomResidueName[1][i].compare("THR")==0) {
          DeltaG[i].push_back (-2.64/1000*std::pow(T,2)+1.591*T-240.267);
        }
        if (AtomResidueName[1][i].compare("TRP")==0) {
          DeltaG[i].push_back (-11.66/1000*std::pow(T,2)+6.916*T-1024.284);
        }
        if (AtomResidueName[1][i].compare("TYR")==0) {
          DeltaG[i].push_back (-7.513/1000*std::pow(T,2)+4.478*T-666.484);
        }
        if (AtomResidueName[1][i].compare("VAL")==0) {
          DeltaG[i].push_back (-4.902/1000*std::pow(T,2)+2.921*T-433.013);
        }
        DeltaG[natoms].push_back (-0.6927/1000*std::pow(T,2)+0.4404*T-68.724);
      }
    }
  }

  Ti = T;

  if (firstStepFlag ==0) {
    if (approach!=2 && approach!=3) {
      cout << "Unknown approach " << approach << endl;
    }
    for(unsigned i=0; i<natoms; i++) {
      if (DeltaG[i].size()==0 ) {
        cout << "Delta G value for residue " << AtomResidueName[1][i] << " not found " << endl;
        error ("missing Delta G parameter\n");
      }
    }
  }
}


//calculates neighbor list
void SASA_LCPO::calcNlist() {
  if(!nopbc) makeWhole();

  for(unsigned i = 0; i < natoms; i++) {
    Nlist[i].clear();
  }

  for(unsigned i = 0; i < natoms; i++) {
    if (LCPOparam[i].size()>0) {
      for (unsigned j = 0; j < i; j++) {
        if (LCPOparam[j].size()>0) {
          const Vector Delta_ij_vec = delta( getPosition(i), getPosition(j) );
          double Delta_ij_mod = Delta_ij_vec.modulo()*10;
          double overlapD = LCPOparam[i][0]+LCPOparam[j][0];
          if (Delta_ij_mod < overlapD) {
            Nlist.at(i).push_back (j);
            Nlist.at(j).push_back (i);
          }
        }
      }
    }
  }
}


//calculates SASA according to LCPO algorithm
void SASA_LCPO::calculate() {
  if(!nopbc) makeWhole();

  if(getExchangeStep()) nl_update = 0;
  if (firstStepFlag ==0) {
    readPDB();
    readLCPOparam();
  }
  if (nl_update == 0) {
    calcNlist();
  }



  double S1, Aij, Ajk, Aijk, Aijt, Ajkt, Aikt;
  double dAdd;
  double sasa = 0;
  vector<Vector> derivatives( natoms );
  Tensor virial;
  vector <double> dAijdc_2t(3);
  vector <double> dSASA_2_neigh_dc(3);
  vector <double> ddij_di(3);
  vector <double> ddik_di(3);

  if( sasa_type==TOTAL ) {
    for(unsigned i = 0; i < natoms; i++) {
      derivatives[i][0] = 0.;
      derivatives[i][1] = 0.;
      derivatives[i][2] = 0.;
      if ( LCPOparam[i].size()>1) {
        if (LCPOparam[i][1]>0.0) {
          Aij = 0.0;
          Aijk = 0.0;
          Ajk = 0.0;
          double ri = LCPOparam[i][0];
          S1 = 4*M_PI*ri*ri;
          vector <double> dAijdc_2(3, 0);
          vector <double> dAijdc_4(3, 0);


          for (unsigned j = 0; j < Nlist[i].size(); j++) {
            const Vector d_ij_vec = delta( getPosition(i), getPosition(Nlist[i][j]) );
            double d_ij = d_ij_vec.modulo()*10;

            double rj = LCPOparam[Nlist[i][j]][0];
            Aijt = (2*M_PI*ri*(ri-d_ij/2-((ri*ri-rj*rj)/(2*d_ij))));
            double sji = (2*M_PI*rj*(rj-d_ij/2+((ri*ri-rj*rj)/(2*d_ij))));

            dAdd = M_PI*rj*(-(ri*ri-rj*rj)/(d_ij*d_ij)-1);

            ddij_di[0] = -10*(getPosition(Nlist[i][j])[0]-getPosition(i)[0])/d_ij;
            ddij_di[1] = -10*(getPosition(Nlist[i][j])[1]-getPosition(i)[1])/d_ij;
            ddij_di[2] = -10*(getPosition(Nlist[i][j])[2]-getPosition(i)[2])/d_ij;

            Ajkt = 0.0;
            Aikt = 0.0;

            vector <double> dSASA_3_neigh_dc(3, 0.0);
            vector <double> dSASA_4_neigh_dc(3, 0.0);
            vector <double> dSASA_3_neigh_dc2(3, 0.0);
            vector <double> dSASA_4_neigh_dc2(3, 0.0);

            dSASA_2_neigh_dc[0] = dAdd * ddij_di[0];
            dSASA_2_neigh_dc[1] = dAdd * ddij_di[1];
            dSASA_2_neigh_dc[2] = dAdd * ddij_di[2];

            dAdd = M_PI*ri*((ri*ri-rj*rj)/(d_ij*d_ij)-1);


            dAijdc_2t[0] = dAdd * ddij_di[0];
            dAijdc_2t[1] = dAdd * ddij_di[1];
            dAijdc_2t[2] = dAdd * ddij_di[2];

            for (unsigned k = 0; k < Nlist[Nlist[i][j]].size(); k++) {
              if (std::find (Nlist[i].begin(), Nlist[i].end(), Nlist[Nlist[i][j]][k]) !=  Nlist[i].end()) {
                const Vector d_jk_vec = delta( getPosition(Nlist[i][j]), getPosition(Nlist[Nlist[i][j]][k]) );
                const Vector d_ik_vec = delta( getPosition(i), getPosition(Nlist[Nlist[i][j]][k]) );

                double d_jk = d_jk_vec.modulo()*10;
                double d_ik = d_ik_vec.modulo()*10;

                double rk = LCPOparam[Nlist[Nlist[i][j]][k]][0];
                double sjk =  (2*M_PI*rj*(rj-d_jk/2-((rj*rj-rk*rk)/(2*d_jk))));
                Ajkt += sjk;
                Aikt += (2*M_PI*ri*(ri-d_ik/2-((ri*ri-rk*rk)/(2*d_ik))));

                dAdd = M_PI*ri*((ri*ri-rk*rk)/(d_ik*d_ik)-1);

                ddik_di[0] = -10*(getPosition(Nlist[Nlist[i][j]][k])[0]-getPosition(i)[0])/d_ik;
                ddik_di[1] = -10*(getPosition(Nlist[Nlist[i][j]][k])[1]-getPosition(i)[1])/d_ik;
                ddik_di[2] = -10*(getPosition(Nlist[Nlist[i][j]][k])[2]-getPosition(i)[2])/d_ik;


                dSASA_3_neigh_dc[0] += dAdd*ddik_di[0];
                dSASA_3_neigh_dc[1] += dAdd*ddik_di[1];
                dSASA_3_neigh_dc[2] += dAdd*ddik_di[2];

                dAdd = M_PI*rk*(-(ri*ri-rk*rk)/(d_ik*d_ik)-1);

                dSASA_3_neigh_dc2[0] += dAdd*ddik_di[0];
                dSASA_3_neigh_dc2[1] += dAdd*ddik_di[1];
                dSASA_3_neigh_dc2[2] += dAdd*ddik_di[2];

                dSASA_4_neigh_dc2[0] += sjk*dAdd*ddik_di[0];
                dSASA_4_neigh_dc2[1] += sjk*dAdd*ddik_di[1];
                dSASA_4_neigh_dc2[2] += sjk*dAdd*ddik_di[2];

              }
            }
            dSASA_4_neigh_dc[0] = sji*dSASA_3_neigh_dc[0] + dSASA_4_neigh_dc2[0];
            dSASA_4_neigh_dc[1] = sji*dSASA_3_neigh_dc[1] + dSASA_4_neigh_dc2[1];
            dSASA_4_neigh_dc[2] = sji*dSASA_3_neigh_dc[2] + dSASA_4_neigh_dc2[2];

            dSASA_3_neigh_dc[0] += dSASA_3_neigh_dc2[0];
            dSASA_3_neigh_dc[1] += dSASA_3_neigh_dc2[1];
            dSASA_3_neigh_dc[2] += dSASA_3_neigh_dc2[2];

            dSASA_4_neigh_dc[0] += dSASA_2_neigh_dc[0] * Aikt;
            dSASA_4_neigh_dc[1] += dSASA_2_neigh_dc[1] * Aikt;
            dSASA_4_neigh_dc[2] += dSASA_2_neigh_dc[2] * Aikt;


            derivatives[i][0] += (dSASA_2_neigh_dc[0]*LCPOparam[Nlist[i][j]][2] + dSASA_3_neigh_dc[0]*LCPOparam[Nlist[i][j]][3]+dSASA_4_neigh_dc[0]*LCPOparam[Nlist[i][j]][4])/10;
            derivatives[i][1] += (dSASA_2_neigh_dc[1]*LCPOparam[Nlist[i][j]][2] + dSASA_3_neigh_dc[1]*LCPOparam[Nlist[i][j]][3]+dSASA_4_neigh_dc[1]*LCPOparam[Nlist[i][j]][4])/10;
            derivatives[i][2] += (dSASA_2_neigh_dc[2]*LCPOparam[Nlist[i][j]][2] + dSASA_3_neigh_dc[2]*LCPOparam[Nlist[i][j]][3]+dSASA_4_neigh_dc[2]*LCPOparam[Nlist[i][j]][4])/10;


            Aijk += (Aijt * Ajkt);
            Aij += Aijt;
            Ajk += Ajkt;

            dAijdc_2[0] += dAijdc_2t[0];
            dAijdc_2[1] += dAijdc_2t[1];
            dAijdc_2[2] += dAijdc_2t[2];


            dAijdc_4[0] += Ajkt*dAijdc_2t[0];
            dAijdc_4[1] += Ajkt*dAijdc_2t[1];
            dAijdc_4[2] += Ajkt*dAijdc_2t[2];


          }
          double sasai = (LCPOparam[i][1]*S1+LCPOparam[i][2]*Aij+LCPOparam[i][3]*Ajk+LCPOparam[i][4]*Aijk);
          if (sasai > 0 ) sasa += sasai/100;
          derivatives[i][0] += (dAijdc_2[0]*LCPOparam[i][2]+dAijdc_4[0]*LCPOparam[i][4])/10;
          derivatives[i][1] += (dAijdc_2[1]*LCPOparam[i][2]+dAijdc_4[1]*LCPOparam[i][4])/10;
          derivatives[i][2] += (dAijdc_2[2]*LCPOparam[i][2]+dAijdc_4[2]*LCPOparam[i][4])/10;
        }
      }
      virial -= Tensor(getPosition(i),derivatives[i]);
    }
  }



  if( sasa_type==TRANSFER ) {

    if (firstStepFlag ==0) {
      readMaxSurf();
    }

    if (firstStepFlag ==0 && DeltaGValues != "absent") {
      readDeltaG();
    }

    if (DeltaGValues == "absent") {
      computeDeltaG();
    }


    for(unsigned i = 0; i < natoms; i++) {



      derivatives[i][0] = 0.;
      derivatives[i][1] = 0.;
      derivatives[i][2] = 0.;

      if ( LCPOparam[i].size()>1) {
        if (LCPOparam[i][1]>0.0) {
          Aij = 0.0;
          Aijk = 0.0;
          Ajk = 0.0;
          double ri = LCPOparam[i][0];
          S1 = 4*M_PI*ri*ri;
          vector <double> dAijdc_2(3, 0);
          vector <double> dAijdc_4(3, 0);


          for (unsigned j = 0; j < Nlist[i].size(); j++) {
            const Vector d_ij_vec = delta( getPosition(i), getPosition(Nlist[i][j]) );
            double d_ij = d_ij_vec.modulo()*10;

            double rj = LCPOparam[Nlist[i][j]][0];
            Aijt = (2*M_PI*ri*(ri-d_ij/2-((ri*ri-rj*rj)/(2*d_ij))));
            double sji = (2*M_PI*rj*(rj-d_ij/2+((ri*ri-rj*rj)/(2*d_ij))));

            dAdd = M_PI*rj*(-(ri*ri-rj*rj)/(d_ij*d_ij)-1);
            ddij_di[0] = -10*(getPosition(Nlist[i][j])[0]-getPosition(i)[0])/d_ij;
            ddij_di[1] = -10*(getPosition(Nlist[i][j])[1]-getPosition(i)[1])/d_ij;
            ddij_di[2] = -10*(getPosition(Nlist[i][j])[2]-getPosition(i)[2])/d_ij;

            Ajkt = 0.0;
            Aikt = 0.0;

            vector <double> dSASA_3_neigh_dc(3, 0.0);
            vector <double> dSASA_4_neigh_dc(3, 0.0);
            vector <double> dSASA_3_neigh_dc2(3, 0.0);
            vector <double> dSASA_4_neigh_dc2(3, 0.0);

            dSASA_2_neigh_dc[0] = dAdd * ddij_di[0];
            dSASA_2_neigh_dc[1] = dAdd * ddij_di[1];
            dSASA_2_neigh_dc[2] = dAdd * ddij_di[2];

            dAdd = M_PI*ri*((ri*ri-rj*rj)/(d_ij*d_ij)-1);

            dAijdc_2t[0] = dAdd * ddij_di[0];
            dAijdc_2t[1] = dAdd * ddij_di[1];
            dAijdc_2t[2] = dAdd * ddij_di[2];

            for (unsigned k = 0; k < Nlist[Nlist[i][j]].size(); k++) {
              if (std::find (Nlist[i].begin(), Nlist[i].end(), Nlist[Nlist[i][j]][k]) !=  Nlist[i].end()) {
                const Vector d_jk_vec = delta( getPosition(Nlist[i][j]), getPosition(Nlist[Nlist[i][j]][k]) );
                const Vector d_ik_vec = delta( getPosition(i), getPosition(Nlist[Nlist[i][j]][k]) );

                double d_jk = d_jk_vec.modulo()*10;
                double d_ik = d_ik_vec.modulo()*10;

                double rk = LCPOparam[Nlist[Nlist[i][j]][k]][0];
                double sjk =  (2*M_PI*rj*(rj-d_jk/2-((rj*rj-rk*rk)/(2*d_jk))));
                Ajkt += sjk;
                Aikt += (2*M_PI*ri*(ri-d_ik/2-((ri*ri-rk*rk)/(2*d_ik))));

                dAdd = M_PI*ri*((ri*ri-rk*rk)/(d_ik*d_ik)-1);

                ddik_di[0] = -10*(getPosition(Nlist[Nlist[i][j]][k])[0]-getPosition(i)[0])/d_ik;
                ddik_di[1] = -10*(getPosition(Nlist[Nlist[i][j]][k])[1]-getPosition(i)[1])/d_ik;
                ddik_di[2] = -10*(getPosition(Nlist[Nlist[i][j]][k])[2]-getPosition(i)[2])/d_ik;


                dSASA_3_neigh_dc[0] += dAdd*ddik_di[0];
                dSASA_3_neigh_dc[1] += dAdd*ddik_di[1];
                dSASA_3_neigh_dc[2] += dAdd*ddik_di[2];

                dAdd = M_PI*rk*(-(ri*ri-rk*rk)/(d_ik*d_ik)-1);

                dSASA_3_neigh_dc2[0] += dAdd*ddik_di[0];
                dSASA_3_neigh_dc2[1] += dAdd*ddik_di[1];
                dSASA_3_neigh_dc2[2] += dAdd*ddik_di[2];

                dSASA_4_neigh_dc2[0] += sjk*dAdd*ddik_di[0];
                dSASA_4_neigh_dc2[1] += sjk*dAdd*ddik_di[1];
                dSASA_4_neigh_dc2[2] += sjk*dAdd*ddik_di[2];

              }
            }
            dSASA_4_neigh_dc[0] = sji*dSASA_3_neigh_dc[0] + dSASA_4_neigh_dc2[0];
            dSASA_4_neigh_dc[1] = sji*dSASA_3_neigh_dc[1] + dSASA_4_neigh_dc2[1];
            dSASA_4_neigh_dc[2] = sji*dSASA_3_neigh_dc[2] + dSASA_4_neigh_dc2[2];

            dSASA_3_neigh_dc[0] += dSASA_3_neigh_dc2[0];
            dSASA_3_neigh_dc[1] += dSASA_3_neigh_dc2[1];
            dSASA_3_neigh_dc[2] += dSASA_3_neigh_dc2[2];

            dSASA_4_neigh_dc[0] += dSASA_2_neigh_dc[0] * Aikt;
            dSASA_4_neigh_dc[1] += dSASA_2_neigh_dc[1] * Aikt;
            dSASA_4_neigh_dc[2] += dSASA_2_neigh_dc[2] * Aikt;

            if (AtomResidueName[0][Nlist[i][j]] == "N" || AtomResidueName[0][Nlist[i][j]] == "CA"  || AtomResidueName[0][Nlist[i][j]] == "C" || AtomResidueName[0][Nlist[i][j]] == "O") {
              derivatives[i][0] += ((dSASA_2_neigh_dc[0]*LCPOparam[Nlist[i][j]][2] + dSASA_3_neigh_dc[0]*LCPOparam[Nlist[i][j]][3]+dSASA_4_neigh_dc[0]*LCPOparam[Nlist[i][j]][4])/MaxSurf[Nlist[i][j]][0]*DeltaG[natoms][0])*10;
              derivatives[i][1] += ((dSASA_2_neigh_dc[1]*LCPOparam[Nlist[i][j]][2] + dSASA_3_neigh_dc[1]*LCPOparam[Nlist[i][j]][3]+dSASA_4_neigh_dc[1]*LCPOparam[Nlist[i][j]][4])/MaxSurf[Nlist[i][j]][0]*DeltaG[natoms][0])*10;
              derivatives[i][2] += ((dSASA_2_neigh_dc[2]*LCPOparam[Nlist[i][j]][2] + dSASA_3_neigh_dc[2]*LCPOparam[Nlist[i][j]][3]+dSASA_4_neigh_dc[2]*LCPOparam[Nlist[i][j]][4])/MaxSurf[Nlist[i][j]][0]*DeltaG[natoms][0])*10;
            }

            if (AtomResidueName[0][Nlist[i][j]] != "N" && AtomResidueName[0][Nlist[i][j]] != "CA"  && AtomResidueName[0][Nlist[i][j]] != "C" && AtomResidueName[0][Nlist[i][j]] != "O") {
              derivatives[i][0] += ((dSASA_2_neigh_dc[0]*LCPOparam[Nlist[i][j]][2] + dSASA_3_neigh_dc[0]*LCPOparam[Nlist[i][j]][3]+dSASA_4_neigh_dc[0]*LCPOparam[Nlist[i][j]][4])/MaxSurf[Nlist[i][j]][1]*DeltaG[Nlist[i][j]][0])*10;
              derivatives[i][1] += ((dSASA_2_neigh_dc[1]*LCPOparam[Nlist[i][j]][2] + dSASA_3_neigh_dc[1]*LCPOparam[Nlist[i][j]][3]+dSASA_4_neigh_dc[1]*LCPOparam[Nlist[i][j]][4])/MaxSurf[Nlist[i][j]][1]*DeltaG[Nlist[i][j]][0])*10;
              derivatives[i][2] += ((dSASA_2_neigh_dc[2]*LCPOparam[Nlist[i][j]][2] + dSASA_3_neigh_dc[2]*LCPOparam[Nlist[i][j]][3]+dSASA_4_neigh_dc[2]*LCPOparam[Nlist[i][j]][4])/MaxSurf[Nlist[i][j]][1]*DeltaG[Nlist[i][j]][0])*10;
            }

            Aijk += (Aijt * Ajkt);
            Aij += Aijt;
            Ajk += Ajkt;

            dAijdc_2[0] += dAijdc_2t[0];
            dAijdc_2[1] += dAijdc_2t[1];
            dAijdc_2[2] += dAijdc_2t[2];

            dAijdc_4[0] += Ajkt*dAijdc_2t[0];
            dAijdc_4[1] += Ajkt*dAijdc_2t[1];
            dAijdc_4[2] += Ajkt*dAijdc_2t[2];

          }
          double sasai = (LCPOparam[i][1]*S1+LCPOparam[i][2]*Aij+LCPOparam[i][3]*Ajk+LCPOparam[i][4]*Aijk);

          if (AtomResidueName[0][i] == "N" || AtomResidueName[0][i] == "CA"  || AtomResidueName[0][i] == "C" || AtomResidueName[0][i] == "O") {
            if (sasai > 0 ) sasa += (sasai/MaxSurf[i][0]*DeltaG[natoms][0]);
            derivatives[i][0] += ((dAijdc_2[0]*LCPOparam[i][2]+dAijdc_4[0]*LCPOparam[i][4])/MaxSurf[i][0]*DeltaG[natoms][0])*10;
            derivatives[i][1] += ((dAijdc_2[1]*LCPOparam[i][2]+dAijdc_4[1]*LCPOparam[i][4])/MaxSurf[i][0]*DeltaG[natoms][0])*10;
            derivatives[i][2] += ((dAijdc_2[2]*LCPOparam[i][2]+dAijdc_4[2]*LCPOparam[i][4])/MaxSurf[i][0]*DeltaG[natoms][0])*10;
          }

          if (AtomResidueName[0][i] != "N" && AtomResidueName[0][i] != "CA"  && AtomResidueName[0][i] != "C" && AtomResidueName[0][i] != "O") {
            if (sasai > 0. ) sasa += (sasai/MaxSurf[i][1]*DeltaG[i][0]);
            derivatives[i][0] += ((dAijdc_2[0]*LCPOparam[i][2]+dAijdc_4[0]*LCPOparam[i][4])/MaxSurf[i][1]*DeltaG[i][0])*10;
            derivatives[i][1] += ((dAijdc_2[1]*LCPOparam[i][2]+dAijdc_4[1]*LCPOparam[i][4])/MaxSurf[i][1]*DeltaG[i][0])*10;
            derivatives[i][2] += ((dAijdc_2[2]*LCPOparam[i][2]+dAijdc_4[2]*LCPOparam[i][4])/MaxSurf[i][1]*DeltaG[i][0])*10;
          }
        }
      }
      virial -= Tensor(getPosition(i),derivatives[i]);
    }
  }


  for(unsigned i=0; i<natoms; i++) { setAtomsDerivatives(i,derivatives[i]);}
  setBoxDerivatives(virial);
  setValue(sasa);
  firstStepFlag = 1;
  ++nl_update;
  if (nl_update == stride) {
    nl_update = 0;
  }
// setBoxDerivativesNoPbc();
}

}//namespace sasa
}//namespace PLMD
