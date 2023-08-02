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

//+PLUMEDOC SASAMOD_COLVAR SASA_HASEL
/*
Calculates the solvent accessible surface area (SASA) of a protein molecule, or other properties related to it.

The atoms for which the SASA is desired should be indicated with the keyword ATOMS, and a pdb file of the protein must be provided in input with the MOLINFO keyword. The algorithm described in \cite Hasel1988 is used for the calculation. The radius of the solvent is assumed to be 0.14 nm, which is the radius of water molecules. Using the keyword NL_STRIDE it is also possible to specify the frequency with which the neighbor list for the calculation of SASA is updated (the default is every 10 steps).

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

The SASA may also be computed using the SASA_LCPO collective variable, which makes use of the LCPO algorithm \cite Weiser1999. SASA_LCPO is more accurate then SASA_HASEL, but the computation is slower.


\par Examples

The following input tells plumed to print the total SASA for atoms 10 to 20 in a protein chain.
\plumedfile
SASA_HASEL TYPE=TOTAL ATOMS=10-20 NL_STRIDE=10 LABEL=sasa
PRINT ARG=sasa STRIDE=1 FILE=colvar
\endplumedfile


The following input tells plumed to compute the transfer free energy for the protein chain containing atoms 10 to 20. Such transfer free energy is then used as a bias in the simulation (e.g., implicit solvent simulations). The free energy of transfer values are read from a file called DeltaG.dat.

\plumedfile
SASA_HASEL TYPE=TRANSFER ATOMS=10-20 NL_STRIDE=10 DELTAGFILE=DeltaG.dat LABEL=sasa

bias: BIASVALUE ARG=sasa

PRINT ARG=sasa,bias.* STRIDE=1 FILE=colvar
\endplumedfile

The following input tells plumed to compute the transfer free energy for the protein chain containing atoms 10 to 20. Such transfer free energy is then used as a bias in the simulation (e.g., implicit solvent simulations). The free energy of transfer values are computed according to \cite Arsiccio-T-SASA-2021, and take into account the effect of temperature using approach 2 as described in the paper.

\plumedfile
SASA_HASEL TYPE=TRANSFER ATOMS=10-20 NL_STRIDE=10 APPROACH=2 LABEL=sasa

bias: BIASVALUE ARG=sasa

PRINT ARG=sasa,bias.* STRIDE=1 FILE=colvar
\endplumedfile

*/
//+ENDPLUMEDOC

class SASA_HASEL : public Colvar {
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
  vector < vector < double > > SASAparam;
  vector < vector < std::string > > CONNECTparam;
  unsigned natoms;
  vector < vector < double > > MaxSurf;
  vector < vector < double > > DeltaG;
  vector < vector < int > > Nlist;
public:
  static void registerKeywords(Keywords& keys);
  explicit SASA_HASEL(const ActionOptions&);
  void readPDB();
  map<string, vector<std::string> > setupHASELparam();
  void readSASAparam();
  void calcNlist();
  map<string, vector<double> > setupMaxSurfMap();
  void readMaxSurf();
  void readDeltaG();
  void computeDeltaG();
  virtual void calculate();
};

PLUMED_REGISTER_ACTION(SASA_HASEL,"SASA_HASEL")

void SASA_HASEL::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords(keys);
  keys.add("atoms","ATOMS","the group of atoms that you are calculating the SASA for");
  keys.add("compulsory","TYPE","TOTAL","The type of calculation you want to perform. Can be TOTAL or TRANSFER");
  keys.add("compulsory", "NL_STRIDE", "The frequency with which the neighbor list for the calculation of SASA is updated.");
  keys.add("optional","DELTAGFILE","a file containing the free energy of transfer values for backbone and sidechains atoms. Necessary only if TYPE = TRANSFER. A Python script for the computation of free energy of transfer values to describe the effect of osmolyte concentration, temperature and pressure is freely available at https://github.com/andrea-arsiccio/DeltaG-calculation. The script automatically outputs a DeltaG.dat file compatible with this SASA module. If TYPE = TRANSFER and no DELTAGFILE is provided, the free energy values are those describing the effect of temperature, and are computed using the temperature value passed by the MD engine");
  keys.add("optional","APPROACH","either approach 2 or 3. Necessary only if TYPE = TRANSFER and no DELTAGFILE is provided. If TYPE = TRANSFER and no DELTAGFILE is provided, the free energy values are those describing the effect of temperature, and the program must know if approach 2 or 3 (as described in Arsiccio and Shea, Protein Cold Denaturation in Implicit Solvent Simulations: A Transfer Free Energy Approach, J. Phys. Chem. B, 2021) needs to be used to compute them");
}


SASA_HASEL::SASA_HASEL(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  nopbc(false),
  stride(10),
  nl_update(0),
  DeltaGValues("absent"),
  Ti(0),
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
  SASAparam.resize(natoms);
  CONNECTparam.resize(natoms);
  MaxSurf.resize(natoms);
  DeltaG.resize(natoms+1);
  Nlist.resize(natoms);


}


//splits strings into tokens. Used to read into SASA parameters file and into reference pdb file
template <class Container>
void split(const std::string& str, Container& cont)
{
  std::istringstream iss(str);
  std::copy(std::istream_iterator<std::string>(iss),
            std::istream_iterator<std::string>(),
            std::back_inserter(cont));
}


//reads input PDB file provided with MOLINFO, and assigns atom and residue names to each atom involved in the CV

void SASA_HASEL::readPDB() {
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


//Hasel et al. parameters database
map<string, vector<std::string> > SASA_HASEL::setupHASELparam() {
  map<string, vector<std::string> > haselmap;
  haselmap["ALA_N"] = { "1.55",  "1.028",  "H",  "CA",  "Z",  "Z", };
  haselmap["ALA_CA"] = { "1.7",  "2.149",  "N",  "C",  "CB",  "Z", };
  haselmap["ALA_C"] = { "1.72",  "1.554",  "CA",  "O",  "Z",  "Z", };
  haselmap["ALA_O"] = { "1.5",  "0.926",  "C",  "Z",  "Z",  "Z", };
  haselmap["ALA_H"] = { "1.1",  "1.128",  "N",  "Z",  "Z",  "Z", };
  haselmap["ALA_CB"] = { "2.0",  "0.88",  "CA",  "Z",  "Z",  "Z", };
  haselmap["ASP_N"] = { "1.55",  "1.028",  "H",  "CA",  "Z",  "Z", };
  haselmap["ASP_CA"] = { "1.7",  "2.149",  "N",  "C",  "CB",  "Z", };
  haselmap["ASP_C"] = { "1.72",  "1.554",  "CA",  "O",  "Z",  "Z", };
  haselmap["ASP_O"] = { "1.5",  "0.926",  "C",  "Z",  "Z",  "Z", };
  haselmap["ASP_H"] = { "1.1",  "1.128",  "N",  "Z",  "Z",  "Z", };
  haselmap["ASP_CB"] = { "1.9",  "1.045",  "CA",  "CG",  "Z",  "Z", };
  haselmap["ASP_CG"] = { "1.72",  "1.554",  "CB",  "OD1",  "OD2",  "Z", };
  haselmap["ASP_OD1"] = { "1.5",  "0.926",  "CG",  "Z",  "Z",  "Z", };
  haselmap["ASP_OD2"] = { "1.7",  "0.922",  "CG",  "Z",  "Z",  "Z", };
  haselmap["ASN_N"] = { "1.55",  "1.028",  "H",  "CA",  "Z",  "Z", };
  haselmap["ASN_CA"] = { "1.7",  "2.149",  "N",  "C",  "CB",  "Z", };
  haselmap["ASN_C"] = { "1.72",  "1.554",  "CA",  "O",  "Z",  "Z", };
  haselmap["ASN_O"] = { "1.5",  "0.926",  "C",  "Z",  "Z",  "Z", };
  haselmap["ASN_H"] = { "1.1",  "1.128",  "N",  "Z",  "Z",  "Z", };
  haselmap["ASN_CB"] = { "1.9",  "1.045",  "CA",  "CG",  "Z",  "Z", };
  haselmap["ASN_CG"] = { "1.7",  "2.149",  "CB",  "OD1",  "ND2",  "Z", };
  haselmap["ASN_OD1"] = { "1.5",  "0.926",  "CG",  "Z",  "Z",  "Z", };
  haselmap["ASN_ND2"] = { "1.6",  "1.215",  "CG",  "1HD2",  "1HD2",  "Z", };
  haselmap["ASN_1HD2"] = { "1.1",  "1.128",  "ND2",  "Z",  "Z",  "Z", };
  haselmap["ASN_2HD2"] = { "1.1",  "1.128",  "ND2",  "Z",  "Z",  "Z", };
  haselmap["ARG_N"] = { "1.55",  "1.028",  "H",  "CA",  "Z",  "Z", };
  haselmap["ARG_CA"] = { "1.7",  "2.149",  "N",  "C",  "CB",  "Z", };
  haselmap["ARG_C"] = { "1.72",  "1.554",  "CA",  "O",  "Z",  "Z", };
  haselmap["ARG_O"] = { "1.5",  "0.926",  "C",  "Z",  "Z",  "Z", };
  haselmap["ARG_H"] = { "1.1",  "1.128",  "N",  "Z",  "Z",  "Z", };
  haselmap["ARG_CB"] = { "1.9",  "1.045",  "CA",  "CG",  "Z",  "Z", };
  haselmap["ARG_CG"] = { "1.9",  "1.045",  "CB",  "CD",  "Z",  "Z", };
  haselmap["ARG_CD"] = { "1.9",  "1.045",  "CG",  "NE",  "Z",  "Z", };
  haselmap["ARG_NE"] = { "1.55",  "1.028",  "CD",  "HE",  "CZ",  "Z", };
  haselmap["ARG_NH1"] = { "1.55",  "1.028",  "CZ",  "1HH1",  "2HH1",  "Z", };
  haselmap["ARG_NH2"] = { "1.55",  "1.028",  "CZ",  "1HH2",  "2HH2",  "Z", };
  haselmap["ARG_CZ"] = { "1.72",  "1.554",  "NE",  "NH1",  "NH2",  "Z", };
  haselmap["ARG_HE"] = { "1.1",  "1.128",  "NE",  "Z",  "Z",  "Z", };
  haselmap["ARG_1HH2"] = { "1.1",  "1.128",  "NH2",  "Z",  "Z",  "Z", };
  haselmap["ARG_2HH2"] = { "1.1",  "1.128",  "NH2",  "Z",  "Z",  "Z", };
  haselmap["ARG_2HH1"] = { "1.1",  "1.128",  "NH1",  "Z",  "Z",  "Z", };
  haselmap["ARG_1HH1"] = { "1.1",  "1.128",  "NH1",  "Z",  "Z",  "Z", };
  haselmap["CYS_N"] = { "1.55",  "1.028",  "H",  "CA",  "Z",  "Z", };
  haselmap["CYS_CA"] = { "1.7",  "2.149",  "N",  "C",  "CB",  "Z", };
  haselmap["CYS_C"] = { "1.72",  "1.554",  "CA",  "O",  "Z",  "Z", };
  haselmap["CYS_O"] = { "1.5",  "0.926",  "C",  "Z",  "Z",  "Z", };
  haselmap["CYS_H"] = { "1.1",  "1.128",  "N",  "Z",  "Z",  "Z", };
  haselmap["CYS_CB"] = { "1.9",  "1.045",  "CA",  "SG",  "Z",  "Z", };
  haselmap["CYS_SG"] = { "1.8",  "1.121",  "CB",  "HG",  "Z",  "Z", };
  haselmap["CYS_HG"] = { "1.2",  "0.928",  "SG",  "Z",  "Z",  "Z", };
  haselmap["GLU_N"] = { "1.55",  "1.028",  "H",  "CA",  "Z",  "Z", };
  haselmap["GLU_CA"] = { "1.7",  "2.149",  "N",  "C",  "CB",  "Z", };
  haselmap["GLU_C"] = { "1.72",  "1.554",  "CA",  "O",  "Z",  "Z", };
  haselmap["GLU_O"] = { "1.5",  "0.926",  "C",  "Z",  "Z",  "Z", };
  haselmap["GLU_H"] = { "1.1",  "1.128",  "N",  "Z",  "Z",  "Z", };
  haselmap["GLU_CB"] = { "1.9",  "1.045",  "CA",  "CG",  "Z",  "Z", };
  haselmap["GLU_CG"] = { "1.9",  "1.045",  "CB",  "CD",  "Z",  "Z", };
  haselmap["GLU_CD"] = { "1.72",  "1.554",  "CG",  "OE1",  "OE2",  "Z", };
  haselmap["GLU_OE1"] = { "1.5",  "0.926",  "CD",  "Z",  "Z",  "Z", };
  haselmap["GLU_OE2"] = { "1.7",  "0.922",  "CD",  "Z",  "Z",  "Z", };
  haselmap["GLN_N"] = { "1.55",  "1.028",  "H",  "CA",  "Z",  "Z", };
  haselmap["GLN_CA"] = { "1.7",  "2.149",  "N",  "C",  "CB",  "Z", };
  haselmap["GLN_C"] = { "1.72",  "1.554",  "CA",  "O",  "Z",  "Z", };
  haselmap["GLN_O"] = { "1.5",  "0.926",  "C",  "Z",  "Z",  "Z", };
  haselmap["GLN_H"] = { "1.1",  "1.128",  "N",  "Z",  "Z",  "Z", };
  haselmap["GLN_CB"] = { "1.9",  "1.045",  "CA",  "CG",  "Z",  "Z", };
  haselmap["GLN_CG"] = { "1.9",  "1.045",  "CB",  "CD",  "Z",  "Z", };
  haselmap["GLN_CD"] = { "1.72",  "1.554",  "CG",  "OE1",  "NE2",  "Z", };
  haselmap["GLN_OE1"] = { "1.5",  "0.926",  "CD",  "Z",  "Z",  "Z", };
  haselmap["GLN_NE2"] = { "1.6",  "1.215",  "CD",  "2HE2",  "1HE2",  "Z", };
  haselmap["GLN_2HE2"] = { "1.1",  "1.128",  "NE2",  "Z",  "Z",  "Z", };
  haselmap["GLN_1HE2"] = { "1.1",  "1.128",  "NE2",  "Z",  "Z",  "Z", };
  haselmap["GLY_N"] = { "1.55",  "1.028",  "H",  "CA",  "Z",  "Z", };
  haselmap["GLY_CA"] = { "1.7",  "2.149",  "N",  "C",  "Z",  "Z", };
  haselmap["GLY_C"] = { "1.72",  "1.554",  "CA",  "O",  "Z",  "Z", };
  haselmap["GLY_O"] = { "1.5",  "0.926",  "C",  "Z",  "Z",  "Z", };
  haselmap["GLY_H"] = { "1.1",  "1.128",  "N",  "Z",  "Z",  "Z", };
  haselmap["HIS_N"] = { "1.55",  "1.028",  "H",  "CA",  "Z",  "Z", };
  haselmap["HIS_CA"] = { "1.7",  "2.149",  "N",  "C",  "CB",  "Z", };
  haselmap["HIS_C"] = { "1.72",  "1.554",  "CA",  "O",  "Z",  "Z", };
  haselmap["HIS_O"] = { "1.5",  "0.926",  "C",  "Z",  "Z",  "Z", };
  haselmap["HIS_H"] = { "1.1",  "1.128",  "N",  "Z",  "Z",  "Z", };
  haselmap["HIS_CB"] = { "1.9",  "1.045",  "CA",  "CG",  "Z",  "Z", };
  haselmap["HIS_CG"] = { "1.72",  "1.554",  "CB",  "CD2",  "ND1",  "Z", };
  haselmap["HIS_ND1"] = { "1.55",  "1.028",  "CG",  "CE1",  "Z",  "Z", };
  haselmap["HIS_CE1"] = { "1.8",  "1.073",  "ND1",  "NE2",  "Z",  "Z", };
  haselmap["HIS_NE2"] = { "1.55",  "1.413",  "CE1",  "2HE",  "CD2",  "Z", };
  haselmap["HIS_CD2"] = { "1.8",  "1.073",  "NE2",  "CG",  "Z",  "Z", };
  haselmap["HIS_2HE"] = { "1.1",  "1.128",  "NE2",  "Z",  "Z",  "Z", };
  haselmap["ILE_N"] = { "1.55",  "1.028",  "H",  "CA",  "Z",  "Z", };
  haselmap["ILE_CA"] = { "1.7",  "2.149",  "N",  "C",  "CB",  "Z", };
  haselmap["ILE_C"] = { "1.72",  "1.554",  "CA",  "O",  "Z",  "Z", };
  haselmap["ILE_O"] = { "1.5",  "0.926",  "C",  "Z",  "Z",  "Z", };
  haselmap["ILE_H"] = { "1.1",  "1.128",  "N",  "Z",  "Z",  "Z", };
  haselmap["ILE_CB"] = { "1.8",  "1.276",  "CA",  "CG2",  "CG1",  "Z", };
  haselmap["ILE_CG2"] = { "2.0",  "0.88",  "CB",  "Z",  "Z",  "Z", };
  haselmap["ILE_CG1"] = { "1.9",  "1.045",  "CB",  "CD1",  "Z",  "Z", };
  haselmap["ILE_CD1"] = { "2.0",  "0.88",  "CG1",  "Z",  "Z",  "Z", };
  haselmap["LEU_N"] = { "1.55",  "1.028",  "H",  "CA",  "Z",  "Z", };
  haselmap["LEU_CA"] = { "1.7",  "2.149",  "N",  "C",  "CB",  "Z", };
  haselmap["LEU_C"] = { "1.72",  "1.554",  "CA",  "O",  "Z",  "Z", };
  haselmap["LEU_O"] = { "1.5",  "0.926",  "C",  "Z",  "Z",  "Z", };
  haselmap["LEU_H"] = { "1.1",  "1.128",  "N",  "Z",  "Z",  "Z", };
  haselmap["LEU_CB"] = { "1.9",  "1.045",  "CA",  "CG",  "Z",  "Z", };
  haselmap["LEU_CG"] = { "1.8",  "1.276",  "CB",  "CD1",  "CD2",  "Z", };
  haselmap["LEU_CD1"] = { "2.0",  "0.88",  "CG",  "Z",  "Z",  "Z", };
  haselmap["LEU_CD2"] = { "2.0",  "0.88",  "CG",  "Z",  "Z",  "Z", };
  haselmap["LYS_N"] = { "1.55",  "1.028",  "H",  "CA",  "Z",  "Z", };
  haselmap["LYS_CA"] = { "1.7",  "2.149",  "N",  "C",  "CB",  "Z", };
  haselmap["LYS_C"] = { "1.72",  "1.554",  "CA",  "O",  "Z",  "Z", };
  haselmap["LYS_O"] = { "1.5",  "0.926",  "C",  "Z",  "Z",  "Z", };
  haselmap["LYS_H"] = { "1.1",  "1.128",  "N",  "Z",  "Z",  "Z", };
  haselmap["LYS_CB"] = { "1.9",  "1.045",  "CA",  "CG",  "Z",  "Z", };
  haselmap["LYS_CG"] = { "1.9",  "1.045",  "CB",  "CD",  "Z",  "Z", };
  haselmap["LYS_CD"] = { "1.9",  "1.045",  "CG",  "CE",  "Z",  "Z", };
  haselmap["LYS_CE"] = { "1.9",  "1.045",  "CD",  "NZ",  "Z",  "Z", };
  haselmap["LYS_NZ"] = { "1.6",  "1.215",  "CE",  "1HZ",  "2HZ",  "3HZ", };
  haselmap["LYS_1HZ"] = { "1.1",  "1.128",  "NZ",  "Z",  "Z",  "Z", };
  haselmap["LYS_2HZ"] = { "1.1",  "1.128",  "NZ",  "Z",  "Z",  "Z", };
  haselmap["LYS_3HZ"] = { "1.1",  "1.128",  "NZ",  "Z",  "Z",  "Z", };
  haselmap["MET_N"] = { "1.55",  "1.028",  "H",  "CA",  "Z",  "Z", };
  haselmap["MET_CA"] = { "1.7",  "2.149",  "N",  "C",  "CB",  "Z", };
  haselmap["MET_C"] = { "1.72",  "1.554",  "CA",  "O",  "Z",  "Z", };
  haselmap["MET_O"] = { "1.5",  "0.926",  "C",  "Z",  "Z",  "Z", };
  haselmap["MET_H"] = { "1.1",  "1.128",  "N",  "Z",  "Z",  "Z", };
  haselmap["MET_CB"] = { "1.9",  "1.045",  "CA",  "CG",  "Z",  "Z", };
  haselmap["MET_CG"] = { "1.9",  "1.045",  "CB",  "SD",  "Z",  "Z", };
  haselmap["MET_SD"] = { "1.8",  "1.121",  "CG",  "CE",  "Z",  "Z", };
  haselmap["MET_CE"] = { "2.0",  "0.88",  "SD",  "Z",  "Z",  "Z", };
  haselmap["PHE_N"] = { "1.55",  "1.028",  "H",  "CA",  "Z",  "Z", };
  haselmap["PHE_CA"] = { "1.7",  "2.149",  "N",  "C",  "CB",  "Z", };
  haselmap["PHE_C"] = { "1.72",  "1.554",  "CA",  "O",  "Z",  "Z", };
  haselmap["PHE_O"] = { "1.5",  "0.926",  "C",  "Z",  "Z",  "Z", };
  haselmap["PHE_H"] = { "1.1",  "1.128",  "N",  "Z",  "Z",  "Z", };
  haselmap["PHE_CB"] = { "1.9",  "1.045",  "CA",  "CG",  "Z",  "Z", };
  haselmap["PHE_CG"] = { "1.72",  "1.554",  "CB",  "CD1",  "CD2",  "Z", };
  haselmap["PHE_CD1"] = { "1.8",  "1.073",  "CG",  "CE1",  "Z",  "Z", };
  haselmap["PHE_CE1"] = { "1.8",  "1.073",  "CD1",  "CZ",  "Z",  "Z", };
  haselmap["PHE_CZ"] = { "1.8",  "1.073",  "CE1",  "CE2",  "Z",  "Z", };
  haselmap["PHE_CE2"] = { "1.8",  "1.073",  "CZ",  "CD2",  "Z",  "Z", };
  haselmap["PHE_CD2"] = { "1.8",  "1.073",  "CE2",  "CG",  "Z",  "Z", };
  haselmap["PRO_N"] = { "1.55",  "1.028",  "CD",  "CA",  "Z",  "Z", };
  haselmap["PRO_CA"] = { "1.7",  "2.149",  "N",  "C",  "CB",  "Z", };
  haselmap["PRO_C"] = { "1.72",  "1.554",  "CA",  "O",  "Z",  "Z", };
  haselmap["PRO_O"] = { "1.5",  "0.926",  "C",  "Z",  "Z",  "Z", };
  haselmap["PRO_CB"] = { "1.9",  "1.045",  "CA",  "CG",  "Z",  "Z", };
  haselmap["PRO_CG"] = { "1.9",  "1.045",  "CB",  "CD",  "Z",  "Z", };
  haselmap["PRO_CD"] = { "1.9",  "1.045",  "CG",  "N",  "Z",  "Z", };
  haselmap["SER_N"] = { "1.55",  "1.028",  "H",  "CA",  "Z",  "Z", };
  haselmap["SER_CA"] = { "1.7",  "2.149",  "N",  "C",  "CB",  "Z", };
  haselmap["SER_C"] = { "1.72",  "1.554",  "CA",  "O",  "Z",  "Z", };
  haselmap["SER_O"] = { "1.5",  "0.926",  "C",  "Z",  "Z",  "Z", };
  haselmap["SER_H"] = { "1.1",  "1.128",  "N",  "Z",  "Z",  "Z", };
  haselmap["SER_CB"] = { "1.9",  "1.045",  "CA",  "OG",  "Z",  "Z", };
  haselmap["SER_OG"] = { "1.52",  "1.08",  "CB",  "HG",  "Z",  "Z", };
  haselmap["SER_HG"] = { "1.0",  "0.944",  "OG",  "Z",  "Z",  "Z", };
  haselmap["THR_N"] = { "1.55",  "1.028",  "H",  "CA",  "Z",  "Z", };
  haselmap["THR_CA"] = { "1.7",  "2.149",  "N",  "C",  "CB",  "Z", };
  haselmap["THR_C"] = { "1.72",  "1.554",  "CA",  "O",  "Z",  "Z", };
  haselmap["THR_O"] = { "1.5",  "0.926",  "C",  "Z",  "Z",  "Z", };
  haselmap["THR_H"] = { "1.1",  "1.128",  "N",  "Z",  "Z",  "Z", };
  haselmap["THR_CB"] = { "1.8",  "1.276",  "CA",  "CG2",  "OG1",  "Z", };
  haselmap["THR_CG2"] = { "2.0",  "0.88",  "CB",  "Z",  "Z",  "Z", };
  haselmap["THR_OG1"] = { "1.52",  "1.08",  "1HG",  "CB",  "Z",  "Z", };
  haselmap["THR_1HG"] = { "1.0",  "0.944",  "OG1",  "Z",  "Z",  "Z", };
  haselmap["TRP_N"] = { "1.55",  "1.028",  "H",  "CA",  "Z",  "Z", };
  haselmap["TRP_CA"] = { "1.7",  "2.149",  "N",  "C",  "CB",  "Z", };
  haselmap["TRP_C"] = { "1.72",  "1.554",  "CA",  "O",  "Z",  "Z", };
  haselmap["TRP_O"] = { "1.5",  "0.926",  "C",  "Z",  "Z",  "Z", };
  haselmap["TRP_H"] = { "1.1",  "1.128",  "N",  "Z",  "Z",  "Z", };
  haselmap["TRP_CB"] = { "1.9",  "1.045",  "CA",  "CG",  "Z",  "Z", };
  haselmap["TRP_CG"] = { "1.72",  "1.554",  "CB",  "CD2",  "CD1",  "Z", };
  haselmap["TRP_CD1"] = { "1.8",  "1.073",  "CG",  "NE1",  "Z",  "Z", };
  haselmap["TRP_NE1"] = { "1.55",  "1.413",  "CD1",  "CE2",  "1HE",  "Z", };
  haselmap["TRP_CE2"] = { "1.72",  "1.554",  "NE1",  "CD2",  "CZ2",  "Z", };
  haselmap["TRP_CZ2"] = { "1.8",  "1.073",  "CE2",  "CH2",  "Z",  "Z", };
  haselmap["TRP_CH2"] = { "1.8",  "1.073",  "CZ2",  "CZ3",  "Z",  "Z", };
  haselmap["TRP_CZ3"] = { "1.8",  "1.073",  "CH2",  "CE3",  "Z",  "Z", };
  haselmap["TRP_CE3"] = { "1.8",  "1.073",  "CZ3",  "CD2",  "Z",  "Z", };
  haselmap["TRP_CD2"] = { "1.72",  "1.554",  "CE3",  "CE2",  "CG",  "Z", };
  haselmap["TRP_1HE"] = { "1.1",  "1.128",  "NE1",  "Z",  "Z",  "Z", };
  haselmap["TYR_N"] = { "1.55",  "1.028",  "H",  "CA",  "Z",  "Z", };
  haselmap["TYR_CA"] = { "1.7",  "2.149",  "N",  "C",  "CB",  "Z", };
  haselmap["TYR_C"] = { "1.72",  "1.554",  "CA",  "O",  "Z",  "Z", };
  haselmap["TYR_O"] = { "1.5",  "0.926",  "C",  "Z",  "Z",  "Z", };
  haselmap["TYR_H"] = { "1.1",  "1.128",  "N",  "Z",  "Z",  "Z", };
  haselmap["TYR_CB"] = { "1.9",  "1.045",  "CA",  "CG",  "Z",  "Z", };
  haselmap["TYR_CG"] = { "1.72",  "1.554",  "CB",  "CD1",  "CD2",  "Z", };
  haselmap["TYR_CD1"] = { "1.8",  "1.073",  "CG",  "CE1",  "Z",  "Z", };
  haselmap["TYR_CE1"] = { "1.8",  "1.073",  "CD1",  "CZ",  "Z",  "Z", };
  haselmap["TYR_CZ"] = { "1.72",  "1.554",  "CE1",  "OH",  "CE2",  "Z", };
  haselmap["TYR_OH"] = { "1.52",  "1.08",  "CZ",  "HH",  "Z",  "Z", };
  haselmap["TYR_HH"] = { "1.0",  "0.944",  "OH",  "Z",  "Z",  "Z", };
  haselmap["TYR_CE2"] = { "1.8",  "1.073",  "CZ",  "CD2",  "Z",  "Z", };
  haselmap["TYR_CD2"] = { "1.8",  "1.073",  "CE2",  "CG",  "Z",  "Z", };
  haselmap["VAL_N"] = { "1.55",  "1.028",  "H",  "CA",  "Z",  "Z", };
  haselmap["VAL_CA"] = { "1.7",  "2.149",  "N",  "C",  "CB",  "Z", };
  haselmap["VAL_C"] = { "1.72",  "1.554",  "CA",  "O",  "Z",  "Z", };
  haselmap["VAL_O"] = { "1.5",  "0.926",  "C",  "Z",  "Z",  "Z", };
  haselmap["VAL_H"] = { "1.1",  "1.128",  "N",  "Z",  "Z",  "Z", };
  haselmap["VAL_CB"] = { "1.8",  "1.276",  "CA",  "CG1",  "CG2",  "Z", };
  haselmap["VAL_CG1"] = { "2.0",  "0.88",  "CB",  "Z",  "Z",  "Z", };
  haselmap["VAL_CG2"] = { "2.0",  "0.88",  "CB",  "Z",  "Z",  "Z", };
  return haselmap;
}

//assigns SASA parameters to each atom reading from HASEL parameter database
void SASA_HASEL::readSASAparam() {

  for(unsigned i=0; i<natoms; i++) {
    SASAparam[i].clear();
    CONNECTparam[i].clear();
  }

  map<string, vector<std::string> > haselmap;
  haselmap = setupHASELparam();
  vector<std::string> HASELparamVector;
  string identifier;


  for(unsigned i=0; i<natoms; i++) {
    identifier = AtomResidueName[1][i]+"_"+AtomResidueName[0][i];
    if (haselmap.find(identifier)!=haselmap.end()) {
      HASELparamVector = haselmap.at(identifier);
      SASAparam[i].push_back (std::atof(HASELparamVector[0].c_str())+rs*10);
      SASAparam[i].push_back (std::atof(HASELparamVector[1].c_str()));
      CONNECTparam[i].push_back (HASELparamVector[2].c_str());
      CONNECTparam[i].push_back (HASELparamVector[3].c_str());
      CONNECTparam[i].push_back (HASELparamVector[4].c_str());
      CONNECTparam[i].push_back (HASELparamVector[5].c_str());
    }
  }


  for(unsigned i=0; i<natoms; i++) {
    if (SASAparam[i].size()==0 ) {
      if ((AtomResidueName[0][i][0]=='O') || (AtomResidueName[0][i][0]=='N') || (AtomResidueName[0][i][0]=='C') || (AtomResidueName[0][i][0]=='S')) {
        cout << "Could not find SASA paramaters for atom " << AtomResidueName[0][i] << " of residue " << AtomResidueName[1][i] << endl;
        error ("missing SASA parameters\n");
      }
    }
  }


}



//Max Surf values, used only if TYPE=TRANSFER
map<string, vector<double> > SASA_HASEL::setupMaxSurfMap() {
  // Max Surface Area for backbone and sidechain, in nm2
  map<string, vector<double> > valuemap;
  valuemap ["ALA"]= {0.56425,0.584851,};
  valuemap ["ARG"]= {0.498656,1.808093,};
  valuemap ["ASN"]= {0.473409,0.818394,};
  valuemap ["ASP"]= {0.477057,0.977303,};
  valuemap ["CYS"]= {0.507512,0.791483,};
  valuemap ["GLN"]= {0.485859,1.281534,};
  valuemap ["GLU"]= {0.495054,1.464718,};
  valuemap ["GLY"]= {0.658632,0,};
  valuemap ["HIS"]= {0.48194,1.118851,};
  valuemap ["ILE"]= {0.461283,1.450569,};
  valuemap ["LEU"]= {0.476315,1.498843,};
  valuemap ["LYS"]= {0.493533,1.619731,};
  valuemap ["MET"]= {0.507019,1.631904,};
  valuemap ["PHE"]= {0.457462, 1.275125,};
  valuemap ["PRO"]= {0.315865,0.859456,};
  valuemap ["SER"]= {0.48636,0.627233,};
  valuemap ["THR"]= {0.45064,0.91088,};
  valuemap ["TRP"]= {0.45762,1.366369,};
  valuemap ["TYR"]= {0.461826,1.425822,};
  valuemap ["VAL"]= {0.477054,1.149101,};
  return valuemap;
}



//reads maximum surface values per residue type and assigns values to each atom, only used if sasa_type = TRANSFER

void SASA_HASEL::readMaxSurf() {
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

void SASA_HASEL::readDeltaG() {

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

void SASA_HASEL::computeDeltaG() {

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
void SASA_HASEL::calcNlist() {
  if(!nopbc) makeWhole();

  for(unsigned i = 0; i < natoms; i++) {
    Nlist[i].clear();
  }

  for(unsigned i = 0; i < natoms; i++) {
    if (SASAparam[i].size()>0) {
      for (unsigned j = 0; j < i; j++) {
        if (SASAparam[j].size()>0) {
          const Vector Delta_ij_vec = delta( getPosition(i), getPosition(j) );
          double Delta_ij_mod = Delta_ij_vec.modulo()*10;
          double overlapD = SASAparam[i][0]+SASAparam[j][0];
          if (Delta_ij_mod < overlapD) {
            Nlist.at(i).push_back (j);
            Nlist.at(j).push_back (i);
          }
        }
      }
    }
  }
}


//calculates SASA according to Hasel et al., Tetrahedron Computer Methodology Vol. 1, No. 2, pp. 103-116, 1988
void SASA_HASEL::calculate() {
  if(!nopbc) makeWhole();

  if(getExchangeStep()) nl_update = 0;
  if (firstStepFlag ==0) {
    readPDB();
    readSASAparam();
  }
  if (nl_update == 0) {
    calcNlist();
  }


  auto* moldat = plumed.getActionSet().selectLatest<GenericMolInfo*>(this);
  if( ! moldat ) error("Unable to find MOLINFO in input");
  double Si, sasai, bij;
  double sasa = 0;
  vector<Vector> derivatives( natoms );
  for(unsigned i = 0; i < natoms; i++) {
    derivatives[i][0] = 0.;
    derivatives[i][1] = 0.;
    derivatives[i][2] = 0.;
  }

  Tensor virial;
  vector <double> ddij_di(3);
  vector <double> dbij_di(3);
  vector <double> dAijt_di(3);

  if( sasa_type==TOTAL ) {
    for(unsigned i = 0; i < natoms; i++) {
      if(SASAparam[i].size() > 0) {
        double ri = SASAparam[i][0];
        Si = 4*M_PI*ri*ri;
        sasai = 1.0;

        vector <vector <double> > derTerm( Nlist[i].size(), vector <double>(3));

        dAijt_di[0] = 0;
        dAijt_di[1] = 0;
        dAijt_di[2] = 0;
        int NumRes_i = moldat->getResidueNumber(atoms[i]);

        for (unsigned j = 0; j < Nlist[i].size(); j++) {
          double pij = 0.3516;

          int NumRes_j = moldat->getResidueNumber(atoms[Nlist[i][j]]);
          if (NumRes_i==NumRes_j) {
            if (CONNECTparam[i][0].compare(AtomResidueName[0][Nlist[i][j]])==0 || CONNECTparam[i][1].compare(AtomResidueName[0][Nlist[i][j]])==0 || CONNECTparam[i][2].compare(AtomResidueName[0][Nlist[i][j]])==0 || CONNECTparam[i][3].compare(AtomResidueName[0][Nlist[i][j]])==0) {
              pij = 0.8875;
            }
          }
          if ( abs(NumRes_i-NumRes_j) == 1 ) {
            if ((AtomResidueName[0][i] == "N"  && AtomResidueName[0][Nlist[i][j]]== "CA") || (AtomResidueName[0][Nlist[i][j]] == "N"  && AtomResidueName[0][i]== "CA")) {
              pij = 0.8875;
            }
          }

          const Vector d_ij_vec = delta( getPosition(i), getPosition(Nlist[i][j]) );
          double d_ij = d_ij_vec.modulo()*10;

          double rj = SASAparam[Nlist[i][j]][0];
          bij = M_PI*ri*(ri+rj-d_ij)*(1+(rj-ri)/d_ij); //Angstrom2

          sasai = sasai*(1-SASAparam[i][1]*pij*bij/Si); //nondimensional

          ddij_di[0] = -10*(getPosition(Nlist[i][j])[0]-getPosition(i)[0])/d_ij; //nondimensional
          ddij_di[1] = -10*(getPosition(Nlist[i][j])[1]-getPosition(i)[1])/d_ij;
          ddij_di[2] = -10*(getPosition(Nlist[i][j])[2]-getPosition(i)[2])/d_ij;

          dbij_di[0] = -M_PI*ri*ddij_di[0]*(1+(ri+rj)*(rj-ri)/(d_ij*d_ij)); //Angstrom
          dbij_di[1] = -M_PI*ri*ddij_di[1]*(1+(ri+rj)*(rj-ri)/(d_ij*d_ij));
          dbij_di[2] = -M_PI*ri*ddij_di[2]*(1+(ri+rj)*(rj-ri)/(d_ij*d_ij));

          dAijt_di[0] += -1/(Si/(SASAparam[i][1]*pij)-bij)*dbij_di[0]; //Angstrom-1
          dAijt_di[1] += -1/(Si/(SASAparam[i][1]*pij)-bij)*dbij_di[1];
          dAijt_di[2] += -1/(Si/(SASAparam[i][1]*pij)-bij)*dbij_di[2];

          derTerm[j][0] = 1/(Si/(SASAparam[i][1]*pij)-bij)*dbij_di[0]; //Angstrom-1
          derTerm[j][1] = 1/(Si/(SASAparam[i][1]*pij)-bij)*dbij_di[1];
          derTerm[j][2] = 1/(Si/(SASAparam[i][1]*pij)-bij)*dbij_di[2];

        }

        sasa += Si*sasai/100; //nm2

        derivatives[i][0] += Si*sasai/10*dAijt_di[0]; //nm
        derivatives[i][1] += Si*sasai/10*dAijt_di[1];
        derivatives[i][2] += Si*sasai/10*dAijt_di[2];

        for (unsigned j = 0; j < Nlist[i].size(); j++) {
          derivatives[Nlist[i][j]][0] += Si*sasai/10*derTerm[j][0]; //nm
          derivatives[Nlist[i][j]][1] += Si*sasai/10*derTerm[j][1];
          derivatives[Nlist[i][j]][2] += Si*sasai/10*derTerm[j][2];
        }
      }
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



      if(SASAparam[i].size() > 0) {
        double ri = SASAparam[i][0];
        Si = 4*M_PI*ri*ri;
        sasai = 1.0;

        vector <vector <double> > derTerm( Nlist[i].size(), vector <double>(3));

        dAijt_di[0] = 0;
        dAijt_di[1] = 0;
        dAijt_di[2] = 0;
        int NumRes_i = moldat->getResidueNumber(atoms[i]);

        for (unsigned j = 0; j < Nlist[i].size(); j++) {
          double pij = 0.3516;

          int NumRes_j = moldat->getResidueNumber(atoms[Nlist[i][j]]);
          if (NumRes_i==NumRes_j) {
            if (CONNECTparam[i][0].compare(AtomResidueName[0][Nlist[i][j]])==0 || CONNECTparam[i][1].compare(AtomResidueName[0][Nlist[i][j]])==0 || CONNECTparam[i][2].compare(AtomResidueName[0][Nlist[i][j]])==0 || CONNECTparam[i][3].compare(AtomResidueName[0][Nlist[i][j]])==0) {
              pij = 0.8875;
            }
          }
          if ( abs(NumRes_i-NumRes_j) == 1 ) {
            if ((AtomResidueName[0][i] == "N"  && AtomResidueName[0][Nlist[i][j]]== "CA") || (AtomResidueName[0][Nlist[i][j]] == "N"  && AtomResidueName[0][i]== "CA")) {
              pij = 0.8875;
            }
          }

          const Vector d_ij_vec = delta( getPosition(i), getPosition(Nlist[i][j]) );
          double d_ij = d_ij_vec.modulo()*10;

          double rj = SASAparam[Nlist[i][j]][0];
          bij = M_PI*ri*(ri+rj-d_ij)*(1+(rj-ri)/d_ij); //Angstrom2

          sasai = sasai*(1-SASAparam[i][1]*pij*bij/Si); //nondimensional

          ddij_di[0] = -10*(getPosition(Nlist[i][j])[0]-getPosition(i)[0])/d_ij; //nondimensional
          ddij_di[1] = -10*(getPosition(Nlist[i][j])[1]-getPosition(i)[1])/d_ij;
          ddij_di[2] = -10*(getPosition(Nlist[i][j])[2]-getPosition(i)[2])/d_ij;

          dbij_di[0] = -M_PI*ri*ddij_di[0]*(1+(ri+rj)*(rj-ri)/(d_ij*d_ij)); //Angstrom
          dbij_di[1] = -M_PI*ri*ddij_di[1]*(1+(ri+rj)*(rj-ri)/(d_ij*d_ij));
          dbij_di[2] = -M_PI*ri*ddij_di[2]*(1+(ri+rj)*(rj-ri)/(d_ij*d_ij));

          dAijt_di[0] += -1/(Si/(SASAparam[i][1]*pij)-bij)*dbij_di[0]; //Angstrom-1
          dAijt_di[1] += -1/(Si/(SASAparam[i][1]*pij)-bij)*dbij_di[1];
          dAijt_di[2] += -1/(Si/(SASAparam[i][1]*pij)-bij)*dbij_di[2];

          derTerm[j][0] = 1/(Si/(SASAparam[i][1]*pij)-bij)*dbij_di[0]; //Angstrom-1
          derTerm[j][1] = 1/(Si/(SASAparam[i][1]*pij)-bij)*dbij_di[1];
          derTerm[j][2] = 1/(Si/(SASAparam[i][1]*pij)-bij)*dbij_di[2];

        }

        if (AtomResidueName[0][i] == "N" || AtomResidueName[0][i] == "CA"  || AtomResidueName[0][i] == "C" || AtomResidueName[0][i] == "O" || AtomResidueName[0][i] == "H") {

          sasa += Si*sasai/MaxSurf[i][0]*DeltaG[natoms][0]; //kJ/mol


          derivatives[i][0] += Si*sasai*dAijt_di[0]/MaxSurf[i][0]*DeltaG[natoms][0]*10; //kJ/mol/nm
          derivatives[i][1] += Si*sasai*dAijt_di[1]/MaxSurf[i][0]*DeltaG[natoms][0]*10;
          derivatives[i][2] += Si*sasai*dAijt_di[2]/MaxSurf[i][0]*DeltaG[natoms][0]*10;
        }

        if (AtomResidueName[0][i] != "N" && AtomResidueName[0][i] != "CA"  && AtomResidueName[0][i] != "C" && AtomResidueName[0][i] != "O" && AtomResidueName[0][i] != "H") {
          sasa += Si*sasai/MaxSurf[i][1]*DeltaG[i][0]; //kJ/mol

          derivatives[i][0] += Si*sasai*dAijt_di[0]/MaxSurf[i][1]*DeltaG[i][0]*10; //kJ/mol/nm
          derivatives[i][1] += Si*sasai*dAijt_di[1]/MaxSurf[i][1]*DeltaG[i][0]*10;
          derivatives[i][2] += Si*sasai*dAijt_di[2]/MaxSurf[i][1]*DeltaG[i][0]*10;
        }


        for (unsigned j = 0; j < Nlist[i].size(); j++) {
          if (AtomResidueName[0][i] == "N" || AtomResidueName[0][i] == "CA"  || AtomResidueName[0][i] == "C" || AtomResidueName[0][i] == "O" || AtomResidueName[0][i] == "H") {
            derivatives[Nlist[i][j]][0] += Si*sasai*10*derTerm[j][0]/MaxSurf[i][0]*DeltaG[natoms][0]; //kJ/mol/nm
            derivatives[Nlist[i][j]][1] += Si*sasai*10*derTerm[j][1]/MaxSurf[i][0]*DeltaG[natoms][0];
            derivatives[Nlist[i][j]][2] += Si*sasai*10*derTerm[j][2]/MaxSurf[i][0]*DeltaG[natoms][0];
          }

          if (AtomResidueName[0][i] != "N" && AtomResidueName[0][i] != "CA"  && AtomResidueName[0][i] != "C" && AtomResidueName[0][i] != "O" && AtomResidueName[0][i] != "H") {
            derivatives[Nlist[i][j]][0] += Si*sasai*10*derTerm[j][0]/MaxSurf[i][1]*DeltaG[i][0]; //kJ/mol/nm
            derivatives[Nlist[i][j]][1] += Si*sasai*10*derTerm[j][1]/MaxSurf[i][1]*DeltaG[i][0];
            derivatives[Nlist[i][j]][2] += Si*sasai*10*derTerm[j][2]/MaxSurf[i][1]*DeltaG[i][0];
          }
        }
      }
    }
  }


  for(unsigned i=0; i<natoms; i++) {
    setAtomsDerivatives(i,derivatives[i]);
    virial -= Tensor(getPosition(i),derivatives[i]);
  }

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
