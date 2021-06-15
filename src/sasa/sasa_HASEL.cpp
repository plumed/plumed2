/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2019 The plumed team
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

#include "Sasa.h"
#include "ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/SetupMolInfo.h"
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
namespace Sasa {

//+PLUMEDOC COLVAR SASA_HASEL
/*
Calculates the solvent accessible surface area (SASA) of a protein molecule, or other properties related to it. The atoms for which the SASA is desired should be indicated with the keyword ATOMS, and a pdb file of the protein must be provided in input with the MOLINFO keyword. The algorithm described in (Hasel et al., Tetrahedron Computer Methodology Vol. 1, No. 2, pp. 103-116, 1988) is used for the calculation. The radius of the solvent is assumed to be 0.14 nm, which is the radius of water molecules. Using the keyword NL_STRIDE it is also possible to specify the frequency with which the neighbor list is updated (the default is every 10 steps).

Different properties can be calculated and selected using the TYPE keyword:
the total SASA (TOTAL); 
the free energy of transfer for the protein according to the transfer model (TRANSFER. This keyword can be used to compute the transfer of a protein to different temperatures, as detailed in Arsiccio and Shea, Protein Cold Denaturation in Implicit Solvent Simulations: A Transfer Free Energy Approach, JPCB, 2021). 



When the TRANSFER keyword is used, a file with the free energy of transfer values for the sidechains and backbone atoms should be provided (using the keyword DELTAGFILE). If the DELTAGFILE is not provided, the program computes the free energy of transfer values as if they had to take into account the effect of temperature according to approaches 2 or 3 in the paper: Arsiccio and Shea, Protein Cold Denaturation in Implicit Solvent Simulations: A Transfer Free Energy Approach, JPCB, 2021. Please read and cite this paper if using the transfer model for computing the effect of temperature in implicit solvent simulations. For this purpose, the keyword APPROACH should be added, and set to either 2 or 3. 

The SASA usually makes sense when atoms used for the calculation are all part of the same molecule. When running with periodic boundary conditions, the atoms should be in the proper periodic image. This is done automatically since PLUMED 2.2, by considering the ordered list of atoms and rebuilding the broken entities using a procedure that is equivalent to that done in \ref WHOLEMOLECULES. Notice that rebuilding is local to this action. This is different from \ref WHOLEMOLECULES which actually modifies the coordinates stored in PLUMED.

In case you want to recover the old behavior you should use the NOPBC flag.
In that case you need to take care that atoms are in the correct periodic image.

The SASA may also be computed using the SASA_LCPO collective variable, which makes us of the LCPO algorithm (Weiser J, Shenkin PS and Still WC. J. Comput. Chem. 1999 (20), 217-230). SASA_LCPO is more accurate then SASA_HASEL, but the computation is slower.


\par Examples

The following input tells plumed to print the total SASA for atoms 10 to 20 in a protein chain.
\plumedfile
SASA_HASEL TYPE=TOTAL ATOMS=10-20 NL_STRIDE=10 LABEL=sasa
PRINT ARG=sasa STRIDE=1 FILE=colvar
\endplumedfile


The following input tells plumed to compute the transfer free energy for the protein chain containing atoms 10 to 20. Such transfer free energy is then used as a bias in the simulation (e.g., implicit solvent simulations). The free energy of transfer values are read from a file called DeltaG.dat, that has the FOLLOWING FORMAT.
\plumedfile
SASA_HASEL TYPE=TRANSFER ATOMS=10-20 NL_STRIDE=10 DELTAGFILE=DeltaG.dat LABEL=sasa 

bias: BIASVALUE ARG=sasa

PRINT ARG=sasa,bias.* STRIDE=1 FILE=colvar
\endplumedfile

The following input tells plumed to compute the transfer free energy for the protein chain containing atoms 10 to 20. Such transfer free energy is then used as a bias in the simulation (e.g., implicit solvent simulations). The free energy of transfer values are computed according to "Arsiccio and Shea, Protein Cold Denaturation in Implicit Solvent Simulations: A Transfer Free Energy Approach, JPCB, 2021", and take into account the effect of temperature using approach 2 as described in the paper.
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
  keys.add("optional","DELTAGFILE","a file containing the free energy of transfer values for backbone and sidechains atoms. Necessary only if TYPE = TRANSFER. If TYPE = TRANSFER and no DELTAGFILE is provided, the free energy values are those describing the effect of temperature, and are computed using the temperature value passed by the MD engine");
  keys.add("optional","APPROACH","either approach 2 or 3. Necessary only if TYPE = TRANSFER and no DELTAGFILE is provided. If TYPE = TRANSFER and no DELTAGFILE is provided, the free energy values are those describing the effect of temperature, and the program must know if approach 2 or 3 (as described in Arsiccio and Shea, Protein Cold Denaturation in Implicit Solvent Simulations: A Transfer Free Energy Approach, JPCB, 2021) needs to be used to compute them");
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

void SASA_HASEL::readPDB(){
 vector<SetupMolInfo*> moldat = plumed.getActionSet().select<SetupMolInfo*>();
 AtomResidueName[0].clear(); 
 AtomResidueName[1].clear();

 for(unsigned i=0; i<natoms; i++) { 
   string Aname = moldat[0]->getAtomName(atoms[i]);
   string Rname = moldat[0]->getResidueName(atoms[i]);
   AtomResidueName[0].push_back (Aname);
   AtomResidueName[1].push_back (Rname);
}

}


//assigns SASA parameters to each atom reading from SASA_HASEL parameter file (hasel.dat, in HASEL folder)
	void SASA_HASEL::readSASAparam() {  

                                     for(unsigned i=0; i<natoms; i++) {  
                                      SASAparam[i].clear(); 
                                      CONNECTparam[i].clear();    }      
               
                    fstream SASAFile;
                    string SASAline; 
                    string path;
                    path = PLMD::config::getPlumedRoot();
                    SASAFile.open(path+"src/sasa/HASEL/hasel.dat");
                    if (SASAFile){
                    while(getline(SASAFile, SASAline)) {
                         if (!SASAline.empty()) {
                           	   std::vector<std::string> SASAtoken;
                                   split(SASAline, SASAtoken);
                                 if(SASAtoken.size() > 1) {
                                      for(unsigned i=0; i<natoms; i++) {              
                                     
                                      if(SASAtoken[0].compare(AtomResidueName[1][i])==0 && SASAtoken[1].compare(AtomResidueName[0][i])==0) {                           
                                SASAparam[i].push_back (std::atof(SASAtoken[3].c_str())+rs*10);
				 SASAparam[i].push_back (std::atof(SASAtoken[4].c_str()));
				 CONNECTparam[i].push_back (SASAtoken[5].c_str());
				 CONNECTparam[i].push_back (SASAtoken[6].c_str());
				 CONNECTparam[i].push_back (SASAtoken[7].c_str());
				 CONNECTparam[i].push_back (SASAtoken[8].c_str());
                                }}}
                      }
}}
             else error ("Unable to open SASA parameters file\n");

for(unsigned i=0; i<natoms; i++) {
      if (SASAparam[i].size()==0 ){
      if ((AtomResidueName[0][i][0]=='O') || (AtomResidueName[0][i][0]=='N') || (AtomResidueName[0][i][0]=='C') || (AtomResidueName[0][i][0]=='S')) {
          cout << "Could not find SASA paramaters for atom " << AtomResidueName[0][i] << " of residue " << AtomResidueName[1][i] << endl;
          error ("missing SASA parameters\n");}}}


}



//Max Surf values, used only if TYPE=TRANSFER
map<string, vector<double> > SASA_HASEL::setupMaxSurfMap() {
  // Max Surface Area for backbone and sidechain, in nm2
  map<string, vector<double> > valuemap;
  valuemap = {
    { "ALA", {
        0.56425,
        0.584851,
      }
    },
    { "ARG", {
        0.498656,
        1.808093,
      }
    },
    { "ASN", {
        0.473409,
        0.818394,
      }
    },
    { "ASP", {
        0.477057,
        0.977303,
      }
    },
    { "CYS", {
        0.507512,
        0.791483,
      }
    },
    { "GLN", {
        0.485859,
        1.281534,
      }
    },
    { "GLU", {
        0.495054,
        1.464718,
      }
    },
    { "GLY", {
        0.658632,
        0,
      }
    },
    { "HIS", {
        0.48194,
        1.118851,
      }
    },
    { "ILE", {
        0.461283,
        1.450569,
      }
    },
    { "LEU", {
        0.476315,
        1.498843,
      }
    },
    { "LYS", {
        0.493533,
        1.619731,
      }
    },
    { "MET", {
        0.507019,
        1.631904,
      }
    },
    { "PHE", {
        0.457462,
        1.275125,
      }
    },
    { "PRO", {
        0.315865,
        0.859456,
      }
    },
    { "SER", {
        0.48636,
        0.627233,
      }
    },
    { "THR", {
        0.45064,
        0.91088,
      }
    },
    { "TRP", {
        0.45762,
        1.366369,
      }
    },
    { "TYR", {
        0.461826,
        1.425822,
      }
    },
    { "VAL", {
        0.477054,
        1.149101,
      }
    }
   };
return valuemap;
}



//reads maximum surface values per residue type and assigns values to each atom, only used if sasa_type = TRANSFER

 void SASA_HASEL::readMaxSurf(){
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

void SASA_HASEL::readDeltaG(){
 
      for(unsigned i=0; i<natoms; i++) { 
              DeltaG[i].clear();   }
 
 string DeltaGline;
 fstream DeltaGFile;
                DeltaGFile.open(DeltaGValues);
                if (DeltaGFile){                
                int backboneflag = 0;
                    while(getline(DeltaGFile, DeltaGline)) {
                             if(!DeltaGline.empty()) { 
                           	std::vector<std::string> DeltaGtoken;
                                split(DeltaGline, DeltaGtoken);  
                                for(unsigned i=0; i<natoms; i++) {                  
                                      if (DeltaGtoken[0].compare(AtomResidueName[1][i])==0 ) {                           
                                         DeltaG[i].push_back (std::atof(DeltaGtoken[1].c_str()));
                                           }}
                                if (DeltaGtoken[0].compare("BACKBONE")==0 ) { 
                                         backboneflag = 1;                                
                                         DeltaG[natoms].push_back (std::atof(DeltaGtoken[1].c_str()));
                                           }}}
                                        if ( backboneflag == 0) error("Cannot find backbone value in Delta G parameters file\n");}
                                 else error("Cannot open DeltaG file");

for(unsigned i=0; i<natoms; i++) {
    if (DeltaG[i].size()==0 ){
          cout << "Delta G value for residue " << AtomResidueName[1][i] << " not found " << endl;
          error ("missing Delta G parameter\n");}}

}

//computes free energy values for the sidechains and for the backbone, and assigns values to each atom. Only used if sasa_type = TRANSFER, and if no DELTAGFILE is provided. In this case, the free energy values are those describing the effect of temperature, and the program must know if approach 2 or 3 (as described in Arsiccio and Shea, Protein Cold Denaturation in Implicit Solvent Simulations: A Transfer Free Energy Approach, JPCB, 2021) needs to be used to compute them.

void SASA_HASEL::computeDeltaG(){
 
      for(unsigned i=0; i<natoms; i++) { 
              DeltaG[i].clear();   }

double T;
T = plumed.getAtoms().getKbT()/plumed.getAtoms().getKBoltzmann();

if (T != Ti) {
for(unsigned i=0; i<natoms; i++) { 
   if (approach==2) {
      if (AtomResidueName[1][i].compare("ALA")==0) {
         DeltaG[i].push_back (-2.995/1000*std::pow(T,2)+1.808*T-272.895);}
      if (AtomResidueName[1][i].compare("ARG")==0) {
         DeltaG[i].push_back (-3.182/1000*std::pow(T,2)+1.894*T-282.032);}   
      if (AtomResidueName[1][i].compare("ASN")==0) {
         DeltaG[i].push_back (-1.047/1000*std::pow(T,2)+0.6068*T-87.846);} 
      if (AtomResidueName[1][i].compare("ASP")==0) {
         DeltaG[i].push_back (-0.1794/1000*std::pow(T,2)+0.1091*T-16.526);}         
      if (AtomResidueName[1][i].compare("CYS")==0) {
         DeltaG[i].push_back (-3.09/1000*std::pow(T,2)+1.835*T-272.26);}  
      if (AtomResidueName[1][i].compare("GLN")==0) {
         DeltaG[i].push_back (-2.23/1000*std::pow(T,2)+1.335*T-199.707);}  
      if (AtomResidueName[1][i].compare("GLU")==0) {
         DeltaG[i].push_back (-1.511/1000*std::pow(T,2)+0.8904*T-131.168);}
      if (AtomResidueName[1][i].compare("GLY")==0) {
         DeltaG[i].push_back (0);} 
      if (AtomResidueName[1][i].compare("HIS")==0) {
         DeltaG[i].push_back (-3.482/1000*std::pow(T,2)+2.084*T-311.694);}
      if (AtomResidueName[1][i].compare("ILE")==0) {
         DeltaG[i].push_back (-6.364/1000*std::pow(T,2)+3.8*T-567.444);}
      if (AtomResidueName[1][i].compare("LEU")==0) {
         DeltaG[i].push_back (-7.466/1000*std::pow(T,2)+4.417*T-653.394);}  
      if (AtomResidueName[1][i].compare("LYS")==0) {
         DeltaG[i].push_back (-2.091/1000*std::pow(T,2)+1.2458*T-185.549);}  
      if (AtomResidueName[1][i].compare("MET")==0) {
         DeltaG[i].push_back (-3.807/1000*std::pow(T,2)+2.272*T-339.007);}  
      if (AtomResidueName[1][i].compare("PHE")==0) {
         DeltaG[i].push_back (-7.828/1000*std::pow(T,2)+4.644*T-688.874);} 
      if (AtomResidueName[1][i].compare("PRO")==0) {
         DeltaG[i].push_back (-2.347/1000*std::pow(T,2)+1.411*T-212.059);} 
      if (AtomResidueName[1][i].compare("SER")==0) {
         DeltaG[i].push_back (1.813/1000*std::pow(T,2)-1.05*T+151.957);}  
      if (AtomResidueName[1][i].compare("THR")==0) {
         DeltaG[i].push_back (-2.64/1000*std::pow(T,2)+1.591*T-239.516);}  
      if (AtomResidueName[1][i].compare("TRP")==0) {
         DeltaG[i].push_back (-11.66/1000*std::pow(T,2)+6.916*T-1025.293);}  
      if (AtomResidueName[1][i].compare("TYR")==0) {
         DeltaG[i].push_back (-7.513/1000*std::pow(T,2)+4.478*T-667.261);} 
      if (AtomResidueName[1][i].compare("VAL")==0) {
         DeltaG[i].push_back (-4.902/1000*std::pow(T,2)+2.921*T-435.309);} 
         DeltaG[natoms].push_back (-0.6962/1000*std::pow(T,2)+0.4426*T-70.015);      
         }
    if (approach==3) {
      if (AtomResidueName[1][i].compare("ALA")==0) {
         DeltaG[i].push_back (-2.995/1000*std::pow(T,2)+1.808*T-272.105);}
      if (AtomResidueName[1][i].compare("ARG")==0) {
         DeltaG[i].push_back (-3.182/1000*std::pow(T,2)+1.894*T-284.086);}   
      if (AtomResidueName[1][i].compare("ASN")==0) {
         DeltaG[i].push_back (-1.047/1000*std::pow(T,2)+0.6068*T-90.597);} 
      if (AtomResidueName[1][i].compare("ASP")==0) {
         DeltaG[i].push_back (-0.1794/1000*std::pow(T,2)+0.1091*T-19.143);}         
      if (AtomResidueName[1][i].compare("CYS")==0) {
         DeltaG[i].push_back (-3.09/1000*std::pow(T,2)+1.835*T-268.527);}  
      if (AtomResidueName[1][i].compare("GLN")==0) {
         DeltaG[i].push_back (-2.23/1000*std::pow(T,2)+1.335*T-201.559);}  
      if (AtomResidueName[1][i].compare("GLU")==0) {
         DeltaG[i].push_back (-1.511/1000*std::pow(T,2)+0.8904*T-133.543);}
      if (AtomResidueName[1][i].compare("GLY")==0) {
         DeltaG[i].push_back (0);} 
      if (AtomResidueName[1][i].compare("HIS")==0) {
         DeltaG[i].push_back (-3.482/1000*std::pow(T,2)+2.084*T-315.398);}
      if (AtomResidueName[1][i].compare("ILE")==0) {
         DeltaG[i].push_back (-6.364/1000*std::pow(T,2)+3.8*T-564.825);}
      if (AtomResidueName[1][i].compare("LEU")==0) {
         DeltaG[i].push_back (-7.466/1000*std::pow(T,2)+4.417*T-651.483);}  
      if (AtomResidueName[1][i].compare("LYS")==0) {
         DeltaG[i].push_back (-2.091/1000*std::pow(T,2)+1.2458*T-187.485);}  
      if (AtomResidueName[1][i].compare("MET")==0) {
         DeltaG[i].push_back (-3.807/1000*std::pow(T,2)+2.272*T-339.242);}  
      if (AtomResidueName[1][i].compare("PHE")==0) {
         DeltaG[i].push_back (-7.828/1000*std::pow(T,2)+4.644*T-687.134);} 
      if (AtomResidueName[1][i].compare("PRO")==0) {
         DeltaG[i].push_back (-2.347/1000*std::pow(T,2)+1.411*T-214.211);} 
      if (AtomResidueName[1][i].compare("SER")==0) {
         DeltaG[i].push_back (1.813/1000*std::pow(T,2)-1.05*T+150.289);}  
      if (AtomResidueName[1][i].compare("THR")==0) {
         DeltaG[i].push_back (-2.64/1000*std::pow(T,2)+1.591*T-240.267);}  
      if (AtomResidueName[1][i].compare("TRP")==0) {
         DeltaG[i].push_back (-11.66/1000*std::pow(T,2)+6.916*T-1024.284);}  
      if (AtomResidueName[1][i].compare("TYR")==0) {
         DeltaG[i].push_back (-7.513/1000*std::pow(T,2)+4.478*T-666.484);} 
      if (AtomResidueName[1][i].compare("VAL")==0) {
         DeltaG[i].push_back (-4.902/1000*std::pow(T,2)+2.921*T-433.013);} 
         DeltaG[natoms].push_back (-0.6927/1000*std::pow(T,2)+0.4404*T-68.724);      
         }        
         }
}

Ti = T;

 if (firstStepFlag ==0) { 
 if (approach!=2 && approach!=3){
         cout << "Unknown approach " << approach << endl;}
for(unsigned i=0; i<natoms; i++) {
    if (DeltaG[i].size()==0 ){
          cout << "Delta G value for residue " << AtomResidueName[1][i] << " not found " << endl;
          error ("missing Delta G parameter\n");}}}
}


//calculates neighbor list
void SASA_HASEL::calcNlist(){
  if(!nopbc) makeWhole();

  for(unsigned i = 0; i < natoms; i++) {
    Nlist[i].clear();}

  for(unsigned i = 0; i < natoms; i++) {
      if (SASAparam[i].size()>0){
        for (unsigned j = 0; j < i; j++) {
              if (SASAparam[j].size()>0){
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
    calcNlist();    }

  
  vector<SetupMolInfo*> moldat = plumed.getActionSet().select<SetupMolInfo*>();
  double Si, sasai, bij;
  double sasa = 0;
  vector<Vector> derivatives( natoms );
  for(unsigned i = 0; i < natoms; i++) {
    derivatives[i][0] = 0.;
    derivatives[i][1] = 0.;
    derivatives[i][2] = 0.;}

  Tensor virial;
  vector <double> ddij_di(3);
  vector <double> dbij_di(3);
  vector <double> dAijt_di(3);

  if( sasa_type==TOTAL ) {
    for(unsigned i = 0; i < natoms; i++) {
        if(SASAparam[i].size() > 0){
         double ri = SASAparam[i][0];
         Si = 4*M_PI*ri*ri;
         sasai = 1.0;
         
         vector <vector <double> > derTerm( Nlist[i].size(), vector <double>(3));
         
         dAijt_di[0] = 0; 
         dAijt_di[1] = 0; 
         dAijt_di[2] = 0;        
         int NumRes_i = moldat[0]->getResidueNumber(atoms[i]);
            
         for (unsigned j = 0; j < Nlist[i].size(); j++) {
             double pij = 0.3516;
 
             int NumRes_j = moldat[0]->getResidueNumber(atoms[Nlist[i][j]]);
             if (NumRes_i==NumRes_j){
                 if (CONNECTparam[i][0].compare(AtomResidueName[0][Nlist[i][j]])==0 || CONNECTparam[i][1].compare(AtomResidueName[0][Nlist[i][j]])==0 || CONNECTparam[i][2].compare(AtomResidueName[0][Nlist[i][j]])==0 || CONNECTparam[i][3].compare(AtomResidueName[0][Nlist[i][j]])==0) {
                     pij = 0.8875;}}
             if ( abs(NumRes_i-NumRes_j) == 1 ){
                if ((AtomResidueName[0][i] == "N"  && AtomResidueName[0][Nlist[i][j]]== "CA") || (AtomResidueName[0][Nlist[i][j]] == "N"  && AtomResidueName[0][i]== "CA")) {
                  pij = 0.8875;
                  }} 
                                           
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
         
            
            
        if(SASAparam[i].size() > 0){
         double ri = SASAparam[i][0];
         Si = 4*M_PI*ri*ri;
         sasai = 1.0;
         
         vector <vector <double> > derTerm( Nlist[i].size(), vector <double>(3));

         dAijt_di[0] = 0; 
         dAijt_di[1] = 0; 
         dAijt_di[2] = 0;        
         int NumRes_i = moldat[0]->getResidueNumber(atoms[i]);
         
         for (unsigned j = 0; j < Nlist[i].size(); j++) {
             double pij = 0.3516;
 
             int NumRes_j = moldat[0]->getResidueNumber(atoms[Nlist[i][j]]);
             if (NumRes_i==NumRes_j){
                 if (CONNECTparam[i][0].compare(AtomResidueName[0][Nlist[i][j]])==0 || CONNECTparam[i][1].compare(AtomResidueName[0][Nlist[i][j]])==0 || CONNECTparam[i][2].compare(AtomResidueName[0][Nlist[i][j]])==0 || CONNECTparam[i][3].compare(AtomResidueName[0][Nlist[i][j]])==0) {
                     pij = 0.8875;}}
             if ( abs(NumRes_i-NumRes_j) == 1 ){
                if ((AtomResidueName[0][i] == "N"  && AtomResidueName[0][Nlist[i][j]]== "CA") || (AtomResidueName[0][Nlist[i][j]] == "N"  && AtomResidueName[0][i]== "CA")) {
                  pij = 0.8875;
                  }} 

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
             
             if (AtomResidueName[0][i] == "N" || AtomResidueName[0][i] == "CA"  || AtomResidueName[0][i] == "C" || AtomResidueName[0][i] == "O" || AtomResidueName[0][i] == "H"){
             
             sasa += Si*sasai/MaxSurf[i][0]*DeltaG[natoms][0]; //kJ/mol
                         

             derivatives[i][0] += Si*sasai*dAijt_di[0]/MaxSurf[i][0]*DeltaG[natoms][0]*10; //kJ/mol/nm
             derivatives[i][1] += Si*sasai*dAijt_di[1]/MaxSurf[i][0]*DeltaG[natoms][0]*10; 
             derivatives[i][2] += Si*sasai*dAijt_di[2]/MaxSurf[i][0]*DeltaG[natoms][0]*10;}
             
             if (AtomResidueName[0][i] != "N" && AtomResidueName[0][i] != "CA"  && AtomResidueName[0][i] != "C" && AtomResidueName[0][i] != "O" && AtomResidueName[0][i] != "H"){
             sasa += Si*sasai/MaxSurf[i][1]*DeltaG[i][0]; //kJ/mol            

             derivatives[i][0] += Si*sasai*dAijt_di[0]/MaxSurf[i][1]*DeltaG[i][0]*10; //kJ/mol/nm
             derivatives[i][1] += Si*sasai*dAijt_di[1]/MaxSurf[i][1]*DeltaG[i][0]*10; 
             derivatives[i][2] += Si*sasai*dAijt_di[2]/MaxSurf[i][1]*DeltaG[i][0]*10;}
            
       
       for (unsigned j = 0; j < Nlist[i].size(); j++) {
             if (AtomResidueName[0][i] == "N" || AtomResidueName[0][i] == "CA"  || AtomResidueName[0][i] == "C" || AtomResidueName[0][i] == "O" || AtomResidueName[0][i] == "H"){          
             derivatives[Nlist[i][j]][0] += Si*sasai*10*derTerm[j][0]/MaxSurf[i][0]*DeltaG[natoms][0]; //kJ/mol/nm
             derivatives[Nlist[i][j]][1] += Si*sasai*10*derTerm[j][1]/MaxSurf[i][0]*DeltaG[natoms][0]; 
             derivatives[Nlist[i][j]][2] += Si*sasai*10*derTerm[j][2]/MaxSurf[i][0]*DeltaG[natoms][0];}
             
             if (AtomResidueName[0][i] != "N" && AtomResidueName[0][i] != "CA"  && AtomResidueName[0][i] != "C" && AtomResidueName[0][i] != "O" && AtomResidueName[0][i] != "H"){          
             derivatives[Nlist[i][j]][0] += Si*sasai*10*derTerm[j][0]/MaxSurf[i][1]*DeltaG[i][0]; //kJ/mol/nm
             derivatives[Nlist[i][j]][1] += Si*sasai*10*derTerm[j][1]/MaxSurf[i][1]*DeltaG[i][0]; 
             derivatives[Nlist[i][j]][2] += Si*sasai*10*derTerm[j][2]/MaxSurf[i][1]*DeltaG[i][0];}
             }                        
}
}
}


  for(unsigned i=0;i<natoms;i++){ 
   setAtomsDerivatives(i,derivatives[i]);
   virial -= Tensor(getPosition(i),derivatives[i]);   }
   
  setBoxDerivatives(virial);
  setValue(sasa);    
  firstStepFlag = 1; 
  ++nl_update;
  if (nl_update == stride) {
    nl_update = 0;
  }
 // setBoxDerivativesNoPbc();
}

}//namespace PLMD
}//namespace Sasa
