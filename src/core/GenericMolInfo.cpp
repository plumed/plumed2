/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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
#include "GenericMolInfo.h"
#include "ActionRegister.h"
#include "ActionSet.h"
#include "PlumedMain.h"
#include "tools/MolDataClass.h"
#include "tools/PDB.h"
#include "tools/Communicator.h"
#include "config/Config.h"
#include <algorithm>
#include <map>

namespace PLMD {


/*
This action is defined in core/ as it is used by other actions.
Anyway, it is registered in setup/, so that excluding that module from
compilation will exclude it from plumed.
*/


void GenericMolInfo::registerKeywords( Keywords& keys ) {
  ActionSetup::registerKeywords(keys);
  keys.add("compulsory","STRUCTURE","a file in pdb format containing a reference structure. "
           "This is used to defines the atoms in the various residues, chains, etc . "
           "For more details on the PDB file format visit http://www.wwpdb.org/docs.html");
  keys.add("compulsory","MOLTYPE","protein","what kind of molecule is contained in the pdb file - usually not needed since protein/RNA/DNA are compatible");
  keys.add("compulsory","PYTHON_BIN","default","python interpreter");
  keys.add("atoms","CHAIN","(for masochists ( mostly Davide Branduardi ) ) The atoms involved in each of the chains of interest in the structure.");
  keys.add("hidden","STRIDE","frequency for resetting the python interpreter. Should be 1.");
  keys.addFlag("WHOLE", false, "The reference structure is whole, i.e. not broken by PBC");
}

GenericMolInfo::~GenericMolInfo() {
// empty destructor to delete unique_ptr
}

GenericMolInfo::GenericMolInfo( const ActionOptions&ao ):
  Action(ao),
  ActionAnyorder(ao),
  ActionPilot(ao),
  ActionAtomistic(ao),
  iswhole_(false) {
  plumed_assert(getStride()==1);
  // Read what is contained in the pdb file
  parse("MOLTYPE",mytype);

  // check if whole
  parseFlag("WHOLE", iswhole_);

  auto* latestmoldat=plumed.getActionSet().selectLatest<GenericMolInfo*>(this);
  if( latestmoldat ) {
    log<<"  overriding last MOLINFO with label " << latestmoldat->getLabel()<<"\n";
  }

  std::vector<AtomNumber> backbone;
  parseAtomList("CHAIN",backbone);
  if( read_backbone.size()==0 ) {
    for(unsigned i=1;; ++i) {
      parseAtomList("CHAIN",i,backbone);
      if( backbone.size()==0 ) {
        break;
      }
      read_backbone.push_back(backbone);
      backbone.resize(0);
    }
  } else {
    read_backbone.push_back(backbone);
  }
  if( read_backbone.size()==0 ) {
    parse("STRUCTURE",reference);

    if( ! pdb.read(reference,usingNaturalUnits(),0.1/getUnits().getLength())) {
      plumed_merror("missing input file " + reference );
    }

    std::vector<std::string> chains;
    pdb.getChainNames( chains );
    log.printf("  pdb file named %s contains %u chains \n",reference.c_str(), static_cast<unsigned>(chains.size()));
    for(unsigned i=0; i<chains.size(); ++i) {
      int start,end;
      std::string errmsg;
      pdb.getResidueRange( chains[i], start, end, errmsg );
      if( errmsg.length()!=0 ) {
        error( errmsg );
      }
      AtomNumber astart,aend;
      pdb.getAtomRange( chains[i], astart, aend, errmsg );
      if( errmsg.length()!=0 ) {
        error( errmsg );
      }
      log.printf("  chain named %s contains residues %d to %d and atoms %u to %u \n",chains[i].c_str(),start,end,astart.serial(),aend.serial());
    }

    std::string python_bin;
    parse("PYTHON_BIN",python_bin);
    if(python_bin=="no") {
      log<<"  python interpreter disabled\n";
    } else {
      pythonCmd=config::getEnvCommand();
      if(python_bin!="default") {
        log<<"  forcing python interpreter: "<<python_bin<<"\n";
        pythonCmd+=" env PLUMED_PYTHON_BIN="+python_bin;
      }
      bool sorted=true;
      const auto & at=pdb.getAtomNumbers();
      for(unsigned i=0; i<at.size(); i++) {
        if(at[i].index()!=i) {
          sorted=false;
        }
      }
      if(!sorted) {
        log<<"  PDB is not sorted, python interpreter will be disabled\n";
      } else if(!Subprocess::available()) {
        log<<"  subprocess is not available, python interpreter will be disabled\n";
      } else {
        enablePythonInterpreter=true;
      }
    }
  }
}

std::map<std::string,std::string> GenericMolInfo::getSpecialKeywords() {
  std::map<std::string,std::string> mapkeys;
  mapkeys.insert( std::pair<std::string,std::string>("@mda:","atom selection built using syntax for <a href=\\\"https://docs.mdanalysis.org/stable/documentation_pages/selections.html\\\">MDAnalysis</a>") );
  mapkeys.insert( std::pair<std::string,std::string>("@mdt:","atom selection built using syntax for <a href=\\\"https://www.mdtraj.org/1.9.8.dev0/index.html\\\">mdtraj</a>") );
  mapkeys.insert( std::pair<std::string,std::string>("@vmdexec:","atom selection built using syntax for <a href=\\\"https://www.ks.uiuc.edu/Research/vmd/\\\">VMD</a>") );
  mapkeys.insert( std::pair<std::string,std::string>("@vmd:","atom selection built using syntax for <a href=\\\"https://www.ks.uiuc.edu/Research/vmd/\\\">VMD</a> python module") );
  mapkeys.insert( std::pair<std::string,std::string>("@nucleic","all atoms that are part of a DNA or RNA molecule") );
  mapkeys.insert( std::pair<std::string,std::string>("@protein","all atoms that are part of a protein") );
  mapkeys.insert( std::pair<std::string,std::string>("@water","all water molecules") );
  mapkeys.insert( std::pair<std::string,std::string>("@ions","all the ions") );
  mapkeys.insert( std::pair<std::string,std::string>("@hydrogens","all hydrogen atoms") );
  mapkeys.insert( std::pair<std::string,std::string>("@nonhydrogens","all non hydrogen atoms") );
  mapkeys.insert( std::pair<std::string,std::string>("@phi-","the four atoms that are required to calculate the phi dihedral for") );
  mapkeys.insert( std::pair<std::string,std::string>("@psi-","the four atoms that are required to calculate the psi dihedral for") );
  mapkeys.insert( std::pair<std::string,std::string>("@omega-","the four atoms that are required to calculate the omega dihedral for") );
  mapkeys.insert( std::pair<std::string,std::string>("@chi1-","the four atoms that are required to calculate the chi1 dihedral for") );
  mapkeys.insert( std::pair<std::string,std::string>("@chi2-","the four atoms that are required to calculate the chi2 dihedral for") );
  mapkeys.insert( std::pair<std::string,std::string>("@chi3-","the four atoms that are required to calculate the chi3 dihedral for") );
  mapkeys.insert( std::pair<std::string,std::string>("@chi4-","the four atoms that are required to calculate the chi4 dihedral for") );
  mapkeys.insert( std::pair<std::string,std::string>("@chi5-","the four atoms that are required to calculate the chi5 dihedral for") );
  mapkeys.insert( std::pair<std::string,std::string>("@sidechain-","the protein sidechain atoms in") );
  mapkeys.insert( std::pair<std::string,std::string>("@back-","the protein/dna/rna backbone atoms in") );
  mapkeys.insert( std::pair<std::string,std::string>("@alpha-","the four atoms that are required to calculate the alpha backbone dihedral for") );
  mapkeys.insert( std::pair<std::string,std::string>("@beta-","the four atoms that are required to calculate the beta backbone dihedral for") );
  mapkeys.insert( std::pair<std::string,std::string>("@gamma-","the four atoms that are required to calculate the gamma backbone dihedral for") );
  mapkeys.insert( std::pair<std::string,std::string>("@delta-","the four atoms that are required to calculate the delta backbone dihedral for") );
  mapkeys.insert( std::pair<std::string,std::string>("@epsilon-","the four atoms that are required to calculate the backbone epsilon dihedral for") );
  mapkeys.insert( std::pair<std::string,std::string>("@zeta-","the four atoms that are required to calculate the zeta backbone dihedral for") );
  mapkeys.insert( std::pair<std::string,std::string>("@v0-","the four atoms that are required to calculate the v0 sugar dihedral for") );
  mapkeys.insert( std::pair<std::string,std::string>("@v1-","the four atoms that are required to calculate the v1 sugar dihedral for") );
  mapkeys.insert( std::pair<std::string,std::string>("@v2-","the four atoms that are required to calculate the v2 sugar dihedral for") );
  mapkeys.insert( std::pair<std::string,std::string>("@v3-","the four atoms that are required to calculate the v3 sugar dihedral for") );
  mapkeys.insert( std::pair<std::string,std::string>("@v4-","the four atoms that are required to calculate the v4 sugar dihedral for") );
  mapkeys.insert( std::pair<std::string,std::string>("@chi-","the four atoms that are required to alcullate the chi dihedral for") );
  mapkeys.insert( std::pair<std::string,std::string>("@sugar-","the heavy atoms of the sugar in") );
  mapkeys.insert( std::pair<std::string,std::string>("@base-","the heavy atoms of the base in") );
  mapkeys.insert( std::pair<std::string,std::string>("@lcs-","an ordered triplet of atoms on the 6-membered ring of the nucleobase in") );
  return mapkeys;
}

void GenericMolInfo::getBackbone( std::vector<std::string>& restrings, const std::string& fortype, std::vector< std::vector<AtomNumber> >& backbone ) {
  if( fortype!=mytype ) {
    error("cannot calculate a variable designed for " + fortype + " molecules for molecule type " + mytype );
  }
  if( MolDataClass::numberOfAtomsPerResidueInBackbone( mytype )==0 ) {
    error("backbone is not defined for molecule type " + mytype );
  }

  if( read_backbone.size()!=0 ) {
    if( restrings.size()!=1 ) {
      error("cannot interpret anything other than all for residues when using CHAIN keywords");
    }
    if( restrings[0]!="all" ) {
      error("cannot interpret anything other than all for residues when using CHAIN keywords");
    }
    backbone.resize( read_backbone.size() );
    for(unsigned i=0; i<read_backbone.size(); ++i) {
      backbone[i].resize( read_backbone[i].size() );
      for(unsigned j=0; j<read_backbone[i].size(); ++j) {
        backbone[i][j]=read_backbone[i][j];
      }
    }
  } else {
    bool useter=false; // This is used to deal with terminal groups in WHOLEMOLECULES
    // residues are identified by (chain, residue number), so that chains that use
    // the same residue numbers (e.g. all starting from 1) are not mixed up
    std::vector< std::pair<std::string,int> > residues;
    std::vector<std::string> chains;
    {
      std::vector<std::string> allchains;
      pdb.getChainNames( allchains );
      // a chain can appear in more than one block of the pdb file: list it once
      for(const auto & c : allchains) {
        if( std::find(chains.begin(),chains.end(),c)==chains.end() ) {
          chains.push_back(c);
        }
      }
    }
    if( restrings.size()==1 ) {
      useter=( restrings[0].find("ter")!=std::string::npos );
    }
    if( restrings.size()==1 && restrings[0].find("all")!=std::string::npos ) {
      // list explicitly the residues that are present in each chain:
      // this works also with negative residue numbers and with gaps in the numbering
      for(const auto & c : chains) {
        const std::vector<int> res = pdb.getResiduesInChain( c );
        if( res.empty() ) {
          continue;
        }
        std::size_t first=0, last=res.size();
        if( !useter ) {
          if( MolDataClass::isTerminalGroup( mytype, pdb.getResidueName( res.front(), c ) ) ) {
            first++;
          }
          if( last>first && MolDataClass::isTerminalGroup( mytype, pdb.getResidueName( res.back(), c ) ) ) {
            last--;
          }
        }
        for(std::size_t k=first; k<last; ++k) {
          residues.emplace_back( c, res[k] );
        }
      }
    } else {
      Tools::interpretRanges(restrings);
      // chains containing each residue number
      std::map< int, std::vector<std::string> > chains_of_residue;
      for(const auto & c : chains) {
        for(const int r : pdb.getResiduesInChain( c )) {
          chains_of_residue[r].push_back( c );
        }
      }
      std::string ambiguous;
      for(const auto & rs : restrings) {
        int r;
        if( !Tools::convertNoexcept( rs, r ) ) {
          error("cannot interpret " + rs + " as a residue number");
        }
        const auto p=chains_of_residue.find( r );
        if( p==chains_of_residue.end() ) {
          error("residue " + rs + " is not in the reference pdb file");
        }
        // a residue number present in more than one chain is taken from the first chain, as in previous versions
        if( p->second.size()>1 ) {
          ambiguous += " " + rs;
        }
        residues.emplace_back( p->second[0], r );
      }
      if( ambiguous.length()>0 ) {
        warning("residues" + ambiguous + " are present in more than one chain of the reference pdb file: they are taken from the first chain that contains them. Use RESIDUES=all to include all the chains");
      }
    }
    if( residues.empty() ) {
      error("no residues found in the reference pdb file");
    }

    // Convert the list of involved residues into a list of segments of consecutive residues of the same chain
    std::vector< std::vector< std::pair<std::string,int> > > segments;
    std::vector< std::pair<std::string,int> > thissegment;
    thissegment.push_back( residues[0] );
    for(unsigned i=1; i<residues.size(); ++i) {
      if( residues[i-1].second+1!=residues[i].second || residues[i-1].first!=residues[i].first ) {
        segments.push_back(thissegment);
        thissegment.resize(0);
      }
      thissegment.push_back( residues[i] );
    }
    segments.push_back( thissegment );

    // And now get the backbone atoms from each segment
    backbone.resize( segments.size() );
    std::vector<AtomNumber> atomnumbers;
    for(unsigned i=0; i<segments.size(); ++i) {
      for(unsigned j=0; j<segments[i].size(); ++j) {
        const std::string & chain=segments[i][j].first;
        const int resnum=segments[i][j].second;
        std::string num;
        Tools::convert( resnum, num );
        num += " of chain " + chain;
        std::string resname=pdb.getResidueName( resnum, chain );
        if( !MolDataClass::allowedResidue(mytype, resname) ) {
          error("residue " + num + " is not recognized for moltype " + mytype );
        }
        if( !useter && MolDataClass::isTerminalGroup( mytype, resname ) ) {
          error("residue " + num + " appears to be a terminal group");
        }
        if( resname=="GLY" ) {
          warning("GLY residues are achiral - assuming HA1 atom is in CB position");
        }
        MolDataClass::getBackboneForResidue( mytype, resnum, chain, pdb, atomnumbers );
        if( atomnumbers.size()==0 ) {
          error("Could not find required backbone atom in residue number " + num );
        } else {
          for(unsigned k=0; k<atomnumbers.size(); ++k) {
            backbone[i].push_back( atomnumbers[k] );
          }
        }
        atomnumbers.resize(0);
      }
    }
  }
}

void GenericMolInfo::interpretSymbol( const std::string& symbol, std::vector<AtomNumber>& atoms ) {
  if(Tools::startWith(symbol,"mdt:") || Tools::startWith(symbol,"mda:") || Tools::startWith(symbol,"vmd:") || Tools::startWith(symbol,"vmdexec:")) {

    plumed_assert(enablePythonInterpreter);

    log<<"  symbol " + symbol + " will be sent to python interpreter\n";
    if(!selector_running) {
      log<<"  MOLINFO "<<getLabel()<<": starting python interpreter\n";
      if(comm.Get_rank()==0) {
        selector=Tools::make_unique<Subprocess>(pythonCmd+" \""+config::getPlumedRoot()+"\"/scripts/selector.sh --pdb " + reference);
        selector->stop();
      }
      selector_running=true;
    }

    atoms.resize(0);

    if(comm.Get_rank()==0) {
      int ok=0;
      std::string error_msg;
      // this is a complicated way to have the exception propagated with MPI.
      // It is necessary since only one process calls the selector.
      // Probably I should recycle the exception propagation at library boundaries
      // to allow transferring the exception to other processes.
      try {
        plumed_assert(selector) << "Python interpreter is disabled, selection " + symbol + " cannot be interpreted";
        auto h=selector->contStop(); // stops again when it goes out of scope
        (*selector) << symbol << "\n";
        selector->flush();
        std::string res;
        std::vector<std::string> words;
        while(true) {
          selector->getline(res);
          words=Tools::getWords(res);
          if(!words.empty() && words[0]=="Error") {
            plumed_error()<<res;
          }
          if(!words.empty() && words[0]=="Selection:") {
            break;
          }
        }
        words.erase(words.begin());
        for(const auto & w : words) {
          int n;
          if(w.empty()) {
            continue;
          }
          Tools::convert(w,n);
          atoms.push_back(AtomNumber::serial(n));
        }
        ok=1;
      } catch (const Exception & e) {
        error_msg=e.what();
      }
      comm.Bcast(ok,0);
      if(!ok) {
        size_t s=error_msg.length();
        comm.Bcast(s,0);
        comm.Bcast(error_msg,0);
        throw Exception()<<error_msg;
      }
      size_t nat=atoms.size();
      comm.Bcast(nat,0);
      comm.Bcast(atoms,0);
    } else {
      int ok=0;
      std::string error_msg;
      comm.Bcast(ok,0);
      if(!ok) {
        size_t s;
        comm.Bcast(s,0);
        error_msg.resize(s);
        comm.Bcast(error_msg,0);
        throw Exception()<<error_msg;
      }
      size_t nat=0;
      comm.Bcast(nat,0);
      atoms.resize(nat);
      comm.Bcast(atoms,0);
    }
    log<<"  selection interpreted using ";
    if(Tools::startWith(symbol,"mdt:")) {
      log<<"mdtraj "<<cite("McGibbon et al, Biophys. J., 109, 1528 (2015)")<<"\n";
    }
    if(Tools::startWith(symbol,"mda:")) {
      log<<"MDAnalysis "<<cite("Gowers et al, Proceedings of the 15th Python in Science Conference, doi:10.25080/majora-629e541a-00e (2016)")<<"\n";
    }
    if(Tools::startWith(symbol,"vmdexec:")) {
      log<<"VMD "<<cite("Humphrey, Dalke, and Schulten, K., J. Molec. Graphics, 14, 33 (1996)")<<"\n";
    }
    if(Tools::startWith(symbol,"vmd:")) {
      log<<"VMD (github.com/Eigenstate/vmd-python) "<<cite("Humphrey, Dalke, and Schulten, K., J. Molec. Graphics, 14, 33 (1996)")<<"\n";
    }
    return;
  }
  MolDataClass::specialSymbol( mytype, symbol, pdb, atoms );
  if(atoms.empty()) {
    error(symbol + " not found in your MOLINFO structure");
  }
}

std::string GenericMolInfo::getAtomName(AtomNumber a)const {
  return pdb.getAtomName(a);
}

bool GenericMolInfo::checkForAtom(AtomNumber a)const {
  return pdb.checkForAtom(a);
}

int GenericMolInfo::getResidueNumber(AtomNumber a)const {
  return pdb.getResidueNumber(a);
}

unsigned GenericMolInfo::getPDBsize()const {
  return pdb.size();
}

std::string GenericMolInfo::getResidueName(AtomNumber a)const {
  return pdb.getResidueName(a);
}

std::string GenericMolInfo::getChainID(AtomNumber a)const {
  return pdb.getChainID(a);
}

Vector GenericMolInfo::getPosition(AtomNumber a)const {
  return pdb.getPosition(a);
}

bool GenericMolInfo::isWhole() const {
  return iswhole_;
}

void GenericMolInfo::prepare() {
  if(selector_running) {
    log<<"  MOLINFO "<<getLabel()<<": killing python interpreter\n";
    selector.reset();
    selector_running=false;
  }
}

}
