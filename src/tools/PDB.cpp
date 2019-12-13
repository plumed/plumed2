/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2019 The plumed team
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
#include "PDB.h"
#include "Tools.h"
#include "Log.h"
#include "h36.h"
#include <cstdio>
#include <iostream>
#include "core/SetupMolInfo.h"
#include "Tensor.h"

using namespace std;

//+PLUMEDOC INTERNAL pdbreader
/*
PLUMED can use the PDB format in several places

- To read molecular structure (\ref MOLINFO).
- To read reference conformations (\ref RMSD, but also many other methods in \ref dists, \ref FIT_TO_TEMPLATE, etc).

The implemented PDB reader expects a file formatted correctly according to the
[PDB standard](http://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html).
In particular, the following columns are read from ATOM records
\verbatim
columns | content
1-6     | record name (ATOM or HETATM)
7-11    | serial number of the atom (starting from 1)
13-16   | atom name
18-20   | residue name
22      | chain id
23-26   | residue number
31-38   | x coordinate
39-46   | y coordinate
47-54   | z coordinate
55-60   | occupancy
61-66   | beta factor
\endverbatim
PLUMED parser is slightly more permissive than the official PDB format
in the fact that the format of real numbers is not fixed. In other words,
any real number that can be parsed is OK and the dot can be placed anywhere. However,
__columns are interpret strictly__. A sample PDB should look like the following
\verbatim
ATOM      2  CH3 ACE     1      12.932 -14.718  -6.016  1.00  1.00
ATOM      5  C   ACE     1      21.312  -9.928  -5.946  1.00  1.00
ATOM      9  CA  ALA     2      19.462 -11.088  -8.986  1.00  1.00
\endverbatim

Notice that serial numbers need not to be consecutive. In the three-line example above,
only the coordinates of three atoms are provided. This is perfectly legal and indicates PLUMED
that information about these atoms only is available. This could be both for structural
information in \ref MOLINFO, where the other atoms would have no name assigned, and for
reference structures used in \ref RMSD, where only the provided atoms would be used to compute RMSD.

\par Occupancy and beta factors

PLUMED reads also occupancy and beta factors that however are given a very special meaning.
In cases where the PDB structure is used as a reference for an alignment (that's the case
for instance in \ref RMSD and in \ref FIT_TO_TEMPLATE), the occupancy column is used
to provide the weight of each atom in the alignment. In cases where, perhaps after alignment,
the displacement between running coordinates and the provided PDB is computed, the beta factors
are used as weight for the displacement.
Since setting the weights to zero is the same as __not__ including an atom in the alignment or
displacement calculation, the two following reference files would be equivalent when used in an \ref RMSD
calculation. First file:
\verbatim
ATOM      2  CH3 ACE     1      12.932 -14.718  -6.016  1.00  1.00
ATOM      5  C   ACE     1      21.312  -9.928  -5.946  1.00  1.00
ATOM      9  CA  ALA     2      19.462 -11.088  -8.986  0.00  0.00
\endverbatim
Second file:
\verbatim
ATOM      2  CH3 ACE     1      12.932 -14.718  -6.016  1.00  1.00
ATOM      5  C   ACE     1      21.312  -9.928  -5.946  1.00  1.00
\endverbatim
However notice that many extra atoms with zero weight might slow down the calculation, so
removing lines is better than setting their weights to zero.
In addition, weights for alignment need not to be equivalent to weights for displacement.

\par Systems with more than 100k atoms

Notice that it very likely does not make any sense to compute the \ref RMSD or any other structural
deviation __using__ so many atoms. However, if the protein for which you want to compute \ref RMSD
has atoms with large serial numbers (e.g. because it is located __after__ solvent in the sorted list of atoms)
you might end up with troubles with the limitations of the PDB format. Indeed, since there are 5
columns available for atom serial number, this number cannot be larger than 99999.
In addition, providing \ref MOLINFO with names associated to atoms with a serial larger than 99999 would be impossible.

Since PLUMED 2.4 we allow [hybrid 36](http://cci.lbl.gov/hybrid_36/) format to be used to specify atom numbers.
This format is not particularly widespread, but has the nice feature that it provides a one-to-one mapping
between numbers up to approximately 80 millions and strings with 5 characters, plus it is backward compatible
for numbers smaller than 100000. This is not true for notations like the hex notation exported by VMD.
Using the hybrid 36 format, the ATOM records for atom ranging from 99997 to 100002 would read like these:
\verbatim
ATOM  99997  Ar      X   1      45.349  38.631  15.116  1.00  1.00
ATOM  99998  Ar      X   1      46.189  38.631  15.956  1.00  1.00
ATOM  99999  Ar      X   1      46.189  39.471  15.116  1.00  1.00
ATOM  A0000  Ar      X   1      45.349  39.471  15.956  1.00  1.00
ATOM  A0000  Ar      X   1      45.349  38.631  16.796  1.00  1.00
ATOM  A0001  Ar      X   1      46.189  38.631  17.636  1.00  1.00
\endverbatim
There are tools that can be found to translate from integers to strings and back using hybrid 36 format
(a simple python script can be found [here](https://sourceforge.net/p/cctbx/code/HEAD/tree/trunk/iotbx/pdb/hybrid_36.py)).
In addition, as of PLUMED 2.5, we provide a \ref pdbrenumber "command line tool" that can be used to renumber atoms in a PDB file.

*/
//+ENDPLUMEDOC


namespace PLMD {

void PDB::setAtomNumbers( const std::vector<AtomNumber>& atoms ) {
  positions.resize( atoms.size() ); occupancy.resize( atoms.size() );
  beta.resize( atoms.size() ); numbers.resize( atoms.size() );
  for(unsigned i=0; i<atoms.size(); ++i) { numbers[i]=atoms[i]; beta[i]=1.0; occupancy[i]=1.0; }
}

void PDB::setArgumentNames( const std::vector<std::string>& argument_names ) {
  argnames.resize( argument_names.size() );
  for(unsigned i=0; i<argument_names.size(); ++i) {
    argnames[i]=argument_names[i];
    arg_data.insert( std::pair<std::string,double>( argnames[i], 0.0 ) );
  }
}

bool PDB::getArgumentValue( const std::string& name, double& value ) const {
  std::map<std::string,double>::const_iterator it = arg_data.find(name);
  if( it!=arg_data.end() ) { value = it->second; return true; }
  return false;
}

void PDB::setAtomPositions( const std::vector<Vector>& pos ) {
  plumed_assert( pos.size()==positions.size() );
  for(unsigned i=0; i<positions.size(); ++i) positions[i]=pos[i];
}

void PDB::setArgumentValue( const std::string& argname, const double& val ) {
  // First set the value of the value of the argument in the map
  arg_data.find(argname)->second = val;
}

// bool PDB::hasRequiredProperties( const std::vector<std::string>& inproperties ){
//   bool hasprop=false;
//   for(unsigned i=0;i<remark.size();++i){
//       if( remark[i].find("PROPERTIES=")!=std::string::npos){ hasprop=true; break; }
//   }
//   if( !hasprop ){
//       std::string mypropstr="PROPERTIES=" + inproperties[0];
//       for(unsigned i=1;i<inproperties.size();++i) mypropstr += "," + inproperties[i];
//       remark.push_back( mypropstr );
//   }
//   // Now check that all required properties are there
//   for(unsigned i=0;i<inproperties.size();++i){
//       hasprop=false;
//       for(unsigned j=0;j<remark.size();++j){
//           if( remark[j].find(inproperties[i]+"=")!=std::string::npos){ hasprop=true; break; }
//       }
//       if( !hasprop ) return false;
//   }
//   return true;
// }

void PDB::addBlockEnd( const unsigned& end ) {
  block_ends.push_back( end );
}

unsigned PDB::getNumberOfAtomBlocks()const {
  return block_ends.size();
}

const std::vector<unsigned> & PDB::getAtomBlockEnds()const {
  return block_ends;
}

const std::vector<Vector> & PDB::getPositions()const {
  return positions;
}

void PDB::setPositions(const std::vector<Vector> &v ) {
  plumed_assert( v.size()==positions.size() );
  positions=v;
}

const std::vector<double> & PDB::getOccupancy()const {
  return occupancy;
}

const std::vector<double> & PDB::getBeta()const {
  return beta;
}

void PDB::addRemark( std::vector<std::string>& v1 ) {
  Tools::parse(v1,"TYPE",mtype);
  Tools::parseVector(v1,"ARG",argnames);
  for(unsigned i=0; i<v1.size(); ++i) {
    if( v1[i].find("=")!=std::string::npos ) {
      std::size_t eq=v1[i].find_first_of('=');
      std::string name=v1[i].substr(0,eq);
      std::string sval=v1[i].substr(eq+1);
      double val; Tools::convert( sval, val );
      arg_data.insert( std::pair<std::string,double>( name, val ) );
    } else {
      flags.push_back(v1[i]);
    }
  }
}

bool PDB::hasFlag( const std::string& fname ) const {
  for(unsigned i=0; i<flags.size(); ++i) {
    if( flags[i]==fname ) return true;
  }
  return false;
}


const std::vector<AtomNumber> & PDB::getAtomNumbers()const {
  return numbers;
}

const Vector & PDB::getBoxAxs()const {
  return BoxXYZ;
}

const Vector & PDB::getBoxAng()const {
  return BoxABG;
}

const Tensor & PDB::getBoxVec()const {
  return Box;
}

std::string PDB::getAtomName(AtomNumber a)const {
  const auto p=number2index.find(a);
  if(p==number2index.end()) {
    std::string num; Tools::convert( a.serial(), num );
    plumed_merror("Name of atom " + num + " not found" );
    return "";
  } else return atomsymb[p->second];
}

unsigned PDB::getResidueNumber(AtomNumber a)const {
  const auto p=number2index.find(a);
  if(p==number2index.end()) {
    std::string num; Tools::convert( a.serial(), num );
    plumed_merror("Residue for atom " + num + " not found" );
    return 0;
  } else return residue[p->second];
}

std::string PDB::getResidueName(AtomNumber a) const {
  const auto p=number2index.find(a);
  if(p==number2index.end()) {
    std::string num; Tools::convert( a.serial(), num );
    plumed_merror("Residue for atom " + num + " not found" );
    return "";
  } else return residuenames[p->second];
}

unsigned PDB::size()const {
  return positions.size();
}

bool PDB::readFromFilepointer(FILE *fp,bool naturalUnits,double scale) {
  //cerr<<file<<endl;
  bool file_is_alive=false;
  if(naturalUnits) scale=1.0;
  string line;
  fpos_t pos; bool between_ters=true;
  while(Tools::getline(fp,line)) {
    //cerr<<line<<"\n";
    fgetpos (fp,&pos);
    while(line.length()<80) line.push_back(' ');
    string record=line.substr(0,6);
    string serial=line.substr(6,5);
    string atomname=line.substr(12,4);
    string residuename=line.substr(17,3);
    string chainID=line.substr(21,1);
    string resnum=line.substr(22,4);
    string x=line.substr(30,8);
    string y=line.substr(38,8);
    string z=line.substr(46,8);
    string occ=line.substr(54,6);
    string bet=line.substr(60,6);
    string BoxX=line.substr(6,9);
    string BoxY=line.substr(15,9);
    string BoxZ=line.substr(24,9);
    string BoxA=line.substr(33,7);
    string BoxB=line.substr(40,7);
    string BoxG=line.substr(47,7);
    Tools::trim(record);
    if(record=="TER") { between_ters=false; block_ends.push_back( positions.size() ); }
    if(record=="END") { file_is_alive=true;  break;}
    if(record=="ENDMDL") { file_is_alive=true;  break;}
    if(record=="REMARK") {
      vector<string> v1;  v1=Tools::getWords(line.substr(6));
      addRemark( v1 );
    }
    if(record=="CRYST1") {
      Tools::convert(BoxX,BoxXYZ[0]);
      Tools::convert(BoxY,BoxXYZ[1]);
      Tools::convert(BoxZ,BoxXYZ[2]);
      Tools::convert(BoxA,BoxABG[0]);
      Tools::convert(BoxB,BoxABG[1]);
      Tools::convert(BoxG,BoxABG[2]);
      BoxXYZ*=scale;
      double cosA=cos(BoxABG[0]*pi/180.);
      double cosB=cos(BoxABG[1]*pi/180.);
      double cosG=cos(BoxABG[2]*pi/180.);
      double sinG=sin(BoxABG[2]*pi/180.);
      for (unsigned i=0; i<3; i++) {Box[i][0]=0.; Box[i][1]=0.; Box[i][2]=0.;}
      Box[0][0]=BoxXYZ[0];
      Box[1][0]=BoxXYZ[1]*cosG;
      Box[1][1]=BoxXYZ[1]*sinG;
      Box[2][0]=BoxXYZ[2]*cosB;
      Box[2][1]=(BoxXYZ[2]*BoxXYZ[1]*cosA-Box[2][0]*Box[1][0])/Box[1][1];
      Box[2][2]=sqrt(BoxXYZ[2]*BoxXYZ[2]-Box[2][0]*Box[2][0]-Box[2][1]*Box[2][1]);
    }
    if(record=="ATOM" || record=="HETATM") {
      between_ters=true;
      AtomNumber a; unsigned resno;
      double o,b;
      Vector p;
      {
        int result;
        auto trimmed=serial;
        Tools::trim(trimmed);
        const char* errmsg = h36::hy36decode(5, trimmed.c_str(),trimmed.length(), &result);
        if(errmsg) {
          std::string msg(errmsg);
          plumed_merror(msg);
        }
        a.setSerial(result);
      }

      Tools::convert(resnum,resno);
      Tools::convert(occ,o);
      Tools::convert(bet,b);
      Tools::convert(x,p[0]);
      Tools::convert(y,p[1]);
      Tools::convert(z,p[2]);
      // scale into nm
      p*=scale;
      numbers.push_back(a);
      number2index[a]=positions.size();
      std::size_t startpos=atomname.find_first_not_of(" \t");
      std::size_t endpos=atomname.find_last_not_of(" \t");
      atomsymb.push_back( atomname.substr(startpos, endpos-startpos+1) );
      residue.push_back(resno);
      chain.push_back(chainID);
      occupancy.push_back(o);
      beta.push_back(b);
      positions.push_back(p);
      residuenames.push_back(residuename);
    }
  }
  if( between_ters ) block_ends.push_back( positions.size() );
  return file_is_alive;
}

bool PDB::read(const std::string&file,bool naturalUnits,double scale) {
  FILE* fp=fopen(file.c_str(),"r");
  if(!fp) return false;
  readFromFilepointer(fp,naturalUnits,scale);
  fclose(fp);
  return true;
}

void PDB::getChainNames( std::vector<std::string>& chains ) const {
  chains.resize(0);
  chains.push_back( chain[0] );
  for(unsigned i=1; i<size(); ++i) {
    if( chains[chains.size()-1]!=chain[i] ) chains.push_back( chain[i] );
  }
}

void PDB::getResidueRange( const std::string& chainname, unsigned& res_start, unsigned& res_end, std::string& errmsg ) const {
  bool inres=false, foundchain=false;
  for(unsigned i=0; i<size(); ++i) {
    if( chain[i]==chainname ) {
      if(!inres) {
        if(foundchain) errmsg="found second start of chain named " + chainname;
        res_start=residue[i];
      }
      inres=true; foundchain=true;
    } else if( inres && chain[i]!=chainname ) {
      inres=false;
      res_end=residue[i-1];
    }
  }
  if(inres) res_end=residue[size()-1];
}

void PDB::getAtomRange( const std::string& chainname, AtomNumber& a_start, AtomNumber& a_end, std::string& errmsg ) const {
  bool inres=false, foundchain=false;
  for(unsigned i=0; i<size(); ++i) {
    if( chain[i]==chainname ) {
      if(!inres) {
        if(foundchain) errmsg="found second start of chain named " + chainname;
        a_start=numbers[i];
      }
      inres=true; foundchain=true;
    } else if( inres && chain[i]!=chainname ) {
      inres=false;
      a_end=numbers[i-1];
    }
  }
  if(inres) a_end=numbers[size()-1];
}

std::string PDB::getResidueName( const unsigned& resnum ) const {
  for(unsigned i=0; i<size(); ++i) {
    if( residue[i]==resnum ) return residuenames[i];
  }
  std::string num; Tools::convert( resnum, num );
  plumed_merror("residue " + num + " not found" );
  return "";
}

std::string PDB::getResidueName(const unsigned& resnum,const std::string& chainid ) const {
  for(unsigned i=0; i<size(); ++i) {
    if( residue[i]==resnum && ( chainid=="*" || chain[i]==chainid) ) return residuenames[i];
  }
  std::string num; Tools::convert( resnum, num );
  plumed_merror("residue " + num + " not found in chain " + chainid );
  return "";
}


AtomNumber PDB::getNamedAtomFromResidue( const std::string& aname, const unsigned& resnum ) const {
  for(unsigned i=0; i<size(); ++i) {
    if( residue[i]==resnum && atomsymb[i]==aname ) return numbers[i];
  }
  std::string num; Tools::convert( resnum, num );
  plumed_merror("residue " + num + " does not contain an atom named " + aname );
  return numbers[0]; // This is to stop compiler errors
}

AtomNumber PDB::getNamedAtomFromResidueAndChain( const std::string& aname, const unsigned& resnum, const std::string& chainid ) const {
  for(unsigned i=0; i<size(); ++i) {
    if( residue[i]==resnum && atomsymb[i]==aname && ( chainid=="*" || chain[i]==chainid) ) return numbers[i];
  }
  std::string num; Tools::convert( resnum, num );
  plumed_merror("residue " + num + " from chain " + chainid + " does not contain an atom named " + aname );
  return numbers[0]; // This is to stop compiler errors
}

std::vector<AtomNumber> PDB::getAtomsInResidue(const unsigned& resnum,const std::string& chainid)const {
  std::vector<AtomNumber> tmp;
  for(unsigned i=0; i<size(); ++i) {
    if( residue[i]==resnum && ( chainid=="*" || chain[i]==chainid) ) tmp.push_back(numbers[i]);
  }
  if(tmp.size()==0) {
    std::string num; Tools::convert( resnum, num );
    plumed_merror("Cannot find residue " + num + " from chain " + chainid  );
  }
  return tmp;
}

std::vector<AtomNumber> PDB::getAtomsInChain(const std::string& chainid)const {
  std::vector<AtomNumber> tmp;
  for(unsigned i=0; i<size(); ++i) {
    if( chainid=="*" || chain[i]==chainid ) tmp.push_back(numbers[i]);
  }
  if(tmp.size()==0) {
    plumed_merror("Cannot find atoms from chain " + chainid  );
  }
  return tmp;
}

std::string PDB::getChainID(const unsigned& resnumber) const {
  for(unsigned i=0; i<size(); ++i) {
    if(resnumber==residue[i]) return chain[i];
  }
  plumed_merror("Not enough residues in pdb input file");
}

bool PDB::checkForResidue( const std::string& name ) const {
  for(unsigned i=0; i<size(); ++i) {
    if( residuenames[i]==name ) return true;
  }
  return false;
}

bool PDB::checkForAtom( const std::string& name ) const {
  for(unsigned i=0; i<size(); ++i) {
    if( atomsymb[i]==name ) return true;
  }
  return false;
}

Log& operator<<(Log& ostr, const PDB&  pdb) {
  char buffer[1000];
  for(unsigned i=0; i<pdb.positions.size(); i++) {
    sprintf(buffer,"ATOM %3d %8.3f %8.3f %8.3f\n",pdb.numbers[i].serial(),pdb.positions[i][0],pdb.positions[i][1],pdb.positions[i][2]);
    ostr<<buffer;
  }
  return ostr;
}

Vector PDB::getPosition(AtomNumber a)const {
  const auto p=number2index.find(a);
  if(p==number2index.end()) plumed_merror("atom not available");
  else return positions[p->second];
}

std::vector<std::string> PDB::getArgumentNames()const {
  return argnames;
}

std::string PDB::getMtype() const {
  return mtype;
}

void PDB::print( const double& lunits, SetupMolInfo* mymoldat, OFile& ofile, const std::string& fmt ) {
  if( argnames.size()>0 ) {
    ofile.printf("REMARK ARG=%s", argnames[0].c_str() );
    for(unsigned i=1; i<argnames.size(); ++i) ofile.printf(",%s",argnames[i].c_str() );
    ofile.printf("\n"); ofile.printf("REMARK ");
  }
  std::string descr2;
  if(fmt.find("-")!=std::string::npos) {
    descr2="%s=" + fmt + " ";
  } else {
    // This ensures numbers are left justified (i.e. next to the equals sign
    std::size_t psign=fmt.find("%");
    plumed_assert( psign!=std::string::npos );
    descr2="%s=%-" + fmt.substr(psign+1) + " ";
  }
  for(std::map<std::string,double>::iterator it=arg_data.begin(); it!=arg_data.end(); ++it) ofile.printf( descr2.c_str(),it->first.c_str(), it->second );
  if( argnames.size()>0 ) ofile.printf("\n");
  if( !mymoldat ) {
    for(unsigned i=0; i<positions.size(); ++i) {
      std::array<char,6> at;
      const char* msg = h36::hy36encode(5,numbers[i].serial(),&at[0]);
      plumed_assert(msg==nullptr) << msg;
      at[5]=0;
      ofile.printf("ATOM  %s  X   RES  %4u    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                   &at[0], i,
                   lunits*positions[i][0], lunits*positions[i][1], lunits*positions[i][2],
                   occupancy[i], beta[i] );
    }
  } else {
    for(unsigned i=0; i<positions.size(); ++i) {
      std::array<char,6> at;
      const char* msg = h36::hy36encode(5,numbers[i].serial(),&at[0]);
      plumed_assert(msg==nullptr) << msg;
      at[5]=0;
      ofile.printf("ATOM  %s %-4s %3s  %4u    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                   &at[0], mymoldat->getAtomName(numbers[i]).c_str(),
                   mymoldat->getResidueName(numbers[i]).c_str(), mymoldat->getResidueNumber(numbers[i]),
                   lunits*positions[i][0], lunits*positions[i][1], lunits*positions[i][2],
                   occupancy[i], beta[i] );
    }
  }
  ofile.printf("END\n");
}


}

