/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2016 The plumed team
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
#include <cstdio>
#include <iostream>
#include "core/SetupMolInfo.h"

using namespace std;

namespace PLMD{

void PDB::setAtomNumbers( const std::vector<AtomNumber>& atoms ){
  positions.resize( atoms.size() ); occupancy.resize( atoms.size() ); 
  beta.resize( atoms.size() ); numbers.resize( atoms.size() );
  for(unsigned i=0;i<atoms.size();++i){ numbers[i]=atoms[i]; beta[i]=1.0; occupancy[i]=1.0; }
}

void PDB::setArgumentNames( const std::vector<std::string>& argument_names ){
  argnames.resize( argument_names.size() );
  for(unsigned i=0;i<argument_names.size();++i){
      argnames[i]=argument_names[i]; 
      arg_data.insert( std::pair<std::string,double>( argnames[i], 0.0 ) ); 
  }
}

bool PDB::getArgumentValue( const std::string& name, double& value ) const {
  std::map<std::string,double>::const_iterator it = arg_data.find(name);
  if( it!=arg_data.end() ){ value = it->second; return true; }
  return false;
}

void PDB::setAtomPositions( const std::vector<Vector>& pos ){
  plumed_assert( pos.size()==positions.size() );
  for(unsigned i=0;i<positions.size();++i) positions[i]=pos[i];
}

void PDB::setArgumentValue( const std::string& argname, const double& val ){
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

void PDB::addBlockEnd( const unsigned& end ){
  block_ends.push_back( end );
}

unsigned PDB::getNumberOfAtomBlocks()const{
  return block_ends.size();
}

const std::vector<unsigned> & PDB::getAtomBlockEnds()const{
  return block_ends;
}

const std::vector<Vector> & PDB::getPositions()const{
  return positions;
}

void PDB::setPositions(const std::vector<Vector> &v ){
	  plumed_assert( v.size()==positions.size() );	
	  positions=v;
}

const std::vector<double> & PDB::getOccupancy()const{
  return occupancy;
}

const std::vector<double> & PDB::getBeta()const{
  return beta;
}

void PDB::addRemark( std::vector<std::string>& v1 ){
  Tools::parse(v1,"TYPE",mtype); 
  Tools::parseVector(v1,"ARG",argnames);
  for(unsigned i=0;i<v1.size();++i){
      if( v1[i].find("=")!=std::string::npos ){
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
  for(unsigned i=0;i<flags.size();++i){
      if( flags[i]==fname ) return true;
  }
  return false;
}


const std::vector<AtomNumber> & PDB::getAtomNumbers()const{
  return numbers;
}

std::string PDB::getAtomName(AtomNumber a)const{
  std::map<AtomNumber,unsigned>::const_iterator p;
  p=number2index.find(a);
  if(p==number2index.end()) return "";
  else return atomsymb[p->second];
}

unsigned PDB::getResidueNumber(AtomNumber a)const{
  std::map<AtomNumber,unsigned>::const_iterator p;
  p=number2index.find(a);
  if(p==number2index.end()) return 0;
  else return residue[p->second];
}

std::string PDB::getResidueName(AtomNumber a) const{
  std::map<AtomNumber,unsigned>::const_iterator p;
  p=number2index.find(a);
  if(p==number2index.end()) return "";
  else return residuenames[p->second];
}

unsigned PDB::size()const{
  return positions.size();
}

bool PDB::readFromFilepointer(FILE *fp,bool naturalUnits,double scale){
  //cerr<<file<<endl;
  bool file_is_alive=false;
  if(naturalUnits) scale=1.0;
  string line;
  fpos_t pos; bool between_ters=true;
  while(Tools::getline(fp,line)){
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
    Tools::trim(record);
    if(record=="TER"){ between_ters=false; block_ends.push_back( positions.size() ); }
    if(record=="END"){ file_is_alive=true;  break;}
    if(record=="ENDMDL"){ file_is_alive=true;  break;}
    if(record=="REMARK"){
         vector<string> v1;  v1=Tools::getWords(line.substr(6));  
         addRemark( v1 );
    }
    if(record=="ATOM" || record=="HETATM"){
      between_ters=true;
      AtomNumber a; unsigned resno;
      double o,b;
      Vector p;
      Tools::convert(serial,a);
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

bool PDB::read(const std::string&file,bool naturalUnits,double scale){
  FILE* fp=fopen(file.c_str(),"r");
  if(!fp) return false;
  readFromFilepointer(fp,naturalUnits,scale);
  fclose(fp);
  return true;
}

void PDB::getChainNames( std::vector<std::string>& chains ) const {
  chains.resize(0);
  chains.push_back( chain[0] );
  for(unsigned i=1;i<size();++i){
     if( chains[chains.size()-1]!=chain[i] ) chains.push_back( chain[i] );
  }
} 

void PDB::getResidueRange( const std::string& chainname, unsigned& res_start, unsigned& res_end, std::string& errmsg ) const {
  bool inres=false, foundchain=false;
  for(unsigned i=0;i<size();++i){
     if( chain[i]==chainname ){
         if(!inres){
           if(foundchain) errmsg="found second start of chain named " + chainname;
           res_start=residue[i];
         }
         inres=true; foundchain=true;
     } else if( inres && chain[i]!=chainname ){
         inres=false;
         res_end=residue[i-1];
     }
  }
  if(inres) res_end=residue[size()-1];
}

void PDB::getAtomRange( const std::string& chainname, AtomNumber& a_start, AtomNumber& a_end, std::string& errmsg ) const {
  bool inres=false, foundchain=false;
  for(unsigned i=0;i<size();++i){
     if( chain[i]==chainname ){
         if(!inres){
           if(foundchain) errmsg="found second start of chain named " + chainname;
           a_start=numbers[i];
         } 
         inres=true; foundchain=true;
     } else if( inres && chain[i]!=chainname ){
         inres=false;
         a_end=numbers[i-1];
     }
  }       
  if(inres) a_end=numbers[size()-1];
} 

std::string PDB::getResidueName( const unsigned& resnum ) const {
  for(unsigned i=0;i<size();++i){
     if( residue[i]==resnum ) return residuenames[i];
  }
  return "";
}

std::string PDB::getResidueName(const unsigned& resnum,const std::string& chainid ) const {
  for(unsigned i=0;i<size();++i){
     if( residue[i]==resnum && ( chainid=="*" || chain[i]==chainid) ) return residuenames[i];
  }
  return "";
}


AtomNumber PDB::getNamedAtomFromResidue( const std::string& aname, const unsigned& resnum ) const {
  for(unsigned i=0;i<size();++i){
     if( residue[i]==resnum && atomsymb[i]==aname ) return numbers[i];
  }
  std::string num; Tools::convert( resnum, num );
  plumed_merror("residue " + num + " does not contain an atom named " + aname );
  return numbers[0]; // This is to stop compiler errors
}

AtomNumber PDB::getNamedAtomFromResidueAndChain( const std::string& aname, const unsigned& resnum, const std::string& chainid ) const{
  for(unsigned i=0;i<size();++i){
     if( residue[i]==resnum && atomsymb[i]==aname && ( chainid=="*" || chain[i]==chainid) ) return numbers[i];
  }
  std::string num; Tools::convert( resnum, num );
  plumed_merror("residue " + num + " from chain " + chainid + " does not contain an atom named " + aname );
  return numbers[0]; // This is to stop compiler errors
}

std::vector<AtomNumber> PDB::getAtomsInResidue(const unsigned& resnum,const std::string& chainid)const {
  std::vector<AtomNumber> tmp;
  for(unsigned i=0;i<size();++i){
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
  for(unsigned i=0;i<size();++i){
     if( chainid=="*" || chain[i]==chainid ) tmp.push_back(numbers[i]);
  }
  if(tmp.size()==0) {
    plumed_merror("Cannot find atoms from chain " + chainid  );
  }
  return tmp;
}

std::string PDB::getChainID(const unsigned& resnumber) const {
  for(unsigned i=0;i<size();++i){
     if(resnumber==residue[i]) return chain[i];
  }
  plumed_merror("Not enough residues in pdb input file");
}

bool PDB::checkForResidue( const std::string& name ) const {
  for(unsigned i=0;i<size();++i){
      if( residuenames[i]==name ) return true;
  }
  return false;
}

bool PDB::checkForAtom( const std::string& name ) const {
  for(unsigned i=0;i<size();++i){
     if( atomsymb[i]==name ) return true;
  }
  return false;
}

Log& operator<<(Log& ostr, const PDB&  pdb){
   char buffer[1000];
   for(unsigned i=0;i<pdb.positions.size();i++){ 
      sprintf(buffer,"ATOM %3d %8.3f %8.3f %8.3f\n",pdb.numbers[i].serial(),pdb.positions[i][0],pdb.positions[i][1],pdb.positions[i][2]);
      ostr<<buffer;
   }
   return ostr;
}

Vector PDB::getPosition(AtomNumber a)const{
     std::map<AtomNumber,unsigned>::const_iterator p;
     p=number2index.find(a);
     if(p==number2index.end()) plumed_merror("atom not available");
     else return positions[p->second];
}

std::vector<std::string> PDB::getArgumentNames()const {
    return argnames;
}

std::string PDB::getMtype() const {
  return mtype;
}

void PDB::print( const double& lunits, SetupMolInfo* mymoldat, OFile& ofile, const std::string& fmt ){
  if( argnames.size()>0 ){
      ofile.printf("REMARK ARG=%s", argnames[0].c_str() ); 
      for(unsigned i=1;i<argnames.size();++i) ofile.printf(",%s",argnames[i].c_str() ); 
      ofile.printf("\n"); ofile.printf("REMARK ");
  }
  std::string descr2;
  if(fmt.find("-")!=std::string::npos){
      descr2="%s=" + fmt + " ";
  } else {
     // This ensures numbers are left justified (i.e. next to the equals sign
     std::size_t psign=fmt.find("%");
     plumed_assert( psign!=std::string::npos );
     descr2="%s=%-" + fmt.substr(psign+1) + " ";
  }
  for(std::map<std::string,double>::iterator it=arg_data.begin(); it!=arg_data.end();it++) ofile.printf( descr2.c_str(),it->first.c_str(), it->second ); 
  if( argnames.size()>0 ) ofile.printf("\n");
  if( !mymoldat ){
      for(unsigned i=0;i<positions.size();++i){
          ofile.printf("ATOM  %4d  X    RES  %4u  %8.3f%8.3f%8.3f%6.2f%6.2f\n",
            numbers[i].serial(), i,
            lunits*positions[i][0], lunits*positions[i][1], lunits*positions[i][2],
            occupancy[i], beta[i] );
      }
  } else {
      for(unsigned i=0;i<positions.size();++i){
          ofile.printf("ATOM  %5d %-4s %3s  %4u    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
            numbers[i].serial(), mymoldat->getAtomName(numbers[i]).c_str(),
            mymoldat->getResidueName(numbers[i]).c_str(), mymoldat->getResidueNumber(numbers[i]),
            lunits*positions[i][0], lunits*positions[i][1], lunits*positions[i][2],
            occupancy[i], beta[i] );
      }
  }
  ofile.printf("END\n");
}


}

