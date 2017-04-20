/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2017 The plumed team
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
#include <cstdio>
#include <iostream>

using namespace std;

namespace PLMD{

/// Tiny namespace for hybrid36 format.
/// This namespace includes freely available tools for h36 format.
/// I place them here for usage within PDB class. In case we need them
/// in other places they might be better encapsulated in a c++ class
/// and placed in a separate file.
namespace h36{


/*! C port of the hy36encode() and hy36decode() functions in the
    hybrid_36.py Python prototype/reference implementation.
    See the Python script for more information.

    This file has no external dependencies, NOT even standard C headers.
    Optionally, use hybrid_36_c.h, or simply copy the declarations
    into your code.

    This file is unrestricted Open Source (cctbx.sf.net).
    Please send corrections and enhancements to cctbx@cci.lbl.gov .

    See also: http://cci.lbl.gov/hybrid_36/

    Ralf W. Grosse-Kunstleve, Feb 2007.
 */

/* The following #include may be commented out.
   It is here only to enforce consistency of the declarations
   and the definitions.
 */
// #include <iotbx/pdb/hybrid_36_c.h>

/* All static functions below are implementation details
   (and not accessible from other translation units).
 */

static
const char*
digits_upper() { return "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"; }

static
const char*
digits_lower() { return "0123456789abcdefghijklmnopqrstuvwxyz"; }

static
const char*
value_out_of_range() { return "value out of range."; }

static
const char* invalid_number_literal() { return "invalid number literal."; }

static
const char* unsupported_width() { return "unsupported width."; }

static
void
fill_with_stars(unsigned width, char* result)
{
  while (width) {
    *result++ = '*';
    width--;
  }
  *result = '\0';
}

static
void
encode_pure(
  const char* digits,
  unsigned digits_size,
  unsigned width,
  int value,
  char* result)
{
  char buf[16];
  int rest;
  unsigned i, j;
  i = 0;
  j = 0;
  if (value < 0) {
    j = 1;
    value = -value;
  }
  while (1) {
    rest = value / digits_size;
    buf[i++] = digits[value - rest * digits_size];
    if (rest == 0) break;
    value = rest;
  }
  if (j) buf[i++] = '-';
  for(j=i;j<width;j++) *result++ = ' ';
  while (i != 0) *result++ = buf[--i];
  *result = '\0';
}

static
const char*
decode_pure(
  const int* digits_values,
  unsigned digits_size,
  const char* s,
  unsigned s_size,
  int* result)
{
  int si, dv;
  int have_minus = 0;
  int have_non_blank = 0;
  int value = 0;
  unsigned i = 0;
  for(;i<s_size;i++) {
    si = s[i];
    if (si < 0 || si > 127) {
      *result = 0;
      return invalid_number_literal();
    }
    if (si == ' ') {
      if (!have_non_blank) continue;
      value *= digits_size;
    }
    else if (si == '-') {
      if (have_non_blank) {
        *result = 0;
        return invalid_number_literal();
      }
      have_non_blank = 1;
      have_minus = 1;
      continue;
    }
    else {
      have_non_blank = 1;
      dv = digits_values[si];
      if (dv < 0 || dv >= digits_size) {
        *result = 0;
        return invalid_number_literal();
      }
      value *= digits_size;
      value += dv;
    }
  }
  if (have_minus) value = -value;
  *result = value;
  return 0;
}

/*! hybrid-36 encoder: converts integer value to string result

      width: must be 4 (e.g. for residue sequence numbers)
                  or 5 (e.g. for atom serial numbers)

      value: integer value to be converted

      result: pointer to char array of size width+1 or greater
              on return result is null-terminated

      return value: pointer to error message, if any,
                    or 0 on success

    Example usage (from C++):
      char result[4+1];
      const char* errmsg = hy36encode(4, 12345, result);
      if (errmsg) throw std::runtime_error(errmsg);
 */
const char*
hy36encode(unsigned width, int value, char* result)
{
  int i = value;
  if (width == 4U) {
    if (i >= -999) {
      if (i < 10000) {
        encode_pure(digits_upper(), 10U, 4U, i, result);
        return 0;
      }
      i -= 10000;
      if (i < 1213056 /* 26*36**3 */) {
        i += 466560 /* 10*36**3 */;
        encode_pure(digits_upper(), 36U, 0U, i, result);
        return 0;
      }
      i -= 1213056;
      if (i < 1213056) {
        i += 466560;
        encode_pure(digits_lower(), 36U, 0U, i, result);
        return 0;
      }
    }
  }
  else if (width == 5U) {
    if (i >= -9999) {
      if (i < 100000) {
        encode_pure(digits_upper(), 10U, 5U, i, result);
        return 0;
      }
      i -= 100000;
      if (i < 43670016 /* 26*36**4 */) {
        i += 16796160 /* 10*36**4 */;
        encode_pure(digits_upper(), 36U, 0U, i, result);
        return 0;
      }
      i -= 43670016;
      if (i < 43670016) {
        i += 16796160;
        encode_pure(digits_lower(), 36U, 0U, i, result);
        return 0;
      }
    }
  }
  else {
    fill_with_stars(width, result);
    return unsupported_width();
  }
  fill_with_stars(width, result);
  return value_out_of_range();
}

/*! hybrid-36 decoder: converts string s to integer result

      width: must be 4 (e.g. for residue sequence numbers)
                  or 5 (e.g. for atom serial numbers)

      s: string to be converted
         does not have to be null-terminated

      s_size: size of s
              must be equal to width, or an error message is
              returned otherwise

      result: integer holding the conversion result

      return value: pointer to error message, if any,
                    or 0 on success

    Example usage (from C++):
      int result;
      const char* errmsg = hy36decode(width, "A1T5", 4, &result);
      if (errmsg) throw std::runtime_error(errmsg);
 */
const char*
hy36decode(unsigned width, const char* s, unsigned s_size, int* result)
{
  static int first_call = 1;
  static int digits_values_upper[128U];
  static int digits_values_lower[128U];
  static const char*
    ie_range = "internal error hy36decode: integer value out of range.";
  unsigned i;
  int di;
  const char* errmsg;
  if (first_call) {
    first_call = 0;
    for(i=0;i<128U;i++) digits_values_upper[i] = -1;
    for(i=0;i<128U;i++) digits_values_lower[i] = -1;
    for(i=0;i<36U;i++) {
      di = digits_upper()[i];
      if (di < 0 || di > 127) {
        *result = 0;
        return ie_range;
      }
      digits_values_upper[di] = i;
    }
    for(i=0;i<36U;i++) {
      di = digits_lower()[i];
      if (di < 0 || di > 127) {
        *result = 0;
        return ie_range;
      }
      digits_values_lower[di] = i;
    }
  }
  if (s_size == width) {
    di = s[0];
    if (di >= 0 && di <= 127) {
      if (digits_values_upper[di] >= 10) {
        errmsg = decode_pure(digits_values_upper, 36U, s, s_size, result);
        if (errmsg == 0) {
          /* result - 10*36**(width-1) + 10**width */
          if      (width == 4U) (*result) -= 456560;
          else if (width == 5U) (*result) -= 16696160;
          else {
            *result = 0;
            return unsupported_width();
          }
          return 0;
        }
      }
      else if (digits_values_lower[di] >= 10) {
        errmsg = decode_pure(digits_values_lower, 36U, s, s_size, result);
        if (errmsg == 0) {
          /* result + 16*36**(width-1) + 10**width */
          if      (width == 4U) (*result) += 756496;
          else if (width == 5U) (*result) += 26973856;
          else {
            *result = 0;
            return unsupported_width();
          }
          return 0;
        }
      }
      else {
        errmsg = decode_pure(digits_values_upper, 10U, s, s_size, result);
        if (errmsg) return errmsg;
        if (!(width == 4U || width == 5U)) {
          *result = 0;
          return unsupported_width();
        }
        return 0;
      }
    }
  }
  *result = 0;
  return invalid_number_literal();
}

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

const std::vector<std::string> & PDB::getRemark()const{
  return remark;
}

void PDB::addRemark( const std::vector<std::string>& v1 ){
  remark.insert(remark.begin(),v1.begin(),v1.end());
}

const std::vector<AtomNumber> & PDB::getAtomNumbers()const{
  return numbers;
}

std::string PDB::getAtomName(AtomNumber a)const{
  const auto p=number2index.find(a);
  if(p==number2index.end()) return "";
  else return atomsymb[p->second];
}

unsigned PDB::getResidueNumber(AtomNumber a)const{
  const auto p=number2index.find(a);
  if(p==number2index.end()) return 0;
  else return residue[p->second];
}

std::string PDB::getResidueName(AtomNumber a) const{
  const auto p=number2index.find(a);
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

      {
        int result;
        const char* errmsg = h36::hy36decode(5, serial.c_str(),5, &result);
        if(errmsg){
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

void PDB::setArgKeyword( const std::string& new_args ){
  bool replaced=false;
  for(unsigned i=0;i<remark.size();++i){
      if( remark[i].find("ARG=")!=std::string::npos){
          remark[i]=new_args; replaced=true;
      }
  }
  plumed_assert( replaced );
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
     const auto p=number2index.find(a);
     if(p==number2index.end()) plumed_merror("atom not available");
     else return positions[p->second];
}



}

