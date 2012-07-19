/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The PLUMED team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of PLUMED, version 2.0.

   PLUMED is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   PLUMED is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with PLUMED.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "PDB.h"
#include "Tools.h"
#include <cstdio>
#include <iostream>

using namespace std;

namespace PLMD{

std::string PDB::documentation(){
  std::ostringstream ostr;
  ostr<<"In PDB files the atomic coordinates and box lengths should be in Angstroms unless you are working with natural units. ";
  ostr<<"If you are working with natural units then the coordinates should be in your natural length unit. In the PDB files used ";
  ostr<<"by plumed the beta column is used to specify the charges on the atoms and the occupancy column is used to specify the atomic masses. ";
  ostr<<"For more details on the PDB file format visit http://www.wwpdb.org/docs.html";
  return ostr.str();
}

const std::vector<Vector> & PDB::getPositions()const{
  return positions;
}

const std::vector<double> & PDB::getOccupancy()const{
  return occupancy;
}

const std::vector<double> & PDB::getBeta()const{
  return beta;
}

const std::vector<AtomNumber> & PDB::getAtomNumbers()const{
  return numbers;
}

unsigned PDB::size()const{
  return positions.size();
}

void PDB::read(const std::string&file,bool naturalUnits,double scale){
  if(naturalUnits) scale=1.0;
  FILE* fp=fopen(file.c_str(),"r");
  //cerr<<file<<endl;
  string line;
  while(Tools::getline(fp,line)){
    while(line.length()<80) line.push_back(' ');
    string record=line.substr(0,6);
    string serial=line.substr(6,5);
    string atomname=line.substr(12,4);
    string chainID=line.substr(21,1);
    string resnum=line.substr(22,4);
    string x=line.substr(30,8);
    string y=line.substr(38,8);
    string z=line.substr(46,8);
    string occ=line.substr(54,6);
    string bet=line.substr(60,6);
    Tools::trim(record);
    if(record=="TER") break;
    if(record=="END") break;
    if(record=="ATOM" || record=="HETATM"){
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
      std::size_t startpos=atomname.find_first_not_of(" \t");
      std::size_t endpos=atomname.find_last_not_of(" \t");
      atomsymb.push_back( atomname.substr(startpos, endpos-startpos+1) );
      residue.push_back(resno);
      chain.push_back(chainID);
      occupancy.push_back(o);
      beta.push_back(b);
      positions.push_back(p);
    }
  }
  fclose(fp);
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

bool PDB::getBackbone( const unsigned& resnumber, const std::vector<std::string>& backnames, std::vector<AtomNumber>& backnumbers ) const {
  plumed_assert( backnames.size()==backnumbers.size() ); 
  bool foundresidue=false; std::vector<bool> found( backnames.size(), false );
  for(unsigned i=0;i<size();++i){
     if( residue[i]==resnumber ){
         foundresidue=true;
         for(unsigned bno=0;bno<backnames.size();++bno){
            if( backnames[bno]==atomsymb[i] ){ backnumbers[bno]=numbers[i]; found[bno]=true; } 
         }
     }
  }
  if(!foundresidue){ return false; }
  for(unsigned i=0;i<found.size();++i){
      if( !found[i] ) return false;
  }
  return true;
}

std::string PDB::getChainID(const unsigned& resnumber) const {
  for(unsigned i=0;i<size();++i){
     if(resnumber==residue[i]) return chain[i];
  }
  plumed_massert(0,"Not enough residues in pdb input file");
}

void PDB::renameAtoms( const std::string& old_name, const std::string& new_name ){
  for(unsigned i=0;i<size();++i){
      if( atomsymb[i]==old_name ) atomsymb[i]=new_name;
  } 
}

}

