#include "PDB.h"
#include "Tools.h"
#include <cstdio>

using namespace std;

namespace PLMD{

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

const std::vector<std::string> & PDB::getAtomNames()const{
  return atomNames;
}

const std::vector<AtomNumber> & PDB::getResidueNumbers()const{
  return residues;
}

unsigned PDB::size()const{
  return positions.size();
}

void PDB::read(const std::string&file,double scale){
  FILE* fp=fopen(file.c_str(),"r");
  string line;
  while(Tools::getline(fp,line)){
    while(line.length()<80) line.push_back(' ');
    string record=line.substr(0,6);
    string atnam=line.substr(12,5);
    string serial=line.substr(6,5);
    string resser=line.substr(22,4);
    string x=line.substr(30,8);
    string y=line.substr(38,8);
    string z=line.substr(46,8);
    string occ=line.substr(54,6);
    string bet=line.substr(60,6);
    Tools::trim(record);
    if(record=="TER") break;
    if(record=="END") break;
    if(record=="ATOM" || record=="HETATM"){
      AtomNumber a,r;
      double o,b;
      Vector p;
      Tools::convert(serial,a);
      Tools::convert(resser,r);
      Tools::convert(occ,o);
      Tools::convert(bet,b);
      Tools::convert(x,p[0]);
      Tools::convert(y,p[1]);
      Tools::convert(z,p[2]);
      p.scale(scale);
      numbers.push_back(a);
      std::vector<std::string> atnam2=Tools::getWords(atnam);
      assert( atnam2.size()==1 );
      atomNames.push_back( atnam2[0] );
      residues.push_back(r);
      occupancy.push_back(o);
      beta.push_back(b);
      positions.push_back(p);
    }
  }
  fclose(fp);
}


}

