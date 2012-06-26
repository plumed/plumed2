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
    string x=line.substr(30,8);
    string y=line.substr(38,8);
    string z=line.substr(46,8);
    string occ=line.substr(54,6);
    string bet=line.substr(60,6);
    Tools::trim(record);
    if(record=="TER") break;
    if(record=="END") break;
    if(record=="ATOM" || record=="HETATM"){
      AtomNumber a;
      double o,b;
      Vector p;
      Tools::convert(serial,a);
      Tools::convert(occ,o);
      Tools::convert(bet,b);
      Tools::convert(x,p[0]);
      Tools::convert(y,p[1]);
      Tools::convert(z,p[2]);
      // scale into nm
      p*=scale;
      numbers.push_back(a);
      occupancy.push_back(o);
      beta.push_back(b);
      positions.push_back(p);
    }
  }
  fclose(fp);
}


}

