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
#include "Units.h"
#include "Tools.h"

using namespace std;

namespace PLMD{

Units::Units():
  energy(1.0),
  energyString("kj/mol"),
  length(1.0),
  lengthString("nm"),
  time(1.0),
  timeString("ps"),
  charge(1.0),
  chargeString("e"),
  mass(1.0),
  massString("amu")
{
}

void Units::setEnergy(const std::string &s){
  energyString=s;
  if(s=="kj/mol"){
    energy=1.0;
  } else if(s=="kcal/mol"){
    energy=4.184;
  } else if(s=="j/mol"){
    energy=0.001;
  } else if(s=="eV"){
    energy=96.48530749925792;
  } else {
    energy=-1.0;
    energyString="";
    Tools::convert(s,energy);
    plumed_massert(energy>0.0,"energy units should be positive");
  }
}

void Units::setLength(const std::string &s){
  lengthString=s;
  if(s=="nm"){
    length=1.0;
  } else if(s=="A"){
    length=0.1;
  } else if(s=="um"){
    length=1000.0;
  } else {
    length=-1.0;
    lengthString="";
    Tools::convert(s,length);
    plumed_massert(length>0.0,"length units should be positive");
  }
}

void Units::setTime(const std::string &s){
  timeString=s;
  if(s=="ps"){
    time=1.0;
  } else if(s=="ns"){
    time=1000.0;
  } else if(s=="fs"){
    time=0.001;
  } else {
    time=-1.0;
    timeString="";
    Tools::convert(s,time);
    plumed_massert(time>0.0,"time units should be positive");
  }
}

void Units::setCharge(const std::string &s){
  chargeString=s;
  if(s=="e"){
    charge=1.0;
  } else {
    charge=-1.0;
    chargeString="";
    Tools::convert(s,charge);
    plumed_massert(charge>0.0,"charge units should be positive");
  }
}

void Units::setMass(const std::string &s){
  massString=s;
  if(s=="amu"){
    mass=1.0;
  } else {
    mass=-1.0;
    massString="";
    Tools::convert(s,mass);
    plumed_massert(mass>0.0,"mass units should be positive");
  }
}

void Units::setEnergy(const double s){
  energyString="";
  energy=s;
}

void Units::setLength(const double s){
  lengthString="";
  length=s;
}

void Units::setTime(const double s){
  timeString="";
  time=s;
}

void Units::setCharge(const double s){
  chargeString="";
  charge=s;
}

void Units::setMass(const double s){
  massString="";
  mass=s;
}



}

