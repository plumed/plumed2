/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2016 The plumed team
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

#include "Stopwatch.h"
#include "Exception.h"

#include <cstdio>
#include <iostream>

/*
Different clocks can be used

gettimeofday (default):
  seems portable
clock_gettime (#define __PLUMED_HAS_CLOCK_GETTIME)
  requires linking -lrt (on linux)
*/

#ifdef __PLUMED_HAS_CLOCK_GETTIME
#include <time.h>
#elif __PLUMED_HAS_GETTIMEOFDAY
// this is the default
#include <sys/time.h>
#endif

using namespace std;

namespace PLMD{

// this is needed for friend operators
std::ostream& operator<<(std::ostream&os,const Stopwatch&sw){
  return sw.log(os);
}

Stopwatch::Time::operator double()const{
  return sec+0.000000001*nsec;
}

Stopwatch::Time Stopwatch::Time::get(){
  Time t;
#ifdef __PLUMED_HAS_CLOCK_GETTIME
  timespec ts;
  clock_gettime(CLOCK_REALTIME,&ts);
  t.sec=ts.tv_sec;
  t.nsec=ts.tv_nsec;
#elif __PLUMED_HAS_GETTIMEOFDAY
  timeval tv;
  gettimeofday(&tv,NULL);
  t.sec=tv.tv_sec;
  t.nsec=1000*tv.tv_usec;
#else
  t.sec=0;
  t.nsec=0;
#endif
  return t;
}

void Stopwatch::Time::reset(){
  sec=0;
  nsec=0;
}

Stopwatch::Time::Time():
  sec(0),nsec(0) { }

Stopwatch::Time Stopwatch::Time::operator-(const Time&t2)const{
  Time t(*this);
  if(t.nsec<t2.nsec){
    t.sec--;
    t.nsec+=1000000000;
  }
  plumed_assert(t.nsec>=t2.nsec);
  t.nsec-=t2.nsec;
  t.sec-=t2.sec;
  return t;
}

const Stopwatch::Time & Stopwatch::Time::operator+=(const Time&t2){
  Time &t(*this);
  t.nsec+=t2.nsec;
  if(t.nsec>1000000000){
    t.nsec-=1000000000;
    t.sec++;
  }
  t.sec+=t2.sec;
  return t;
}

Stopwatch::Watch::Watch():
  cycles(0),running(0) { }

void Stopwatch::Watch::start(){
  running++;
  lastStart=Time::get();
}

void Stopwatch::Watch::stop(){
  pause();
  cycles++;
  total+=lap;
  if(lap>max)max=lap;
  if(min>lap || cycles==1)min=lap;
  lap.reset();
} 

void Stopwatch::Watch::pause(){
  plumed_assert(running>0);
  running--;
  if(running!=0) return;
  lap+=Time::get()-lastStart;
} 

void Stopwatch::start(const std::string & name){
  watches[name].start();
}

void Stopwatch::stop(const std::string & name){
  watches[name].stop();
}

void Stopwatch::pause(const std::string & name){
  watches[name].pause();
}


std::ostream& Stopwatch::log(std::ostream&os)const{
  char buffer[1000];
  buffer[0]=0;
  for(unsigned i=0;i<40;i++) os<<" ";
  os<<"      Cycles        Total      Average      Minumum      Maximum\n";
  for(map<string,Watch>::const_iterator it=watches.begin();it!=watches.end();++it){
    const Watch&t((*it).second);
    std::string name((*it).first);
    os<<name;
    for(unsigned i=name.length();i<40;i++) os<<" ";
    std::sprintf(buffer,"%12u %12.6f %12.6f %12.6f %12.6f\n", t.cycles, double(t.total), double(t.total/t.cycles), double(t.min),double(t.max));
    os<<buffer;
  }
  return os;
}

}




