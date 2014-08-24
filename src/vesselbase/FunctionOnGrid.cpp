/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#include "FunctionOnGrid.h"
#include "VesselRegister.h"

namespace PLMD{
namespace vesselbase{

PLUMED_REGISTER_VESSEL(FunctionOnGrid,"GRID_NOSPLINE")

FunctionOnGrid* FunctionOnGrid::spawn( const GridVesselBase* myg, const std::vector<std::string>& args, const std::vector<unsigned>& nbins, ActionWithVessel* aa ){
  
  std::string ogmin, ogmax, ogbin, ognam, ogper, sbin; unsigned n=0;
  std::vector<std::string> gmin( myg->getMin() ), gmax( myg->getMax() );
  for(unsigned i=0;i<args.size();++i){
      bool found=false;
      for(unsigned j=0;j<myg->getDimension();++j){
          found=true;
          if( args[i]==myg->getQuantityDescription(j) ){
              Tools::convert( nbins[j], sbin ); n++;
              if(n==1){
                 ogmin="MIN=" + gmin[j];
                 ogmax="MAX=" + gmax[j]; 
                 ogbin="NBIN=" + sbin;
                 ognam="ARGS=" + myg->getQuantityDescription(j);
                 if( myg->pbc[j] ) ogper="PERIODIC=yes"; 
                 else ogper="PERIODIC=no";
              } else {
                 ogmin += "," + gmin[j]; 
                 ogmax += "," + gmax[j]; 
                 ogbin += "," + sbin;
                 ognam += "," + myg->getQuantityDescription(j);
                 if( myg->pbc[j] ) ogper += ",yes";
                 else ogper += ",no";
              }
          }
      }
      if(!found) aa->error( args[i] + " is not in the input grid" );
  }
  if( n!=args.size() ) aa->error("problems finding required arguments in input grid");

  std::string grid_input = ogmin + " " + ogmax + " " + ogbin + " " + ognam + " " + ogper;
  VesselOptions da( "GRID_NOSPLINE", "", 0, grid_input, aa );
  Keywords mykeys; FunctionOnGrid::registerKeywords( mykeys );
  VesselOptions ba( da, mykeys );
  FunctionOnGrid* outgrid = new FunctionOnGrid( ba );
  return outgrid;
}

void FunctionOnGrid::reserveKeyword( Keywords& keys ){
  keys.reserve("optional","GRID_NOSPLINE","create a grid to store a function");
} 

void FunctionOnGrid::registerKeywords( Keywords& keys ){
  GridVesselBase::registerKeywords( keys );
  keys.add("compulsory","ARGS","names of arguments for grid dimensions");
  keys.add("compulsory","PERIODIC","is this input variable periodic");
}

FunctionOnGrid::FunctionOnGrid( const VesselOptions& da ):
GridVesselBase(da)
{
  std::vector<std::string> names; parseVector("ARGS",names);  
  names.push_back( getAction()->getLabel() );

  std::vector<std::string> periodic; 
  parseVector("PERIODIC",periodic);
  std::vector<bool> mypbc( periodic.size() );
  for(unsigned i=0;i<periodic.size();++i){
      if( periodic[i]=="yes" ) mypbc[i]=true;
      else mypbc[i]=false;
  }
  finishSetup( 1, mypbc, names );
}

std::string FunctionOnGrid::description(){
  return getGridDescription();
}

bool FunctionOnGrid::calculate(){
  plumed_merror("This should not be called");
  return true;
}

void FunctionOnGrid::finish(){
  plumed_merror("This should not be called");
}

bool FunctionOnGrid::applyForce( std::vector<double>& forces ){
  plumed_merror("This should not be called");
  return false;
}

}
}
