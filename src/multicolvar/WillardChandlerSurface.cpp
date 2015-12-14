/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2015 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#include "core/ActionAtomistic.h"
#include "core/ActionPilot.h"
#include "core/ActionRegister.h"
#include "tools/Pbc.h"
#include "tools/File.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/Units.h"
#include <cstdio>
#include "core/SetupMolInfo.h"
#include "core/ActionSet.h"
#include "MultiColvarBase.h"
#include "tools/Grid.h"
#include "tools/OFile.h"
#include "tools/KernelFunctions.h"
#include "vesselbase/ActionWithInputVessel.h"
#include "vesselbase/StoreDataVessel.h"

using namespace std;

namespace PLMD {
namespace multicolvar {

//+PLUMEDOC ANALYSIS WCSURFACE
/*
Calculate the Willard Chanlder dividing surface

\par Examples


*/
//+ENDPLUMEDOC

class WillardChandlerSurface :
  public ActionPilot,
  public ActionAtomistic,
  public vesselbase::ActionWithInputVessel
{
  std::string kerneltype;
  unsigned rstride;
  std::string filename;
  unsigned dir;
  double contour;
  OFile of;
  double lenunit;
  std::string fmt_xyz;
  MultiColvarBase* mycolv; 
  std::vector<unsigned> nbins;
  std::vector<double> bw;
public:
  explicit WillardChandlerSurface(const ActionOptions&);
  static void registerKeywords( Keywords& keys );
  void calculate(){}
  void calculateNumericalDerivatives( ActionWithValue* a=NULL ){ plumed_error(); }
  void apply(){}
  void update();
};

PLUMED_REGISTER_ACTION(WillardChandlerSurface,"WCSURFACE")

void WillardChandlerSurface::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  ActionWithInputVessel::registerKeywords( keys );
  keys.add("compulsory","STRIDE","1","the frequency with which the data should be collected and added to the grid");
  keys.add("atoms","ORIGIN","we will use the position of this atom as the origin");
  keys.add("compulsory","NBINS","the number of bins to use to represent the density profile");
  keys.add("compulsory","BANDWIDTH","the bandwidths for kernel density esimtation");
  keys.add("compulsory","KERNEL","gaussian","the kernel function you are using.  More details on  the kernels available "
                                            "in plumed plumed can be found in \\ref kernelfunctions.");
  keys.add("compulsory","CONTOUR","the value we would like to draw the contour at in the space");
// We want a better way of doing this bit
  keys.add("compulsory", "FILE", "file on which to output coordinates");
  keys.add("compulsory", "UNITS","PLUMED","the units in which to print out the coordinates. PLUMED means internal PLUMED units");
  keys.add("optional", "PRECISION","The number of digits in trajectory file");
// Better way for this bit 
}

WillardChandlerSurface::WillardChandlerSurface(const ActionOptions&ao):
  Action(ao),
  ActionPilot(ao),
  ActionAtomistic(ao),
  ActionWithInputVessel(ao)
{
  std::vector<AtomNumber> atom; parseAtomList("ORIGIN",atom);
  if( atom.size()!=1 ) error("should only be one atom specified");
  log.printf("  origin is at position of atom : %d\n",atom[0].serial() );

  readArgument("store");
  mycolv = dynamic_cast<MultiColvarBase*>( getDependencies()[0] );
  plumed_assert( getDependencies().size()==1 ); 
  if(!mycolv) error("action labeled " + mycolv->getLabel() + " is not a multicolvar");

  parse("KERNEL",kerneltype); parseVector("BANDWIDTH",bw); parseVector("NBINS",nbins); parse("CONTOUR",contour);
  if( bw.size()!=3 || nbins.size()!=3 ) error("wrong size found for bandwidth vector or number of bins");
  log.printf("  calculating dividing surface where average equals %f along ", contour);
  log.printf(" for colvars calculated by action %s \n",mycolv->getLabel().c_str() );

  // START OF BIT TO IMPROVE
  std::string file; parse("FILE",file);
  if(file.length()==0) error("name out output file was not specified");
  std::string type=Tools::extension(file);
  log<<"  file name "<<file<<"\n";
  if(type!="xyz") error("can only print xyz file type with DUMPMULTICOLVAR");

  fmt_xyz="%f";
  std::string precision; parse("PRECISION",precision);
  if(precision.length()>0){
    int p; Tools::convert(precision,p);
    log<<"  with precision "<<p<<"\n";
    string a,b;
    Tools::convert(p+5,a);
    Tools::convert(p,b);
    fmt_xyz="%"+a+"."+b+"f";
  }

  std::string unitname; parse("UNITS",unitname);
  if(unitname!="PLUMED"){
    Units myunit; myunit.setLength(unitname);
    lenunit=plumed.getAtoms().getUnits().getLength()/myunit.getLength();
  }
  else lenunit=1.0;
  // END OF BIT TO IMPROVE

  of.link(*this);
  of.open(file);

  checkRead(); requestAtoms(atom); 
  // Stupid dependencies cleared by requestAtoms - why GBussi why? That's got me so many times
  addDependency( mycolv );
}

void WillardChandlerSurface::update(){
  std::vector<bool> pbc(nbins.size()); std::vector<double> min(nbins.size()), max(nbins.size());
  std::vector<std::string> args(nbins.size()), gmin(nbins.size()), gmax(nbins.size());;
  args[0]="x"; args[1]="y"; args[2]="z";
  for(unsigned i=0;i<3;++i){ min[i]=-0.5; max[i]=0.5; pbc[i]=true; }

  if( !mycolv->getPbc().isOrthorombic() ){
      error("I think that Willard-Chandler surfaces with non-orthorhombic cells don't work.  If you want it have a look and see if you can work it out");
  }

  // Get the box (this is going to give us the extent of the grid)
  for(unsigned i=0;i<3;++i){ min[i]*=mycolv->getBox()(i,i); max[i]*=mycolv->getBox()(i,i); }

  // This will be the only thing we keep eventually
  for(unsigned i=0;i<3;++i){ Tools::convert(min[i],gmin[i]); Tools::convert(max[i],gmax[i]); }

  // Setup the grid
  std::string funcl=mycolv->getLabel() + ".dens";
  Grid gg(funcl,args,gmin,gmax,nbins,true,true,true,pbc,gmin,gmax);

  vesselbase::StoreDataVessel* stash=dynamic_cast<vesselbase::StoreDataVessel*>( getPntrToArgument() );
  plumed_dbg_assert( stash ); std::vector<double> cvals( mycolv->getNumberOfQuantities() );
  Vector origin = getPosition(0); std::vector<double> pp( 3 ); Vector fpos;

  unsigned rank=comm.Get_rank(), size=comm.Get_size();

  for(unsigned i=rank;i<mycolv->getFullNumberOfTasks();i+=size){
      stash->retrieveValue( i, false, cvals );
      Vector apos = pbcDistance( mycolv->getCentralAtomPos( mycolv->getTaskCode(i) ), origin );

      // fpos = getPbc().realToScaled( apos ); Ideally want to do with scaled coordinates eventually GAT
      for(unsigned j=0;j<3;++j) pp[j]=apos[j];

      KernelFunctions kernel( pp, bw, kerneltype, false, cvals[0]*cvals[1], true );
      gg.addKernel( kernel ); 
  }
  gg.mpiSumValuesAndDerivatives( comm );

  unsigned npoints=0; std::vector<std::vector<double> > contour_points;
  gg.findSetOfPointsOnContour( contour, npoints, contour_points );
  if(npoints==0 ) warning("found no points on Willard-Chandler surface try changing the CONTOUR parameter"); 

  of.printf("%u\n",npoints);
  const Tensor & t(mycolv->getPbc().getBox());
  if(mycolv->getPbc().isOrthorombic()){
    of.printf((" "+fmt_xyz+" "+fmt_xyz+" "+fmt_xyz+"\n").c_str(),lenunit*t(0,0),lenunit*t(1,1),lenunit*t(2,2));
  }else{
    of.printf((" "+fmt_xyz+" "+fmt_xyz+" "+fmt_xyz+" "+fmt_xyz+" "+fmt_xyz+" "+fmt_xyz+" "+fmt_xyz+" "+fmt_xyz+" "+fmt_xyz+"\n").c_str(),
                 lenunit*t(0,0),lenunit*t(0,1),lenunit*t(0,2),
                 lenunit*t(1,0),lenunit*t(1,1),lenunit*t(1,2),
                 lenunit*t(2,0),lenunit*t(2,1),lenunit*t(2,2)
           );
  }
  const char* nn="X"; Vector cpos;
  for(unsigned i=0;i<npoints;++i){
     for(unsigned j=0;j<3;++j) cpos[j]=contour_points[i][j];
     // cpos=mycolv->getPbc().scaledToReal(fpos);
     of.printf(("%s "+fmt_xyz+" "+fmt_xyz+" "+fmt_xyz).c_str(),nn,lenunit*cpos[0], lenunit*cpos[1], lenunit*cpos[2] );
     of.printf("\n");
  }
}

}
}
