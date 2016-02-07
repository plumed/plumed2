/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2015 The plumed team
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
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "ActionWithInputGrid.h"
#include "tools/OFile.h"
#include "tools/RootFindingBase.h"
#include "tools/Random.h"

namespace PLMD {
namespace gridtools {

class FindSphericalContour : public ActionWithInputGrid {
private:
  OFile of;
  double lenunit;
  std::string fmt_xyz;
  int rnd;
  unsigned mycomp;
  double contour, offset, increment;
  double min, max;
  unsigned npoints;
public:
  static void registerKeywords( Keywords& keys );
  explicit FindSphericalContour(const ActionOptions&ao);
  void performOperationsWithGrid( const bool& from_update );
  double getDifferenceFromContour( const std::vector<double>& x, std::vector<double>& der );
  unsigned getNumberOfDerivatives(){ return 0; }
  void performTask( const unsigned& task_index, const unsigned& current, MultiValue& myvals ) const {}
  bool isPeriodic(){ return false; }
};

PLUMED_REGISTER_ACTION(FindSphericalContour,"FIND_SPHERICAL_CONTOUR")

void FindSphericalContour::registerKeywords( Keywords& keys ){
  ActionWithInputGrid::registerKeywords( keys );
  keys.add("compulsory","NPOINTS","the number of points for which we are looking for the contour");
  keys.add("compulsory","CONTOUR","the value we would like to draw the contour at in the space");
  keys.add("compulsory","INNER_RADIUS","the minimum radius on which to look for the contour");
  keys.add("compulsory","OUTER_RADIUS","the outer radius on which to look for the contour");
// We want a better way of doing this bit
  keys.add("compulsory", "FILE", "file on which to output coordinates");
  keys.add("compulsory", "UNITS","PLUMED","the units in which to print out the coordinates. PLUMED means internal PLUMED units");
  keys.add("optional","COMPONENT","if your input is a vector field use this to specifiy the component of the input vector field for which you wish to find the contour");
  keys.add("optional", "PRECISION","The number of digits in trajectory file");  
}

FindSphericalContour::FindSphericalContour(const ActionOptions&ao):
Action(ao),
ActionWithInputGrid(ao)
{
  if( mygrid->getDimension()!=3 ) error("input grid must be three dimensional");

  if( mygrid->getNumberOfComponents()==1 ){ 
     mycomp=0; 
  } else {
     int tcomp=-1; parse("COMPONENT",tcomp);
     if( tcomp<0 ) error("component of vector field was not specified - use COMPONENT keyword");
     mycomp=tcomp;
  }
  if( mygrid->noDerivatives() ) error("cannot find contours if input grid has no derivatives");

  parse("NPOINTS",npoints);
  log.printf("  searching for %d points on dividing surface \n");
  parse("CONTOUR",contour); 
  log.printf("  calculating dividing surface along which function equals %f \n", contour);
  parse("INNER_RADIUS",min); parse("OUTER_RADIUS",max);
  log.printf("  expecting to find dividing surface at radii between %f and %f \n",min,max);

  // Set this here so the same set of grid points are used on every turn
  Random random; rnd = std::floor( npoints*random.RandU01() );
  offset = 2 / static_cast<double>( npoints );
  increment = pi*( 3 - sqrt(5) );

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
    std::string a,b;
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
  checkRead();
}

double FindSphericalContour::getDifferenceFromContour( const std::vector<double>& x, std::vector<double>& der ){
  return mygrid->getValueAndDerivatives( x, mycomp, der ) - contour;
}

void FindSphericalContour::performOperationsWithGrid( const bool& from_update ){
  std::vector< std::vector<double> > contour_points( npoints );

  of.printf("%u\n",npoints); const char* nn="X";
  of.printf("Points found on isocontour\n");
  RootFindingBase<FindSphericalContour> mymin( this );
  std::vector<double> tmp( mygrid->getDimension() ), der( mygrid->getDimension() );
  std::vector<double> contour_point( mygrid->getDimension() );
  std::vector<double> direction( mygrid->getDimension() );
  for(unsigned i=0;i<npoints;++i){
      // Generate contour point on inner sphere
      contour_point[1] = ((i*offset) - 1) + (offset/2);   
      double r = sqrt( 1 - pow(contour_point[1],2) );  
      double phi = ((i + rnd)%npoints)*increment;
      contour_point[0] = r*cos(phi);
      contour_point[2] = r*sin(phi);

      // normalize direction vector
      double norm=0; 
      for(unsigned j=0;j<3;++j) norm+=contour_point[j]*contour_point[j];
      norm = sqrt( norm );
      for(unsigned j=0;j<3;++j) direction[j] = contour_point[j] / norm;

      // Now set up direction as vector from inner sphere to outer sphere
      for(unsigned j=0;j<3;++j){
         contour_point[j] = min*direction[j];
         direction[j] = (max-min)*direction[j];
         tmp[j] = contour_point[j] + direction[j]; 
      }

      double val1 = getDifferenceFromContour( contour_point, der );
      double val2 = getDifferenceFromContour( tmp, der );
      if( val1*val2>=0 ) error("range does not bracket the dividing surface");
 
      mymin.linesearch( direction, contour_point, &FindSphericalContour::getDifferenceFromContour );
      of.printf(("%s "+fmt_xyz+" "+fmt_xyz+" "+fmt_xyz).c_str(),nn,lenunit*contour_point[0], lenunit*contour_point[1], lenunit*contour_point[2] );
      of.printf("\n");
  }

  // Clear the grid ready for next time
  if( from_update ) mygrid->reset();
}

}
}
