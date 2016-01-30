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

namespace PLMD {
namespace gridtools {

class FindContour : public ActionWithInputGrid {
private:
  OFile of;
  double lenunit;
  unsigned gbuffer;
  std::string fmt_xyz;
  unsigned mycomp;
  double contour;
  std::vector<bool> nosearch_dirs;
public:
  static void registerKeywords( Keywords& keys );
  explicit FindContour(const ActionOptions&ao);
  bool checkAllActive() const { return gbuffer==0; }
  void performOperationsWithGrid( const bool& from_update );
  double getDifferenceFromContour( const std::vector<double>& x, std::vector<double>& der );
  unsigned getNumberOfDerivatives(){ return 0; }
  void performTask( const unsigned& task_index, const unsigned& current, MultiValue& myvals ) const {}
  bool isPeriodic(){ return false; }
};

PLUMED_REGISTER_ACTION(FindContour,"FIND_CONTOUR")

void FindContour::registerKeywords( Keywords& keys ){
  ActionWithInputGrid::registerKeywords( keys );
  keys.add("compulsory","CONTOUR","the value we would like to draw the contour at in the space");
// We want a better way of doing this bit
  keys.add("compulsory", "FILE", "file on which to output coordinates");
  keys.add("compulsory", "UNITS","PLUMED","the units in which to print out the coordinates. PLUMED means internal PLUMED units");
  keys.add("compulsory","SEARCHDIR","all","In which directions do you wish to search for the contour.  By default the code searches "
                                          "for points on the dividing surface in all directions.  You may whish, however to search along one direction "
                                          "for points on the dividing surface.  You can specify the a subset of directions as numbers so if your grid "
                                          "is thre dimensional specifying SEARCHDIR=1,2 means that only points found along grid lines parallel to the "
                                          "first and second axis of the grid will be looked for.  The code will not search for the dividing surface "
                                          "along the third axis of the grid");
  keys.add("compulsory","BUFFER","0","number of buffer grid points around location where grid was found on last step.  If this is zero the full grid is calculated on each step");
  keys.add("optional","COMPONENT","if your input is a vector field use this to specifiy the component of the input vector field for which you wish to find the contour");
  keys.add("optional", "PRECISION","The number of digits in trajectory file");  
}

FindContour::FindContour(const ActionOptions&ao):
Action(ao),
ActionWithInputGrid(ao),
nosearch_dirs( mygrid->getDimension() )
{
  if( mygrid->getNumberOfComponents()==1 ){ 
     mycomp=0; 
  } else {
     int tcomp=-1; parse("COMPONENT",tcomp);
     if( tcomp<0 ) error("component of vector field was not specified - use COMPONENT keyword");
     mycomp=tcomp;
  }
  if( mygrid->noDerivatives() ) error("cannot find contours if input grid has no derivatives");

  parse("CONTOUR",contour); parse("BUFFER",gbuffer);
  log.printf("  calculating dividing surface along which function equals %f \n", contour);
  if( gbuffer>0 ) log.printf("  after first step a subset of only %d grid points around where the countour was found will be checked\n",gbuffer);

  std::string searchdir_str; parse("SEARCHDIR",searchdir_str);
  if( searchdir_str=="all" ){
     nosearch_dirs.assign( nosearch_dirs.size(), false );
  } else {
     std::vector<std::string> searchdirs = Tools::getWords(searchdir_str,"\t\n ,");
     if( searchdirs.size()>mygrid->getDimension() ) error("number of search directions is larger than number of input dimensions");
     nosearch_dirs.assign( nosearch_dirs.size(), true );
     unsigned nn; Tools::convert(searchdirs[0],nn); nosearch_dirs[nn-1]=false;
     log.printf("  searching for contour along %d", nn);
     for(unsigned i=1;i<searchdirs.size()-1;++i){
         Tools::convert(searchdirs[i],nn); 
         nosearch_dirs[nn-1]=false; log.printf(", %d",nn);
     }
     if( searchdirs.size()>1 ){ 
         Tools::convert(searchdirs[searchdirs.size()-1],nn); 
         nosearch_dirs[nn-1]=false; log.printf(" and %d directions \n",nn);
     } else log.printf(" direction \n");
  }
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

double FindContour::getDifferenceFromContour( const std::vector<double>& x, std::vector<double>& der ){
  return mygrid->getValueAndDerivatives( x, mycomp, der ) - contour;
}

void FindContour::performOperationsWithGrid( const bool& from_update ){
  std::vector< std::vector<double> > contour_points( mygrid->getDimension()*mygrid->getNumberOfPoints() );

  // Two points for search
  std::vector<unsigned> ind( mygrid->getDimension() );
  std::vector<double> dx( mygrid->getGridSpacing() );
  std::vector<double> direction( mygrid->getDimension(), 0 );
  // Retrieve nbin from grid (remember this returns nbin for restart)
  std::vector<unsigned> nbin( mygrid->getNbin() );
  for(unsigned i=0;i<nbin.size();++i){
     if( !mygrid->isPeriodic(i) ) nbin[i]+=1;
  }

  // Run over whole grid
  std::vector<bool> active( mygrid->getNumberOfPoints(), false );
  std::vector<unsigned> neighbours; unsigned num_neighbours;
  std::vector<unsigned> gbuffer_vec( mygrid->getDimension(), gbuffer );
  std::vector<unsigned> ugrid_indices( mygrid->getDimension() );
  unsigned npoints=0; RootFindingBase<FindContour> mymin( this );
  for(unsigned i=0;i<mygrid->getNumberOfPoints();++i){
     // Get the index of the current grid point
     mygrid->getIndices( i, ind ); 
     // Ensure inactive grid points are ignored
     bool cycle=false;
     if( mygrid->inactive(i) ) continue ;
     for(unsigned j=0;j<mygrid->getDimension();++j) ugrid_indices[j]=ind[j];
     for(unsigned j=0;j<mygrid->getDimension();++j){
        if( mygrid->isPeriodic(j) ) ugrid_indices[j]=(ugrid_indices[j]+1)%nbin[j];
        else if( (ugrid_indices[j]+1)==nbin[j] ) continue ;
        else ugrid_indices[j]+=1;
        // Now check the grid is active
        if( mygrid->inactive( mygrid->getIndex( ugrid_indices ) ) ){ cycle=true; break; }
        // And reset ugrid_indices
        ugrid_indices[j]=ind[j];
     }
     if( cycle ) continue ;

     // Get the value of a point on the grid
     double val1=getGridElement( i, mycomp*(mygrid->getDimension()+1) ) - contour;

     bool edge=false;
     for(unsigned j=0;j<mygrid->getDimension();++j){
         if( nosearch_dirs[j] ) continue ;
         // Make sure we don't search at the edge of the grid
         if( !mygrid->isPeriodic(j) && (ind[j]+1)==nbin[j] ) continue;
         else if( (ind[j]+1)==nbin[j] ){ edge=true; ind[j]=0; }
         else ind[j]+=1; 
         double val2=getGridElement(ind,mycomp*(mygrid->getDimension()+1)) - contour;
         if( val1*val2<0 ){
             // Use initial point location as first guess for search
             contour_points[npoints].resize( mygrid->getDimension() );  
             mygrid->getGridPointCoordinates( i, contour_points[npoints] );  
             // Setup direction vector
             direction[j]=0.999999999*dx[j];
             // And do proper search for contour point
             mymin.linesearch( direction, contour_points[npoints], &FindContour::getDifferenceFromContour );
             direction[j]=0.0; 
             // For next run through activate buffer region around this grid point
             if( gbuffer>0 ){
                 mygrid->getNeighbors( contour_points[npoints], gbuffer_vec, num_neighbours, neighbours );
                 for(unsigned n=0;n<num_neighbours;++n) active[ neighbours[n] ]=true; 
             }
             npoints++;
         }
         if( mygrid->isPeriodic(j) && edge ){ edge=false; ind[j]=nbin[j]-1; }
         else ind[j]-=1;
     } 
   
  }
  if( gbuffer>0 ) mygrid->activateThesePoints( active );
  // Clear the grid ready for next time
  if( from_update ) mygrid->reset();

  of.printf("%u\n",npoints);
  of.printf("Points found on isocontour\n");
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
