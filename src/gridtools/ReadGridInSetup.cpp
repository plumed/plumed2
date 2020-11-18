/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2018 The plumed team
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
#include "CreateGridInSetup.h"
#include "core/ActionRegister.h"
#include "tools/IFile.h"

using namespace std;

namespace PLMD {
namespace gridtools {

class ReadGridInSetup : public CreateGridInSetup {
public:
  static void registerKeywords( Keywords& keys );
  explicit ReadGridInSetup(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(ReadGridInSetup,"REFERENCE_GRID")

void ReadGridInSetup::registerKeywords( Keywords& keys ) {
  CreateGridInSetup::registerKeywords(keys);
  keys.add("compulsory","FILE","the name of the file that contains the reference data");
  keys.add("compulsory","VALUE","the name of the value that should be read from the grid");
}

ReadGridInSetup::ReadGridInSetup(const ActionOptions&ao):
  Action(ao),
  CreateGridInSetup(ao)
{
   std::string valuestr; parse("VALUE",valuestr);
   std::string tstyle, filen; parse("FILE",filen);
   if( filen.length()>0 ) {
       std::size_t dot=filen.find_first_of(".");
       if( dot!=std::string::npos ) tstyle=filen.substr(dot+1);
       if( tstyle!="grid" ) error("can only read in grids using read value in setup");
       log.printf("  reading function %s on grid from file %s \n", valuestr.c_str(), filen.c_str() );
   }
   IFile ifile; ifile.open(filen);
   if( !ifile.FieldExist( valuestr ) ) error("could not find grid value in input file"); 
   std::vector<std::string> fieldnames; ifile.scanFieldList( fieldnames );
 
   // Retrieve the names of the variables the grid is computed over
   bool flatgrid=false;
   for(unsigned i=0;i<fieldnames.size();++i) {
       if( fieldnames[i].find("min_")!=std::string::npos ) flatgrid=true;
       std::size_t dot = fieldnames[i].find_first_of("d" + valuestr + "_" );
       if( fieldnames[i].find("d" + valuestr + "_")!=std::string::npos ) labels.push_back( fieldnames[i].substr(dot+2+valuestr.length()) );
   }
   if( flatgrid && labels.size()==0 ) error("could not find any derivatives for value " + valuestr + " in input file.  Header should contain at least columns with a name starting d" + valuestr + "_");
   // Now get all the header data for the grid
   std::vector<std::string> gmin( labels.size() ), gmax( labels.size() ); std::string pstring;
   int gbin1; std::vector<unsigned> gbin( labels.size() ); std::vector<bool> ipbc( labels.size() );
   if( !flatgrid ) {
     ifile.scanField( "nbins", gbin1); gbin[0]=gbin1; 
     createGridAndValue( "fibonacci", ipbc, gbin[0], gmin, gmax, gbin );
   } else {
     for(unsigned i=0; i<labels.size(); ++i) {
       ifile.scanField( "min_" + labels[i], gmin[i]);
       ifile.scanField( "max_" + labels[i], gmax[i]);
       ifile.scanField( "periodic_" + labels[i], pstring );
       ifile.scanField( "nbins_" + labels[i], gbin1); gbin[i]=gbin1;
       if( pstring=="true" ) {
           log.printf("   for periodic coordinate %s minimum is %s maximum is %s and number of bins is %d \n",labels[i].c_str(),gmin[i].c_str(),gmax[i].c_str(),gbin[i]);
           ipbc[i]=true; 
       } else if( pstring=="false" ) {
           log.printf("   for coordinate %s minimum is %s maximum is %s and number of bins is %d \n",labels[i].c_str(),gmin[i].c_str(),gmax[i].c_str(),gbin[i]);
           ipbc[i]=false; 
       } else error("do not understand periodidicy of " + labels[i] );
 
       bool hasder=ifile.FieldExist( "d" + valuestr + "_" + labels[i] );
       if( !hasder ) plumed_merror("missing derivatives from grid file");
     }
     createGridAndValue( "flat", ipbc, 0, gmin, gmax, gbin );
   }
   // And finally read all the grid points
   Value* valout=getPntrToOutput(0); 
   std::vector<double> dder( labels.size() ), xx( labels.size() );
   for(unsigned i=0;i<valout->getNumberOfValues( getLabel() );++i) {
     double x, val, norm; ifile.scanField( valuestr, val ); ifile.scanField( "normalisation", norm );
     for(unsigned j=0; j<labels.size(); ++j) {
       ifile.scanField(labels[j],x); 
       if( !flatgrid ) {
          ifile.scanField("nbins", gbin1);
       } else {
          xx[j]=x+gridobject.getGridSpacing()[j]/2.0;
          ifile.scanField( "min_" + labels[j], gmin[j]);
          ifile.scanField( "max_" + labels[j], gmax[j]);
          ifile.scanField( "nbins_" + labels[j], gbin1);
          ifile.scanField( "periodic_" + labels[j], pstring );
       }
     }
     for(unsigned j=0; j<labels.size(); ++j) ifile.scanField( "d" + valuestr + "_" + labels[j], dder[j] );  

     unsigned index=gridobject.getIndex(xx); if( !flatgrid ) index=i;
     valout->setNorm( norm ); valout->add( index*(labels.size()+1), norm*val );
     for(unsigned j=0; j<labels.size(); ++j) valout->add( index*(labels.size()+1)+1+j, norm*dder[j] );  
     ifile.scanField();
   }
   ifile.close();
}

}
}

