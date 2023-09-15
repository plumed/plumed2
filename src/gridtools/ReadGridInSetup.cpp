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
#include "ActionWithGrid.h"
#include "core/ActionRegister.h"
#include "tools/IFile.h"

using namespace std;

namespace PLMD {
namespace gridtools {

class ReadGridInSetup : public ActionWithGrid {
private:
  GridCoordinatesObject gridobject;
  std::vector<std::string> dernames;
  void createGridAndValue( const std::string& gtype, const std::vector<bool>& ipbc, const unsigned& nfermi,
                           const std::vector<std::string>& gmin, const std::vector<std::string>& gmax,
                           const std::vector<unsigned>& gbin );
public:
  static void registerKeywords( Keywords& keys );
  explicit ReadGridInSetup(const ActionOptions&ao);
  unsigned getNumberOfDerivatives() override ;
  void setupOnFirstStep() override { plumed_merror("should not be in ReadGridInSetup setupOnFirstStep" ); }
  void performTask( const unsigned& current, MultiValue& myvals ) const override { plumed_merror("should not be in readGridInSetup performTask"); }
  std::vector<std::string> getGridCoordinateNames() const override ;
  const GridCoordinatesObject& getGridCoordinatesObject() const override ;
};

PLUMED_REGISTER_ACTION(ReadGridInSetup,"REFERENCE_GRID")

void ReadGridInSetup::registerKeywords( Keywords& keys ) {
  ActionWithGrid::registerKeywords(keys); keys.remove("SERIAL");
  keys.add("optional","FUNC","the function to compute on the grid");
  keys.add("compulsory","FILE","the name of the file that contains the reference data");
  keys.add("compulsory","VALUE","the name of the value that should be read from the grid");
}

ReadGridInSetup::ReadGridInSetup(const ActionOptions&ao):
  Action(ao),
  ActionWithGrid(ao)
{
   std::string func; parse("FUNC",func);
   if( func.length()>0 ) {

   } else {
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
           if( fieldnames[i].find("d" + valuestr + "_")!=std::string::npos ) dernames.push_back( fieldnames[i].substr(dot+2+valuestr.length()) );
       }
       if( flatgrid && dernames.size()==0 ) error("could not find any derivatives for value " + valuestr + " in input file.  Header should contain at least columns with a name starting d" + valuestr + "_");
       // Now get all the header data for the grid
       std::vector<std::string> gmin( dernames.size() ), gmax( dernames.size() ); std::string pstring;
       int gbin1; std::vector<unsigned> gbin( dernames.size() ); std::vector<bool> ipbc( dernames.size() );
       if( !flatgrid ) {
         ifile.scanField( "nbins", gbin1); gbin[0]=gbin1; 
         createGridAndValue( "fibonacci", ipbc, gbin[0], gmin, gmax, gbin );
       } else {
         for(unsigned i=0; i<dernames.size(); ++i) {
           ifile.scanField( "min_" + dernames[i], gmin[i]);
           ifile.scanField( "max_" + dernames[i], gmax[i]);
           ifile.scanField( "periodic_" + dernames[i], pstring );
           ifile.scanField( "nbins_" + dernames[i], gbin1); gbin[i]=gbin1;
           if( pstring=="true" ) {
               log.printf("   for periodic coordinate %s minimum is %s maximum is %s and number of bins is %d \n",dernames[i].c_str(),gmin[i].c_str(),gmax[i].c_str(),gbin[i]);
               ipbc[i]=true; 
           } else if( pstring=="false" ) {
               log.printf("   for coordinate %s minimum is %s maximum is %s and number of bins is %d \n",dernames[i].c_str(),gmin[i].c_str(),gmax[i].c_str(),gbin[i]);
               ipbc[i]=false; 
           } else error("do not understand periodidicy of " + dernames[i] );
 
           bool hasder=ifile.FieldExist( "d" + valuestr + "_" + dernames[i] );
           if( !hasder ) plumed_merror("missing derivatives from grid file");
         }
         createGridAndValue( "flat", ipbc, 0, gmin, gmax, gbin );
       }
       // And finally read all the grid points
       Value* valout=getPntrToComponent(0); 
       std::vector<double> dder( dernames.size() ), xx( dernames.size() );
       for(unsigned i=0;i<valout->getNumberOfValues();++i) {
         double x, val; ifile.scanField( valuestr, val ); 
         for(unsigned j=0; j<dernames.size(); ++j) {
           ifile.scanField(dernames[j],x); 
           if( !flatgrid ) {
              ifile.scanField("nbins", gbin1);
           } else {
              xx[j]=x+gridobject.getGridSpacing()[j]/2.0;
              ifile.scanField( "min_" + dernames[j], gmin[j]);
              ifile.scanField( "max_" + dernames[j], gmax[j]);
              ifile.scanField( "nbins_" + dernames[j], gbin1);
              ifile.scanField( "periodic_" + dernames[j], pstring );
           }
         }
         for(unsigned j=0; j<dernames.size(); ++j) ifile.scanField( "d" + valuestr + "_" + dernames[j], dder[j] );  

         unsigned index=gridobject.getIndex(xx); if( !flatgrid ) index=i;
         valout->add( index*(dernames.size()+1), val );
         for(unsigned j=0; j<dernames.size(); ++j) valout->add( index*(dernames.size()+1)+1+j, dder[j] );  
         ifile.scanField();
       }
       ifile.close();
   }
}

void ReadGridInSetup::createGridAndValue( const std::string& gtype, const std::vector<bool>& ipbc, const unsigned& nfermi,
                                          const std::vector<std::string>& gmin, const std::vector<std::string>& gmax,
                                          const std::vector<unsigned>& gbin ) {
  gridobject.setup( gtype, ipbc, nfermi, 0.0 ); std::vector<double> gspacing;
  if( gtype=="flat" ) {
      gridobject.setBounds( gmin, gmax, gbin, gspacing );
      // Now create the value
      std::vector<unsigned> shape( gridobject.getNbin(true) );
      ActionWithValue::addValueWithDerivatives( shape ); setNotPeriodic();
  } else {
      std::vector<unsigned> shape( 3 ); shape[0]=gbin[0]; shape[1]=shape[2]=1;
      ActionWithValue::addValueWithDerivatives( shape ); setNotPeriodic();
  }
  for(unsigned i=0;i<getNumberOfComponents();++i) {
      getPntrToComponent(i)->setConstant();
      getPntrToComponent(i)->setDerivativeIsZeroWhenValueIsZero();
  }
  // This ensures we set the flag to never active the action.  We can say we have atoms here as we don't need them 
  // to calculate the CV
  setupConstantValues( true );
}

unsigned ReadGridInSetup::getNumberOfDerivatives() {
   return dernames.size();
}

std::vector<std::string> ReadGridInSetup::getGridCoordinateNames() const {
  return dernames;
}             
          
const GridCoordinatesObject& ReadGridInSetup::getGridCoordinatesObject() const {
  return gridobject;
}   

}
}

