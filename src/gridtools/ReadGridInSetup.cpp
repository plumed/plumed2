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
#include "core/ActionSetup.h"
#include "core/ActionWithValue.h"
#include "core/ActionRegister.h"
#include "tools/IFile.h"
#include "GridCoordinatesObject.h"

using namespace std;

namespace PLMD {
namespace gridtools {

class ReadGridInSetup :
  public virtual ActionSetup,
  public virtual ActionWithValue
{
private:
  std::vector<std::string> labels;
  GridCoordinatesObject gridobject;
public:
  static void registerKeywords( Keywords& keys );
  explicit ReadGridInSetup(const ActionOptions&ao);
  unsigned getNumberOfDerivatives() const ;
  void clearDerivatives( const bool& force=false ) {}
  void getInfoForGridHeader( std::string& gtype, std::vector<std::string>& argn, std::vector<std::string>& min,
                             std::vector<std::string>& max, std::vector<unsigned>& out_nbin,
                             std::vector<double>& spacing, std::vector<bool>& pbc, const bool& dumpcube ) const ;
  void getGridPointIndicesAndCoordinates( const unsigned& ind, std::vector<unsigned>& indices, std::vector<double>& coords ) const ; 
  void getGridPointAsCoordinate( const unsigned& ind, const bool& setlength, std::vector<double>& coords ) const ;
};

PLUMED_REGISTER_ACTION(ReadGridInSetup,"REFERENCE_GRID")

void ReadGridInSetup::registerKeywords( Keywords& keys ) {
  ActionSetup::registerKeywords(keys); ActionWithValue::registerKeywords(keys);
  keys.remove("NUMERICAL_DERIVATIVES"); keys.remove("SERIAL"); keys.remove("TIMINGS"); keys.add("hidden","LABEL","");
  keys.add("compulsory","FILE","the name of the file that contains the reference data");
  keys.add("compulsory","VALUE","the name of the value that should be read from the grid");
}

ReadGridInSetup::ReadGridInSetup(const ActionOptions&ao):
  Action(ao),
  ActionSetup(ao),
  ActionWithValue(ao)
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
  for(unsigned i=0;i<fieldnames.size();++i) {
      std::size_t dot = fieldnames[i].find_first_of("min_");
      if( fieldnames[i].find("min_")!=std::string::npos ) labels.push_back( fieldnames[i].substr(dot+4) );
  }
  // Now get all the header data for the grid
  std::vector<std::string> gmin( labels.size() ), gmax( labels.size() ); std::string pstring;
  int gbin1; std::vector<unsigned> gbin( labels.size() ); std::vector<bool> ipbc( labels.size() );
  for(unsigned i=0; i<labels.size(); ++i) {
    ifile.scanField( "min_" + labels[i], gmin[i]);
    ifile.scanField( "max_" + labels[i], gmax[i]);
    ifile.scanField( "periodic_" + labels[i], pstring );
    ifile.scanField( "nbins_" + labels[i], gbin1); gbin[i]=gbin1;
    if( pstring=="true" ) ipbc[i]=true; 
    else if( pstring=="false" ) ipbc[i]=false; 
    else error("do not understand periodidicy of " + labels[i] );

    bool hasder=ifile.FieldExist( "d" + valuestr + "_" + labels[i] );
    if( !hasder ) plumed_merror("missing derivatives from grid file");
  }
  gridobject.setup( "flat", ipbc, 0, 0.0 ); std::vector<double> gspacing;
  gridobject.setBounds( gmin, gmax, gbin, gspacing );
  // Now create the value
  std::vector<unsigned> shape( gridobject.getNbin(true) );
  ActionWithValue::addValueWithDerivatives( shape ); setNotPeriodic();
  // And finally read all the grid points
  Value* valout=getPntrToOutput(0); 
  std::vector<double> dder( labels.size() ), xx( labels.size() );
  for(unsigned i=0;i<valout->getNumberOfValues( getLabel() );++i) {
    double x, val, norm; ifile.scanField( valuestr, val ); ifile.scanField( "normalisation", norm );
    for(unsigned j=0; j<labels.size(); ++j) {
      ifile.scanField(labels[j],x); xx[j]=x+gridobject.getGridSpacing()[j]/2.0;
      ifile.scanField( "min_" + labels[j], gmin[j]);
      ifile.scanField( "max_" + labels[j], gmax[j]);
      ifile.scanField( "nbins_" + labels[j], gbin1);
      ifile.scanField( "periodic_" + labels[j], pstring );
    }
    for(unsigned j=0; j<labels.size(); ++j) ifile.scanField( "d" + valuestr + "_" + labels[j], dder[j] );  

    unsigned index=gridobject.getIndex(xx);
    valout->setNorm( norm ); valout->add( index*(labels.size()+1), norm*val );
    for(unsigned j=0; j<labels.size(); ++j) valout->add( index*(labels.size()+1)+1+j, norm*dder[j] );  
    ifile.scanField();
  }
  ifile.close();
}

unsigned ReadGridInSetup::getNumberOfDerivatives() const {
  return labels.size();
}

void ReadGridInSetup::getInfoForGridHeader( std::string& gtype, std::vector<std::string>& argn, std::vector<std::string>& min,
                                            std::vector<std::string>& max, std::vector<unsigned>& out_nbin,
                                            std::vector<double>& spacing, std::vector<bool>& pbc, const bool& dumpcube ) const {
  plumed_assert( !dumpcube ); gtype="flat"; std::vector<unsigned> nbin( gridobject.getNbin(false) );
  for(unsigned i=0; i<getPntrToOutput(0)->getRank(); ++i) {
    argn[i]=labels[i]; double gmin, gmax;
    if( gridobject.getMin().size()>0 ) {
      Tools::convert( gridobject.getMin()[i], gmin ); Tools::convert( gmin, min[i] );
      Tools::convert( gridobject.getMax()[i], gmax ); Tools::convert( gmax, max[i] );
    }
    if( nbin.size()>0 ) out_nbin[i]=nbin[i];
    if( spacing.size()>0 ) spacing[i]=gridobject.getGridSpacing()[i];
    pbc[i]=gridobject.isPeriodic(i);
  }
}

void ReadGridInSetup::getGridPointIndicesAndCoordinates( const unsigned& ind, std::vector<unsigned>& indices, std::vector<double>& coords ) const {
  gridobject.getGridPointCoordinates( ind, indices, coords );
}

void ReadGridInSetup::getGridPointAsCoordinate( const unsigned& ind, const bool& setlength, std::vector<double>& coords ) const {
  if( setlength ) gridobject.putCoordinateAtValue( ind, getPntrToOutput(0)->get(ind), coords );
  else  gridobject.putCoordinateAtValue( ind, 1.0, coords );
}

}
}

