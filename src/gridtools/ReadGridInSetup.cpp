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
#include "lepton/Lepton.h"
#include "tools/IFile.h"

//+PLUMEDOC GRIDCALC REFERENCE_GRID
/*
Setup a constant grid by either reading values from a file or definining a function in input

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace gridtools {

static std::map<std::string, double> leptonConstants= {
  {"e", std::exp(1.0)},
  {"log2e", 1.0/std::log(2.0)},
  {"log10e", 1.0/std::log(10.0)},
  {"ln2", std::log(2.0)},
  {"ln10", std::log(10.0)},
  {"pi", pi},
  {"pi_2", pi*0.5},
  {"pi_4", pi*0.25},
//  {"1_pi", 1.0/pi},
//  {"2_pi", 2.0/pi},
//  {"2_sqrtpi", 2.0/std::sqrt(pi)},
  {"sqrt2", std::sqrt(2.0)},
  {"sqrt1_2", std::sqrt(0.5)}
};

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
  void setupOnFirstStep( const bool incalc ) override {
    plumed_merror("should not be in ReadGridInSetup setupOnFirstStep" );
  }
  void performTask( const unsigned& current, MultiValue& myvals ) const override {
    plumed_merror("should not be in readGridInSetup performTask");
  }
  std::vector<std::string> getGridCoordinateNames() const override ;
  const GridCoordinatesObject& getGridCoordinatesObject() const override ;
};

PLUMED_REGISTER_ACTION(ReadGridInSetup,"REFERENCE_GRID")

void ReadGridInSetup::registerKeywords( Keywords& keys ) {
  ActionWithGrid::registerKeywords(keys);
  keys.remove("SERIAL");
  keys.add("optional","FUNC","the function to compute on the grid");
  keys.add("compulsory","GRID_MIN","auto","the lower bounds for the grid");
  keys.add("compulsory","GRID_MAX","auto","the upper bounds for the grid");
  keys.add("compulsory","PERIODIC","are the grid directions periodic");
  keys.add("optional","GRID_BIN","the number of bins for the grid");
  keys.add("optional","GRID_SPACING","the approximate grid spacing (to be used as an alternative or together with GRID_BIN)");
  keys.add("optional","VAR","the names to give each of the grid directions in the function.  If you have up to three grid coordinates in your function you can use x, y and z to refer to them.  Otherwise you must use this flag to give your variables names.");
  keys.add("compulsory","FILE","the name of the file that contains the reference data");
  keys.add("compulsory","VALUE","the name of the value that should be read from the grid");
  keys.setValueDescription("the constant function on the grid that was specified in input");
}

ReadGridInSetup::ReadGridInSetup(const ActionOptions&ao):
  Action(ao),
  ActionWithGrid(ao) {
  std::string func;
  parse("FUNC",func);
  if( func.length()>0 ) {
    // Read in stuff for grid
    std::vector<std::string> gmin;
    parseVector("GRID_MIN",gmin);
    std::vector<std::string> gmax(gmin.size());
    parseVector("GRID_MAX",gmax);
    std::vector<unsigned> gbin(gmin.size());
    parseVector("GRID_BIN",gbin);
    std::vector<std::string> pbc(gmin.size());
    parseVector("PERIODIC",pbc);
    std::vector<bool> ipbc( pbc.size() );
    for(unsigned i=0; i<ipbc.size(); ++i) {
      if( pbc[i]=="YES" ) {
        ipbc[i]=true;
      } else if( pbc[i]=="NO" ) {
        ipbc[i]=false;
      } else {
        error( pbc[i] + " is not a valid instruction to the PERIODIC keyword");
      }
    }

    // Read in the variables
    parseVector("VAR",dernames);
    if(dernames.size()==0) {
      dernames.resize(gmin.size());
      if(gmin.size()>3) {
        error("Using more than 3 arguments you should explicitly write their names with VAR");
      }
      if(dernames.size()>0) {
        dernames[0]="x";
      }
      if(dernames.size()>1) {
        dernames[1]="y";
      }
      if(dernames.size()>2) {
        dernames[2]="z";
      }
    }
    if(dernames.size()!=gmin.size()) {
      error("Size of VAR array should be the same as number of grid dimensions");
    }

    // Create the grid and the value of the grid
    createGridAndValue( "flat", ipbc, 0, gmin, gmax, gbin );

    // Read in stuff for function
    log.printf("  evaluating function : %s\n",func.c_str());
    log.printf("  with variables :");
    for(unsigned i=0; i<dernames.size(); i++) {
      log.printf(" %s",dernames[i].c_str());
    }
    log.printf("\n");
    log.printf("  on %d", gbin[0]);
    for(unsigned i=1; i<gbin.size(); ++i) {
      log.printf(" by %d \n", gbin[i]);
    }
    log.printf(" grid of points between (%s", gmin[0].c_str() );
    for(unsigned i=1; i<gmin.size(); ++i) {
      log.printf(", %s", gmin[i].c_str() );
    }
    log.printf(") and (%s", gmax[0].c_str() );
    for(unsigned i=1; i<gmax.size(); ++i) {
      log.printf(", %s", gmax[i].c_str() );
    }
    log.printf(")\n");

    lepton::ParsedExpression pe=lepton::Parser::parse(func).optimize(leptonConstants);
    log<<"  function as parsed by lepton: "<<pe<<"\n";
    lepton::CompiledExpression expression=pe.createCompiledExpression();
    for(auto &p: expression.getVariables()) {
      if(std::find(dernames.begin(),dernames.end(),p)==dernames.end()) {
        error("variable " + p + " is not defined");
      }
    }
    log<<"  derivatives as computed by lepton:\n";
    std::vector<lepton::CompiledExpression> expression_deriv( dernames.size() );
    for(unsigned i=0; i<dernames.size(); i++) {
      lepton::ParsedExpression pe=lepton::Parser::parse(func).differentiate(dernames[i]).optimize(leptonConstants);
      log<<"    "<<pe<<"\n";
      expression_deriv[i]=pe.createCompiledExpression();
    }
    // And finally calculate all the grid points
    std::vector<double> dder( dernames.size() ), xx( dernames.size() );
    Value* valout=getPntrToComponent(0);
    for(unsigned index=0; index<valout->getNumberOfValues(); ++index) {
      gridobject.getGridPointCoordinates( index, xx );
      for(unsigned j=0; j<xx.size(); ++j) {
        try {
          expression.getVariableReference(dernames[j])=xx[j];
        } catch(PLMD::lepton::Exception& exc) {
// this is necessary since in some cases lepton things a variable is not present even though it is present
// e.g. func=0*x
        }
      }
      valout->set( index, expression.evaluate() );
      for(unsigned k=0; k<xx.size(); ++k) {
        for(unsigned j=0; j<xx.size(); ++j) {
          try {
            expression_deriv[k].getVariableReference(dernames[j])=xx[j];
          } catch(PLMD::lepton::Exception& exc) {
// this is necessary since in some cases lepton things a variable is not present even though it is present
// e.g. func=0*x
          }
        }
        valout->addGridDerivatives( index, k, expression_deriv[k].evaluate() );
      }
    }
  } else {
    std::string valuestr;
    parse("VALUE",valuestr);
    std::string tstyle, filen;
    parse("FILE",filen);
    if( filen.length()>0 ) {
      std::size_t dot=filen.find_first_of(".");
      if( dot!=std::string::npos ) {
        tstyle=filen.substr(dot+1);
      }
      if( tstyle!="grid" ) {
        error("can only read in grids using read value in setup");
      }
      log.printf("  reading function %s on grid from file %s \n", valuestr.c_str(), filen.c_str() );
    }
    IFile ifile;
    ifile.open(filen);
    if( !ifile.FieldExist( valuestr ) ) {
      error("could not find grid value in input file");
    }
    std::vector<std::string> fieldnames;
    ifile.scanFieldList( fieldnames );

    // Retrieve the names of the variables the grid is computed over
    bool flatgrid=false;
    for(unsigned i=0; i<fieldnames.size(); ++i) {
      if( fieldnames[i].find("min_")!=std::string::npos ) {
        flatgrid=true;
      }
      std::size_t dot = fieldnames[i].find_first_of("d" + valuestr + "_" );
      if( fieldnames[i].find("d" + valuestr + "_")!=std::string::npos ) {
        dernames.push_back( fieldnames[i].substr(dot+2+valuestr.length()) );
      }
    }
    if( flatgrid && dernames.size()==0 ) {
      error("could not find any derivatives for value " + valuestr + " in input file.  Header should contain at least columns with a name starting d" + valuestr + "_");
    }
    // Now get all the header data for the grid
    std::vector<std::string> gmin( dernames.size() ), gmax( dernames.size() );
    std::string pstring;
    int gbin1;
    std::vector<unsigned> gbin( dernames.size() );
    std::vector<bool> ipbc( dernames.size() );
    if( !flatgrid ) {
      ifile.scanField( "nbins", gbin1);
      gbin[0]=gbin1;
      createGridAndValue( "fibonacci", ipbc, gbin[0], gmin, gmax, gbin );
    } else {
      for(unsigned i=0; i<dernames.size(); ++i) {
        ifile.scanField( "min_" + dernames[i], gmin[i]);
        ifile.scanField( "max_" + dernames[i], gmax[i]);
        ifile.scanField( "periodic_" + dernames[i], pstring );
        ifile.scanField( "nbins_" + dernames[i], gbin1);
        gbin[i]=gbin1;
        if( pstring=="true" ) {
          log.printf("   for periodic coordinate %s minimum is %s maximum is %s and number of bins is %d \n",dernames[i].c_str(),gmin[i].c_str(),gmax[i].c_str(),gbin[i]);
          ipbc[i]=true;
        } else if( pstring=="false" ) {
          log.printf("   for coordinate %s minimum is %s maximum is %s and number of bins is %d \n",dernames[i].c_str(),gmin[i].c_str(),gmax[i].c_str(),gbin[i]);
          ipbc[i]=false;
        } else {
          error("do not understand periodidicy of " + dernames[i] );
        }

        bool hasder=ifile.FieldExist( "d" + valuestr + "_" + dernames[i] );
        if( !hasder ) {
          plumed_merror("missing derivatives from grid file");
        }
      }
      createGridAndValue( "flat", ipbc, 0, gmin, gmax, gbin );
    }
    // And finally read all the grid points
    Value* valout=getPntrToComponent(0);
    std::vector<double> dder( dernames.size() ), xx( dernames.size() );
    for(unsigned i=0; i<valout->getNumberOfValues(); ++i) {
      double x, val;
      ifile.scanField( valuestr, val );
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
      for(unsigned j=0; j<dernames.size(); ++j) {
        ifile.scanField( "d" + valuestr + "_" + dernames[j], dder[j] );
      }

      unsigned index=gridobject.getIndex(xx);
      if( !flatgrid ) {
        index=i;
      }
      valout->set( index, val );
      for(unsigned j=0; j<dernames.size(); ++j) {
        valout->addGridDerivatives( index, j, dder[j] );
      }
      ifile.scanField();
    }
    ifile.close();
  }
}

void ReadGridInSetup::createGridAndValue( const std::string& gtype, const std::vector<bool>& ipbc, const unsigned& nfermi,
    const std::vector<std::string>& gmin, const std::vector<std::string>& gmax,
    const std::vector<unsigned>& gbin ) {
  gridobject.setup( gtype, ipbc, nfermi, 0.0 );
  std::vector<double> gspacing;
  if( gtype=="flat" ) {
    gridobject.setBounds( gmin, gmax, gbin, gspacing );
    // Now create the value
    std::vector<unsigned> shape( gridobject.getNbin(true) );
    ActionWithValue::addValueWithDerivatives( shape );
    setNotPeriodic();
  } else {
    std::vector<unsigned> shape( 3 );
    shape[0]=gbin[0];
    shape[1]=shape[2]=1;
    ActionWithValue::addValueWithDerivatives( shape );
    setNotPeriodic();
  }
  for(unsigned i=0; i<getNumberOfComponents(); ++i) {
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

