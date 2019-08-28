/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2019 The plumed team
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
#include "core/ActionRegister.h"
#include "SketchMapBase.h"
#include "reference/ReferenceConfiguration.h"
#include "reference/MetricRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/PDB.h"

//+PLUMEDOC DIMRED SKETCHMAP_READ
/*
Read in a sketch-map projection from an input file

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace dimred {

class SketchMapRead : public SketchMapBase {
private:
  std::string mtype;
  PDB mypdb;
  std::vector<Value*> all_values;
  std::vector<double> weights;
/// The list of properties in the property map
  std::map<std::string,std::vector<double> > property;
/// The data collection objects we have
  std::vector<analysis::DataCollectionObject> data;
/// The frames that we are using to calculate distances
  std::vector<std::unique_ptr<ReferenceConfiguration> > myframes;
public:
  static void registerKeywords( Keywords& keys );
  explicit SketchMapRead( const ActionOptions& ao );
  void minimise( Matrix<double>& ) override;
  analysis::DataCollectionObject& getStoredData( const unsigned& idata, const bool& calcdist ) override;
  unsigned getNumberOfDataPoints() const override;
  std::vector<Value*> getArgumentList() override;
  unsigned getDataPointIndexInBase( const unsigned& idata ) const override;
  double getDissimilarity( const unsigned& i, const unsigned& j ) override;
  double getWeight( const unsigned& idata ) override;
};

PLUMED_REGISTER_ACTION(SketchMapRead,"SKETCHMAP_READ")

void SketchMapRead::registerKeywords( Keywords& keys ) {
  SketchMapBase::registerKeywords( keys ); keys.remove("USE_OUTPUT_DATA_FROM");
  keys.add("compulsory","TYPE","OPTIMAL-FAST","the manner in which distances are calculated. More information on the different "
           "metrics that are available in PLUMED can be found in the section of the manual on "
           "\\ref dists");
  keys.add("compulsory","REFERENCE","the file containing the sketch-map projection");
  keys.add("compulsory","PROPERTY","the property to be used in the index. This should be in the REMARK of the reference");
  keys.addFlag("DISABLE_CHECKS",false,"disable checks on reference input structures.");
  useCustomisableComponents(keys);
}

SketchMapRead::SketchMapRead( const ActionOptions& ao ):
  Action(ao),
  SketchMapBase(ao)
{
  std::vector<std::string> propnames; parseVector("PROPERTY",propnames);
  if(propnames.size()==0) error("no properties were specified");
  nlow=propnames.size();
  for(unsigned i=0; i<nlow; ++i) {
    std::size_t dot=propnames[i].find_first_of( getLabel() + "." ); std::string substr=propnames[i].c_str();
    if( dot!=std::string::npos ) { substr.erase(dot,getLabel().length()+1); }
    log.printf(",%s", propnames[i].c_str() ); addComponent( substr ); componentIsNotPeriodic( substr );
    property.insert( std::pair<std::string,std::vector<double> >( propnames[i], std::vector<double>() ) );
  }
  log.printf("  mapped properties are %s ",propnames[0].c_str() );
  for(unsigned i=1; i<nlow; ++i) log.printf(",%s", propnames[i].c_str() );
  log.printf("\n");

  parse("TYPE",mtype); bool skipchecks; parseFlag("DISABLE_CHECKS",skipchecks);
  std::string ifilename; parse("REFERENCE",ifilename);
  FILE* fp=fopen(ifilename.c_str(),"r");
  if(!fp) error("could not open reference file " + ifilename );

  // Read in the embedding
  bool do_read=true; double val, ww, wnorm=0, prop; unsigned nfram=0;
  while (do_read) {
    PDB inpdb;
    // Read the pdb file
    do_read=inpdb.readFromFilepointer(fp,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength());
    // Break if we are done
    if( !do_read ) break ;
    // Check for required properties
    for(std::map<std::string,std::vector<double> >::iterator it=property.begin(); it!=property.end(); ++it) {
      if( !inpdb.getArgumentValue( it->first, prop ) ) error("pdb input does not have contain property named " + it->first );
      it->second.push_back(prop);
    }
    // And read the frame ( create a measure )
    myframes.emplace_back( metricRegister().create<ReferenceConfiguration>( mtype, inpdb ) );
    if( !inpdb.getArgumentValue( "WEIGHT", ww ) ) error("could not find weights in input pdb");
    weights.push_back( ww ); wnorm += ww; nfram++;
    // And create a data collection object
    analysis::DataCollectionObject new_data; new_data.setAtomNumbersAndArgumentNames( getLabel(), inpdb.getAtomNumbers(), inpdb.getArgumentNames() );
    new_data.setAtomPositions( inpdb.getPositions() );
    for(unsigned i=0; i<inpdb.getArgumentNames().size(); ++i) {
      std::string aname = inpdb.getArgumentNames()[i];
      if( !inpdb.getArgumentValue( aname, val ) ) error("failed to find argument named " + aname);
      new_data.setArgument( aname, val );
    }
    data.push_back( new_data );
  }
  fclose(fp);
  // Finish the setup of the object by getting the arguments and atoms that are required
  std::vector<AtomNumber> atoms; std::vector<std::string> args;
  for(unsigned i=0; i<myframes.size(); ++i) { weights[i] /= wnorm; myframes[i]->getAtomRequests( atoms, skipchecks ); myframes[i]->getArgumentRequests( args, skipchecks ); }
  requestAtoms( atoms ); std::vector<Value*> req_args; std::vector<std::string> fargs;
  for(unsigned i=0; i<args.size(); ++i) {
    bool found=false;
    for(std::map<std::string,std::vector<double> >::iterator it=property.begin(); it!=property.end(); ++it) {
      if( args[i]==it->first ) { found=true; break; }
    }
    if( !found ) { fargs.push_back( args[i] ); }
  }
  interpretArgumentList( fargs, req_args ); mypdb.setArgumentNames( fargs ); requestArguments( req_args );

  if(nfram==0 ) error("no reference configurations were specified");
  log.printf(" found %u configurations in file %s\n",nfram,ifilename.c_str() );
}

void SketchMapRead::minimise( Matrix<double>& projections ) {
  unsigned j=0;
  for(std::map<std::string,std::vector<double> >::iterator it=property.begin(); it!=property.end(); ++it) {
    for(unsigned i=0; i<myframes.size(); ++i) projections(i,j) = it->second[i];
    j++;
  }
}

analysis::DataCollectionObject & SketchMapRead::getStoredData( const unsigned& idata, const bool& calcdist ) {
  return data[idata];
}

unsigned SketchMapRead::getNumberOfDataPoints() const {
  return myframes.size();
}

unsigned SketchMapRead::getDataPointIndexInBase( const unsigned& idata ) const {
  error("cannot use read in sketch-map to out of sample project data");
  return idata;
}

std::vector<Value*> SketchMapRead::getArgumentList() {
  std::vector<Value*> arglist( ActionWithArguments::getArguments() );
  for(unsigned i=0; i<nlow; ++i) arglist.push_back( getPntrToComponent(i) );
  return arglist;
}

// Highly unsatisfactory solution to problem GAT
double SketchMapRead::getDissimilarity( const unsigned& i, const unsigned& j ) {
  plumed_assert( i<myframes.size() && j<myframes.size() );
  if( i!=j ) {
    double dd;
    getStoredData( i, true ).transferDataToPDB( mypdb );
    auto myref1=metricRegister().create<ReferenceConfiguration>(mtype, mypdb);
    getStoredData( j, true ).transferDataToPDB( mypdb );
    auto myref2=metricRegister().create<ReferenceConfiguration>(mtype, mypdb);
    dd=distance( getPbc(), ActionWithArguments::getArguments(), myref1.get(), myref2.get(), true );
    return dd;
  }
  return 0.0;
}

double SketchMapRead::getWeight( const unsigned& idata ) {
  plumed_assert( idata<weights.size() );
  return weights[idata];
}

}
}
