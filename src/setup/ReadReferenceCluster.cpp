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
#include "ReadReferenceCluster.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "tools/IFile.h"

namespace PLMD {
namespace setup {

PLUMED_REGISTER_ACTION(ReadReferenceCluster,"READ_VECTOR")
PLUMED_REGISTER_ACTION(ReadReferenceCluster,"READ_CLUSTER")

void ReadReferenceCluster::registerKeywords( Keywords& keys ) {
  SetupReferenceBase::registerKeywords( keys ); 
  keys.add("compulsory","CENTER","the position of the center of the cluster");
  keys.add("optional","SIGMA","square root of variance of the cluster");
  keys.add("compulsory","COVAR","the covariance of the cluster");
  keys.add("optional","REFERENCE","a file containing information on the reference cluster.");
  keys.add("compulsory","NUMBER","1","if there are multiple clusters in the input file which structure would you like to read in here");
  keys.add("hidden","READ_ARG","this is used by pathtool to get the arguments that must be read in");
  keys.addOutputComponent("center","default","the position of the center of the cluster in CV space");
  keys.addOutputComponent("variance","SIGMA","the vector of variances for the CVs that describes the extent of the cluster");
  keys.addOutputComponent("covariance","COVAR","the covariance matrix for the CVs that describes the extent of the cluster");
}

ReadReferenceCluster::ReadReferenceCluster(const ActionOptions&ao):
Action(ao),
SetupReferenceBase(ao)
{
  if( getNumberOfArguments()==0 ) parseVector("READ_ARG",read_args);
  std::string reference; parse("REFERENCE",reference); 
  if( reference.length()==0 ) {
      // Read in the position of the center
      std::vector<double> center; parseVector("CENTER",center); std::vector<unsigned> shape( 1 ); shape[0] = getNumberOfArguments();
      if( getNumberOfArguments()==0 && read_args.size()==0 ) { 
          shape[0]=center.size();
      } else if( getNumberOfArguments()==0 ) { 
          shape[0]=read_args.size(); 
      } else if( !numberedkeys ) {
          shape[0]=0; for(unsigned i=0;i<getNumberOfArguments();++i) shape[0] += getPntrToArgument(i)->getNumberOfValues();
      }
      log.printf("  read in center of cluster in space of dimension %d \n", shape[0] );
      if( getName()=="READ_VECTOR" ) {
          addValue( shape ); setNotPeriodic(); setCenterFromVector( center );
      } else { 
          addComponent( "center", shape ); componentIsNotPeriodic("center"); setCenterFromVector( center );
          // Read in the covariance
          std::vector<double> sigma; parseVector("SIGMA",sigma);
          if( sigma.size()==shape[0] ) {
              addComponent("variance", shape); componentIsNotPeriodic("variance"); getPntrToComponent(1)->alwaysStoreValues();
              for(unsigned i=0;i<shape[0];++i) getPntrToComponent(1)->set( i, sigma[i]*sigma[i] );
          } else if( sigma.size()==0 ) {
              std::vector<double> covar; parseVector("COVAR",covar); 
              if( covar.size()==shape[0] ) {
                  addComponent("variance", shape); componentIsNotPeriodic("variance"); getPntrToComponent(1)->alwaysStoreValues();
                  for(unsigned i=0;i<shape[0];++i) getPntrToComponent(1)->set( i, covar[i] );
              } else if( covar.size()==shape[0]*shape[0] ) {
                  shape.push_back( shape[0] ); addComponent("covariance", shape); componentIsNotPeriodic("covariance"); getPntrToComponent(1)->alwaysStoreValues();
                  for(unsigned i=0;i<covar.size();++i) getPntrToComponent(1)->set( i, covar[i] ); 
              } else error("covariance has the wrong shape");
          } else error("sigma has the wrong shape");
      }
  } else {
      unsigned number; parse("NUMBER",number); std::vector<std::string> names; 
      for(unsigned i=0;i<getNumberOfArguments();++i) {
          if( numberedkeys ) {
              names.push_back( getPntrToArgument(i)->getName() );
          } else if( getPntrToArgument(i)->getRank()==0 ) {
              names.push_back( getPntrToArgument(i)->getName() );
          } else if( getPntrToArgument(i)->getRank()==1 ) {
              for(unsigned j=0;j<getPntrToArgument(i)->getShape()[0];++j) {
                  std::string num; Tools::convert( j+1, num );
                  names.push_back( getPntrToArgument(i)->getName() + "." + num );
              }
          } else if( getPntrToArgument(i)->getRank()==2 ) {
              for(unsigned j=0;j<getPntrToArgument(i)->getShape()[0];++j) {
                  std::string jnum; Tools::convert( j+1, jnum );
                  for(unsigned k=0;k<getPntrToArgument(k)->getShape()[2];++k) {
                      std::string knum; Tools::convert( k+1, knum );
                      names.push_back( getPntrToArgument(i)->getName() + "." + jnum + "." + knum );
                  }
              }
          } else error("cannot deal with objects with ranks greater than 2");
      }
      log.printf("  reading %dth reference structure from file %s \n", number, reference.c_str());
      log.printf("  which contains %d arguments \n", names.size() );
      log.printf("  labels of arguments are : ");
      for(unsigned i=0;i<names.size();++i) log.printf("%s ", names[i].c_str() );
      log.printf("\n"); std::string input; 
      if( getNumberOfArguments()==0 ) {
          input = "READ_ARGS=" + read_args[0]; for(unsigned i=1;i<read_args.size();++i) input += "," + read_args[i];
          input += " " + convertFileToLine( reference, number, read_args );
      } else {
          input = "ARG=" + getPntrToArgument(0)->getName(); for(unsigned i=1;i<getNumberOfArguments();++i) input += "," + getPntrToArgument(i)->getName();
          input += " " + convertFileToLine( reference, number, names );
      }
      std::string slab = getLabel(); plumed.readInputLine( slab + ": READ_CLUSTER " + input );
  }
}

void ReadReferenceCluster::setCenterFromVector( const std::vector<double>& center ) {
  getPntrToComponent(0)->alwaysStoreValues();
  std::vector<unsigned> shape( getPntrToComponent(0)->getShape() );
  if( center.size()!=shape[0] ) error("size of center does not match number of arguments");
  for(unsigned i=0;i<shape[0];++i) getPntrToComponent(0)->set( i, center[i] );
}

std::string ReadReferenceCluster::convertFileToLine( const std::string& reference, const unsigned& number, const std::vector<std::string>& names ) {
  bool readline=false; IFile ifile; ifile.open(reference); ifile.allowIgnoredFields(); std::string input;
  for(unsigned line=0;line<number;++line) {
    // Read in the position of the center of the cluster
    std::string val; input="CENTER=";
    for(unsigned i=0;i<names.size();++i) { ifile.scanField(names[i], val); if(i==0) input += val; else input += "," + val; }
    if( ifile.FieldExist("sigma_" + names[0]) ) {
        input += " COVAR=";
        for(unsigned i=0;i<names.size();++i) { ifile.scanField("sigma_" + names[i], val); if( i==0 ) input += val; else input += "," + val; } 
    } else {
        input += " COVAR="; 
        for(unsigned i=0; i<names.size(); ++i) {
            for(unsigned j=0; j<names.size(); j++) {
                ifile.scanField("sigma_" +names[i] + "_" + names[j], val ); 
                if(i==0 && j==0 ) input += val; else input += "," + val;
            }
        }
    }
    if( line==number-1 ) { readline=true; break; }
    ifile.scanField(); 
  }
  if( !readline ) plumed_merror("could not read reference configuration");
  ifile.scanField(); ifile.close();
  return input;
}

std::string ReadReferenceCluster::getArgName( const unsigned& k ) const {
  if( read_args.size()>0 ) return read_args[k];
  return SetupReferenceBase::getArgName(k);
}

}
}
