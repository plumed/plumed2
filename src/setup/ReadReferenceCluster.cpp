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
#include "SetupReferenceBase.h"
#include "core/ActionRegister.h"
#include "tools/IFile.h"

namespace PLMD {
namespace setup {

class ReadReferenceCluster: public SetupReferenceBase {
private:
  std::vector<std::string> read_args;
public: 
  static void registerKeywords( Keywords& keys );
  explicit ReadReferenceCluster(const ActionOptions&ao);
  std::string getArgName( const unsigned& k ) const ;
};

PLUMED_REGISTER_ACTION(ReadReferenceCluster,"READ_CLUSTER")

void ReadReferenceCluster::registerKeywords( Keywords& keys ) {
  SetupReferenceBase::registerKeywords( keys ); 
  keys.add("compulsory","CENTER","the position of the center of the cluster");
  keys.add("optional","SIGMA","square root of variance of the cluster");
  keys.add("compulsory","COVAR","the covariance of the cluster");
  keys.add("optional","REFERENCE","a file containing information on the reference cluster.");
  keys.add("compulsory","NUMBER","1","if there are multiple clusters in the input file which structure would you like to read in here");
  keys.add("hidden","READ_ARG","this is used by pathtool to get the arguments that must be read in");
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
      if( !numberedkeys ) {
          shape[0]=0; for(unsigned i=0;i<getNumberOfArguments();++i) shape[0] += getPntrToArgument(i)->getNumberOfValues( getLabel() );
      }
      addComponent( "center", shape ); componentIsNotPeriodic("center"); getPntrToComponent(0)->buildDataStore( getLabel() );
      if( center.size()!=shape[0] ) error("size of center does not match number of arguments");
      for(unsigned i=0;i<shape[0];++i) getPntrToComponent(0)->set( i, center[i] );
      // Read in the covariance
      std::vector<double> sigma; parseVector("SIGMA",sigma);
      if( sigma.size()==shape[0] ) {
          addComponent("variance", shape); componentIsNotPeriodic("variance"); getPntrToComponent(1)->buildDataStore( getLabel() );
          for(unsigned i=0;i<shape[0];++i) getPntrToComponent(1)->set( i, sigma[i]*sigma[i] );
      } else if( sigma.size()==0 ) {
          std::vector<double> covar; parseVector("COVAR",covar); 
          if( covar.size()==shape[0] ) {
              addComponent("variance", shape); componentIsNotPeriodic("variance"); getPntrToComponent(1)->buildDataStore( getLabel() );
              for(unsigned i=0;i<shape[0];++i) getPntrToComponent(1)->set( i, covar[i] );
          } else if( covar.size()==shape[0]*shape[0] ) {
              shape.push_back( shape[0] ); addComponent("covariance", shape); componentIsNotPeriodic("covariance"); getPntrToComponent(1)->buildDataStore( getLabel() );
              for(unsigned i=0;i<covar.size();++i) getPntrToComponent(1)->set( i, covar[i] ); 
          } else error("covariance has the wrong shape");
      } else error("sigma has the wrong shape");
  } else {
      IFile ifile; ifile.open(reference); ifile.allowIgnoredFields();
      unsigned number; parse("NUMBER",number); bool readline=false;
      for(unsigned line=0;line<number;++line) {
        // Read in the position of the center of the cluster
        std::vector<std::string> names; std::vector<double> values; double val;
        for(unsigned i=0;i<getNumberOfArguments();++i) {
            if( numberedkeys ) {
                ifile.scanField(getPntrToArgument(i)->getName(), val); 
                values.push_back(val); names.push_back( getPntrToArgument(i)->getName() );
            } else if( getPntrToArgument(i)->getRank()==0 ) {
                ifile.scanField(getPntrToArgument(i)->getName(), val); 
                values.push_back( val ); names.push_back( getPntrToArgument(i)->getName() );
            } else if ( !getPntrToArgument(i)->usingAllVals( getLabel() ) ) {
                for(unsigned j=0;j<getPntrToArgument(i)->getNumberOfValues( getLabel() );++j) {  
                    std::string argname = getPntrToArgument(i)->getOutputDescription( getLabel(), j );
                    ifile.scanField(argname,val); values.push_back( val ); names.push_back( argname );
                }
            } else if( getPntrToArgument(i)->getRank()==1 ) {
                for(unsigned j=0;j<getPntrToArgument(i)->getShape()[0];++j) {
                    std::string num; Tools::convert( j+1, num );
                    ifile.scanField(getPntrToArgument(i)->getName() + "." + num, val);
                    values.push_back( val ); names.push_back( getPntrToArgument(i)->getName() + "." + num );
                }
            } else if( getPntrToArgument(i)->getRank()==2 ) {
                for(unsigned j=0;j<getPntrToArgument(i)->getShape()[0];++j) {
                    std::string jnum; Tools::convert( j+1, jnum );
                    for(unsigned k=0;k<getPntrToArgument(k)->getShape()[2];++k) {
                        std::string knum; Tools::convert( k+1, knum );
                        ifile.scanField(getPntrToArgument(i)->getName() + "." + jnum + "." + knum, val);
                        values.push_back( val ); names.push_back( getPntrToArgument(i)->getName() + "." + jnum + "." + knum );
                    }
                }
            } else error("cannot deal with objects with ranks greater than 2");
        }
        if(line==number-1) {
             readline=true;
             log.printf("  reading %dth reference structure from file %s \n", number, reference.c_str());
             log.printf("  which contains");
             if( read_args.size()>0 ) log.printf(" %d arguments \n", read_args.size() );
             else log.printf(" %d arguments \n", getNumberOfArguments() );
             if( getNumberOfArguments()>0 ) {
                 log.printf("  labels of arguments are : ");
                 for(unsigned i=0;i<getNumberOfArguments();++i) log.printf("%s ", getPntrToArgument(i)->getName().c_str() );
                 log.printf("\n");
             } 
             if( read_args.size()>0 ) {
                 log.printf("  labels of arguments are : ");
                 for(unsigned i=0;i<read_args.size();++i) log.printf("%s ", read_args[i].c_str() );
                 log.printf("\n");
             }
             
             if( getNumberOfArguments()>0 ) {
                 // Create the component that will hold the position of the center of the cluster
                 std::vector<unsigned> shape( 1 ); shape[0] = getNumberOfArguments();  
                 if( !numberedkeys ) {
                     shape[0]=0; for(unsigned i=0;i<getNumberOfArguments();++i) shape[0] += getPntrToArgument(i)->getNumberOfValues( getLabel() );
                 }
                 addComponent( "center", shape ); componentIsNotPeriodic("center"); getPntrToComponent(0)->buildDataStore( getLabel() );

                 // Read in the position of the center of the cluster
                 unsigned n=0;
                 for(unsigned i=0;i<getNumberOfArguments();++i) {
                     if( numberedkeys ) {
                         getPntrToComponent(0)->set( i, values[n] ); n++;
                     } else if( getPntrToArgument(i)->getRank()==0 ) {
                         getPntrToComponent(0)->set( i, values[n] ); n++; 
                     } else if ( !getPntrToArgument(i)->usingAllVals( getLabel() ) ) {
                         for(unsigned j=0;j<getPntrToArgument(i)->getNumberOfValues( getLabel() );++j) {
                             getPntrToComponent(0)->set( n, values[n] ); n++;
                         }
                     } else if( getPntrToArgument(i)->getRank()==1 ) {
                         for(unsigned j=0;j<getPntrToArgument(i)->getShape()[0];++j) {
                             getPntrToComponent(0)->set( n, values[n] ); n++; 
                         }
                     } else if( getPntrToArgument(i)->getRank()==2 ) { 
                         for(unsigned j=0;j<getPntrToArgument(i)->getShape()[0];++j) {
                             std::string jnum; Tools::convert( j+1, jnum );
                             for(unsigned k=0;k<getPntrToArgument(k)->getShape()[2];++k) {
                                 getPntrToComponent(0)->set( n, values[n] ); n++; 
                             }
                         }
                     } else error("cannot deal with objects with ranks greater than 2");
                 }
                 if( ifile.FieldExist("sigma_" + names[0]) ) {
                     addComponent("variance", shape); componentIsNotPeriodic("variance"); getPntrToComponent(1)->buildDataStore( getLabel() );
                     for(unsigned i=0;i<names.size();++i) {
                         ifile.scanField("sigma_" + names[i], val); getPntrToComponent(1)->set( i, val ); 
                     }
                 } else {
                     std::vector<unsigned> nshape(2); nshape[0]=nshape[1]=shape[0];
                     addComponent("covariance", nshape); componentIsNotPeriodic("covariance"); getPntrToComponent(1)->buildDataStore( getLabel() );
                     for(unsigned i=0; i<names.size(); ++i) {
                         for(unsigned j=0; j<names.size()-i; j++) {
                             ifile.scanField("sigma_" +names[j+i] + "_" + names[j], val ); 
                             getPntrToComponent(1)->set( (i+j)*names.size()+j, val );
                             getPntrToComponent(1)->set( j*names.size()+j+i, val );
                         }
                     }
                 }
             }
             if( read_args.size()>0 ) {
                 std::vector<unsigned> shape( 1 ); shape[0] = read_args.size(); double val;
                 addComponent( "center", shape ); componentIsNotPeriodic("center"); getPntrToComponent(0)->buildDataStore( getLabel() );
                 for(unsigned i=0;i<read_args.size();++i) {
                     ifile.scanField( read_args[i], val ); getPntrToComponent(0)->set( i, val );
                 }
                 if( ifile.FieldExist("sigma_" + read_args[0]) ) {
                     addComponent("variance", shape); componentIsNotPeriodic("variance"); getPntrToComponent(1)->buildDataStore( getLabel() );
                     for(unsigned i=0;i<read_args.size();++i) {
                         ifile.scanField("sigma_" + read_args[i], val); getPntrToComponent(1)->set( i, val ); 
                     }
                 } else {
                     std::vector<unsigned> nshape(2); nshape[0]=nshape[1]=shape[0];
                     addComponent("covariance", nshape); componentIsNotPeriodic("covariance"); getPntrToComponent(1)->buildDataStore( getLabel() );
                     for(unsigned i=0; i<read_args.size(); ++i) {
                         for(unsigned j=0; j<read_args.size()-i; j++) {
                             ifile.scanField("sigma_" +read_args[j+i] + "_" + read_args[j], val );
                             getPntrToComponent(1)->set( (i+j)*read_args.size()+j, val );
                             getPntrToComponent(1)->set( j*read_args.size()+i+j, val );
                         }   
                     }       
                 }  
             }
             break;
          }
          ifile.scanField(); 
      }
      if( !readline ) error("could not read reference configuration");
      ifile.scanField(); ifile.close();
  }
}

std::string ReadReferenceCluster::getArgName( const unsigned& k ) const {
  if( read_args.size()>0 ) return read_args[k];
  return SetupReferenceBase::getArgName(k);
}

}
}
