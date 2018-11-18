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
  keys.add("compulsory","REFERENCE","a file containing information on the reference cluster.");
  keys.add("compulsory","NUMBER","1","if there are multiple clusters in the input file which structure would you like to read in here");
  keys.addFlag("READ_VARIANCE",false,"the variance should be read from the input");
  keys.addFlag("READ_COVARIANCE",false,"the covariance should also be read from the input");
  keys.add("hidden","READ_ARG","this is used by pathtool to get the arguments that must be read in");
}

ReadReferenceCluster::ReadReferenceCluster(const ActionOptions&ao):
Action(ao),
SetupReferenceBase(ao)
{
  if( getNumberOfArguments()==0 ) parseVector("READ_ARG",read_args);
  std::string reference; parse("REFERENCE",reference); 
  IFile ifile; ifile.open(reference); ifile.allowIgnoredFields();
  unsigned number; parse("NUMBER",number); bool readline=false;
  bool read_covar; parseFlag("READ_COVARIANCE",read_covar);
  bool read_var; parseFlag("READ_VARIANCE",read_var);
  if( read_covar && read_var ) error("cannot read both the variance and the covariance matrix");
  for(unsigned i=0;i<number;++i) {
    ifile.scanField();  
    if(i==number-1) {
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
             std::vector<unsigned> shape( 1 ); shape[0] = 0; unsigned n=0; double val;
             for(unsigned i=0;i<getNumberOfArguments();++i) shape[0] += getPntrToArgument(i)->getNumberOfValues( getLabel() );
             addComponent( "center", shape ); componentIsNotPeriodic("center"); getPntrToComponent(0)->buildDataStore( getLabel() );
  
             // Read in the position of the center of the cluster
             std::vector<std::string> names;
             for(unsigned i=0;i<getNumberOfArguments();++i) {
                 if( getPntrToArgument(i)->getRank()==0 ) {
                     ifile.scanField(getPntrToArgument(i)->getName(), val);
                     getPntrToComponent(0)->set( n, val ); n++; names.push_back( getPntrToArgument(i)->getName() );
                 } else if( getPntrToArgument(i)->getRank()==1 ) {
                     for(unsigned j=0;j<getPntrToArgument(i)->getShape()[0];++j) {
                         std::string num; Tools::convert( j+1, num );
                         ifile.scanField(getPntrToArgument(i)->getName() + "." + num, val);
                         getPntrToComponent(0)->set( n, val ); n++; names.push_back( getPntrToArgument(i)->getName() + "." + num );
                     }
                 } else if( getPntrToArgument(i)->getRank()==2 ) { 
                     for(unsigned j=0;j<getPntrToArgument(i)->getShape()[0];++j) {
                         std::string jnum; Tools::convert( j+1, jnum );
                         for(unsigned k=0;k<getPntrToArgument(k)->getShape()[2];++k) {
                             std::string knum; Tools::convert( k+1, knum );
                             ifile.scanField(getPntrToArgument(i)->getName() + "." + jnum + "." + knum, val); 
                             getPntrToComponent(0)->set( n, val ); n++; names.push_back( getPntrToArgument(i)->getName() + "." + jnum + "." + knum );
                         }
                     }
                 } else error("cannot deal with objects with ranks greater than 2");
             }
             if( read_var ) {
                 addComponent("variance", shape); componentIsNotPeriodic("variance"); getPntrToComponent(1)->buildDataStore( getLabel() );
                 for(unsigned i=0;i<names.size();++i) {
                     ifile.scanField("sigma_" + names[i], val); getPntrToComponent(1)->set( i, val ); 
                 }
             } else if( read_covar ) {
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
             if( read_var ) {
                 addComponent("variance", shape); componentIsNotPeriodic("variance"); getPntrToComponent(1)->buildDataStore( getLabel() );
                 for(unsigned i=0;i<read_args.size();++i) {
                     ifile.scanField("sigma_" + read_args[i], val); getPntrToComponent(1)->set( i, val ); 
                 }
             } else if( read_covar ) {
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
  }
  if( !readline ) error("could not read reference configuration");
  ifile.scanField(); ifile.close();
}

std::string ReadReferenceCluster::getArgName( const unsigned& k ) const {
  if( read_args.size()>0 ) return read_args[k];
  return SetupReferenceBase::getArgName(k);
}

}
}
