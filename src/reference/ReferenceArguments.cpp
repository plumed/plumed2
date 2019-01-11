/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2019 The plumed team
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
#include "ReferenceArguments.h"
#include "ReferenceAtoms.h"
#include "tools/OFile.h"
#include "core/Value.h"
#include "tools/PDB.h"

namespace PLMD {

ReferenceArguments::ReferenceArguments( const ReferenceConfigurationOptions& ro ):
  ReferenceConfiguration(ro),
  hasweights(false),
  hasmetric(false)
{
}

void ReferenceArguments::readArgumentsFromPDB( const PDB& pdb ) {
  ReferenceAtoms* aref=dynamic_cast<ReferenceAtoms*>( this );
  arg_names.resize( pdb.getArgumentNames().size() );
  for(unsigned i=0; i<arg_names.size(); ++i) arg_names[i]=pdb.getArgumentNames()[i];
  if( !aref && arg_names.size()==0 ) error("no arguments in input PDB file");

  reference_args.resize( arg_names.size() ); arg_der_index.resize( arg_names.size() );
  for(unsigned i=0; i<arg_names.size(); ++i) {
    if( !pdb.getArgumentValue(arg_names[i], reference_args[i]) ) error("argument " + arg_names[i] + " was not set in pdb input");
    arg_der_index[i]=i;
  }

  if( hasweights ) {
    plumed_massert( !hasmetric, "should not have weights if we are using metric");
    weights.resize( arg_names.size() ); sqrtweight.resize( arg_names.size() );
    for(unsigned i=0; i<reference_args.size(); ++i) {
      if( !pdb.getArgumentValue("sigma_" + arg_names[i], weights[i]) ) error("value sigma_" + arg_names[i] + " was not set in pdb input");
      sqrtweight[i] = sqrt( weights[i] );
    }
  } else if( hasmetric ) {
    plumed_massert( !hasweights, "should not have weights if we are using metric");
    double thissig; metric.resize( arg_names.size(), arg_names.size() );
    for(unsigned i=0; i<reference_args.size(); ++i) {
      for(unsigned j=i; j<reference_args.size(); ++j) {
        if( !pdb.getArgumentValue("sigma_" + arg_names[i] + "_" + arg_names[j], thissig) ) {
          error("value sigma_" + arg_names[i] + "_" + arg_names[j] + " was not set in pdb input");
        }
        metric(i,j)=metric(j,i)=thissig;
      }
    }
  } else {
    weights.resize( arg_names.size() ); sqrtweight.resize( arg_names.size() );
    for(unsigned i=0; i<weights.size(); ++i) sqrtweight[i]=weights[i]=1.0;
  }
}

void ReferenceArguments::setReferenceArguments( const std::vector<double>& arg_vals, const std::vector<double>& sigma ) {
  moveReferenceArguments( arg_vals );

  if( hasmetric ) {
    unsigned k=0;
    for(unsigned i=0; i<reference_args.size(); ++i) {
      for(unsigned j=i; j<reference_args.size(); ++j) {
        metric(i,j)=metric(j,i)=sigma[k]; k++;
      }
    }
    plumed_assert( k==sigma.size() );
  } else {
    plumed_assert( reference_args.size()==sigma.size() );
    for(unsigned i=0; i<reference_args.size(); ++i) weights[i]=sigma[i];
  }
}

void ReferenceArguments::moveReferenceArguments( const std::vector<double>& arg_vals ) {
  plumed_dbg_assert( reference_args.size()==arg_vals.size() );
  for(unsigned i=0; i<arg_vals.size(); ++i) reference_args[i]=arg_vals[i];
}

void ReferenceArguments::getArgumentRequests( std::vector<std::string>& argout, bool disable_checks ) {
  arg_der_index.resize( arg_names.size() );

  if( argout.size()==0 ) {
    for(unsigned i=0; i<arg_names.size(); ++i) {
      argout.push_back( arg_names[i] );
      arg_der_index[i]=i;
    }
  } else {
    if(!disable_checks) {
      if( arg_names.size()!=argout.size() ) error("mismatched numbers of arguments in pdb frames");
    }
    for(unsigned i=0; i<arg_names.size(); ++i) {
      if(!disable_checks) {
        if( argout[i]!=arg_names[i] ) error("found mismatched arguments in pdb frames");
        arg_der_index[i]=i;
      } else {
        bool found=false;
        for(unsigned j=0; j<arg_names.size(); ++j) {
          if( argout[j]==arg_names[i] ) { found=true; arg_der_index[i]=j; break; }
        }
        if( !found ) {
          arg_der_index[i]=argout.size(); argout.push_back( arg_names[i] );
        }
      }
    }
  }
}

const std::vector<double>& ReferenceArguments::getReferenceMetric() {
  if( hasmetric ) {
    unsigned ntot=(reference_args.size() / 2 )*(reference_args.size()+1);
    if( trig_metric.size()!=ntot ) trig_metric.resize( ntot );
    unsigned k=0;
    for(unsigned i=0; i<reference_args.size(); ++i) {
      for(unsigned j=i; j<reference_args.size(); ++j) {
        plumed_dbg_assert( fabs( metric(i,j)-metric(j,i) ) < epsilon );
        trig_metric[k]=metric(i,j); k++;
      }
    }
  } else {
    if( trig_metric.size()!=reference_args.size() ) trig_metric.resize( reference_args.size() );
    for(unsigned i=0; i<reference_args.size(); ++i) trig_metric[i]=weights[i];
  }
  return trig_metric;
}

double ReferenceArguments::calculateArgumentDistance( const std::vector<Value*> & vals, const std::vector<double>& arg,
    ReferenceValuePack& myder, const bool& squared ) const {
  double r=0; std::vector<double> arg_ders( vals.size() );
  if( hasmetric ) {
    for(unsigned i=0; i<reference_args.size(); ++i) {
      unsigned ik=arg_der_index[i]; arg_ders[ ik ]=0;
      double dp_i=vals[ik]->difference( reference_args[i], arg[ik] );
      for(unsigned j=0; j<reference_args.size(); ++j) {
        double dp_j;
        unsigned jk=arg_der_index[j];
        if(i==j) dp_j=dp_i;
        else dp_j=vals[jk]->difference( reference_args[j], arg[jk] );

        arg_ders[ ik ]+=2.0*metric(i,j)*dp_j;    // Factor of two for off diagonal terms as you have terms from ij and ji
        r+=dp_i*dp_j*metric(i,j);
      }
    }
  } else {
    for(unsigned i=0; i<reference_args.size(); ++i) {
      unsigned ik=arg_der_index[i];
      double dp_i=vals[ik]->difference( reference_args[i], arg[ik] );
      r+=weights[i]*dp_i*dp_i; arg_ders[ik]=2.0*weights[i]*dp_i;
    }
  }
  if(!squared) {
    r=sqrt(r); double ir=1.0/(2.0*r);
    for(unsigned i=0; i<arg_ders.size(); ++i) myder.setArgumentDerivatives( i, arg_ders[i]*ir );
  } else {
    for(unsigned i=0; i<arg_ders.size(); ++i) myder.setArgumentDerivatives( i, arg_ders[i] );
  }
  return r;
}

void ReferenceArguments::extractArgumentDisplacement( const std::vector<Value*>& vals, const std::vector<double>& arg, std::vector<double>& dirout ) const {
  if( hasmetric ) {
    plumed_error();
  } else {
    for(unsigned j=0; j<reference_args.size(); ++j) {
      unsigned jk=arg_der_index[j]; dirout[jk]=sqrtweight[j]*vals[jk]->difference( reference_args[j], arg[jk] );
    }
  }
}

double ReferenceArguments::projectArgDisplacementOnVector( const std::vector<double>& eigv, const std::vector<Value*>& vals, const std::vector<double>& arg, ReferenceValuePack& mypack ) const {
  if( hasmetric ) {
    plumed_error();
  } else {
    double proj=0;
    for(unsigned j=0; j<reference_args.size(); ++j) {
      unsigned jk=arg_der_index[j];
      proj += eigv[j]*sqrtweight[j]*vals[jk]->difference( reference_args[j], arg[jk] );
      mypack.setArgumentDerivatives( jk, eigv[j]*sqrtweight[j] );
    }
    return proj;
  }
}

void ReferenceArguments::displaceReferenceArguments( const double& weight, const std::vector<double>& displace ) {
  plumed_dbg_assert( displace.size()==getNumberOfReferenceArguments() );
  for(unsigned i=0; i<displace.size(); ++i) reference_args[i] += weight*displace[i];
}

}
