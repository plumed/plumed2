/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#include "RMSDBase.h"
#include "MetricRegister.h"

namespace PLMD{

class SimpleRMSD : public RMSDBase {
private:
  std::vector<double> weights;
public:
  SimpleRMSD( const ReferenceConfigurationOptions& ro );
  void read( const PDB& );
  double calc( const std::vector<Vector>& pos, const bool& squared );
  double simpleAlignment(const  std::vector<double>  & align,
                         const  std::vector<double>  & displace,
                         const std::vector<Vector> & positions,
                         bool squared);
};

PLUMED_REGISTER_METRIC(SimpleRMSD,"SIMPLE")

SimpleRMSD::SimpleRMSD( const ReferenceConfigurationOptions& ro ):
ReferenceConfiguration( ro ),
RMSDBase( ro )
{
}

void SimpleRMSD::read( const PDB& pdb ){
  readAtomsFromPDB( pdb ); weights.resize( pdb.getAtomNumbers().size() );
  for(unsigned i=0;i<weights.size();++i) weights[i]=getDisplace()[i];
}

double SimpleRMSD::calc( const std::vector<Vector>& pos, const bool& squared ){
  return simpleAlignment( getAlign(), getDisplace(), pos, squared );
}

double SimpleRMSD::simpleAlignment(const  std::vector<double>  & align,
                                   const  std::vector<double>  & displace,
                                   const std::vector<Vector> & positions,
                                   bool squared) {

  double dist(0), anorm(0), dnorm(0);
  unsigned n=getNumberOfReferencePositions();
  // Number of members of align and displace is the number of reference positions
  plumed_dbg_assert( n==align.size() && n==displace.size() );
  // Positions array might contain vectors that are not particularly interesting
  plumed_dbg_assert( positions.size()==getNumberOfAtoms() );
  
  Vector iref, apositions, areference;
  Vector dpositions, dreference;
  
  for(unsigned i=0;i<n;i++){
    unsigned iatom=getAtomIndex(i);
    double aw=align[i];
    double dw=displace[i];
    anorm+=aw;
    dnorm+=dw;
    iref=getReferencePosition(i);
    apositions+=positions[iatom]*aw;
    areference+=iref*aw;
    dpositions+=positions[iatom]*dw;
    dreference+=iref*dw;
  }

  double invdnorm=1.0/dnorm;
  double invanorm=1.0/anorm;

  apositions*=invanorm;
  areference*=invanorm;
  dpositions*=invdnorm;
  dreference*=invdnorm;

  Vector shift=((apositions-areference)-(dpositions-dreference));
  for(unsigned i=0;i<n;i++){
    unsigned iatom=getAtomIndex(i);
    Vector d=(positions[iatom]-apositions)-(getReferencePosition(i)-areference);
    dist+=displace[i]*d.modulo2();
    addAtomicDerivatives( i,2*(invdnorm*displace[i]*d+invanorm*align[i]*shift) );
  }
  dist*=invdnorm;

  if(!squared){
     // sqrt and normalization
     dist=sqrt(dist);
     ///// sqrt and normalization on derivatives
     for(unsigned i=0;i<getNumberOfAtoms();i++){atom_ders[i]*=(0.5/dist);}
  }
  return dist;
}

}
