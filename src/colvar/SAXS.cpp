/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015,2016 The plumed team
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
/* 
 This class was originally written by Alexander Jussupow and
 Carlo Camilloni 
*/

#include "Colvar.h"
#include "ActionRegister.h"
#include "tools/Communicator.h"
#include "tools/OpenMP.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD{
namespace colvar{

//+PLUMEDOC COLVAR SAXS
/*

Calculate SAXS scattered intensity

*/
//+ENDPLUMEDOC
   
class SAXS : public Colvar {
private:
  bool                     pbc;
  bool                     serial;
  unsigned                 numq;
  vector<double>           q_list;
  vector<vector<double> >  FF_value;
  vector<double>           FF_rank;

public:
  static void registerKeywords( Keywords& keys );
  explicit SAXS(const ActionOptions&);
  virtual void calculate();
};

PLUMED_REGISTER_ACTION(SAXS,"SAXS")

void SAXS::registerKeywords(Keywords& keys){
  Colvar::registerKeywords( keys );
  componentsAreNotOptional(keys);
  useCustomisableComponents(keys);
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.add("compulsory","NUMQ","Number of used q values");
  keys.add("compulsory","SCEXP","SCALING value of the experimental data");
  keys.add("atoms","ATOMS","The atoms to be included in the calculation, e.g. the whole protein.");
  keys.add("numbered","qvalue","Used qvalue Keywords like q_value1, q_value2, qvalue_3,... should be listed.");
  keys.add("numbered","parameter","Used parameter Keywords like parameter1, parameter2, parameter3,... should be listed.");
  keys.addFlag("ADDEXPVALUES",false,"Set to TRUE if you want to have fixed components with the experimetnal values.");
  keys.add("numbered","EXPINT","Add an experimental value for each PRE.");
}

SAXS::SAXS(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
pbc(true),
serial(false)
{
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  const unsigned size = atoms.size();

  parseFlag("SERIAL",serial);

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  numq = 0;
  parse("NUMQ",numq);
  if(numq==0) error("NUMQ must be set");  

  double scexp = 0;
  parse("SCEXP",scexp);
  if(scexp==0) error("SCEXP must be set");  

  q_list.resize( numq );
  unsigned ntarget=0;
  for(unsigned i=0;i<numq;++i){
    if( !parseNumbered( "qvalue", i+1, q_list[i]) ) break;
    ntarget++;
  }
  if( ntarget!=numq ) error("found wrong number of qvalue values");

  for(unsigned i=0;i<numq;i++) {
    log.printf("  my q: %lf \n",q_list[i]);    
    std::string num; Tools::convert(i,num);
    addComponentWithDerivatives("q_"+num);
    componentIsNotPeriodic("q_"+num);
  }

  //read in parameter vector
  vector<vector<long double> > parameter;
  parameter.resize(size);
  ntarget=0;
  for(unsigned i=0;i<size;++i){
    if( !parseNumberedVector( "parameter", i+1, parameter[i]) ) break;
    ntarget++;
  }
  if( ntarget!=size ) error("found wrong number of parameter vectors");

  FF_value.resize(numq,vector<double>(size));
  vector<vector<long double> >  FF_tmp;
  FF_tmp.resize(numq,vector<long double>(size));
  for(unsigned i=0;i<size;++i) {
    for(unsigned j=0;j<parameter[i].size();++j) {
      for(unsigned k=0;k<numq;++k){
        FF_tmp[k][i]+= parameter[i][j]*powl(static_cast<long double>(q_list[k]),j);
      }
    }
  }

  // Calculate Rank of FF_matrix
  FF_rank.resize(numq);
  for(unsigned k=0;k<numq;++k){
    for(unsigned i=0;i<size;i++){
       FF_value[k][i] = static_cast<double>(FF_tmp[k][i]);
       FF_rank[k]+=FF_value[k][i]*FF_value[k][i];
    }
  }

  bool exp=false;
  parseFlag("ADDEXPVALUES",exp);
  if(exp) {
    vector<double>   expint;
    expint.resize( numq ); 
    unsigned ntarget=0;
    for(unsigned i=0;i<numq;++i){
       if( !parseNumbered( "EXPINT", i+1, expint[i] ) ) break;
       ntarget++; 
    }
    if( ntarget!=numq ) error("found wrong number of NOEDIST values");

    for(unsigned i=0;i<numq;i++) {
      std::string num; Tools::convert(i,num);
      addComponent("exp_"+num);
      componentIsNotPeriodic("exp_"+num);
      Value* comp=getPntrToComponent("exp_"+num); comp->set(expint[i]*scexp);
    }
  }

  // convert units to nm^-1
  for(unsigned i=0;i<numq;++i){
    q_list[i]=q_list[i]*10.0;    //factor 10 to convert from A^-1 to nm^-1
  } 

  requestAtoms(atoms);
  checkRead();
}

void SAXS::calculate(){
  if(pbc) makeWhole();

  const unsigned size=getNumberOfAtoms();
  
  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  if(serial){
    stride=1;
    rank=0;
  }

  vector<Vector> deriv(numq*size);
  vector<Tensor> deriv_box(numq);
  vector<double> sum(numq);

  #pragma omp parallel for num_threads(OpenMP::getNumThreads())
  for (unsigned k=0; k<numq; k++) {
    const unsigned kdx=k*size;
    for (unsigned i=rank; i<size-1; i+=stride) {
      const unsigned kdxi=kdx+i;
      const double FF=2.*FF_value[k][i];
      const Vector posi=getPosition(i);
      for (unsigned j=i+1; j<size ; j++) {
        const Vector c_distances = delta(posi,getPosition(j));
        const double m_distances = c_distances.modulo();
        const double qdist     = q_list[k]*m_distances;
        const double FFF = FF*FF_value[k][j];
        const double tsq = FFF*sin(qdist)/qdist;
        const double tcq = FFF*cos(qdist);
        const double tmp = (tcq-tsq)/(m_distances*m_distances);
        const Vector dd  = c_distances*tmp;
        sum[k] += tsq;
        deriv_box[k] += Tensor(c_distances,dd);
        deriv[kdxi ] -= dd;
        deriv[kdx+j] += dd;
      }
    }
  }

  if(!serial) {
    comm.Sum(&deriv[0][0], 3*deriv.size());
    comm.Sum(&deriv_box[0][0][0], numq*9);
    comm.Sum(&sum[0], numq);
  }

  #pragma omp parallel for num_threads(OpenMP::getNumThreads())
  for (unsigned k=0; k<numq; k++) {
    const unsigned kdx=k*size;
    Value* val=getPntrToComponent(k);
    for(unsigned i=0; i<size; i++) setAtomsDerivatives(val, i, deriv[kdx+i]);
    sum[k]+=FF_rank[k];
    setBoxDerivatives(val, -deriv_box[k]);
    val->set(sum[k]);
  }

}

}
}
