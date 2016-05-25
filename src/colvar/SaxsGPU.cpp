/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2016 The plumed team
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
#include "Colvar.h"
#include "ActionRegister.h"
#include "tools/Communicator.h"
#include "tools/OpenMP.h"

#include <string>
#include <cmath>

#ifdef __PLUMED_HAS_ARRAYFIRE
#include <arrayfire.h>
#include <af/util.h>
#endif

namespace PLMD{
namespace colvar{

//+PLUMEDOC COLVAR SAXSGPU
/*

Calculate SAXSGPU scattered intensity

*/
//+ENDPLUMEDOC
   
class SAXSGPU : public Colvar {
private:
  bool                pbc;
  bool                serial;
  unsigned            numq;
  unsigned            splitb;
  std::vector<double> q_list;
  int                 total_device;
  af::array           allFFa;
  af::array          *sum_device;
  af::array          *box_device;
  af::array          *deriv_device;

public:
  static void registerKeywords( Keywords& keys );
  explicit SAXSGPU(const ActionOptions&);
  ~SAXSGPU();
  virtual void calculate();
};

PLUMED_REGISTER_ACTION(SAXSGPU,"SAXSGPU")

void SAXSGPU::registerKeywords(Keywords& keys){
  Colvar::registerKeywords( keys );
  componentsAreNotOptional(keys);
  useCustomisableComponents(keys);
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.add("compulsory","NUMQ","Number of used q values");
  keys.add("compulsory","SCEXP","SCALING value of the experimental data");
  keys.add("compulsory","SPLITB","Spliting the length of the atom array");
  keys.add("atoms","ATOMS","The atoms to be included in the calculation, e.g. the whole protein.");
  keys.add("numbered","qvalue","Used qvalue Keywords like q_value1, q_value2, qvalue_3,... should be listed.");
  keys.add("numbered","parameter","Used parameter Keywords like parameter1, parameter2, parameter3,... should be listed.");
  keys.addFlag("ADDEXPVALUES",false,"Set to TRUE if you want to have fixed components with the experimetnal values.");
  keys.add("numbered","EXPINT","Add an experimental value for each PRE.");
  keys.addFlag("MULTIGPU",false,"Set to TRUE if you want to use multiple GPU");
}

SAXSGPU::SAXSGPU(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
pbc(true),
serial(false)
{
#ifndef __PLUMED_HAS_ARRAYFIRE
  error("SAXSGPU can only be used if ARRAYFIRE is installed");
#else
  std::vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  const unsigned size = atoms.size();

  parseFlag("SERIAL",serial);

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  numq = 0;
  parse("NUMQ",numq);
  if(numq==0) error("NUMQ must be set");  

  splitb = 0;
  parse("SPLITB",splitb);
  if(splitb==0) error("SPLITB must be set");  

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
    std::string num; Tools::convert(i,num);
    addComponentWithDerivatives("q_"+num);
    componentIsNotPeriodic("q_"+num);
  }

  //read in parameter vector
  std::vector<std::vector<double> > parameter;
  parameter.resize(size);
  ntarget=0;
  for(unsigned i=0;i<size;++i){
    if( !parseNumberedVector( "parameter", i+1, parameter[i]) ) break;
    ntarget++;
  }
  if( ntarget!=size ) error("found wrong number of parameter vectors");

  // Calculate Rank of FF_matrix
  float *FF_new = new float[numq*size];  
  for(unsigned i=0;i<size;++i) {
    for(unsigned j=0;j<parameter[i].size();++j) {
      for(unsigned k=0;k<numq;++k){
        FF_new[k+i*numq]+=parameter[i][j]*pow(q_list[k],j);
      }
    }
  }
  allFFa= af::array(numq, size, FF_new);
  delete[] FF_new;

  bool exp=false;
  parseFlag("ADDEXPVALUES",exp);
  if(exp) {
    std::vector<double>   expint;
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

  bool multi=false;
  parseFlag("MULTIGPU",multi);
  if(multi) {
    total_device = af::getDeviceCount();
  } else {
    total_device = 1;
  }

  sum_device = new af::array[total_device*numq];
  box_device = new af::array[total_device*numq];
  deriv_device = new af::array[total_device*numq];

  requestAtoms(atoms);
  checkRead();
#endif
}

SAXSGPU::~SAXSGPU(){
  delete[] sum_device;
  delete[] box_device;
  delete[] deriv_device;
}

void SAXSGPU::calculate(){
#ifdef __PLUMED_HAS_ARRAYFIRE
  if(pbc) makeWhole();

  const unsigned size=getNumberOfAtoms();
  
  float* posi;
  posi = new float[3*size];
  for (unsigned i=0; i<size; i++) {
    const Vector tmp = getPosition(i);
    posi[i*3+0] = tmp[0];
    posi[i*3+1] = tmp[1];
    posi[i*3+2] = tmp[2];
  }

  for(unsigned i=0;i<total_device; i++) {
     af::setDevice(i);
     sum_device[i] = af::constant(0, numq, f32);
     box_device[i] = af::constant(0, numq, 6, f32);
     deriv_device[i] = af::constant(0, numq, size, 3, f32);
  }

  for (unsigned i=0; i<size; i=i+splitb) {
    //multiple device
    int dnumber=(i/splitb) % total_device;
    af::setDevice(dnumber);
    //first step calculate the short size of the matrix
    int sizeb; 
    if (size-i > splitb) { 
      sizeb=splitb;
    } else {
      sizeb=size-i;
    }

    af::seq seqb(i, i+sizeb-1);
    af::seq seqa(0, size-1);

    // create array a and b containing atomic coordinates
    af::array a = af::array(3, size, posi);
    af::array b = a(af::span, seqb);
    b += 0.000001; // crapy solution

    // calculate distance matrix
    af::array b_mod = af::moddims(b, 3, 1, sizeb);
    // xyz_dist is 3,size,sizeb
    af::array xyz_dist = af::tile(a, 1, 1, sizeb) -  af::tile(b_mod, 1, size, 1);
    // square is 1,size,sizeb
    af::array square = af::sum(xyz_dist*xyz_dist);
    // dist_sqrt is 1,size,sizeb
    af::array dist_sqrt = af::sqrt(square);
    //this can be alternative to the above addition, but it is a bit slower
    //dist_sqrt(dist_sqrt == 0.) = 0.000001;
    

    // calculate distance vectors
    // xyz_dist is now size,sizeb,3
    xyz_dist = af::reorder(xyz_dist, 1, 2, 0);
    af::array atom_box = af::moddims(xyz_dist, size*sizeb, 3);
    af::array allFFb = allFFa(af::span, seqb);

    for (unsigned k=0; k<numq; k++) {
      // calculate FF matrix
      af::array FFdist_mod = (af::tile(af::moddims(allFFa.row(k), size, 1), 1, sizeb) * 
                              af::tile(af::moddims(allFFb.row(k), 1, sizeb), size, 1));

      // get q*dist and sin
      float qvalue = q_list[k];
      af::array dist_q = qvalue*dist_sqrt;
      af::array dist_sin = (af::sin(dist_q)/dist_q);
      
      // flat it and get the intensity
      sum_device[dnumber](k) += af::sum(af::flat(dist_sin)*af::flat(FFdist_mod));

      // array get cos and tmp
      // tmp is 1,size,sizeb
      af::array tmp = (dist_sin - af::cos(dist_q))/square;
      // tmp to size,sizeb
      tmp = af::moddims(tmp, size, sizeb);
      tmp = tmp*FFdist_mod;

      // increase the tmp size and calculate dd
      af::array tmp_tile = af::tile(tmp, 1, 1, 3);
      af::array dd_all = (tmp_tile*xyz_dist);
      af::array dd = af::sum(dd_all, 0);
 
      deriv_device[dnumber](k, seqb, af::span) = dd(0, af::span, af::span);

      dd_all = af::moddims(dd_all, size*sizeb, 3);
      af::array box_pre = af::array(size*sizeb, 6, f32);
      box_pre(af::span,0) = atom_box(af::span,0)*dd_all(af::span,0);
      box_pre(af::span,1) = atom_box(af::span,0)*dd_all(af::span,1);
      box_pre(af::span,2) = atom_box(af::span,0)*dd_all(af::span,2);
      box_pre(af::span,3) = atom_box(af::span,1)*dd_all(af::span,1);
      box_pre(af::span,4) = atom_box(af::span,1)*dd_all(af::span,2);
      box_pre(af::span,5) = atom_box(af::span,2)*dd_all(af::span,2);

      af::array box_sum = af::sum(box_pre);
      box_device[dnumber](k,af::span) += box_sum(0,af::span);
    }    
  }

  // read out results
  double* box;
  box = new double[numq*6];

  double* deriv;
  deriv = new double[size*3*numq];

  double* inten;
  inten = new double[numq];

  for (unsigned j=0; j<numq; j++) {
    inten[j] = 0;
  }

  for (unsigned j=0; j<6*numq; j++) {
    box[j] = 0;
  }

  for (unsigned j=0; j<size*3*numq; j++) {
    deriv[j] = 0;
  }

  for (unsigned i=0; i< total_device; i++) {
    af::setDevice(i);
    float* tmp_inten;
    tmp_inten = new float[numq];
    sum_device[i].host(tmp_inten);
    for (unsigned j=0; j<numq; j++) {
      inten[j] += tmp_inten[j];
    }

    float* tmp_box;
    tmp_box = new float[numq*6];
    box_device[i] = af::flat(box_device[i].T());
    box_device[i].host(tmp_box);
    for(unsigned j=0; j<6*numq; j++) {
      box[j] += tmp_box[j];
    }

    float* tmp_deriv;
    tmp_deriv = new float[size*3*numq];
    deriv_device[i] = af::reorder(deriv_device[i], 2, 1, 0);
    deriv_device[i] = af::flat(deriv_device[i]); 
    deriv_device[i].host(tmp_deriv);

    for(unsigned j=0; j<size*3*numq; j++) {
      deriv[j] += tmp_deriv[j];
    }

    delete[] tmp_inten;
    delete[] tmp_box;
    delete[] tmp_deriv;
  }

  #pragma omp parallel for num_threads(OpenMP::getNumThreads())
  for(unsigned k=0; k<numq; k++) {
    Value* val=getPntrToComponent(k);
    val->set(inten[k]);
    // get deriv Tensor
    Tensor deriv_box;

    deriv_box[0][0]=box[k*6+0];
    deriv_box[0][1]=box[k*6+1];
    deriv_box[0][2]=box[k*6+2];

    deriv_box[1][0]=box[k*6+1];
    deriv_box[1][1]=box[k*6+3];
    deriv_box[1][2]=box[k*6+4];

    deriv_box[2][0]=box[k*6+2];
    deriv_box[2][1]=box[k*6+4];
    deriv_box[2][2]=box[k*6+5];

    setBoxDerivatives(val, deriv_box);
    for(unsigned i=0;i<size;i++) {
      Vector dd;
      dd[0] = deriv[k*size*3+i*3+0];
      dd[1] = deriv[k*size*3+i*3+1];
      dd[2] = deriv[k*size*3+i*3+2];

      setAtomsDerivatives(val, i, 2*dd);
    }    
  }

  delete[] inten;
  delete[] box;
  delete[] deriv;
  delete[] posi;
#endif
}

}
}
