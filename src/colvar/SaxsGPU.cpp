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
Calculate SAXS scattered intensity on GPU
*/
//+ENDPLUMEDOC
   
class SAXSGPU : public Colvar {
private:
  bool                pbc;
  bool                serial;
  bool                npt;
  unsigned            numq;
  unsigned            splitb;
  std::vector<double> q_list;
  unsigned            total_device;
#ifdef __PLUMED_HAS_ARRAYFIRE
  af::array          *allFFa;
  af::array          *sum_device;
  af::array          *box_device;
  af::array          *deriv_device;
#endif

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
  keys.addFlag("NPT",false,"Set to TRUE if you want to use NPT");

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
  double *FF_tmp = new double[numq*size];  
  for(unsigned i=0;i<size;++i) {
    for(unsigned j=0;j<parameter[i].size();++j) {
      for(unsigned k=0;k<numq;++k){
        FF_tmp[k+i*numq] += parameter[i][j]*pow(q_list[k],j);
      }
    }
  }
  float *FF_new = new float[numq*size];  
  for(unsigned i=0;i<numq*size;i++) FF_new[i] = FF_tmp[i];
  delete[] FF_tmp;

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

  npt=false;
  parseFlag("NPT",npt);

  sum_device = new af::array[total_device*numq];
  box_device = new af::array[total_device*numq];
  deriv_device = new af::array[total_device*numq];

  allFFa = new af::array[total_device];
  for(unsigned i=0;i<total_device; i++) {
     af::setDevice(i);
     allFFa[i] = af::array(numq, size, FF_new);
  }

  delete[] FF_new;


  requestAtoms(atoms);
  checkRead();
#endif
}

SAXSGPU::~SAXSGPU(){
#ifdef __PLUMED_HAS_ARRAYFIRE
  delete[] sum_device;
  delete[] box_device;
  delete[] deriv_device;
  delete[] allFFa;
#endif
}

void SAXSGPU::calculate(){
#ifdef __PLUMED_HAS_ARRAYFIRE
  if(pbc) makeWhole();

  const unsigned size=getNumberOfAtoms();
  
  float* posi;
  posi = new float[3*size];
  #pragma omp parallel for num_threads(OpenMP::getNumThreads())
  for (unsigned i=0; i<size; i++) {
    const Vector tmp = getPosition(i);
    posi[i]        = tmp[0];
    posi[i+size]   = tmp[1];
    posi[i+2*size] = tmp[2];
  }

  for(unsigned i=0;i<total_device; i++) {
    af::setDevice(i);
    sum_device[i]   = af::constant(0, numq, f32);
    if(npt) box_device[i]   = af::constant(0, numq, 6, f32);
    deriv_device[i] = af::constant(0, numq, size, 3, f32);
  }

  for (unsigned i=0; i<size; i=i+splitb) {
    //multiple device
    const int dnumber=(i/splitb) % total_device;
    af::setDevice(dnumber);

    //first step calculate the short size of the matrix
    unsigned sizeb = size - i; 
    if(sizeb > splitb) sizeb = splitb;
    af::seq seqb(i, i+sizeb-1);
    const unsigned size2=size*sizeb;

    // create array a and b containing atomic coordinates
    af::array a = af::array(size, 3, posi);
    af::array b = a(seqb, af::span);
    a += 0.000001; // crapy solution

    a = af::moddims(a, size, 1, 3);
    b = af::moddims(b, 1, sizeb, 3);
    af::array xyz_dist = af::tile(a, 1, sizeb, 1) - af::tile(b, size, 1, 1);

    // atom_box is now size*sizeb,3
   // if(npt) {
    af::array atom_box;
    if(npt) {
      atom_box = af::moddims(xyz_dist, size2, 3);
    }

    // square size,sizeb,1
    af::array square = af::sum(xyz_dist*xyz_dist,2);
    // dist_sqrt is size,sizeb,1
    af::array dist_sqrt = af::sqrt(square);

    // allFA numq,size
    // allFB numq,sizeb
    af::array allFFb = allFFa[dnumber](af::span, seqb);

    for (unsigned k=0; k<numq; k++) {
      // calculate FF matrix
      // FFdist_mod size,sizeb,1
      af::array FFdist_mod = (af::tile(af::moddims(allFFa[dnumber].row(k), size, 1, 1), 1, sizeb)* 
                              af::tile(af::moddims(allFFb.row(k), 1, sizeb), size, 1, 1));

      // get q*dist and sin
      const float qvalue = q_list[k];
      // distq size,sizeb,1
      af::array dist_q = qvalue*dist_sqrt;
      // dist_sin size,sizeb,1
      af::array dist_sin = af::sin(dist_q)/dist_q;      
      // flat it and get the intensity
      sum_device[dnumber](k) += af::sum(af::flat(dist_sin)*af::flat(FFdist_mod));

      // array get cos and tmp
      // tmp is size,sizeb
      af::array tmp = af::moddims(FFdist_mod*(dist_sin - af::cos(dist_q))/square, size, sizeb);

      // increase the tmp size and calculate dd
      // now is size, sizeb, 3
      af::array dd_all = af::tile(tmp, 1, 1, 3)*xyz_dist;
      deriv_device[dnumber](k, seqb, af::span) = af::sum(dd_all);

      if(npt) {
        // dd_all to size*sizeb, 3
        dd_all = af::moddims(dd_all, size2, 3);
        af::array box_pre = af::array(size2, 6, f32);
        af::seq so0(0,2,1);
        af::seq so1(3,4,1);
        af::seq so2(0,1,1);
        af::seq so3(1,2,1);
        box_pre(af::span,so0) = atom_box*dd_all;
        box_pre(af::span,so1) = atom_box(af::span,so2)*dd_all(af::span,so3);
        box_pre(af::span,5)   = atom_box(af::span,0)*dd_all(af::span,2);

        af::array box_sum = af::sum(box_pre);
        box_device[dnumber](k,af::span) += box_sum(0,af::span);
      }
    }   
  }
  delete[] posi;

  // accumulate the results
  std::vector<double> inten; inten.resize(numq,0);
  std::vector<double> deriv; deriv.resize(numq*size*3,0);
  std::vector<double> box; 
  if(npt) box.resize(numq*6,0);

  // read out results
  for (unsigned i=0; i<total_device; i++) {
    af::setDevice(i);
    float* tmp_inten;
    tmp_inten = new float[numq];
    sum_device[i].host(tmp_inten);

    float* tmp_box;
    if (npt) {
      tmp_box = new float[numq*6];
      box_device[i] = af::flat(box_device[i].T());
      box_device[i].host(tmp_box);
    }

    float* tmp_deriv;
    tmp_deriv = new float[size*3*numq];
    deriv_device[i] = af::reorder(deriv_device[i], 2, 1, 0);
    deriv_device[i] = af::flat(deriv_device[i]); 
    deriv_device[i].host(tmp_deriv);

    #pragma omp parallel num_threads(OpenMP::getNumThreads())
    {
      
      #pragma omp for nowait
      for(unsigned i=0; i<numq; i++) {
        inten[i] += tmp_inten[i];
        if (npt) {
          const unsigned wi = 6*i;
          box[wi+0] += tmp_box[wi+0];
          box[wi+1] += tmp_box[wi+3];
          box[wi+2] += tmp_box[wi+5];
          box[wi+3] += tmp_box[wi+1];
          box[wi+4] += tmp_box[wi+4];
          box[wi+5] += tmp_box[wi+2];
        }
      }
      #pragma omp for nowait
      for(unsigned i=0; i<size*3*numq; i++) {
        deriv[i] += tmp_deriv[i];
      }
    }
 
    delete[] tmp_inten;
    if (npt) {
      delete[] tmp_box;
    }
    delete[] tmp_deriv;
  }

  #pragma omp parallel for num_threads(OpenMP::getNumThreads())
  for(unsigned k=0; k<numq; k++) {
    Value* val=getPntrToComponent(k);
    val->set(inten[k]);
    // get deriv Tensor
    if(npt) {
      const unsigned bi = k*6;
      const Tensor deriv_box(box[bi+0],box[bi+1],box[bi+2],
                             box[bi+1],box[bi+3],box[bi+4],
                             box[bi+2],box[bi+4],box[bi+5]);
   
      setBoxDerivatives(val, deriv_box);
    }
    for(unsigned i=0;i<size;i++) {
      const unsigned di = k*size*3+i*3;
      const Vector dd(deriv[di+0],deriv[di+1],deriv[di+2]);
      setAtomsDerivatives(val, i, 2*dd);
    }    
  }
#endif
}

}
}
