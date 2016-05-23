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

#include "test1.h"

#include <arrayfire.h>
#include <af/util.h>

//#include <algorithm>


using namespace af;
//#include "array_test.h"
//#include "test2.h"

using namespace std;
//using namespace af;

namespace PLMD{
namespace colvar{

//+PLUMEDOC COLVAR SAXSGPU
/*

Calculate SAXSGPU scattered intensity

*/
//+ENDPLUMEDOC
   
class SAXSGPU : public Colvar {
private:
  bool                     pbc;
  bool                     serial;
  unsigned                 numq;
  unsigned                 splitb;
  vector<double>           q_list;
  vector<vector<double> >  FF_value;
  double*                  FF_new;
  vector<double>           FF_rank;
  int                      total_device;
  // test if FFF values can be keep in device memory
  //array                    FF_device;

public:
  static void registerKeywords( Keywords& keys );
  explicit SAXSGPU(const ActionOptions&);
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
    //log.printf("  my q: %lf \n",q_list[i]);    
    std::string num; Tools::convert(i,num);
    addComponentWithDerivatives("q_"+num);
    componentIsNotPeriodic("q_"+num);
  }

  //read in parameter vector
  vector<vector<double> > parameter;
  parameter.resize(size);
  ntarget=0;
  for(unsigned i=0;i<size;++i){
    if( !parseNumberedVector( "parameter", i+1, parameter[i]) ) break;
    ntarget++;
  }
  if( ntarget!=size ) error("found wrong number of parameter vectors");


  FF_value.resize(numq,vector<double>(size));
  for(unsigned i=0;i<size;++i) {
    for(unsigned j=0;j<parameter[i].size();++j) {
      for(unsigned k=0;k<numq;++k){
        FF_value[k][i]+=parameter[i][j]*pow(q_list[k],j);
      }
    }
  }

  FF_new = new double[numq*size];  
  // Calculate Rank of FF_matrix
  FF_rank.resize(numq);
  for(unsigned k=0;k<numq;++k){
    for(unsigned i=0;i<size;i++){
       FF_rank[k]+=FF_value[k][i]*FF_value[k][i];
       FF_new[k+i*numq]=FF_value[k][i];
       // save FF on device
    }
  }
/*
  array a= array(numq, size, FF_new);
  int feat_len = a.dims(0);
  // Shape of a is (feat_len, alen, 1)
  array a_mod = a;
  // Reshape b from (feat_len, blen) to (feat_len, 1, blen)
  array b_mod = moddims(a, feat_len, 1, size);

  // Tile both matrices to be (feat_len, alen, blen)
  array a_tiled = tile(a_mod, 1, 1, size);
  array b_tiled = tile(b_mod, 1, size, 1);

  // Do The sum operation along first dimension
  // Output is of shape (1, alen, blen)
  array dist_mod = (a_tiled * b_tiled);

  // Reshape dist_mat from (1, alen, blen) to (alen, blen)
  array FF_device = moddims(dist_mod, numq, size*size);
//  FF_new.resize(numq*size);
//  for(unsigned k=0;k<numq;++k){
//    for(unsigned i=0;i<size;++i) {
//      FF_new[k*size+i]=FF_value[k][i];
//    }
//  }
*/

  for(unsigned k=0;k<numq;++k){
     //log.printf("  FF_rank: %lf \n",FF_rank[k]);
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

  //log.printf("  size: %lf \n",double(FF_device.dims(0)));
  //log.printf("  size: %lf \n",double(FF_device.dims(1)));

bool multi=false;
  parseFlag("MULTIGPU",multi);
  if(multi) {
    total_device = devicecount();
  } else {
    total_device = 1; // devicecount();
  }

  requestAtoms(atoms);
  checkRead();
}

void SAXSGPU::calculate(){
  if(pbc) makeWhole();

  const unsigned size=getNumberOfAtoms();
  
  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  if(serial){
    stride=1;
    rank=0;
  }

  double* posi;
  posi = new double[3*size];
  //#pragma omp parallel for num_threads(OpenMP::getNumThreads())
  for (unsigned i=0; i<size; i++) {
    const Vector tmp = getPosition(i);
    posi[i*3+0] = tmp[0];
    posi[i*3+1] = tmp[1];
    posi[i*3+2] = tmp[2];
  }

  //af::setDevice(0);
  //af::info();

  double* inten;
  inten = new double[numq];
  for (unsigned k=0; k<numq; k++) {
    inten[k] = 0;
  }
  // transform form factors and atom list
  // array allpos = array(3, size, posi);
  // array allFF = array(numq, sizea, FF_value);
  // array q_device = array(numq, &q_list[0]);

  // perform calculation on multiple device
  //int total_device = devicecount();

//  array *sum_sync = new array[total_device];
//  array *box_sync = new array[total_device];
//  array *deriv_sync = new array[total_device];

/*
  array sum_device = constant(0, numq, f64);
  array box_device = constant(0, numq, 6, f64);
  array deriv_device = constant(0, numq, size, 3, f64);
*/


  //log.printf("  devicecount: %lf \n",double(total_device));

  array *sum_device = new array[total_device*numq];
  array *box_device = new array[total_device*numq];
  array *deriv_device = new array[total_device*numq];


  for (unsigned i=0;i<total_device; i++) {
     af::deviceset(i);
     sum_device[i] = constant(0, numq, f64);
     box_device[i] = constant(0, numq, 6, f64);
     deriv_device[i] = constant(0, numq, size, 3, f64);
     //log.printf("  box_pre: %lf \n",double(sum_device[i].dims(0)));
  }





  //
  for (unsigned i=0; i<size; i=i+splitb) {
    //multiple device
    int dnumber=(i/splitb) % total_device;
    deviceset(dnumber);
    //log.printf("  deviceset: %lf \n",double((i/splitb) % total_device));
    //first step calculate the short size of the matrix
    int sizeb; 
    if (size-i > splitb) { 
      sizeb=splitb;
    } else {
      sizeb=size-i;
    }

    //log.printf("  sizeb: %lf \n", double(sizeb));
    seq seqb(i, i+sizeb-1);
    seq seqa(0, size-1);



    // create array a and b containing atomic coordinates
    array a = array(3, size, posi);
    array b = a(span, seqb);
    a += 0.0000001; // crapy solution
    //double* posib;
    //posib = new double[3*sizeb];
    //std::copy(posi+(3*i), posi+(3*(sizeb+i)), posib);
    //array b = array(3, sizeb, posib);
    //delete[] posib;

    // calculate distance matrix

    //log.printf("  array erfolgreich \n"); 

    array b_mod = moddims(b, 3, 1, sizeb);
    array a_tiled = tile(a, 1, 1, sizeb);
    array b_tiled = tile(b_mod, 1, size, 1);
    array xyz_dist = a_tiled - b_tiled;

    array square = sum(xyz_dist*xyz_dist);
    array dist_sqrt = sqrt(square);

    // calculate distance vectors
    //array a_trans = a.T();
    //array b_trans = b.T();

    //array a_mod = moddims(a_trans, size, 1, 3);
    //b_mod = moddims(b_trans, 1, sizeb, 3);
    //a_tiled = tile(a_mod, 1, sizeb, 1);
    //b_tiled = tile(b_mod, size, 1, 1);
    //xyz_dist = a_tiled - b_tiled;
    xyz_dist = reorder(xyz_dist, 1, 2, 0);
    array atom_box = moddims(xyz_dist, size*sizeb, 3);
    //log.printf("  atom_box: %lf \n",max<double>((atom_box(0,0))));
    //log.printf("  atom_box2: %lf \n",max<double>(sum(atom_box,0)));

    // read in FF_values
    //double* FF_newb;
    //FF_newb = new double[numq*sizeb];
    //std::copy(FF_new+(numq*i), FF_new+(numq*(sizeb+i)), FF_newb);
    array allFFa= array(numq, size, FF_new);
    array allFFb= allFFa(span, seqb);
    //delete[] FF_newb;

    //log.printf("  bis k loop erfolgreich \n"); 

    for (unsigned k=0; k<numq; k++) {

      //deviceset((i/splitb) % total_device);

      // calculate FF matrix
      array FFtmpa = allFFa.row(k);
      array FFtmpb = allFFb.row(k);
      
      array FFa_mod = moddims(FFtmpa, size, 1);
      array FFb_mod = moddims(FFtmpb, 1, sizeb);
      array FFa_tiled = tile(FFa_mod, 1, sizeb);
      array FFb_tiled = tile(FFb_mod, size, 1);
      array FFdist_mod = (FFa_tiled * FFb_tiled);

      // get q*dist and sin
      double qvalue = q_list[k];
      array dist_q = qvalue*dist_sqrt;
      array dist_sin = (sin(dist_q)/dist_q);

      //log.printf("  bis dist_sin erfolgreich \n"); 
      //replace(dist_sin, isNaN(dist_sin), 0);
      //array nans = isNaN(dist_sin);
      //if (anyTrue<bool>(nans)) {
      //   dist_sin(nans) = 1;
      //}
      
      // flat it and get the intensity
      array flat_vec1 = flat(dist_sin);
      array flat_vec2 = flat(FFdist_mod);
      //log.printf("  bis flat_vec erfolgreich \n"); 
      //sum_device(k) += sum(flat_vec1*flat_vec2, 0, 0.0);       
      //sum_device(k) += sum(flat_vec1*flat_vec2);  

      //log.printf(" deviceget: %lf \n", double(deviceget()));     
      sum_device[dnumber](k) += sum(flat_vec1*flat_vec2);       
      //log.printf(" ich lebe noch2");

      // array get cos and tmp
      array tmp = cos(dist_q);
      tmp = (dist_sin - tmp)/square;

      //replace(tmp, isNaN(tmp), 0);

      tmp = moddims(tmp, size, sizeb);
      tmp = tmp*FFdist_mod;

      // increase the tmp size and calculate dd
      array tmp_tile = tile(tmp, 1, 1, 3);
      array dd_all = (tmp_tile*xyz_dist);
      //log.printf("  dd_all: %lf \n",max<double>(dd_all(0,1)));
      array dd = sum(dd_all, 0);
      //log.printf("  dd: %lf \n",max<double>(dd(0,1,1)));
      //dd = moddims(dd, sizeb, 3);
      //deriv_device(k, seqb, span) = dd(0, span, span); 
 
      deriv_device[dnumber](k, seqb, span) = dd(0, span, span);

      //log.printf("  deriv_device: %lf \n",max<double>(deriv_device(k,1,1)));



      //log.printf("  dd: %lf \n",max<double>(dd(0,0)));
      // calculate box values
      //array atom_box = moddims(xyz_dist, size*sizeb, 3);

      dd_all = moddims(dd_all, size*sizeb, 3);
      //log.printf("  dd: %lf \n",max<double>(dd_all(0)));
     // array box_pre = array(size*sizeb, 6);
      array box_pre = array(size*sizeb, 6, f64);
      //box_pre = atom_box * dd_all;
      box_pre(span,0) = atom_box(span,0)*dd_all(span,0);
      box_pre(span,1) = atom_box(span,0)*dd_all(span,1);
      box_pre(span,2) = atom_box(span,0)*dd_all(span,2);
      box_pre(span,3) = atom_box(span,1)*dd_all(span,1);
      box_pre(span,4) = atom_box(span,1)*dd_all(span,2);
      box_pre(span,5) = atom_box(span,2)*dd_all(span,2);
      //gfor (seq j, size*sizeb) {
      //   box_pre(j, 0) = atom_box(j, 0); // * dd_all(j, 0);
     //    box_pre(j, 1) = atom_box(j, 0) ; //* dd_all(j, 1);
     //    box_pre(j, 2) = atom_box(j, 0) ; //* dd_all(j, 2);
    //     box_pre(j, 3) = atom_box(j, 1) ; //* dd_all(j, 1);
    //     box_pre(j, 4) = atom_box(j, 1) ; //* dd_all(j, 2);
    //     box_pre(j, 5) = atom_box(j, 2) ; //* dd_all(j, 2);
   //   }
      //log.printf("  box_pre: %lf \n",double(box_pre.dims(0)));
      //log.printf("  box_pre: %lf \n",double(box_pre.dims(1)));
     // for(unsigned j=0; j<size*sizeb;j++) if(max<double>(box_pre(j,0))!=max<double>(atom_box(j,0))) printf("%lf %lf %lf\n",max<double>(atom_box(j,0)), max<double>(box_pre(j,0)), j);
      //log.printf("  box_pre_last: %i %lf \n",k, max<double>(box_pre(0, 0)));
      //log.printf("  atom_box_last: %i %lf \n",k, max<double>(atom_box(0, 0)));

      //array box_sum = sum(box_pre, 0, 0.0);
     // array box_sum = sum(box_pre);

      array box_sum = sum(box_pre);
      //log.printf("  size: %lf \n",double(box_sum.dims(0)));
      //log.printf("  size: %lf \n",double(box_sum.dims(1)));
/*
      box_device(k,0) += box_sum(0,0);
      box_device(k,1) += box_sum(0,1);
      box_device(k,2) += box_sum(0,2);
      box_device(k,3) += box_sum(0,3);
      box_device(k,4) += box_sum(0,4);
      box_device(k,5) += box_sum(0,5);
*/
      box_device[dnumber](k,0) += box_sum(0,0);
      box_device[dnumber](k,1) += box_sum(0,1);
      box_device[dnumber](k,2) += box_sum(0,2);
      box_device[dnumber](k,3) += box_sum(0,3);
      box_device[dnumber](k,4) += box_sum(0,4);
      box_device[dnumber](k,5) += box_sum(0,5);

      //log.printf("  box_pre: %lf \n",max<double>(box_pre(0,0)));
      //log.printf("  box_pre: %lf \n",max<double>(box_pre(0,1)));
      //log.printf("  box_sum: %i %lf \n",k, max<double>(box_sum));
      //log.printf("  atom_box3: %lf \n",max<double>(sum(atom_box,0)));
      //log.printf("  atom_box: %lf \n",max<double>(atom_box(0)));
      //log.printf("  atom_box2: %lf \n",max<double>(sum(atom_box(0))));
    //log.printf("  atom_box3: %lf \n",max<double>((atom_box(0,1))));
    //log.printf("  atom_box4: %lf \n",max<double>(sum(atom_box,0)));
    }    
  }

  // read out results


  //sum_device.host(inten);

  double* box;
  box = new double[numq*6];
  //box_device = flat(box_device.T()); 
  //box_device = flat(box_device); 
  //box_device.host(box);

  double* deriv;
  deriv = new double[size*3*numq];
  //deriv_device = reorder(deriv_device, 2, 1, 0);
  //deriv_device = flat(deriv_device); 
  //deriv_device.host(deriv);


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
    deviceset(i);
    double* tmp_inten;
    tmp_inten = new double[numq];
    sum_device[i].host(tmp_inten);
    for (unsigned j=0; j<numq; j++) {
      inten[j] += tmp_inten[j];
    }

    double* tmp_box;
    tmp_box = new double[numq*6];
    box_device[i] = flat(box_device[i].T());
    box_device[i].host(tmp_box);
    for (unsigned j=0; j<6*numq; j++) {
      box[j] += tmp_box[j];
    }


    double* tmp_deriv;
    tmp_deriv = new double[size*3*numq];
    deriv_device[i] = reorder(deriv_device[i], 2, 1, 0);
    deriv_device[i] = flat(deriv_device[i]); 
    deriv_device[i].host(tmp_deriv);

    for (unsigned j=0; j<size*3*numq; j++) {
      deriv[j] += tmp_deriv[j];
    }

    delete[] tmp_inten;
    delete[] tmp_box;
    delete[] tmp_deriv;
  }

  #pragma omp parallel for num_threads(OpenMP::getNumThreads())
  for (unsigned k=0; k<numq; k++) {
    Value* val=getPntrToComponent(k);
    //inten[k]+=FF_rank[k];
    val->set(inten[k]);
    // get deriv Tensor
    Tensor deriv_box;
   // deriv_box[0][0]=box[k];
   // deriv_box[0][1]=box[(6*numq)*k+1];
   // log.printf("  box: %lf \n",double(box[k*6]));
   // log.printf("  box: %lf \n",double(box[k*6+1]));

    deriv_box[0][0]=box[k*6+0];
    deriv_box[0][1]=box[k*6+1];
    deriv_box[0][2]=box[k*6+2];

    deriv_box[1][0]=box[k*6+1];
    deriv_box[1][1]=box[k*6+3];
    deriv_box[1][2]=box[k*6+4];

    deriv_box[2][0]=box[k*6+2];
    deriv_box[2][1]=box[k*6+4];
    deriv_box[2][2]=box[k*6+5];

    //log.printf("  box: %lf \n",double(box[k*6]));
    //log.printf("  box: %lf \n",double(box[k*6+1]));
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

  delete[] sum_device;
  delete[] box_device;
  delete[] deriv_device;

/*
  #pragma omp parallel for num_threads(OpenMP::getNumThreads())
  for (unsigned k=0; k<numq; k++) {
    const unsigned kdx=k*size;
    Value* val=getPntrToComponent(k);
    for (unsigned i=0; i<size; i++) {
        setAtomsDerivatives(val, i, deriv[kdx+i]);

    }
    inten[k]+=FF_rank[k];
    setBoxDerivatives(val, -deriv_box[k]);
    val->set(sum[k]);
  }
*/

/*  // arrayfire try
  //int splitb = 256; // as dummy value
  for (unsigned i=0; i<size; i=i+splitb) {
    int sizeb; 
    if (size-i > splitb) { 
      sizeb=splitb;
    } else {
      sizeb=size-i;
    }
    //log.printf("  sizeb: %lf \n",double(sizeb));
    //log.printf("  i: %lf \n",double(i));
    int sizea = size;
    double* posib;
    posib = new double[3*size];
    std::copy(posi+(3*i), posi+(3*(sizeb+i)), posib);

    // create distance matrix
    array a = array(3, sizea, posi);
    array b = array(3, sizeb, posib);



    int feat_len = a.dims(0);
    array a_mod = a;
    array b_mod = moddims(b, feat_len, 1, sizeb);
    array a_tiled = tile(a_mod, 1, 1, sizeb);
    array b_tiled = tile(b_mod, 1, sizea, 1);
    array xyz_dist = a_tiled - b_tiled;
    array dist_mod = sqrt(sum(xyz_dist*xyz_dist));


    //array dist_mat = moddims(dist_mod, sizea*sizeb);
    array dist_mat = moddims(dist_mod, sizea*sizeb, 1);

    // create FF matrix
    double* FF_newb;
    FF_newb = new double[numq*sizeb];
    std::copy(FF_new+(numq*i), FF_new+(numq*(sizeb+i)), FF_newb);
    array FFtmpa= array(numq, sizea, FF_new);
    array FFtmpb= array(numq, sizeb, FF_newb);

    feat_len = numq;
    array FFa = FFtmpa.T();
    array FFb = FFtmpb.T();
    //log.printf("  size: %lf \n",double(FFb.dims(0)));
    //log.printf("  size: %lf \n",double(FFb.dims(1)));

    array FFa_mod = moddims(FFa, sizea, 1, feat_len);
    array FFb_mod = moddims(FFb, 1, sizeb, feat_len);
    array FFa_tiled = tile(FFa_mod, 1, sizeb, 1);
    array FFb_tiled = tile(FFb_mod, sizea, 1, 1);
    array FFdist_mod = (FFa_tiled * FFb_tiled);

    array FF_device = moddims(FFdist_mod, sizea*sizeb, numq);
    // To transpose a matrix cost a lot of time
    //array FF_device = FF_devicetmp.T();

    // calculate q*dist
    array dist_mat_tile = tile(dist_mat, 1, numq);
    array q_device = array(1, numq, &q_list[0]);
    //array q_device = moddims(q_device_tmp, 1, numq);
    array q_tile = tile(q_device, sizea*sizeb, 1);
    array qdist = dist_mat_tile*q_tile;
    array qdistsin = sin(qdist)/qdist*FF_device;

    const int cx = 1;
    const double cy = 0.0;
    //qdist=qdist.T();
    //log.printf("  size: %lf \n",double(qdist.dims(0)));
    //log.printf("  size: %lf \n",double(qdist.dims(1)));
    //printMemInfo();
    array sum_device = sum(qdistsin, 0, cy);
    //printMemInfo();
    double sumtmp[numq];
    sum_device.host(sumtmp);

    for (unsigned k=0; k<numq; k++) {
      inten[k] += sumtmp[k];
    }
    // derivative first step calculate tmp
    qdist = (cos(qdist) - qdistsin)/(dist_mat_tile*dist_mat_tile);
    array tmp = moddims(qdist, sizea, sizeb, 1, numq);
    array tmp_tile = tile(tmp, 1, 1, 3, 1);
    array xyz_dist_tile = tile (xyz_dist, 1, 1, 1, numq); 
    
    //array tmp = moddims(qdist_tile, sizea, sizeb, numq, 3);
    
  }

  #pragma omp parallel for num_threads(OpenMP::getNumThreads())
  for (unsigned k=0; k<numq; k++) {
    Value* val=getPntrToComponent(k);
    inten[k]+=FF_rank[k];
    val->set(inten[k]);
  }
*/

/*
  double inten[numq];
  //#pragma omp parallel for num_threads(OpenMP::getNumThreads())
  for (unsigned k=0; k<numq; k++) {
    
    //array FFtmp= array(size, &FF_value[k][0]);
    //array FFa_mod = moddims(FFtmp, size, 1);
    //array FFb_mod = moddims(FFtmp, 1, size);
    //array FFa_tiled = tile(FFa_mod, 1, size);
    //array FFb_tiled = tile(FFb_mod, size, 1);
    //array FFdist_mod = FFa_tiled*FFb_tiled; 
    //FF_device = moddims(FFdist_mod, size*size);
    
    //int size2 = 808201;
    //gfor (seq i, 0, 1, size*size-1){
    //  dist(i)=sin(q_list[k]*dist(i))/(q_list[k]*dist(i))*FF_device(i);
    //}

    Value* val=getPntrToComponent(k);
    const int cx = -1;
    const double cy = 0.0;

    inten[k] = (af::sum<double>((sin(dist_mat*q_list[k])/(dist_mat*q_list[k])), cy));    
    //inten[k] = (af::sum<double>((sin(dist_mat*q_list[k])/(dist_mat*q_list[k])*FF_device), cy));
    inten[k]+=FF_rank[k];
    val->set(inten[k]);
  }
*/
/*
  // test2 with thrust and cublas
  const unsigned nqsize = numq*size;

  double x_pos[size];
  double y_pos[size];
  double z_pos[size];

  
  for (unsigned i=0; i<size; i++) {
    const Vector posi = getPosition(i);
    x_pos[i] = posi[0];
    y_pos[i] = posi[1];
    z_pos[i] = posi[2];
  }

  //double* sum;
  double sum[numq];
  //sum = all_distance(x_pos, y_pos, z_pos, &FF_new[0], &q_list[0], size, numq);
  //double x = sort_test();
  //sum = test_boxderiv(x_pos, y_pos, z_pos, &FF_new[0], &q_list[0], size, numq);
  test_cublas(x_pos, y_pos, z_pos, &FF_new[0], &q_list[0], size, numq, sum);
  #pragma omp parallel for num_threads(OpenMP::getNumThreads())
  for (unsigned k=0; k<numq; k++) {
    Value* val=getPntrToComponent(k);
    sum[k]+=FF_rank[k];
    val->set(sum[k]);
  }
*/
/*
  vector<double> sum(numq);
  #pragma omp parallel for num_threads(OpenMP::getNumThreads())
  for (unsigned k=0; k<numq; k++) {
    Value* val=getPntrToComponent(k);
    vector<double> FFk = FF_value[k];
    sum[k]=all_distance(x_pos, y_pos, z_pos, &FFk[0], q_list[k], size);
    sum[k]+=FF_rank[k];
    val->set(sum[k]);
  }
*/

// old code
/*
  Vector c_distances;
  vector<double> m_distances(size*(size-1));
  vector<vector<double> >  FF;
  FF.resize(numq,vector<double>(size*(size-1)));


  for (unsigned i=0; i<size-1; i++) {
    const unsigned idx = i*(size-1);
    const Vector posi=getPosition(i);
    for (unsigned j=i+1; j<size ; j++) {
       unsigned index = idx+j;
       c_distances = delta(posi,getPosition(j));
       m_distances[index] = c_distances.modulo();
       for (unsigned k=0; k<numq; k++) {
         FF[k][index]=2.*FF_value[k][i]*FF_value[k][j];
       }
    }
  }

  vector<double> sum(numq);
  #pragma omp parallel for num_threads(OpenMP::getNumThreads())
  for (unsigned k=0; k<numq; k++) {
    Value* val=getPntrToComponent(k);
    vector<double> FFk = FF[k];
    sum[k]=SAXS_one_q(&m_distances[0], &FFk[0], q_list[k], size*(size-1));
    sum[k]+=FF_rank[k];
    val->set(sum[k]);
  }
*/
}

}
}
