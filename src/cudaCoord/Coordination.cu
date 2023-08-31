/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#include "plumed/colvar/CoordinationBase.h"
#include "plumed/tools/SwitchingFunction.h"
#include "plumed/tools/NeighborList.h"
#include "plumed/core/ActionRegister.h"

#include <cusparse_v2.h>
#include "ndReduction.h"
#include "cudaHelpers.cuh"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <numeric>
#include <limits>
#include <iostream>

using std::cerr;

#define vdbg(...) std::cerr << __LINE__ <<":" << #__VA_ARGS__ << " " << (__VA_ARGS__) <<'\n'
//#define vdbg(...)

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR CUDACOORDINATION
/*
Calculate coordination numbers. Like coordination, but on nvdia gpu and with no swithcing function

This keyword can be used to calculate the number of contacts between two groups of atoms
and is defined as
\f[
\sum_{i\in A} \sum_{i\in B} s_{ij}
\f]
where \f$s_{ij}\f$ is 1 if the contact between atoms \f$i\f$ and \f$j\f$ is formed,
zero otherwise.
In actuality, \f$s_{ij}\f$ is replaced with a switching function so as to ensure that the calculated CV has continuous derivatives.
The default switching function is:
\f[
s_{ij} = \frac{ 1 - \left(\frac{{\bf r}_{ij}-d_0}{r_0}\right)^n } { 1 - \left(\frac{{\bf r}_{ij}-d_0}{r_0}\right)^m }
\f]
but it can be changed using the optional SWITCH option.

To make your calculation faster you can use a neighbor list, which makes it that only a
relevant subset of the pairwise distance are calculated at every step.

If GROUPB is empty, it will sum the \f$\frac{N(N-1)}{2}\f$ pairs in GROUPA. This avoids computing
twice permuted indexes (e.g. pair (i,j) and (j,i)) thus running at twice the speed.

Notice that if there are common atoms between GROUPA and GROUPB the switching function should be
equal to one. These "self contacts" are discarded by plumed (since version 2.1),
so that they actually count as "zero".


\par Examples

The following example instructs plumed to calculate the total coordination number of the atoms in group 1-10 with the atoms in group 20-100.  For atoms 1-10 coordination numbers are calculated that count the number of atoms from the second group that are within 0.3 nm of the central atom.  A neighbor list is used to make this calculation faster, this neighbor list is updated every 100 steps.
\plumedfile
COORDINATION GROUPA=1-10 GROUPB=20-100 R_0=0.3 NLIST NL_CUTOFF=0.5 NL_STRIDE=100
\endplumedfile

The following is a dummy example which should compute the value 0 because the self interaction
of atom 1 is skipped. Notice that in plumed 2.0 "self interactions" were not skipped, and the
same calculation should return 1.
\plumedfile
c: COORDINATION GROUPA=1 GROUPB=1 R_0=0.3
PRINT ARG=c STRIDE=10
\endplumedfile

Here's an example that shows what happens when providing COORDINATION with
a single group:
\plumedfile
# define some huge group:
group: GROUP ATOMS=1-1000
# Here's coordination of a group against itself:
c1: COORDINATION GROUPA=group GROUPB=group R_0=0.3
# Here's coordination within a single group:
x: COORDINATION GROUPA=group R_0=0.3
# This is just multiplying times 2 the variable x:
c2: COMBINE ARG=x COEFFICIENTS=2 PERIODIC=NO

# the two variables c1 and c2 should be identical, but the calculation of c2 is twice faster
# since it runs on half of the pairs.
PRINT ARG=c1,c2 STRIDE=10
\endplumedfile



*/
//+ENDPLUMEDOC


//these constant will be used within the kernels
struct rationalSwitchParameters{
  double dmaxSQ=std::numeric_limits<double>::max();
  double invr0_2=1.0;//r0=1
  double stretch=1.0;
  double shift=0.0;
  int nn=6;
  int mm=12;
};

//does not inherit from coordination base because nl is private
class CudaCoordination : public Colvar {
  std::unique_ptr<NeighborList> nl;
  ///the pointer to the coordinates on the GPU
  CUDAHELPERS::memoryHolder<double>  cudaCoords;
  ///the pointer to the nn list on the GPU
  CUDAHELPERS::memoryHolder<unsigned>  cudaPairList;
  CUDAHELPERS::memoryHolder<double> cudaCoordination;
  CUDAHELPERS::memoryHolder<double> cudaDerivatives;
  CUDAHELPERS::memoryHolder<double> cudaDerivatives_sparse;
  CUDAHELPERS::memoryHolder<int64_t> cudaDerivatives_sparserow;
  CUDAHELPERS::memoryHolder<int64_t> cudaDerivatives_sparsecols;
  CUDAHELPERS::memoryHolder<double> cudaVirial;
  CUDAHELPERS::memoryHolder<double> reductionMemoryDerivatives;
  CUDAHELPERS::memoryHolder<double> bufferDerivatives;
  CUDAHELPERS::memoryHolder<double> reductionMemoryVirial;
  CUDAHELPERS::memoryHolder<double> reductionMemoryCoord;
  CUDAHELPERS::memoryHolder<double> ones;
  cusparseHandle_t sparseMDevHandle;
  cudaStream_t streamDerivatives;
  cudaStream_t streamVirial;
  cudaStream_t streamCoordination;
  cusparseDnVecDescr_t outVecDescr;
  unsigned maxNumThreads=512;
  SwitchingFunction switchingFunction;
  rationalSwitchParameters switchingParameters;

  bool pbc{true};
  bool serial{false};
  bool invalidateList{true};
  bool firsttime{true};
  void setUpPermanentGPUMemory();
public:
  explicit CudaCoordination(const ActionOptions&);
  virtual ~CudaCoordination();
// active methods:
  static void registerKeywords( Keywords& keys );
  void prepare() override;
  void calculate() override;
};

PLUMED_REGISTER_ACTION(CudaCoordination,"CUDACOORDINATION")

void CudaCoordination::setUpPermanentGPUMemory(){
  auto nat = getPositions().size();
  //coordinates values are updated at each step
  cudaCoords.resize(3*nat);
  
  //the neighbour list will be updated at each request of prepare
  auto pairList = nl->getClosePairs();
  const unsigned nn=nl->size();
  cudaPairList.resize(2*nn);
  cudaPairList.copyToCuda(pairList.data());
  reductionMemoryDerivatives.resize(3*nat);
  cusparseCreateDnVec(&outVecDescr,
                    3*nat,
                    reductionMemoryDerivatives.pointer(),
                    CUDA_R_64F);
}

void CudaCoordination::prepare() {
  if(nl->getStride()>0) {
    if(firsttime || (getStep()%nl->getStride()==0)) {
      requestAtoms(nl->getFullAtomList());
      setUpPermanentGPUMemory();
      invalidateList=true;
      firsttime=false;
    } else {
      requestAtoms(nl->getReducedAtomList());
      setUpPermanentGPUMemory();
      invalidateList=false;
      if(getExchangeStep()) error("Neighbor lists should be updated on exchange steps - choose a NL_STRIDE which divides the exchange stride!");
    }
    if(getExchangeStep()) firsttime=true;
  }
}
void CudaCoordination::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords(keys);
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.addFlag("PAIR",false,"Pair only 1st element of the 1st group with 1st element in the second, etc");
  keys.addFlag("NLIST",false,"Use a neighbor list to speed up the calculation");
  keys.add("optional","NL_CUTOFF","The cutoff for the neighbor list");
  keys.add("optional","NL_STRIDE","The frequency with which we are updating the atoms in the neighbor list");
  keys.add("optional","THREADS","The upper limit of the number of threads");
  keys.add("atoms","GROUPA","First list of atoms");
  keys.add("atoms","GROUPB","Second list of atoms (if empty, N*(N-1)/2 pairs in GROUPA are counted)");
  keys.add("compulsory","NN","6","The n parameter of the switching function ");
  keys.add("compulsory","MM","0","The m parameter of the switching function; 0 implies 2*NN");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
}

//these constant will be used within the kernels
__constant__ double cu_epsilon;

__device__ double pcuda_fastpow(double base, int expo) {
  if(expo<0) {
    expo=-expo;
    base=1.0/base;
  }
  double result = 1.0;
  while (expo) {
    if (expo & 1) {
      result *= base;
    }
    expo >>= 1;
    base *= base;
  }
  return result;
}

__device__ double pcuda_Rational(const double rdist,int NN, int MM,double&dfunc) {
  double result;
  if(2*NN==MM) {
// if 2*N==M, then (1.0-rdist^N)/(1.0-rdist^M) = 1.0/(1.0+rdist^N)
    double rNdist=pcuda_fastpow(rdist,NN-1);
    double iden=1.0/(1+rNdist*rdist);
    dfunc = -NN*rNdist*iden*iden;
    result = iden;
  } else {
    if(rdist>(1.-100.0*cu_epsilon) && rdist<(1+100.0*cu_epsilon)) {
      result=NN/MM;
      dfunc=0.5*NN*(NN-MM)/MM;
    } else {
      double rNdist=pcuda_fastpow(rdist,NN-1);
      double rMdist=pcuda_fastpow(rdist,MM-1);
      double num = 1.-rNdist*rdist;
      double iden = 1.0/(1.0-rMdist*rdist);
      double func = num*iden;
      result = func;
      dfunc = ((-NN*rNdist*iden)+(func*(iden*MM)*rMdist));
    }
  }
  return result;
}

__global__ void getpcuda_Rational(const double *rdists,const int NN, const int MM,
    double *dfunc,
    double*res) {
  const int i = threadIdx.x + blockIdx.x * blockDim.x;
  if(rdists[i]<=0.) {
    res[i]=1.;
    dfunc[i]=0.0;
  }else{
  res[i]=pcuda_Rational(rdists[i],NN,MM,dfunc[i]);
  }
  //printf("CUDA: %i :: d=%f -> %f, %f\n", i,rdists[i],res[i],dfunc[i]);
}


__global__ void getConst() {
  //printf("Cuda: cu_epsilon = %f\n", cu_epsilon);
}

CudaCoordination::CudaCoordination(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao) {
  parseFlag("SERIAL",serial);

  std::vector<AtomNumber> ga_lista,gb_lista;
  parseAtomList("GROUPA",ga_lista);
  parseAtomList("GROUPB",gb_lista);

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

// pair stuff
  bool dopair=false;
  parseFlag("PAIR",dopair);

// neighbor list stuff
  bool doneigh=false;
  double nl_cut=0.0;
  int nl_st=0;
  parseFlag("NLIST",doneigh);
  if(doneigh) {
    parse("NL_CUTOFF",nl_cut);
    if(nl_cut<=0.0) error("NL_CUTOFF should be explicitly specified and positive");
    parse("NL_STRIDE",nl_st);
    if(nl_st<=0) error("NL_STRIDE should be explicitly specified and positive");
  }
  parse("THREADS",maxNumThreads);
  if(maxNumThreads<=0) error("THREADS should be positive");
  addValueWithDerivatives(); setNotPeriodic();
  if(gb_lista.size()>0) {
    if(doneigh)
    {nl=Tools::make_unique<NeighborList>(ga_lista,gb_lista,serial,dopair,pbc,getPbc(),comm,nl_cut,nl_st);}
    else
    {nl=Tools::make_unique<NeighborList>(ga_lista,gb_lista,serial,dopair,pbc,getPbc(),comm);}
  } else {
    if(doneigh)
    {nl=Tools::make_unique<NeighborList>(ga_lista,serial,pbc,getPbc(),comm,nl_cut,nl_st);}
    else
    {nl=Tools::make_unique<NeighborList>(ga_lista,serial,pbc,getPbc(),comm);}
  }

  requestAtoms(nl->getFullAtomList());

  log.printf("  between two groups of %u and %u atoms\n",static_cast<unsigned>(ga_lista.size()),static_cast<unsigned>(gb_lista.size()));
  log.printf("  first group:\n");
  for(unsigned int i=0; i<ga_lista.size(); ++i) {
    if ( (i+1) % 25 == 0 ) log.printf("  \n");
    log.printf("  %d", ga_lista[i].serial());
  }
  log.printf("  \n  second group:\n");
  for(unsigned int i=0; i<gb_lista.size(); ++i) {
    if ( (i+1) % 25 == 0 ) log.printf("  \n");
    log.printf("  %d", gb_lista[i].serial());
  }
  log.printf("  \n");
  if(pbc) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");
  if(dopair) log.printf("  with PAIR option\n");
  if(doneigh) {
    log.printf("  using neighbor lists with\n");
    log.printf("  update every %d steps and cutoff %f\n",nl_st,nl_cut);
  }
  std::string sw,errors;

  {//loading data to the GPU
    int nn_=6;
    int mm_=0;
    double d0_=0.0;
    double r0_=0.0;
    parse("R_0",r0_);
    if(r0_<=0.0) {
      error("R_0 should be explicitly specified and positive");
    }
    parse("D_0",d0_);
    parse("NN",nn_);
    parse("MM",mm_);
    if(mm_==0) {
      mm_=2*nn_;
      }
      
    switchingParameters.nn=nn_;
    switchingParameters.mm=mm_;
    switchingParameters.stretch=1.0;
    switchingParameters.shift=0.0;
    double dmax=d0_+r0_*std::pow(0.00001,1./(nn_-mm_));
    constexpr bool dostretch=true;
    if (dostretch){
      std::vector<double> inputs = {0.0,dmax};
      double *inputsc,*dummy;
      double *sc;
      cudaMalloc(&inputsc, 2 *sizeof(double));
      cudaMalloc(&dummy, 2*sizeof(double));
      cudaMalloc(&sc, 2*sizeof(double));
      cudaMemcpy(inputsc, inputs.data(), 2* sizeof(double),
                cudaMemcpyHostToDevice);
      getpcuda_Rational<<<1,2>>>(inputsc,nn_,mm_,dummy,sc);
      std::vector<double> s = {0.0,0.0};
      cudaMemcpy(s.data(), sc, 2* sizeof(double),
                cudaMemcpyDeviceToHost);
      cudaFree(inputsc);
      cudaFree(dummy);
      cudaFree(sc);
      switchingParameters.stretch=1.0/(s[0]-s[1]);
      switchingParameters.shift=-s[1]*switchingParameters.stretch;
    }
    
    cudaMemcpyToSymbol(cu_epsilon, &epsilon, sizeof(double));
    switchingParameters.dmaxSQ = dmax* dmax;
    double invr0 = 1.0/r0_;
    switchingParameters.invr0_2 = invr0*=invr0;
  }
  checkRead();
  cudaStreamCreate(&streamDerivatives);
  cudaStreamCreate(&streamVirial);
  cudaStreamCreate(&streamCoordination);
  cusparseCreate(&sparseMDevHandle);
  //cusparseSetPointerMode(sparseMDevHandle, CUSPARSE_POINTER_MODE_HOST or CUSPARSE_POINTER_MODE_DEVICE);
  cusparseSetStream(sparseMDevHandle, streamDerivatives);
  setUpPermanentGPUMemory();
  log<<"  contacts are counted with cutoff "<<switchingFunction.description()<<"\n";
}

CudaCoordination::~CudaCoordination(){
   cudaStreamDestroy(streamDerivatives);
   cudaStreamDestroy(streamVirial);
   cudaStreamDestroy(streamCoordination);
   cusparseDestroyDnVec(outVecDescr);
   cusparseDestroy(sparseMDevHandle);
}

__device__ double calculateSqr(const double distancesq,
    const rationalSwitchParameters switchingParameters,
    double& dfunc) {
  double result=0.0;
  dfunc=0.0;
  if(distancesq<switchingParameters.dmaxSQ) {
    const double rdist_2 = distancesq*switchingParameters.invr0_2;
    result=pcuda_Rational(rdist_2,switchingParameters.nn/2,switchingParameters.mm/2,dfunc);
    // chain rule:
    dfunc*=2*switchingParameters.invr0_2;
    // cu_stretch:
    result=result*switchingParameters.stretch+switchingParameters.shift;
    dfunc*=switchingParameters.stretch;
  }
  //printf("%f\n",result);
  return result;
}

#define X(I) 3*I
#define Y(I) 3*I+1
#define Z(I) 3*I+2

#define Xd(I) i+ (   I*3 )* numOfPairs
#define Yd(I) i+ (1 +I*3 )* numOfPairs
#define Zd(I) i+ (2 +I*3 )* numOfPairs
__global__ void getCoord(
                        const unsigned numOfPairs,
                        const rationalSwitchParameters switchingParameters,
                        const double *coordinates,
                        const unsigned *pairList,
                        double *ncoordOut,
                        double *ddOut,//will be used for the virial
                        double *ddOut_sparse,//will be used for the derivatives
                        int64_t *sparseRows,
                        int64_t *sparseCols,
                        double* ones
                        ) {
  //blockDIm are the number of threads in your block
  const unsigned i = threadIdx.x + blockIdx.x * blockDim.x;
  
  if (i >=numOfPairs) {
    //Safeguard
    return;
  }
  const unsigned i0= pairList[i*2];
  const unsigned i1= pairList[i*2+1]; 
  if (i0 == i1) {
    return;
  }
  double d[3]={
    coordinates[X(i1)] - coordinates[X(i0)],
    coordinates[Y(i1)] - coordinates[Y(i0)],
    coordinates[Z(i1)] - coordinates[Z(i0)]
  };

  double dsq=(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
  double dfunc=0.;
  ncoordOut[i]= calculateSqr(dsq,switchingParameters, dfunc);

  ddOut[i] = d[0];
  ddOut[i + 1 * numOfPairs] = d[1];
  ddOut[i + 2 * numOfPairs] = d[2];
  ddOut[i + 3 * numOfPairs] = dfunc;
  ddOut_sparse[0 + 6 * i] = -d[0] * dfunc;
  ddOut_sparse[1 + 6 * i] = -d[1] * dfunc;
  ddOut_sparse[2 + 6 * i] = -d[2] * dfunc;
  ddOut_sparse[3 + 6 * i] = d[0] * dfunc;
  ddOut_sparse[4 + 6 * i] = d[1] * dfunc;
  ddOut_sparse[5 + 6 * i] = d[2] * dfunc;
  sparseRows[0 + 6 * i] = 3*i0 + 0;
  sparseRows[1 + 6 * i] = 3*i0 + 1;
  sparseRows[2 + 6 * i] = 3*i0 + 2;
  sparseRows[3 + 6 * i] = 3*i1 + 0;
  sparseRows[4 + 6 * i] = 3*i1 + 1;
  sparseRows[5 + 6 * i] = 3*i1 + 2;
  ones[i] = 1.0;
  sparseCols[i]=6*i;
  // ncoord[i]= 1;
  // printf("Cuda: %i,%i %i %i\n", i,threadIdx.x , blockIdx.x,  blockDim.x);
}

void CudaCoordination::calculate() {
  
  auto positions = getPositions();
  auto nat = positions.size();
  if(nl->getStride()>0 && invalidateList) {
    nl->update(getPositions());
  }
  auto pairList = nl->getClosePairs();
  const unsigned nn=nl->size();
  
  constexpr unsigned nthreads=512;
    
  unsigned ngroups=ceil(double(nn)/nthreads);

  /****************allocating the memory on the GPU****************/
  cudaCoords.copyToCuda(&positions[0][0]);
  
  ones.resize(nn);
  cudaCoordination.resize(nn);  
  cudaDerivatives.resize(nn * 4);
  //cudaVirial.resize(nn * 9);//it will be eventually updated in the reduceDVS function

  cudaDerivatives_sparse.resize(nn * 6);
  cudaDerivatives_sparserow.resize(nn * 6);
  cudaDerivatives_sparsecols.resize(nn);
  
  /****************starting the calculations****************/
  getCoord<<<ngroups,nthreads>>> (nn,switchingParameters,
    cudaCoords.pointer(),
    cudaPairList.pointer(),
    cudaCoordination.pointer(),
    cudaDerivatives.pointer(),
    cudaDerivatives_sparse.pointer(),
    cudaDerivatives_sparserow.pointer(),
    cudaDerivatives_sparsecols.pointer(),
    ones.pointer()
    );
  cusparseSpMatDescr_t spMatDescr;
  cusparseDnVecDescr_t dnVecDescr;
  cusparseCreateDnVec(&dnVecDescr,
                  nn,
                  ones.pointer(),
                  CUDA_R_64F);
  
  cusparseCreateCsc(&spMatDescr,
                  getPositions().size() * 3,
                  nn,
                  nn*6,
                  cudaDerivatives_sparsecols.pointer(),
                  cudaDerivatives_sparserow.pointer(),
                  cudaDerivatives_sparse.pointer(),
                  CUSPARSE_INDEX_64I,
                  CUSPARSE_INDEX_64I,
                  CUSPARSE_INDEX_BASE_ZERO,
                  CUDA_R_64F);
  double one=1.0;
  double zero=0.0;
  size_t bufferSize = 0;
  cusparseSpMV_bufferSize( sparseMDevHandle,
                        CUSPARSE_OPERATION_NON_TRANSPOSE,
                        &one,
                        spMatDescr,
                        dnVecDescr,
                        &zero,
                        outVecDescr,
                        CUDA_R_64F,
                        CUSPARSE_SPMV_ALG_DEFAULT,//may be improved
                        &bufferSize);
bufferDerivatives.resize(bufferSize);
cusparseSpMV( sparseMDevHandle,
                        CUSPARSE_OPERATION_NON_TRANSPOSE,
                        &one,
                        spMatDescr,
                        dnVecDescr,
                        &zero,
                        outVecDescr,
                        CUDA_R_64F,
                        CUSPARSE_SPMV_ALG_DEFAULT,//may be improved
                        bufferDerivatives.pointer());
  CUDAHELPERS::DVS ret = CUDAHELPERS::reduceDVS(
    cudaDerivatives,
    cudaVirial,
    cudaCoordination,
    cudaPairList,
    reductionMemoryVirial,
    reductionMemoryCoord,
    streamVirial,
    streamCoordination,
    nn,nat,
    maxNumThreads
  );
  
  //cusparseDnVecGetValues(outVecDescr,**void); returns reductionMemoryDerivatives.pointer()
   
  reductionMemoryDerivatives.copyFromCuda(&ret.deriv[0][0]);
  cusparseDestroySpMat(spMatDescr);
  cusparseDestroyDnVec(dnVecDescr);
  
  for(unsigned i=0; i<ret.deriv.size(); ++i) {
    setAtomsDerivatives(i,ret.deriv[i]);
  }
  
  setValue           (ret.scalar);
  setBoxDerivatives  (ret.virial);
}

} // namespace colvar
} // namespace PLMD
