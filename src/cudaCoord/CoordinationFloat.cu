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

//#define USE_CUDA_SPARSE_PRECISION CUDA_R_64F
#define USE_CUDA_SPARSE_PRECISION CUDA_R_32F
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

using calculateFloat=float;
//does not inherit from coordination base because nl is private
class CudaCoordination_f : public Colvar {
  std::unique_ptr<NeighborList> nl;
  ///the pointer to the coordinates on the GPU
  CUDAHELPERS::memoryHolder<calculateFloat>  cudaCoords;
  ///the pointer to the nn list on the GPU
  CUDAHELPERS::memoryHolder<unsigned>  cudaPairList;
  CUDAHELPERS::memoryHolder<calculateFloat> cudaCoordination;
  CUDAHELPERS::memoryHolder<calculateFloat> cudaDerivatives;
  CUDAHELPERS::memoryHolder<calculateFloat> cudaDerivatives_sparse;
  CUDAHELPERS::memoryHolder<int64_t> cudaDerivatives_sparserows;
  CUDAHELPERS::memoryHolder<int64_t> cudaDerivatives_sparsecols;
  CUDAHELPERS::memoryHolder<calculateFloat> cudaVirial;
  CUDAHELPERS::memoryHolder<calculateFloat> resultDerivatives;
  CUDAHELPERS::memoryHolder<calculateFloat> bufferDerivatives;
  CUDAHELPERS::memoryHolder<calculateFloat> reductionMemoryVirial;
  CUDAHELPERS::memoryHolder<calculateFloat> reductionMemoryCoord;
  
  cusparseHandle_t sparseMDevHandle;
  cusparseDnVecDescr_t outDevDescr;

  cudaStream_t streamDerivatives;
  cudaStream_t streamVirial;
  cudaStream_t streamCoordination;
  unsigned maxNumThreads=512;
  SwitchingFunction switchingFunction;
  rationalSwitchParameters switchingParameters;

  bool pbc{true};
  bool serial{false};
  bool invalidateList{true};
  bool firsttime{true};
  void setUpPermanentGPUMemory();
public:
  explicit CudaCoordination_f(const ActionOptions&);
  virtual ~CudaCoordination_f();
// active methods:
  static void registerKeywords( Keywords& keys );
  void prepare() override;
  void calculate() override;
};

PLUMED_REGISTER_ACTION(CudaCoordination_f,"CUDACOORDINATIONFLOAT")

void CudaCoordination_f::setUpPermanentGPUMemory(){
  auto nat = getPositions().size();
  //coordinates values are updated at each step
  if ((3*nat) != cudaCoords.size()){
    cudaCoords.resize(3*nat);
    if(resultDerivatives.size()!=0)
      cusparseDestroyDnVec(outDevDescr);
    resultDerivatives.resize(3*nat);
    cusparseCreateDnVec(&outDevDescr,
      3*nat,
      resultDerivatives.pointer(),
      USE_CUDA_SPARSE_PRECISION);
  }
  //the neighbour list will be updated at each request of prepare
  auto pairList = nl->getClosePairs();
  const unsigned nn=nl->size();
  cudaPairList.resize(2*nn);
  //streamDerivatives makes the copy async
  cudaPairList.copyToCuda(pairList.data(),streamDerivatives);
  resultDerivatives.resize(3*nat);
  cusparseCreateDnVec(&outDevDescr,
                    3*nat,
                    resultDerivatives.pointer(),
                    USE_CUDA_SPARSE_PRECISION);
}

void CudaCoordination_f::prepare() {
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
      if(getExchangeStep())
        error("Neighbor lists should be updated on exchange steps - choose a "
          "NL_STRIDE which divides the exchange stride!");
    }
    if(getExchangeStep())
      firsttime=true;
  }
}

void CudaCoordination_f::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords(keys);
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.addFlag("PAIR",false,"Pair only 1st element of the 1st group with 1st "
    "element in the second, etc");
  keys.addFlag("NLIST",false,"Use a neighbor list to speed up the calculation");
  keys.add("optional","NL_CUTOFF","The cutoff for the neighbor list");
  keys.add("optional","NL_STRIDE","The frequency with which we are updating the "
    "atoms in the neighbor list");
  keys.add("optional","THREADS","The upper limit of the number of threads");
  keys.add("atoms","GROUPA","First list of atoms");
  keys.add("atoms","GROUPB","Second list of atoms (if empty, N*(N-1)/2 pairs in "
    "GROUPA are counted)");
  keys.add("compulsory","NN","6","The n parameter of the switching function ");
  keys.add("compulsory","MM","0","The m parameter of the switching function; 0 implies 2*NN");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
}

//these constant will be used within the kernels
__constant__ calculateFloat cu_epsilon;

__device__ calculateFloat pcuda_fastpow(calculateFloat base,int expo) {
  if(expo<0) {
    expo=-expo;
    base=1.0/base;
  }
  calculateFloat result = 1.0;
  while (expo) {
    if (expo & 1)
      result *= base;
    expo >>= 1;
    base *= base;
  }
  return result;
}

__device__ calculateFloat pcuda_Rational(const calculateFloat rdist,int NN, int MM,calculateFloat&dfunc) {
  calculateFloat result;
  if(2*NN==MM) {
// if 2*N==M, then (1.0-rdist^N)/(1.0-rdist^M) = 1.0/(1.0+rdist^N)
    calculateFloat rNdist=pcuda_fastpow(rdist,NN-1);
    calculateFloat iden=1.0/(1+rNdist*rdist);
    dfunc = -NN*rNdist*iden*iden;
    result = iden;
  } else {
    if(rdist>(1.-100.0*cu_epsilon) && rdist<(1+100.0*cu_epsilon)) {
      result=NN/MM;
      dfunc=0.5*NN*(NN-MM)/MM;
    } else {
      calculateFloat rNdist=pcuda_fastpow(rdist,NN-1);
      calculateFloat rMdist=pcuda_fastpow(rdist,MM-1);
      calculateFloat num = 1.-rNdist*rdist;
      calculateFloat iden = 1.0/(1.0-rMdist*rdist);
      calculateFloat func = num*iden;
      result = func;
      dfunc = ((-NN*rNdist*iden)+(func*(iden*MM)*rMdist));
    }
  }
  return result;
}

__global__ void getpcuda_Rational(const calculateFloat *rdists,const int NN, const int MM,
    calculateFloat *dfunc,
    calculateFloat*res) {
  const int i = threadIdx.x + blockIdx.x * blockDim.x;
  if(rdists[i]<=0.) {
    res[i]=1.;
    dfunc[i]=0.0;
  }else
    res[i]=pcuda_Rational(rdists[i],NN,MM,dfunc[i]);
}

// __global__ void getConst() {
//   printf("Cuda: cu_epsilon = %f\n", cu_epsilon);
// }

CudaCoordination_f::CudaCoordination_f(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao) {
  parseFlag("SERIAL",serial);

  std::vector<AtomNumber> GroupA,GroupB;
  parseAtomList("GROUPA",GroupA);
  parseAtomList("GROUPB",GroupB);

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
    if(nl_cut<=0.0)
      error("NL_CUTOFF should be explicitly specified and positive");
    parse("NL_STRIDE",nl_st);
    if(nl_st<=0)
      error("NL_STRIDE should be explicitly specified and positive");
  }
  parse("THREADS",maxNumThreads);
  if(maxNumThreads<=0)
    error("THREADS should be positive");
  addValueWithDerivatives();
  setNotPeriodic();
  if(GroupB.size()>0) {
    if(doneigh){
      nl=Tools::make_unique<NeighborList>(
        GroupA,
        GroupB,
        serial,
        dopair,
        pbc,
        getPbc(),
        comm,
        nl_cut,
        nl_st);
    } else {
      nl=Tools::make_unique<NeighborList>(
        GroupA,
        GroupB,
        serial,
        dopair,
        pbc,
        getPbc(),
        comm);}
  } else {
    if(doneigh){
      nl=Tools::make_unique<NeighborList>(
      GroupA,
      serial,
      pbc,
      getPbc(),
      comm,
      nl_cut,
      nl_st);
    }else{
      nl=Tools::make_unique<NeighborList>(GroupA,serial,pbc,getPbc(),comm);
    }
  }

  requestAtoms(nl->getFullAtomList());

  log.printf(
    "  between two groups of %u and %u atoms\n",
    static_cast<unsigned>(GroupA.size()),
    static_cast<unsigned>(GroupB.size()));
  log.printf("  first group:\n");
  for(unsigned int i=0; i<GroupA.size(); ++i) {
    if ( (i+1) % 25 == 0 )
      log.printf("  \n");
    log.printf("  %d", GroupA[i].serial());
  }
  log.printf("  \n  second group:\n");
  for(unsigned int i=0; i<GroupB.size(); ++i) {
    if ( (i+1) % 25 == 0 )
      log.printf("  \n");
    log.printf("  %d", GroupB[i].serial());
  }
  log.printf("  \n");
  if(pbc)
    log.printf("  using periodic boundary conditions\n");
  else
    log.printf("  without periodic boundary conditions\n");
  if(dopair)
    log.printf("  with PAIR option\n");
  if(doneigh) {
    log.printf("  using neighbor lists with\n");
    log.printf("  update every %d steps and cutoff %f\n",nl_st,nl_cut);
  }
  std::string sw,errors;

  {//loading data to the GPU
    int nn_=6;
    int mm_=0;
    calculateFloat d0_=0.0;
    calculateFloat r0_=0.0;
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
    calculateFloat dmax=d0_+r0_*std::pow(0.00001,1./(nn_-mm_));
    constexpr bool dostretch=true;
    if (dostretch){
      std::vector<calculateFloat> inputs = {0.0,dmax};
      calculateFloat *inputsc,*dummy;
      calculateFloat *sc;
      cudaMalloc(&inputsc, 2 *sizeof(calculateFloat));
      cudaMalloc(&dummy, 2*sizeof(calculateFloat));
      cudaMalloc(&sc, 2*sizeof(calculateFloat));
      cudaMemcpy(inputsc, inputs.data(), 2* sizeof(calculateFloat),
                cudaMemcpyHostToDevice);
      getpcuda_Rational<<<1,2>>>(inputsc,nn_,mm_,dummy,sc);
      std::vector<calculateFloat> s = {0.0,0.0};
      cudaMemcpy(s.data(), sc, 2* sizeof(calculateFloat),
                cudaMemcpyDeviceToHost);
      cudaFree(inputsc);
      cudaFree(dummy);
      cudaFree(sc);
      switchingParameters.stretch=1.0/(s[0]-s[1]);
      switchingParameters.shift=-s[1]*switchingParameters.stretch;
    }
    
    cudaMemcpyToSymbol(cu_epsilon, &epsilon, sizeof(calculateFloat));
    switchingParameters.dmaxSQ = dmax* dmax;
    calculateFloat invr0 = 1.0/r0_;
    switchingParameters.invr0_2 = invr0*=invr0;
  }
  checkRead();
  cudaStreamCreate(&streamDerivatives);
  cudaStreamCreate(&streamVirial);
  cudaStreamCreate(&streamCoordination);
  cusparseCreate(&sparseMDevHandle);
  cusparseSetStream(sparseMDevHandle, streamDerivatives);
  setUpPermanentGPUMemory();
  log << "  contacts are counted with cutoff "
      << switchingFunction.description() << "\n";
}

CudaCoordination_f::~CudaCoordination_f(){
  cusparseDestroyDnVec(outDevDescr);
  cusparseDestroy(sparseMDevHandle);
  cudaStreamDestroy(streamDerivatives);
  cudaStreamDestroy(streamVirial);
  cudaStreamDestroy(streamCoordination);
}

__device__ calculateFloat calculateSqr(const calculateFloat distancesq,
    const rationalSwitchParameters switchingParameters,
    calculateFloat& dfunc) {
  calculateFloat result=0.0;
  dfunc=0.0;
  if(distancesq<switchingParameters.dmaxSQ) {
    const calculateFloat rdist_2 = distancesq*switchingParameters.invr0_2;
    result=pcuda_Rational(
      rdist_2,
      switchingParameters.nn/2,
      switchingParameters.mm/2,
      dfunc);
    // chain rule:
    dfunc*=2*switchingParameters.invr0_2;
    // cu_stretch:
    result=result*switchingParameters.stretch+switchingParameters.shift;
    dfunc*=switchingParameters.stretch;
  }
  return result;
}

#define X(I) 3*I
#define Y(I) 3*I+1
#define Z(I) 3*I+2

__global__ void getCoord(
                        const unsigned numOfPairs,
                        const rationalSwitchParameters switchingParameters,
                        const calculateFloat *coordinates,
                        const unsigned *pairList,
                        calculateFloat *ncoordOut,
                        calculateFloat *ddOut,
                        calculateFloat *ddOut_sparse,
                        int64_t *sparseRows,
                        int64_t *sparseCols
                        ) {
  //blockDIm are the number of threads in your block
  const unsigned i = threadIdx.x + blockIdx.x * blockDim.x;

  //Safeguard
  if (i >=numOfPairs)
    return;

  const unsigned i0= pairList[i*2];
  const unsigned i1= pairList[i*2+1]; 
  if (i0 == i1)
    return;

  calculateFloat d[3]={
    coordinates[X(i1)] - coordinates[X(i0)],
    coordinates[Y(i1)] - coordinates[Y(i0)],
    coordinates[Z(i1)] - coordinates[Z(i0)]
  };

  calculateFloat dsq=(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
  calculateFloat dfunc=0.;
  ncoordOut[i]= calculateSqr(dsq,switchingParameters, dfunc);
  const int sparsePlace = 6 * i;
  ddOut[i] = d[0];
  ddOut[i + 1 * numOfPairs] = d[1];
  ddOut[i + 2 * numOfPairs] = d[2];
  ddOut[i + 3 * numOfPairs] = dfunc;
  ddOut_sparse[0 + sparsePlace] = -d[0];// * dfunc;
  ddOut_sparse[1 + sparsePlace] = -d[1];// * dfunc;
  ddOut_sparse[2 + sparsePlace] = -d[2];// * dfunc;
  ddOut_sparse[3 + sparsePlace] = d[0] ;//* dfunc;
  ddOut_sparse[4 + sparsePlace] = d[1] ;//* dfunc;
  ddOut_sparse[5 + sparsePlace] = d[2] ;//* dfunc;
  sparseRows[0 + sparsePlace] = 3 * i0 + 0;
  sparseRows[1 + sparsePlace] = 3 * i0 + 1;
  sparseRows[2 + sparsePlace] = 3 * i0 + 2;
  sparseRows[3 + sparsePlace] = 3 * i1 + 0;
  sparseRows[4 + sparsePlace] = 3 * i1 + 1;
  sparseRows[5 + sparsePlace] = 3 * i1 + 2;
  sparseCols[i]=sparsePlace;
}

//In ndReduction.cu there are the more "general" reductions
  //(vector, tensor and scalar) here we have the adapted one for how the virial
  //information are stored
  template <unsigned numThreads, typename T>
  __global__ void reductionVirial(T *inputArray, T *outputArray, const unsigned int len) {
  //we saved the dd in an array x0 x1 x2..,xn-1,y0 y1 y2..,yn-1,z0 z1 z2..,zn-1
  //virialOut[ii*3+jj]-=d[ii]*d[jj]*dfunc;
  auto sdata = CUDAHELPERS::shared_memory_proxy<T>();
  const unsigned int ii = blockIdx.y;
  const unsigned int jj = blockIdx.z;
  const unsigned int vcoord = ii*3+jj;
  const unsigned int place = threadIdx.x;
      
  // each thread preprocess some element from global to shared memory
  const unsigned int diplacementI = ii*len;
  const unsigned int diplacementJ = jj*len;
  unsigned int i = (numThreads*2)*blockIdx.x + place + diplacementI;
  unsigned int j = (numThreads*2)*blockIdx.x + place + diplacementJ;
  unsigned int dfunc = (numThreads*2)*blockIdx.x + place + 3*len;
  const unsigned int gridSize = (numThreads*2)*gridDim.x;
  //the first element is in blockIdx.y*len, the last element to sum in (blockIdx.y+1)*len-1
  const unsigned int trgt=diplacementI + len;


  sdata[place] = T(0);
  
  while (i+numThreads < trgt) {
    sdata[place] -= inputArray[i]*inputArray[j]*inputArray[dfunc]
     + inputArray[i+numThreads]*inputArray[j+numThreads]*inputArray[dfunc+numThreads];
    i+=gridSize;
    j+=gridSize;
    dfunc+=gridSize;
  }
  while (i < trgt) {
    sdata[place] -= inputArray[i]*inputArray[j]*inputArray[dfunc];
    i+=gridSize;
    j+=gridSize;
    dfunc+=gridSize;
  }
    
  __syncthreads();
  // do the actual reduction, its contained in cudaHelpers.cuh
  CUDAHELPERS::reductor<numThreads>(
    sdata,
    outputArray,
    blockIdx.x + vcoord * gridDim.x
  );
}

template <typename T>
void doReductionVirial (T *inputArray, T *outputArray, const unsigned int len,
 const dim3 blocks, const unsigned nthreads,cudaStream_t stream=0){
  switch (nthreads) {
  case 512:
    reductionVirial<512,T><<<blocks,512,512*sizeof(T),stream>>>(inputArray,outputArray, len);
    break;
  case 256:
    reductionVirial<256,T><<<blocks,256,256*sizeof(T),stream>>>(inputArray,outputArray, len);
    break;
  case 128:
    reductionVirial<128,T><<<blocks,128,128*sizeof(T),stream>>>(inputArray,outputArray, len);
    break;
  case 64:
    reductionVirial<64, T><<<blocks,64,64*sizeof(T),stream>>>(inputArray,outputArray, len);
    break;
  case 32:
    reductionVirial<32, T><<<blocks,32,32*sizeof(T),stream>>>(inputArray,outputArray, len);
    break;
  default:
    plumed_merror("Reduction can be called only with 512, 256, 128, 64 or 32 threads.");
  }
}

void CudaCoordination_f::calculate() {
  auto positions = getPositions();
  auto nat = positions.size();

  Tensor virial;
  double coordination;
  auto deriv = std::vector<Vector>(nat);
  
  if(nl->getStride()>0 && invalidateList)
    nl->update(getPositions());
  
  auto pairList = nl->getClosePairs();
  const unsigned nn=nl->size();
  constexpr unsigned nthreads=512;
    
  unsigned ngroups=ceil(double(nn)/nthreads);

  /**********************allocating the memory on the GPU**********************/
  cudaCoords.copyToCuda(&positions[0][0],streamDerivatives);
  
  cudaCoordination.resize(nn);  
  cudaDerivatives.resize(nn *4);
  cudaDerivatives_sparse.resize(nn *6);
  cudaVirial.resize(nn*9);
  cudaDerivatives_sparserows.resize(nn*6);
  cudaDerivatives_sparsecols.resize(nn);
  /**************************starting the calculations*************************/
  //this initializes the memory to be accumulated
  getCoord<<<ngroups,nthreads,0,streamDerivatives>>> (
    nn,
    switchingParameters,
    cudaCoords.pointer(),
    cudaPairList.pointer(),
    cudaCoordination.pointer(),
    cudaDerivatives.pointer(),
    cudaDerivatives_sparse.pointer(),
    cudaDerivatives_sparserows.pointer(),
    cudaDerivatives_sparsecols.pointer()
  );
  //since the virial and the coordination reduction are on different streams, 
  //this barrier makes sure that the memory is ready for the operations
  
  /**************************accumulating the results**************************/
  cusparseDnVecDescr_t dfuncDense;
  cusparseCreateDnVec(
    &dfuncDense,
    nn,
    //this part of the array contains the nn dfuncs
    cudaDerivatives.pointer()+3*nn,
    USE_CUDA_SPARSE_PRECISION);
  cusparseSpMatDescr_t derivativesSparse;
  cusparseCreateCsc(
    &derivativesSparse,
    getPositions().size() * 3,//rows
    nn,//columns
    nn*6,//non zerp elements
    cudaDerivatives_sparsecols.pointer(),
    cudaDerivatives_sparserows.pointer(),
    cudaDerivatives_sparse.pointer(),
    CUSPARSE_INDEX_64I,
    CUSPARSE_INDEX_64I,
    CUSPARSE_INDEX_BASE_ZERO,
    USE_CUDA_SPARSE_PRECISION);
  calculateFloat one=1.0;
  calculateFloat zero=0.0;
  size_t bufferSize = 0;
  //this computes the buffersize
  cusparseSpMV_bufferSize(
    sparseMDevHandle,
    CUSPARSE_OPERATION_NON_TRANSPOSE,
    &one,
    derivativesSparse,
    dfuncDense,
    &zero,
    outDevDescr,
    USE_CUDA_SPARSE_PRECISION,
    CUSPARSE_SPMV_ALG_DEFAULT,//may be improved
    &bufferSize);
  bufferDerivatives.resize(bufferSize);
  cudaStreamSynchronize ( streamDerivatives );
  cusparseSpMV(
    sparseMDevHandle,
    CUSPARSE_OPERATION_NON_TRANSPOSE,
    &one,
    derivativesSparse,
    dfuncDense,
    &zero,
    outDevDescr,
    USE_CUDA_SPARSE_PRECISION,
    CUSPARSE_SPMV_ALG_DEFAULT,//may be improved
    bufferDerivatives.pointer());
  
  auto N=nn;
  bool first=true;
  while(N>1){
    size_t runningThreads = CUDAHELPERS::threadsPerBlock(N,maxNumThreads);
    unsigned nGroups = CUDAHELPERS::idealGroups(N, runningThreads);
    
    
    reductionMemoryCoord.resize(nGroups);
    reductionMemoryVirial.resize(9* nGroups);

    if (first){
      dim3 ngroupsVirial(nGroups,3,3);
      doReductionVirial (
        cudaDerivatives.pointer(),
        reductionMemoryVirial.pointer(),
        N,
        ngroupsVirial,
        runningThreads,
        streamVirial);
        //putting this here improves the overlap between mem movements and kernels
        resultDerivatives.copyFromCuda(&deriv[0][0],streamDerivatives);
    }else{
      dim3 ngroupsVirial(nGroups,9);
      CUDAHELPERS::doReductionND (
        cudaVirial.pointer(),
        reductionMemoryVirial.pointer(),
        N,
        ngroupsVirial,
        runningThreads,
        streamVirial);
    }

    CUDAHELPERS::doReduction1D (
      cudaCoordination.pointer(),//reduceScalarIn->pointer(),
      reductionMemoryCoord.pointer(),//reduceSOut->pointer(),
      N,
      nGroups,
      runningThreads,
      streamCoordination);
  
    if (nGroups==1){
      reductionMemoryVirial.copyFromCuda(&virial[0][0],streamVirial);
      //reduceSOut->copyFromCuda(&coordination,streamCoordination);
      reductionMemoryCoord.copyFromCuda(&coordination,streamCoordination);
    } else {
      //std::swap(reduceScalarIn,reduceSOut);
      reductionMemoryCoord.swap(cudaCoordination);
      reductionMemoryVirial.swap(cudaVirial);
    }
    first=false;
    N=nGroups;
  }
  //this ensures that the memory is fully in the host ram
  cudaDeviceSynchronize ();
  cusparseDestroySpMat(derivativesSparse);
  cusparseDestroyDnVec(dfuncDense);
  for(unsigned i=0; i<deriv.size(); ++i) {
    setAtomsDerivatives(i,deriv[i]);
  }
  
  setValue           (coordination);
  setBoxDerivatives  (virial);
}

} // namespace colvar
} // namespace PLMD
