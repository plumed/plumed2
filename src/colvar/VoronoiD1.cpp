/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
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
#include "tools/NeighborList.h"
#include "tools/Communicator.h"
#include "tools/OpenMP.h"
#include "Colvar.h"

#include "tools/Matrix.h"
#include "core/ActionRegister.h"
#include <string>
#include <cmath>
#include <iostream>
using namespace std;

namespace PLMD {
namespace colvar {

class VoronoiD1 : public Colvar {
  bool pbc;
  bool serial;
  std::unique_ptr<NeighborList> nl;
  std::vector<PLMD::AtomNumber> list_a,list_b,list_c;
  std::vector<PLMD::AtomNumber> atomsToRequest;
  bool invalidateList;
  bool firsttime;
  double  lambda;
  int  nrx,num_atomsa,num_atomsb,num_atoms,num_atomso;
  double d0, d1, d2, d3, r0, sum_exp;

public:
  explicit VoronoiD1(const ActionOptions&);
  ~VoronoiD1();
// active methods:
  void calculate() override;
  void prepare() override;
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(VoronoiD1,"VORONOID1")

void VoronoiD1::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords(keys);
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.addFlag("PAIR",false,"Pair only 1st element of the 1st group with 1st element in the second, etc");
  keys.addFlag("NLIST",false,"Use a neighbor list to speed up the calculation");
  keys.add("optional","NL_CUTOFF","The cutoff for the neighbor list");
  keys.add("optional","NL_STRIDE","The frequency with which we are updating the atoms in the neighbor list");
  keys.add("atoms","GROUPA","First list of atoms");
  keys.add("atoms","GROUPB","Second list of atoms (if empty, N*(N-1)/2 pairs in GROUPA are counted)");
  keys.add("compulsory","LAMBDA","1","The lambda parameter of the sum_exp function; 0 implies 1");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","D_1","0.0","The d_1 parameter of the switching function");
  keys.add("compulsory","D_2","0.0","The d_2 parameter of the switching function");
  keys.add("compulsory","D_3","0.0","The d_3 parameter of the switching function");
  keys.add("compulsory","NRX","0.0","The number of reactive sites");
  //keys.add("compulsory","NN_THETA","1.0","The number of reactive sites");
}

VoronoiD1::VoronoiD1(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc(true),
  serial(false),
  invalidateList(true),
  firsttime(true)
{

  parseFlag("SERIAL",serial);

  std::vector<AtomNumber> ga_lista,gb_lista;
  parseAtomList("GROUPA",ga_lista);
  parseAtomList("GROUPB",gb_lista);

  list_a = ga_lista;
  list_b = gb_lista;

  num_atomsa = list_a.size();
  num_atomsb = list_b.size();
  num_atoms = num_atomsa + num_atomsb;

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  parse("D_0",d0);
  parse("D_1",d1);
  parse("D_2",d2);
  parse("D_3",d3);
  parse("NRX",nrx);
  //parse("NN_THETA",nn_theta);
  parse("LAMBDA",lambda);
  num_atomso = num_atomsa - nrx;

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

  addValueWithDerivatives(); setNotPeriodic();
  if(gb_lista.size()>0) {
    if(doneigh)  nl=Tools::make_unique<NeighborList>(ga_lista,gb_lista,serial,dopair,pbc,getPbc(),comm,nl_cut,nl_st);
    else         nl=Tools::make_unique<NeighborList>(ga_lista,gb_lista,serial,dopair,pbc,getPbc(),comm);
  } else {
    if(doneigh)  nl=Tools::make_unique<NeighborList>(ga_lista,serial,pbc,getPbc(),comm,nl_cut,nl_st);
    else         nl=Tools::make_unique<NeighborList>(ga_lista,serial,pbc,getPbc(),comm);
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
}

VoronoiD1::~VoronoiD1() {
// destructor required to delete forward declared class
}

void VoronoiD1::prepare() {
  if(nl->getStride()>0) {
    if(firsttime || (getStep()%nl->getStride()==0)) {
      requestAtoms(nl->getFullAtomList());
      invalidateList=true;
      firsttime=false;
    } else {
      //requestAtoms(nl->getReducedAtomList());
      requestAtoms(nl->getFullAtomList()); //I had to take the FULL list all the time, otherwise the different order in the local list of atoms create a lot of trouble later. Slightly more inefficient but unavoidable
      invalidateList=false;
      if(getExchangeStep()) error("Neighbor lists should be updated on exchange steps - choose a NL_STRIDE which divides the exchange stride!");
    }
    if(getExchangeStep()) firsttime=true;
  }
}

// calculator
void VoronoiD1::calculate()
{

  double IonDistance=0.0;
  vector<double> sum_exp(getNumberOfAtoms());
  fill(sum_exp.begin(),sum_exp.end(),0.);

  //getNumberOfAtoms() gives the total number of atoms requested. It depends on Nlist. Its index is LOCAL, not absolute!!! Every stride it is ALL atoms requested and the order is different.
  //the link between local and absolute is getAbsoluteIndexes()[k]. It is a LOCAL list, where the k local element indicates its ABSOLUTE index !!!!!!!!!!
  //nn below is different, it's the list of couples of atoms within range

  Tensor virial;
  vector<Vector> deriv(getNumberOfAtoms());
  Vector zeros;
  zeros.zero();
  fill(deriv.begin(), deriv.end(), zeros);


  if(nl->getStride()>0 && invalidateList) {
    nl->update(getPositions());
  }

  unsigned stride;
  unsigned rank;
  if(serial) {
    stride=1;
    rank=0;
  } else {
    stride=comm.Get_size();
    rank=comm.Get_rank();
  }

  unsigned nt=OpenMP::getNumThreads();
  const unsigned nn=nl->size(); //it is the list of neighbors!
  //first, atom0 in group A and all its neighbours from groupB, then another from A and so on. It is a way to put all couples in a 1D list
  if(nt*stride*10>nn) nt=1;

  //NEW INSTRUCTIONS!!!!!!
  vector<double> nnexp(nn);
  vector<unsigned> nni0(nn);
  vector<unsigned> nni1(nn);
  vector<double> nnexpnorm(getNumberOfAtoms()); //normalization of the colvar, local list
  vector<vector<double>> c(getNumberOfAtoms(),vector<double>(getNumberOfAtoms()));  //store exponentials in local indeces
  vector<double> charge(num_atomsa); //if I request reducedatomlist, it is a lot of trouble here because the local index of atoms in groupa changes in the update step
  vector<vector<double>> distGA(num_atomsa,vector<double>(num_atomsa));  //vector of the distances in groupA
  vector<vector<Vector>> distAB(num_atomsa,vector<Vector>(getNumberOfAtoms())); //I store all distances between groupA and B so not to compute them again in the derivative loop
  vector<vector<double>> distABinvmod(num_atomsa,vector<double>(getNumberOfAtoms())); //the inverse modulus of distAB, it saves time to store it and read it in the double loop

  // --- precomputed quantities to make LOOP4 O(nn) instead of O(nn*num_atomsa^2) ---
  vector<double> A(num_atomsa, 0.0);          // A(m) = sum_{n != m, both water-O} distGA(m,n)*charge(n)
  vector<double> abar(getNumberOfAtoms(), 0.0); // abar(h) = sum_o c[o][h]*A[o]

  #pragma omp parallel num_threads(nt)
  {
    std::vector<Vector> omp_deriv(getPositions().size());  //defined locally so that each thread will modify it independently. Later all will gather into deriv
    Tensor omp_virial;                                      //NEW: thread-local virial accumulator
    vector<double> nnexpnormt(getNumberOfAtoms()); //local version
    vector<double> charget(num_atomsa);  //local version
    //LOOP1 compute the distances, the exponentials and the normalizations over groupb elements
    #pragma omp for
    for(unsigned int i=0; i<nn; i+=1) {

        unsigned i0=nl->getClosePair(i).first;  //local index
        unsigned i1=nl->getClosePair(i).second;

        //if(getAbsoluteIndex(i0)==getAbsoluteIndex(i1)) continue;

        if(pbc) {
          distAB[i0][i1]=pbcDistance(getPosition(i0),getPosition(i1));
        } else {
          distAB[i0][i1]=delta(getPosition(i0),getPosition(i1));
        }

        distABinvmod[i0][i1]=1.0/distAB[i0][i1].modulo();
        nnexp[i]=exp(lambda * distAB[i0][i1].modulo());
        nni0[i]=i0;  //local index of the atom in groupA
        nni1[i]=i1;
        //printf("%d %f %e \n",i,nndist[i],nnexp[i] );
        //printf("%d %d %d \n",i,nn,getNumberOfAtoms() );
        //printf("local %d %d couple %d \n", i0,i1,i);
        //printf("absolute index %d %d couple %d \n", getAbsoluteIndex(i0),getAbsoluteIndex(i1),i);
        //to fill this, I need the LOCAL index in the reduced list, not i1!!
        nnexpnormt[i1]+=nnexp[i]; //build normalization, only on second index
    }
    #pragma omp critical  //EACH thread must take their local version and write it into the global one!
    for(unsigned i=0; i<getNumberOfAtoms(); i++) nnexpnorm[i]+=nnexpnormt[i];
    #pragma omp barrier

    //LOOP2 compute the unshifted charge on atoms A
    #pragma omp for
    for(unsigned int i=0; i<nn; i+=1) {
      c[nni0[i]][nni1[i]]=nnexp[i]/nnexpnorm[nni1[i]];  //the couple of indeces is unique, there is no overwriting in the array
      charget[nni0[i]]+=c[nni0[i]][nni1[i]];
    }
    #pragma omp critical
    for(unsigned i=0; i<num_atomsa; i++) charge[i]+=charget[i];
    #pragma omp barrier

    //LOOP2.5 shift charge
    #pragma omp for
    for(unsigned int j=0; j<num_atomsa; j+=1) {
      if(j==num_atomso)     charge[j]-=d1;       //N-2   equals the charge of N in glycine
      else if(j==num_atomso+1) charge[j]-=d2;    //O-0.5 equals the charge of each O in glycine
      else if(j==num_atomso+2) charge[j]-=d3;    //O-0.5 equals the charge of each O in glycine
      else charge[j]-=d0;                        //O-2   equals the charge of each O in water
    }

    //LOOP3 id CV, compute the distance in groupA
    //      also accumulate A(m) = sum_{n!=m} distGA(m,n)*charge(n), restricted to water oxygens
    vector<double> At(num_atomsa, 0.0); // thread-local version of A
    #pragma omp for reduction(-:IonDistance)
    for(unsigned int j=0;j<num_atomsa;j++) {
      for(unsigned int k=j+1;k<num_atomsa;k++) {

        Vector distance_jk;
        if(pbc){
          distance_jk=pbcDistance(getPosition(j),getPosition(k));
        } else {
          distance_jk=delta(getPosition(j),getPosition(k));
        }

                //if(j==k) distGA[j][k] = nn_theta;
        //else distGA[j][k]=distance_jk.modulo();
                distGA[j][k]=distance_jk.modulo();

        double buffer,buffer2;
                //only consider distance between [O in water] and [O in water]
                if (j<num_atomso && k<num_atomso) {
                        buffer = distGA[j][k]*(charge[j] * charge[k]);
                        buffer2=(charge[j] * charge[k])/distGA[j][k];
                        At[j] += charge[k]*distGA[j][k];   // contributes to A(j) (the "greater" side, since k>j)
                        At[k] += charge[j]*distGA[j][k];   // contributes to A(k) (the "less" side, since j<k)
                } else {
                    buffer  =  0.0;
                        buffer2 =  0.0;
                }
        IonDistance -= buffer; //this is the IonDistance for water, the result!
        Vector dd(buffer2*distance_jk);
        //derivatives on groupA only
        omp_deriv[j]+=dd;
        omp_deriv[k]-=dd;
        //NEW: virial contribution for this pair.
        //Convention (consistent with PLUMED's Distance.cpp): for a pair vector
        //r = pbcDistance(pos(j),pos(k)) with dCV/dR_j = +dd, dCV/dR_k = -dd,
        //the box-derivative contribution is +Tensor(r,dd).
        omp_virial += Tensor(distance_jk,dd);
      }
    }
    #pragma omp critical
    for(unsigned i=0; i<num_atomsa; i++) A[i]+=At[i];
    #pragma omp barrier



    //LOOP4 (rewritten): compute the derivatives on couples AB in O(nn) instead of O(nn*num_atomsa^2)
    //
    // Using: dCV/d(O-H pair i) = lambda * c[ind0][ind1] * ( A(ind0) - abar(ind1) )
    // where  A(m)    = sum_{n!=m, water-O only} distGA(m,n)*charge(n)      [already computed above]
    //        abar(h) = sum_o  c[o][h] * A(o)                              [computed below, O(nn)]
    //
    // This is mathematically identical to the original double loop over (k,j) for every
    // neighbor pair i - see derivation notes - just computed in two O(nn) passes instead of
    // one O(nn * num_atomsa^2) pass.

    vector<double> abart(getNumberOfAtoms(), 0.0); // thread-local version of abar

    // Pass 1: build abar(h) = sum_o c[o][h] * A[o]
    #pragma omp for
    for(unsigned int i=0; i<nn; i+=1) {
      unsigned ind0=nni0[i];
      unsigned ind1=nni1[i];
      abart[ind1] += c[ind0][ind1] * A[ind0];
    }
    #pragma omp critical
    for(unsigned i=0; i<getNumberOfAtoms(); i++) abar[i]+=abart[i];
    #pragma omp barrier

    // Pass 2: compute forces using A(ind0) - abar(ind1)
    #pragma omp for
    for(unsigned int i=0; i<nn; i+=1) {
      unsigned ind0=nni0[i];
      unsigned ind1=nni1[i];

      double buf = lambda * c[ind0][ind1] * ( A[ind0] - abar[ind1] );

      Vector dd(buf*distABinvmod[ind0][ind1]*distAB[ind0][ind1]);
      omp_deriv[ind0]+=dd;
      omp_deriv[ind1]-=dd;
      //NEW: virial contribution, same convention as LOOP3
      //(distAB[ind0][ind1] = pbcDistance(pos(ind0),pos(ind1)), dd added to ind0, subtracted from ind1)
      omp_virial += Tensor(distAB[ind0][ind1],dd);
    }

    #pragma omp critical
    {
      for(unsigned i=0; i<getPositions().size(); i++) deriv[i]+=omp_deriv[i];
      virial += omp_virial;   //NEW: merge thread-local virial into the global one
    }
    #pragma omp barrier

    }

  for(unsigned i=0; i<deriv.size(); ++i) setAtomsDerivatives(i,deriv[i]);
  setValue           (IonDistance);
  setBoxDerivatives  (virial);

}
}
}

