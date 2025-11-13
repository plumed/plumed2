#include "plumed/tools/Communicator.h"
#include "plumed/tools/LinkCells.h"
#include "plumed/tools/Pbc.h"
#include "plumed/tools/Stopwatch.h"
#include "plumed/tools/AtomDistribution.h"

#include "testUtils.h"
#include "mpi.h"

#include <fstream>
#include <iostream>
#include <vector>
#include <numeric>

using namespace PLMD;

void test (PLMD::Communicator &comm);
void testIndexes (PLMD::Communicator &comm);
void bench (PLMD::Communicator &comm);

int main(int argc,char**argv) {
  MPI_Init(&argc,&argv);
  {
    MPI_Comm c;
    MPI_Comm_dup(MPI_COMM_WORLD,&c);
    PLMD::Communicator comm;
    comm.Set_comm(&c);
    //test(comm);
    testIndexes(comm);
    //bench(comm);
  }
  MPI_Finalize();
}

void testIndexes (PLMD::Communicator &comm) {
  tee ofs ("outputIndexes");

  for(Tensor mybox: {
  Tensor{
  10.0,  0.0,  0.0,
  0.0, 10.0,  0.0,
  0.0,  0.0, 10.0
}, {
  // non ortho
  10.0,  10.0,  0.0,
  0.0, 10.0,  0.0,
  0.0,  0.0, 10.0
}, {
  // heavily non ortho
  10.0,  5.0,  3.0,
  5.0, 10.0,  2.0,
  3.0,  2.0, 10.0
},
    }) {
    Pbc pbc;
    pbc.setBox(mybox);
    LinkCells cells(comm);
    //setting the cutoff as the original cell for sake of simplicity
    cells.setCutoff(1.5);
//PBC are definded, but we still need a dummy atom for setting up the cell
    std::vector<Vector> atoms(1);
    cells.setupCells(make_const_view(atoms),pbc);
    ofs <<"Ncells: "<< cells.getNumberOfCells()<<"\n";
    auto maxs = cells.getCellLimits();
    for(unsigned k=0; k<maxs[2]; ++k) {
      for(unsigned j=0; j<maxs[1]; ++j) {
        for(unsigned i=0; i<maxs[0]; ++i) {
          std::array<unsigned,3> cell{i,j,k};
          ofs << cells.convertIndicesToIndex(cell)<< " ("
              << i <<" "<< j << " " << k
              <<")";
          auto cellID = cells.findMyCell(cells.convertIndicesToIndex(cell));
          ofs << " == ("
              << cellID[0] <<" "<< cellID[1] << " " << cellID[2]
              <<") "
              << ((cellID[0]==i && cellID[1] ==j && cellID[2]==k )?"true":"false")
              <<"\n";
        }
      }
    }
  }
}


void test (PLMD::Communicator &comm) {
  tee ofs ("output");
  Random rng;
  const std::array<unsigned,3> replicas= {2,2,1};
  std::vector<Vector> atoms(4*4*4);
  std::vector<double> basebox(9);
  unsigned nat = atoms.size();
  auto rep= std::make_unique<repliedTrajectory>([&]() {
    auto d = AtomDistribution::getAtomDistribution("sphere");
//getting some base informations
    d->frame(atoms,basebox,0,rng);
    return d;
  }
  (),
  replicas[0],
  replicas[1],
  replicas[2],
  nat);
  const unsigned startingnat=nat;
  rep->overrideNat(nat);
  const unsigned nrep = nat/startingnat;
  atoms.resize(nat);
//for sphere the box dimension is slighly bigger than the sphere:
  double radius = basebox[0]/2;
  std::vector<double> box(9);
  ofs <<"Half cell:"<< radius<< "\n";
  rep->frame(atoms,box,0,rng);
  Tensor mybox{box[0], box[1], box[2], box[3], box[4], box[5], box[6], box[7], box[8]};
  Pbc pbc;
  pbc.setBox(mybox);
  LinkCells cells(comm);
  //setting the cutoff as the original cell for sake of simplicity
  cells.setCutoff(radius);
  cells.setupCells(make_const_view(atoms),pbc);
  ofs <<"Ncells: "<< cells.getNumberOfCells()<<"\n";
  ofs <<"atomGroup: cells\n";
  //for how the sphere is made, each cube should have its own group
  std::ofstream f("sphere.xyz");
  f << nat << "\n";
  f << "Lattice=\"" << pbc.getBox() << "\" Properties=species:S:1:pos:R:3:cell:I:1\n";
  std::vector<unsigned int> indices(atoms.size());
  std::iota(indices.begin(),indices.end(),0);
  auto collection = cells.getCollection(make_const_view(atoms),
                                        make_const_view(indices));
  std::vector <unsigned> thecell(atoms.size());
  for(unsigned i=0; i< cells.getNumberOfCells(); ++i) {
    auto thiscell=collection.getCellIndexes(i);
    for (auto x:thiscell) {
      thecell[x]=i;
    }
  }
  bool allequal=true;
  for (unsigned i=0; i< nrep; ++i) {
    ofs << i <<":";
    for (unsigned j=0; j< startingnat; ++j) {
      auto mycell=cells.findCell(atoms[startingnat*i+j]);
      ofs << " "<<mycell <<","<<thecell[startingnat*i+j];
      allequal &= mycell==thecell[startingnat*i+j];
      f << "X " << atoms[startingnat*i+j] << " " << mycell << "\n";
    }
    ofs << "\n";
  }
  ofs << "the cells coincide: " << (allequal?"true":"false") << "\n";
}

void bench (PLMD::Communicator &comm) {
  tee ofs ("output");
  Random rng;
  const std::array<unsigned,3> replicas= {2,2,1};
  std::vector<Vector> atoms(50000);
  std::vector<double> basebox(9);
  unsigned nat = atoms.size();
  auto rep= std::make_unique<repliedTrajectory>([&]() {
    auto d = AtomDistribution::getAtomDistribution("cube");
//getting some base informations
    d->frame(atoms,basebox,0,rng);
    return d;
  }
  (),
  replicas[0],
  replicas[1],
  replicas[2],
  nat);
  const unsigned startingnat=nat;
  rep->overrideNat(nat);
  const unsigned nrep = nat/startingnat;
  atoms.resize(nat);
//for sphere the box dimension is slighly bigger than the sphere:
  double radius = basebox[0]/2;
  std::vector<double> box(9);
  ofs <<"Half cell:"<< radius<< "\n";
  LinkCells cells(comm);
  Pbc pbc;
  for (int i=0; i<5000; ++i) {
    rep->frame(atoms,box,0,rng);
    Tensor mybox{box[0], box[1], box[2], box[3], box[4], box[5], box[6], box[7], box[8]};
    pbc.setBox(mybox);
    //setting the cutoff as the original cell for sake of simplicity
    cells.setCutoff(radius/5);

    std::vector<unsigned int> indices(atoms.size());
    std::iota(indices.begin(),indices.end(),0);
    cells.buildCellLists(make_const_view(atoms),
                         make_const_view(indices),
                         //cells.buildCellLists(atoms,
                         //                      indices,
                         pbc);
  }
  ofs <<"Ncells: "<< cells.getNumberOfCells()<<"\n";
  /*
    cells.setupCells(make_const_view(atoms),pbc);
    ofs <<"atomGroup: cells\n";
    //for how the sphere is made, each cube should have its own group
    std::ofstream f("sphere.xyz");
    f << nat << "\n";
    f << "Lattice=\"" << pbc.getBox() << "\" Properties=species:S:1:pos:R:3:cell:I:1\n";
    for (unsigned i=0; i< nrep; ++i) {
      ofs << i <<":";
      for (unsigned j=0; j< startingnat; ++j) {
        auto mycell=cells.findCell(atoms[startingnat*i+j]);
        ofs << " "<<mycell;
        f << "X " << atoms[startingnat*i+j] << " " << mycell << "\n";
      }
      ofs << "\n";
    }*/
}
