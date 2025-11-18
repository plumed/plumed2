#include "plumed/tools/Communicator.h"
#include "plumed/tools/LinkCells.h"
#include "plumed/tools/Pbc.h"
#include "plumed/tools/Stopwatch.h"
#include "plumed/tools/AtomDistribution.h"

#include "testUtils.h"
#include "mpi.h"

#include <array>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <numeric>

#include <vector>
#include <numeric>

using namespace PLMD;

void test (PLMD::Communicator &comm);
void testIndexes (PLMD::Communicator &comm);
void testRequiredCells (PLMD::Communicator &comm);
void bench (PLMD::Communicator &comm);

int main(int argc,char**argv) {
  MPI_Init(&argc,&argv);
  {
    MPI_Comm c;
    MPI_Comm_dup(MPI_COMM_WORLD,&c);
    PLMD::Communicator comm;
    comm.Set_comm(&c);
    test(comm);
    testIndexes(comm);
    testRequiredCells(comm);
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


void testRequiredCells (PLMD::Communicator &comm) {
  Pbc pbc;
  std::vector<Vector> atoms(1, {0,0,0});
  tee ofs("outputRequiredCells");
  LinkCells cells(comm);
  pbc.setBox(
  Tensor{
    10.0,  0.0,  0.0,
    0.0, 10.0,  0.0,
    0.0,  0.0, 10.0
  });
  std::vector<unsigned> cells_required(27);
  unsigned ncells_required=0;
  for( const auto cutoff: {
         10.0,5.0,3.0,2.0
       }) {
    ofs << "####cutoff (in cell-size): " << cutoff/pbc.getBox()[0][0] << "\n";
    cells.setCutoff(cutoff);
    cells.setupCells(atoms,pbc);
    ofs << "NUMBER of cells: " <<  cells.getNumberOfCells() << "\n";
    PLMDTests::testcases<std::array<unsigned,3>> celPos= {{"vertex",std::array<unsigned,3>{0,0,0}}};
    if (cells.getNumberOfCells() > 8) {
      //the opposite end
      auto last_cell = cells.getCellLimits();
      last_cell[0]--;
      last_cell[1]--;
      last_cell[2]--;
      celPos.push_back({"opposite vertex",last_cell});
//in the middle of the cell
      last_cell = cells.getCellLimits();
      last_cell[0]/=2;
      last_cell[1]/=2;
      last_cell[2]/=2;
      celPos.push_back({"middle",last_cell});
//in the middle of an edge
      last_cell = cells.getCellLimits();
      last_cell[0]=0;
      last_cell[1]/=2;
      last_cell[2]=0;
      celPos.push_back({"edge",last_cell});
//in the middle of an edge
      last_cell = cells.getCellLimits();
      last_cell[0]/=2;
      last_cell[1]=0;
      last_cell[2]--;
      celPos.push_back({"edge opposite",last_cell});
// in the middle of a face
      last_cell = cells.getCellLimits();
      last_cell[0]/=2;
      last_cell[1]=0;
      last_cell[2]/=2;
      celPos.push_back({"face",last_cell});
// in the middle of a face
      last_cell = cells.getCellLimits();
      last_cell[0]--;
      last_cell[1]/=2;
      last_cell[2]/=2;
      celPos.push_back({"face opposite",last_cell});
    }
    for (const auto& [name,cp]: celPos) {
      ncells_required=0;

      auto ss= std::ostringstream() ;
      ss << "[co=" << cutoff << " -"<<name<<"- {"<< cp[0] << ", "<<cp[1] <<", " <<cp[2]<< "}" << "] ";
      //the "header" is useful in undersanding where is the test case that broke the test
      auto head = ss.str();
      cells.addRequiredCells(cp,ncells_required, cells_required);
      ofs <<head<< "Cells Required:" << ncells_required << "\n";
      //ordering for "test stability"
      std::sort(cells_required.begin(),cells_required.begin()+ncells_required);
      ofs <<head<< "Cells:";
      for (unsigned i=0; i< ncells_required; ++i) {
        ofs << " " << cells_required[i];
      }
      ofs << "\n";
    }
    //Now we redo the calculations, but this time with no PBCs
    for (const auto& [name,cp]: celPos) {
      ncells_required=0;

      auto ss= std::ostringstream() ;
      ss << "[co=" << cutoff << " -"<<name<<"- {"<< cp[0] << ", "<<cp[1] <<", " <<cp[2]<< "}" << " noPBC] ";
      //the "header" is useful in undersanding where is the test case that broke the test
      auto head = ss.str();
      cells.addRequiredCells(cp,ncells_required, cells_required, false);
      ofs <<head<< "Cells Required:" << ncells_required << "\n";
      //ordering for "test stability"
      std::sort(cells_required.begin(),cells_required.begin()+ncells_required);
      ofs <<head<< "Cells:";
      for (unsigned i=0; i< ncells_required; ++i) {
        ofs << " " << cells_required[i];
      }
      ofs << "\n";
    }
  }
}

void test (PLMD::Communicator &comm) {
  tee ofs ("output");
  Random rng;
  const std::array<unsigned,3> replicas= {1,1,1};
  //std::vector<Vector> atoms(4*4*4);
  std::vector<Vector> atoms(5*5*5);
  std::vector<double> basebox(9);
  unsigned nat = atoms.size();
  auto rep= std::make_unique<repliedTrajectory>([&]() {
    auto d = AtomDistribution::getAtomDistribution("sc");
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
  //const unsigned startingnat=nat;
  rep->overrideNat(nat);
  //const unsigned nrep = nat/startingnat;
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
