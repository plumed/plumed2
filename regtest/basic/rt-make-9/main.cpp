#include "plumed/tools/Tools.h"
#include <fstream>
#include <iostream>
#include <vector>

using namespace PLMD;

void run(bool queue) {
  std::ofstream ofs("output",std::ios_base::app);
  std::vector<int> v1={1,3,5,10,11};
  std::vector<int> v2={2,9,10,11,12,13};
  std::vector<int> v3={};
  std::vector<int> v4={5};

  {
    std::vector<int> result;
    std::vector<std::vector<int>*> vecs;

    vecs.push_back(&v1);
    vecs.push_back(&v2);
    vecs.push_back(&v3);
    vecs.push_back(&v4);
    Tools::mergeSortedVectors(vecs,result);
    for(const auto &v : result) ofs<<" "<<v; ofs<<"\n";
  }


  {
    // note: result will be appended! So these will be at the beginning of the vector
    std::vector<int> result={999,1000};
    std::vector<std::vector<int>*> vecs;

    // note: pointers can be repeated
    vecs.push_back(&v4);
    vecs.push_back(&v3);
    vecs.push_back(&v2);
    vecs.push_back(&v1);
    vecs.push_back(&v1);
    vecs.push_back(&v2);
    vecs.push_back(&v3);
    vecs.push_back(&v4);
    Tools::mergeSortedVectors(vecs,result);
    for(const auto &v : result) ofs<<" "<<v; ofs<<"\n";
  }
}

int main() {
  run(false);
  run(true);

  return 0;
}
