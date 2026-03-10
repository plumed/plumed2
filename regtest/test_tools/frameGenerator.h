#ifndef __PLUMED_TEST_FRAMEGENERATOR
#define __PLUMED_TEST_FRAMEGENERATOR
#include "plumed/tools/AtomNumber.h"
#include "plumed/tools/Communicator.h"
#include "plumed/tools/NeighborList.h"
#include "plumed/tools/Pbc.h"
#include "plumed/tools/AtomDistribution.h"
#include "plumed/tools/Random.h"
#include <vector>

class frameGenerator {
  PLMD::Random rng;
  std::vector<double> box=std::vector<double>(9,0.0);
  std::vector<PLMD::Vector> atoms;
  std::unique_ptr<PLMD::AtomDistribution> d;
  unsigned step=0;
public:
  frameGenerator(unsigned const nat,std::string_view type);
  std::size_t size() const;
  unsigned getFrame() const;
  void generateFrame();
  PLMD::Tensor getBox()const;
  const std::vector<PLMD::Vector>& getAtoms() const;
  std::vector<PLMD::Vector> requestAtoms (const std::vector<PLMD::AtomNumber>& AtomsList) const;
};
#endif //__PLUMED_TEST_FRAMEGENERATOR
