#include "frameGenerator.h"

frameGenerator::frameGenerator(unsigned const nat,std::string_view type):
  atoms(nat),
  d(PLMD::AtomDistribution::getAtomDistribution(type)) {
  generateFrame();
}
std::size_t frameGenerator::size() const {
  return atoms.size();
}
unsigned frameGenerator::getFrame() const {
  return step-1;
}
void frameGenerator::generateFrame() {
  d->frame(atoms,box,step,rng);
  ++step;
}
PLMD::Tensor frameGenerator::getBox()const {
  return   {box[0], box[1], box[2], box[3], box[4], box[5], box[6], box[7], box[8]};
}
const std::vector<PLMD::Vector>& frameGenerator::getAtoms() const {
  return atoms;
}
std::vector<PLMD::Vector> frameGenerator::requestAtoms (const std::vector<PLMD::AtomNumber>& AtomsList) const {
  std::vector<PLMD::Vector> atoms_indexed(AtomsList.size());
  unsigned i=0;
  for(const auto idx: AtomsList) {
    atoms_indexed[i] = atoms[idx.index()];
    ++i;
  }
  return atoms_indexed;
}

