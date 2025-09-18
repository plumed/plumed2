/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2021, Andrea Arsiccio

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not
   claim that you wrote the original software. If you use this software
   in a product, an acknowledgment in the product documentation would be
   appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_pines_PINES_h
#define __PLUMED_pines_PINES_h

// Core PLUMED functionality
#include "core/ActionWithVirtualAtom.h"
#include "core/ActionWithValue.h"
#include "core/PlumedMain.h"
#include "tools/Stopwatch.h"
#include "tools/PDB.h"
#include "tools/Vector.h"
#include "tools/AtomNumber.h"
#include "tools/SwitchingFunction.h"
#include "tools/Communicator.h"
#include "tools/Units.h"

// Colvar interface
#include "colvar/Colvar.h"
#include "core/ActionRegister.h"

// STL
#include <string>
#include <vector>

#include <unordered_map>
#include <set>
#include <fstream>
#include <unordered_set>
#include <utility>
#include <functional>

namespace PLMD {
namespace pines {

struct AtomNumberLess {
  bool operator()(const AtomNumber& a, const AtomNumber& b) const {
    return a.index() < b.index();
  }
};

struct pair_hash {
  std::size_t operator()(const std::pair<AtomNumber, AtomNumber>& p) const {
    std::size_t h1 = std::hash<int>()(p.first.index());
    std::size_t h2 = std::hash<int>()(p.second.index());
    return h1 ^ (h2 << 1);
  }
};



class PINES      : public PLMD::colvar::Colvar {
private:
  PLMD::Stopwatch timer;
  int N_Blocks;
  int total_PIV_length;
  int inited_step;
  int last_step_latched;
  bool driver_mode;
  bool freeze_selection;
  std::vector<int> steps_since_update;
  std::vector<int> nstride;
  std::string ref_file;
  std::vector<SwitchingFunction> sfs;
  std::vector<std::string> sw;
  std::vector<double> r00;
  std::vector<std:: vector<double> > PIV;
  std::vector<std:: vector<Vector> > ann_deriv;

  std::vector<std:: vector<AtomNumber> > listall;
  std::vector<std:: vector<AtomNumber> > listreduced;
  std::set<AtomNumber, AtomNumberLess> listreducedall;
  std::vector<AtomNumber> listreducedall_vec;
  std::unordered_map<int,int> atom_ind_hashmap;

  std::vector<bool> stale_tolerance;
  PDB mypdb;
  std::vector<std::string> block_params;
  std::vector<std::vector<std::vector<AtomNumber> > > block_groups_atom_list;
  std::vector<int> block_lengths;
  std::vector<int> Buffer_Pairs;
  std::vector<int> tot_num_pairs;
  std::vector<std::vector<std::vector<bool> > > input_filters;
  std::vector<double> delta_pd;
  std::vector<double> r_tolerance;
  std::vector<std::vector<Vector> > PL_atoms_ref_coords;
  std::vector<std::vector<std::pair<AtomNumber,AtomNumber> > > Exclude_Pairs;
  std::vector<std::vector<std::vector<std::string> > > Name_list;
  std::vector<std::vector<std::vector<AtomNumber> > > ID_list;
  std::vector<std::vector<std::vector<int> > > ResID_list;
  std::vector<char> all_g1g2_pairs;
  std::vector<std::vector<std::pair<double, std::pair<AtomNumber,AtomNumber> > > > vecMaxHeapVecs;
  std::vector<std::vector<std::pair<AtomNumber,AtomNumber> > > latched_pairs;
  std::vector<char> preupdated_block;
  std::vector<bool> isFirstBuild;

  bool atomMatchesFilters(int n, int g, AtomNumber ind, int resid, const std::string& atom_name);
  void buildMaxHeapVecBlock(int n, const PDB& mypdb, std::vector<std::pair<double, std::pair<AtomNumber, AtomNumber>>>& heap);
  void updateBlockPairList(int n, std::vector<std::pair<double, std::pair<AtomNumber, AtomNumber>>>& heap);
  double calculateDistance(int n, const AtomNumber& ind0, const AtomNumber& ind1, const PDB& mypdb);
  std::ofstream log;
  bool ensureBlockUpdated(int n);
  void latchFromCurrentHeaps();
  void logMsg(const std::string& msg, const std::string& section);
  void logMsg(const Vector& vec, const std::string& section);
  void resizeAllContainers(int N);

public:
  static void registerKeywords( Keywords& keys );
  explicit PINES(const ActionOptions&);
  //~PINES();
  // active methods:
  struct MinCompareDist {
    bool operator()(const std::pair<double, std::pair<AtomNumber, AtomNumber>>& p1, const std::pair<double, std::pair<AtomNumber, AtomNumber>>& p2) {
      return p1.first < p2.first; // Min heap
    }
  };
  struct MaxCompareDist {
    bool operator()(const std::pair<double, std::pair<AtomNumber, AtomNumber>>& p1, const std::pair<double, std::pair<AtomNumber, AtomNumber>>& p2) {
      return p1.first > p2.first; // Max heap
    }
  };

  virtual void calculate();
  void checkFieldsAllowed() {}
  // -- SD prepare to requestAtoms during simulation
  void prepare() override;
};

}
}

#endif
