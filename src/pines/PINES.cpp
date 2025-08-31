/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2023 of Nicholas Herringer and Siva Dasetty.

The PINES module is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The PINES module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "colvar/Colvar.h"
#include "colvar/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionWithVirtualAtom.h"
#include "tools/NeighborList.h"
#include "tools/SwitchingFunction.h"
#include "tools/PDB.h"
#include "tools/Pbc.h"
#include "tools/Stopwatch.h"

// -- SD header file for PINES
#include "PINES.h"

#include <string>
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <unordered_map>
#include <set>
#include <numeric>
#include <omp.h>
#include <unordered_set>
#include <utility>
#include <functional>
#include <vector>

constexpr int CACHE_LINE_SIZE = 64;

template<typename T>
struct alignas(CACHE_LINE_SIZE) Padded {
    T value;

    // Optional constructors
    Padded() : value() {}
    explicit Padded(const T& val) : value(val) {}
};

using namespace std;

namespace PLMD
{
  namespace PINES
  {

    PLUMED_REGISTER_ACTION(PINES, "PINES")

    void PINES::registerKeywords(Keywords &keys)
    {
      Colvar::registerKeywords(keys);
      keys.add("numbered", "SWITCH", "The switching functions parameter. You must specify a Switching function for all PINES blocks."
                                     "Details of the various switching functions you can use are provided on \\ref switchingfunction.");
      keys.add("numbered", "BLOCK", "Each block of the PIV");
      keys.add("compulsory", "REF_FILE", "PDB file name that contains the information about system connectivity and labels.");
      keys.add("compulsory", "N_BLOCKS", "Number of blocks in PIV");
      keys.add("compulsory", "BUFFER", "Number of additional pairwise distances to include as a buffer for each PIV block");
      keys.add("compulsory", "PL_REFRESH", "Upper limit refresh rate for each block pair list");
      keys.add("compulsory", "SIZE", "Length of each PIV Block");
      keys.add("optional", "ATOMID", "AtomIDs");
      keys.add("optional", "RESID", "ResIDs");
      keys.add("optional", "NAME", "Atom Names");
      keys.add("optional", "EXCLUDE_PAIRS", "Excluded pairs");
      keys.addFlag("DRIVERMODE", false, "Use when post-processing a trajectory. This will force the pair list to be rebuilt at every step (i.e. no pair list)");
      keys.addFlag("FREEZE_SELECTION", false, "Freeze top-K membership during derivative checks");
      componentsAreNotOptional(keys);
      keys.addOutputComponent("ELEMENT", "default", "Elements of the PINES block"); 
      keys.reset_style("SWITCH", "compulsory");
    }

    bool PINES::atomMatchesFilters(int n, int g, AtomNumber ind, int resid, const std::string& atom_name) {
      bool id_check = !input_filters[n][g][0];
      bool res_check = !input_filters[n][g][1];
      bool name_check = !input_filters[n][g][2];
    
      if (input_filters[n][g][0]) {
        for (auto& atomID : ID_list[n][g]) {
          if (ind == atomID) { id_check = true; break; }
        }
      }
      if (input_filters[n][g][1]) {
        for (auto& rID : ResID_list[n][g]) {
          if (resid == rID) { res_check = true; break; }
        }
      }
      if (input_filters[n][g][2]) {
        for (auto& name : Name_list[n][g]) {
          if (atom_name == name) { name_check = true; break; }
        }
      }
      return id_check && res_check && name_check;
    }
    
    void PINES::buildMaxHeapVecBlock(int n, const PDB& mypdb, std::vector<std::pair<double, std::pair<AtomNumber, AtomNumber>>>& heap) {
      //logMsg("Building Pair List for block: " + std::to_string(n), "buildMaxHeapVecBlock");

      int nthreads = omp_get_max_threads();
      std::vector<std::vector<std::pair<double, std::pair<AtomNumber, AtomNumber>>>> thread_heaps(nthreads);

      //logMsg("isFirstBuild[n]? " + std::to_string(isFirstBuild[n]), "buildMaxHeapVecBlock");
      
      // Rough estimate, tweak as needed
      int estimated_pairs_per_thread = block_groups_atom_list[n][0].size() * block_groups_atom_list[n][1].size() / nthreads;

      for (int tid = 0; tid < nthreads; ++tid) {
	thread_heaps[tid].reserve(estimated_pairs_per_thread);
      }

      heap.clear();
      #pragma omp parallel for schedule(dynamic)
      for (int i = 0; i < block_groups_atom_list[n][0].size(); i++) {
        AtomNumber ind0 = block_groups_atom_list[n][0][i];
	Vector Pos0 = isFirstBuild[n] ? mypdb.getPosition(ind0) : getPosition(atom_ind_hashmap[ind0.index()]);
        int tid = omp_get_thread_num();
        auto& local_heap = thread_heaps[tid];
        std::unordered_set<std::pair<AtomNumber, AtomNumber>, pair_hash> local_unique_pairs;
        // logMsg("ind0: " + std::to_string(ind0.index()), "buildMaxHeapVecBlock");
        // logMsg("Num pairs in heap: " + std::to_string(block_groups_atom_list[n][0].size()), "buildMaxHeapVecBlock");
        // logMsg(Pos0, "buildMaxHeapVecBlock");
        for (int j = 0; j < block_groups_atom_list[n][1].size(); j++) {
          AtomNumber ind1 = block_groups_atom_list[n][1][j];
          // logMsg("ind1: " + std::to_string(ind1.index()), "buildMaxHeapVecBlock");
          // logMsg("Num pairs in heap: " + std::to_string(block_groups_atom_list[n][1].size()), "buildMaxHeapVecBlock");
          if (ind1 == ind0) continue;

          auto test_pair = std::make_pair(ind0, ind1);
          auto reverse_pair = std::make_pair(ind1, ind0);
    
          if (std::find(Exclude_Pairs[n].begin(), Exclude_Pairs[n].end(), test_pair) != Exclude_Pairs[n].end() ||
              std::find(Exclude_Pairs[n].begin(), Exclude_Pairs[n].end(), reverse_pair) != Exclude_Pairs[n].end()) {
            continue;
          }
          if (local_unique_pairs.count(reverse_pair)) continue;
          local_unique_pairs.insert(test_pair);

          Vector Pos1 = isFirstBuild[n] ? mypdb.getPosition(ind1) : getPosition(atom_ind_hashmap[ind1.index()]);
          
          //logMsg(Pos1, "buildMaxHeapVecBlock");
          double mag = pbcDistance(Pos0, Pos1).modulo();
          //logMsg("pairdist: " + std::to_string(mag), "buildMaxHeapVecBlock");
          local_heap.emplace_back(mag, std::make_pair(ind0, ind1));
        }
      }

        // Flatten all thread-local heaps into one vector
      std::vector<std::pair<double, std::pair<AtomNumber, AtomNumber>>> all_pairs;
      for (auto& th : thread_heaps) {
        all_pairs.insert(all_pairs.end(), th.begin(), th.end());
      }

      // Build max-heap from top-K shortest pairs
      std::sort(all_pairs.begin(), all_pairs.end(), MinCompareDist());
      size_t K = std::min(static_cast<size_t>(tot_num_pairs[n]), all_pairs.size());
      heap.assign(all_pairs.begin(), all_pairs.begin() + K);
 
      //logMsg("Heap created succesfully", "buildMaxHeapVecBlock");
      std::sort(heap.begin(), heap.end(), MaxCompareDist());
      int chk1 = 0;
      int chk2 = Buffer_Pairs[n];
      //logMsg("heap[chk1].first: " + std::to_string(heap[chk1].first), "buildMaxHeapVecBlock");
      //logMsg("heap[chk2].first: " + std::to_string(heap[chk2].first), "buildMaxHeapVecBlock");
      delta_pd[n] = heap[chk1].first - heap[chk2].first;
      r_tolerance[n] = delta_pd[n]/4;
      listreduced[n].clear();
      std::set<AtomNumber, AtomNumberLess> uniqueAtoms;
      for (const auto& pair : heap) {
        uniqueAtoms.insert(pair.second.first);  // Atom1
        uniqueAtoms.insert(pair.second.second); // Atom2
      }
      for (const auto& atom : uniqueAtoms) {
        listreduced[n].push_back(atom);
      }
      //logMsg("listreduced.size(): " + std::to_string(listreduced[n].size()), "buildMaxHeapVecBlock");
      if(driver_mode || isFirstBuild[n]){
        listreduced[n] = listall[n];
      }
      //logMsg("listreduced.size() after drivermode change: " + std::to_string(listreduced[n].size()), "buildMaxHeapVecBlock");
      //logMsg("listreduced[" + std::to_string(n) + "].size(): " + std::to_string(listreduced[n].size()), "buildMaxHeapVecBlock");
      PL_atoms_ref_coords[n].resize(listreduced[n].size());
      //logMsg("PL_atoms_ref_coords[" + std::to_string(n) + "].size(): " + std::to_string(PL_atoms_ref_coords[n].size()), "buildMaxHeapVecBlock");
      for (int i=0; i<listreduced[n].size(); i++)
      {
        PL_atoms_ref_coords[n][i].zero();
        //logMsg("listreduced[" + std::to_string(n) + "][" + std::to_string(i) + "].index(): " + std::to_string(listreduced[n][i].index()), "buildMaxHeapVecBlock");
        PL_atoms_ref_coords[n][i] = isFirstBuild[n] ? mypdb.getPosition(listreduced[n][i]) : getPosition(atom_ind_hashmap[listreduced[n][i].index()]);
      }
      isFirstBuild[n] = false;

    }
    
    void PINES::updateBlockPairList(int n, std::vector<std::pair<double, std::pair<AtomNumber, AtomNumber>>>& heap) {
      // This is separate from staleness, which triggers a refresh instead of an update.
      // This is simply to update the pair distances in the MaxHeap/Pairlist with the new atom positions
      // and rearrange the ordering in the MaxHeap/Pairlisat

      plumed_massert((int)heap.size() >= tot_num_pairs[n], "heap too small");
      plumed_massert(!listreducedall_vec.empty(), "positions requested before requestAtoms");

      #pragma omp parallel for
      for (int i = 0; i < tot_num_pairs[n]; i++)
      {
        AtomNumber ind0 = heap[i].second.first;
        AtomNumber ind1 = heap[i].second.second;
	auto it0 = atom_ind_hashmap.find(ind0.index());
        auto it1 = atom_ind_hashmap.find(ind1.index());
        plumed_massert(it0 != atom_ind_hashmap.end() && it1 != atom_ind_hashmap.end(),"updateBlockPairList: atom index not requested");
        heap[i].first = calculateDistance(n,ind0,ind1,mypdb);
      }
      std::nth_element(heap.begin(), heap.begin() + block_lengths[n], heap.end(), MinCompareDist());
      std::sort(heap.begin(), heap.begin() + block_lengths[n], MaxCompareDist());
      // for (int i = 0; i < tot_num_pairs[n]; i++)
      // {
      //   logMsg("Pair: " + std::to_string(heap[i].second.first.index()) + ", " + std::to_string(heap[i].second.second.index()) + "; Dist: " + std::to_string(heap[i].first),"updateBlockPairList");
      // }
    }
    
    double PINES::calculateDistance(int n, const AtomNumber& ind0, const AtomNumber& ind1, const PDB& mypdb) {
      Vector Pos0 = isFirstBuild[n] ? mypdb.getPosition(ind0) : getPosition(atom_ind_hashmap[ind0.index()]);
      Vector Pos1 = isFirstBuild[n] ? mypdb.getPosition(ind1) : getPosition(atom_ind_hashmap[ind1.index()]);
      return pbcDistance(Pos0, Pos1).modulo();
    }

    void PINES::resizeAllContainers(int N) {
      // Outer containers
      nstride.resize(N,1);
      steps_since_update.resize(N, 0);
      block_params.resize(N);
      block_groups_atom_list.resize(N);
      block_lengths.resize(N);
      Buffer_Pairs.resize(N);
      tot_num_pairs.resize(N);
      Exclude_Pairs.resize(N);
      all_g1g2_pairs.resize(N);
      vecMaxHeapVecs.resize(N);
      PIV.resize(N);
      listall.resize(N);
      listreduced.resize(N);
      stale_tolerance.resize(N);
      r00.resize(N);
      sw.resize(N);
      sfs.resize(N);
      delta_pd.resize(N, 0.0);
      r_tolerance.resize(N, 0.0);
      PL_atoms_ref_coords.resize(N);
      input_filters.resize(N);
      ID_list.resize(N);
      ResID_list.resize(N);
      Name_list.resize(N);
      atom_ind_hashmap.clear();
      latched_pairs.resize(N);
      isFirstBuild.resize(N,true);
    
      // Per-block inner structures
      for (int n = 0; n < N; n++) {
        block_groups_atom_list[n].resize(2);
        ID_list[n].resize(2);
        ResID_list[n].resize(2);
        Name_list[n].resize(2);
        input_filters[n].resize(2);
        PL_atoms_ref_coords[n].resize(0);  // will be filled per atom if needed
      }
    
      // Final 3D input_filters init
      for (int n = 0; n < N; n++) {
        for (int g = 0; g < 2; g++) {
          input_filters[n][g].resize(3, false);  // [ID, ResID, Name]
        }
      }
    }    

    // member state
    // static long inited_step = -1;
    // static std::vector<char> preupdated_block;  // size N_Blocks, resets each step

    bool PINES::ensureBlockUpdated(int n) {
      const bool do_update = true; //(plumed.getStep() == 0 || driver_mode);
      if (!do_update) return false;          // not an update step
      if ((int)vecMaxHeapVecs[n].size() != tot_num_pairs[n]) {
        buildMaxHeapVecBlock(n, mypdb, vecMaxHeapVecs[n]);  // safe: PDB-based build
      }
      if (preupdated_block[n]) return false; // already updated this step
      updateBlockPairList(n, vecMaxHeapVecs[n]);            // uses getPosition(); safe in calculate()
      preupdated_block[n] = 1;
      return true;
    }

    void PINES::latchFromCurrentHeaps() {
      if ((int)latched_pairs.size() != N_Blocks) latched_pairs.resize(N_Blocks);
      for (int n = 0; n < N_Blocks; ++n) {
        latched_pairs[n].resize(block_lengths[n]);
        const bool do_update = true; // (plumed.getStep() == 0 || steps_since_update[n] != 1 || driver_mode);
        if (do_update) {
          // top-K prefix (already sorted by ensureBlockUpdated if it ran)
          for (int i = 0; i < block_lengths[n]; ++i) {
            latched_pairs[n][i] = vecMaxHeapVecs[n][i].second;
	  }
        } else {
          // buffered slice
          for (int i = 0; i < block_lengths[n]; ++i) {
	    latched_pairs[n][i] = vecMaxHeapVecs[n][i + Buffer_Pairs[n]].second;
	  }
        }
      }
    }

    void PINES::logMsg(const std::string& msg, const std::string& section) {
      log << "[" << plumed.getStep() << "] "
          << "[" << section << "] " << msg << std::endl;
    }

    void PINES::logMsg(const Vector& vec, const std::string& section) {
      log << "[" << plumed.getStep() << "] "
          << "[" << section << "] Vector = (" 
          << vec[0] << ", " << vec[1] << ", " << vec[2] << ")" << std::endl;
    }

    PINES::PINES(const ActionOptions &ao) : PLUMED_COLVAR_INIT(ao),
                                            N_Blocks(1),
                                            total_PIV_length(1),
                                            driver_mode(false),
					    freeze_selection(false),
					    last_step_latched(-1),
					    inited_step(-1),
                                            steps_since_update(std::vector<int>(1)),
                                            nstride(std::vector<int>(1,10)),
                                            ref_file(std::string()),
                                            sfs(),
                                            sw(std::vector<string>(1)),
                                            r00(std::vector<double>(1)),
                                            PIV(std::vector<std::vector<double>>(1)),
                                            ann_deriv(std::vector<std::vector<Vector> >(1)),
                                            listall(std::vector<std::vector<AtomNumber> >()),
                                            listreduced(std::vector<std::vector<AtomNumber> >()),
                                            listreducedall(std::set<AtomNumber, AtomNumberLess>()),
                                            listreducedall_vec(std::vector<AtomNumber>()),
                                            atom_ind_hashmap(),
                                            stale_tolerance(std::vector<bool>(1,false)),
                                            mypdb(),
                                            block_params(std::vector<string>(1)),
                                            block_groups_atom_list(std::vector<std::vector<std::vector<AtomNumber> > >()),
                                            block_lengths(std::vector<int>(1,1)),
                                            Buffer_Pairs(std::vector<int>(1,0)),
                                            tot_num_pairs(std::vector<int>(1,1)),
                                            input_filters(
                                              std::vector<std::vector<std::vector<bool> > >(
                                                1,  // outermost dimension
                                                std::vector<std::vector<bool>>(
                                                  1,  // middle dimension
                                                  std::vector<bool>(1, false)  // innermost dimension with value `false`
                                                )
                                              )
                                            ),
                                            delta_pd(std::vector<double>()),
                                            r_tolerance(std::vector<double>()),
                                            PL_atoms_ref_coords(),
                                            Exclude_Pairs(std::vector<std::vector<std::pair<AtomNumber, AtomNumber> > >()),
                                            Name_list(std::vector<std::vector<std::vector<string> > >()),
                                            ID_list(std::vector<std::vector<std::vector<AtomNumber> > >()),
                                            ResID_list(std::vector<std::vector<std::vector<int> > >()),
                                            all_g1g2_pairs(std::vector<char>()),
                                            vecMaxHeapVecs(),
					    latched_pairs(std::vector<std::vector<std::pair<AtomNumber,AtomNumber> > >()),
					    preupdated_block(std::vector<char>()),
					    isFirstBuild(std::vector<bool>())
    {
      log.open("pines_debug.log", std::ios::out);
      if (!log.is_open()) {
        error("Failed to open debug log file.");
      }

      // Reference PDB file from which atom names, types, ids, and initial positions are determined
      parse("REF_FILE", ref_file);
      FILE *fp = fopen(ref_file.c_str(), "r");
      if (fp != NULL)
      {
        mypdb.readFromFilepointer(fp, plumed.getAtoms().usingNaturalUnits(), 0.1 / atoms.getUnits().getLength());
        fclose(fp);
      }
      else error("Error in reference PDB file");

      // Create variable to get number of blocks
      parse("N_BLOCKS", N_Blocks);
      parseFlag("DRIVERMODE", driver_mode);
      parseFlag("FREEZE_SELECTION", freeze_selection);
      //logMsg("N_Blocks: " + std::to_string(N_Blocks), "Constructor");
      resizeAllContainers(N_Blocks);


      // Check that the correct number of Blocks are specified
      for (unsigned n = 0; n < N_Blocks; n++)
      {
        if (!parseNumbered("BLOCK", n + 1, block_params[n]))
          break;
      }

      // Parse blocks for keywords
      for (int n = 0; n < N_Blocks; n++)
      {
        string block_length;
        std::vector<string> ex_pairs_n;
        std::vector<string> g1_ids;
        std::vector<string> g2_ids;
        std::vector<string> g1_resids;
        std::vector<string> g2_resids;
        std::vector<string> g1_names;
        std::vector<string> g2_names;

        std::vector<string> block_data = Tools::getWords(block_params[n]);


        std::vector<string> G1_data;
        Tools::parseVector(block_data, "G1", G1_data);

        std::vector<string> G2_data;
        Tools::parseVector(block_data, "G2", G2_data);

        Tools::parseVector(G1_data, "ATOMID", g1_ids);

        for (int i = 0; i < g1_ids.size(); i++)
        {
          AtomNumber g1i;
          g1i.setIndex(std::stoi(g1_ids[i]));
          ID_list[n][0].push_back(g1i);
        }

        Tools::parseVector(G2_data, "ATOMID", g2_ids);

        for (int i = 0; i < g2_ids.size(); i++)
        {
          AtomNumber g2i;
          g2i.setIndex(std::stoi(g2_ids[i]));
          ID_list[n][1].push_back(g2i);
        }

        Tools::parseVector(G1_data, "RESID", g1_resids);
        for (int i = 0; i < g1_resids.size(); i++)
        {
          ResID_list[n][0].push_back(std::stoi(g1_resids[i]));
        }
        Tools::parseVector(G2_data, "RESID", g2_resids);
        for (int i = 0; i < g2_resids.size(); i++)
        {
          ResID_list[n][1].push_back(std::stoi(g2_resids[i]));
        }
        Tools::parseVector(G1_data, "NAME", g1_names);
        for (int i = 0; i < g1_names.size(); i++)
        {
          Name_list[n][0].push_back(g1_names[i]);
        }
        Tools::parseVector(G2_data, "NAME", g2_names);
        for (int i = 0; i < g2_names.size(); i++)
        {
          Name_list[n][1].push_back(g2_names[i]);
        }

        Tools::parse(block_data, "SIZE", block_length);
        block_lengths[n] = std::stoi(block_length);

        Tools::parseVector(block_data, "EXCLUDE_PAIRS", ex_pairs_n);
        if (!ex_pairs_n.empty())
        {
          for (int i = 0; i < ex_pairs_n.size() - 1; i+=2)
          {
            AtomNumber atom1; 
            atom1.setIndex(std::stoi(ex_pairs_n[i]));
            AtomNumber atom2; 
            atom2.setIndex(std::stoi(ex_pairs_n[i+1]));

            std::pair<AtomNumber, AtomNumber> excluded_pair;
            excluded_pair = {atom1, atom2};
            Exclude_Pairs[n].push_back(excluded_pair);
          }
        }
        string buffer_pairs;
        string pl_refresh;
        Tools::parse(block_data, "BUFFER", buffer_pairs);
        Tools::parse(block_data, "PL_REFRESH", pl_refresh);

        if (!buffer_pairs.empty())
        {
          Buffer_Pairs[n] = std::stoi(buffer_pairs);
        }
        else
        {
          Buffer_Pairs[n] = 0;
        }

        if (!pl_refresh.empty())
        {
          nstride[n] = std::stoi(pl_refresh);
        }
        else
        {
          nstride[n] = 25;
        }
        tot_num_pairs[n] = block_lengths[n] + Buffer_Pairs[n];
        // logMsg("N: " + std::to_string(n), "Constructor");
        // logMsg("Nstride: " + std::to_string(nstride[n]), "Constructor");
        // logMsg("Buffer_pairs: " + std::to_string(Buffer_Pairs[n]), "Constructor");
        // logMsg("Block_lengths: " + std::to_string(block_lengths[n]), "Constructor");
        // logMsg("Total_num_pairs: " + std::to_string(tot_num_pairs[n]), "Constructor");
      }

      total_PIV_length = std::accumulate(block_lengths.begin(), block_lengths.end(), 0);
      // logMsg("total_PIV_length: " + std::to_string(total_PIV_length), "Constructor");
      for (int n = 0; n < N_Blocks; n++)
      {
        // pseudo-code ish
        for (int g = 0; g < 2; g++)
        {
          if (ID_list[n][g].size() > 0)
          {
            input_filters[n][g][0] = true;
            //logMsg("n: " + std::to_string(n) + "; g: " + std::to_string(g) + "-- ID Keyword", "Constructor");
          }
          if (ResID_list[n][g].size() > 0)
          {
            input_filters[n][g][1] = true;
            //logMsg("n: " + std::to_string(n) + "; g: " + std::to_string(g) + "-- ResID Keyword", "Constructor");
          }
          if (Name_list[n][g].size() > 0)
          {
            input_filters[n][g][2] = true;
            //logMsg("n: " + std::to_string(n) + "; g: " + std::to_string(g) + "-- Name Keyword", "Constructor");
          }
        }
      }

      for (int n=0; n < N_Blocks; n++) listall[n].clear();

      for (int i=0; i < mypdb.getAtomNumbers().size(); i++)
      {
        AtomNumber ind = mypdb.getAtomNumbers()[i];
        int resid = mypdb.getResidueNumber(ind);
        string atom_name = mypdb.getAtomName(ind);
        for (int n = 0; n < N_Blocks; n++)
        {
          bool atom_added = false;
          for (int g = 0; g < 2; g++)
          {
            if (atomMatchesFilters(n, g, ind, resid, atom_name)) {
              block_groups_atom_list[n][g].push_back(ind);
              //logMsg("n: " + std::to_string(n) + "; g: " + std::to_string(g) + ";  block_groups_atom_list[n][g]: " + std::to_string(ind.index()), "Constructor");
              atom_added = true;
            }
          }
          if (atom_added) listall[n].push_back(ind);
        }
      }

      int g1g2_pairs;
      for (int n = 0; n < N_Blocks; n++)
      {
        g1g2_pairs = block_groups_atom_list[n][0].size() * block_groups_atom_list[n][1].size() - Exclude_Pairs[n].size();
        if (g1g2_pairs == tot_num_pairs[n])
        {
          all_g1g2_pairs[n] = 1;
        } 
        else
        {
          all_g1g2_pairs[n] = 0;
        }
      }

      r00.resize(N_Blocks);
      sw.resize(N_Blocks);
      sfs.resize(N_Blocks);

      for (unsigned n = 0; n < N_Blocks; n++) if (!parseNumbered("SWITCH", n + 1, sw[n])) break;
      
      std::string errors;
      for (unsigned n = 0; n < N_Blocks; n++)
      {
        sfs[n].set(sw[n], errors);
        std::string num;
        Tools::convert(n + 1, num);
        if (errors.length() != 0) error("problem reading SWITCH" + num + " keyword : " + errors);
        r00[n] = sfs[n].get_r0();
        //logMsg("n: " + std::to_string(n) + "; r00: " + std::to_string(r00[n]), "Constructor");
      }
      checkRead();

      int total_count = 0;
      for (int n = 0; n < N_Blocks; n++)
      {
        for (int i = 0; i < block_lengths[n]; i++)
        {
          string comp = "ELEMENT-" + to_string(total_count);
          addComponentWithDerivatives(comp);
          componentIsNotPeriodic(comp);
          total_count += 1;
        }
      }
    }

    void PINES::prepare(){
      bool heap_refreshed = false;
      bool stale_refresh_prep = false;
      for (int n = 0; n < N_Blocks; n++){
        if (steps_since_update[n] == 0 || plumed.getStep() == 0 || driver_mode || (int)vecMaxHeapVecs[n].size() != tot_num_pairs[n]){
          buildMaxHeapVecBlock(n, mypdb, vecMaxHeapVecs[n]);
          heap_refreshed = true;
        }
        else if((steps_since_update[n] >= nstride[n] || stale_tolerance[n]) && !all_g1g2_pairs[n]){
          //logMsg("Pair List " + std::to_string(n) + " is stale. All g1/g2 atoms are being requested to rebuild pair list from scratch.", "Prepare");
          //logMsg( "Nstride: " + std::to_string(nstride[n]) + "; Steps Since Update: " + std::to_string(steps_since_update[n]), "Prepare");
          listreduced[n] = listall[n];
          stale_refresh_prep = true;
          steps_since_update[n] = -1;
          stale_tolerance[n] = false;
        }
        steps_since_update[n]+=1;
      }
      if (heap_refreshed || stale_refresh_prep){
        // Collate all atoms from all lists
        listreducedall.clear();
        for (const auto& blockList : listreduced) {        // Loop over each block's vector
          listreducedall.insert(blockList.begin(), blockList.end()); // Insert all atoms into set
        }
        atom_ind_hashmap.clear();
        listreducedall_vec = std::vector<AtomNumber>(listreducedall.begin(), listreducedall.end());
        //logMsg("listreducedall_vec size: " + std::to_string(listreducedall_vec.size()), "Prepare");
        for (int i=0; i < listreducedall_vec.size(); i++) atom_ind_hashmap[listreducedall_vec[i].index()] = i;
        requestAtoms(listreducedall_vec);
        ann_deriv.resize(listreducedall_vec.size());

        for (int i=0; i < ann_deriv.size(); i++) ann_deriv[i].resize(total_PIV_length);
        
      }
    }
      
    void PINES::calculate()
    {
      const long step = plumed.getStep();
      if (inited_step == -1) {
	for (int n =0; n < N_Blocks; n++) buildMaxHeapVecBlock(n, mypdb, vecMaxHeapVecs[n]);
      }
      if (step != inited_step) {
        if ((int)preupdated_block.size() != N_Blocks) preupdated_block.assign(N_Blocks, 0);
        else std::fill(preupdated_block.begin(), preupdated_block.end(), 0);

        // One pass: ensure any block that should update does so now (once per step)
        for (int n = 0; n < N_Blocks; n++) ensureBlockUpdated(n);

        // Optional: latch here (membership/order frozen for this step)
        if (freeze_selection) latchFromCurrentHeaps();

	for (int n = 0; n < N_Blocks; n++){
	  logMsg("_________________N: "+std::to_string(n)+"____________________\n","Calculate");
	  for (int p = 0; p < vecMaxHeapVecs[n].size(); p++) logMsg("P: "+std::to_string(vecMaxHeapVecs[n][p].first)+"\n","Calculate");
	}
        inited_step = step;
      }

      Vector ref_xyz, step_xyz;
      AtomNumber aID;
      double delta_r;

      //timer.start("toleranceCheck");
      
      for(int n = 0; n < N_Blocks; n++){
        if(steps_since_update[n] > 0 && steps_since_update[n] < nstride[n] && r_tolerance[n] > 0.0 && !all_g1g2_pairs[n]){
          for(int i = 0; i < listreduced[n].size(); i++){
            aID = listreduced[n][i];
            ref_xyz = PL_atoms_ref_coords[n][i];
            step_xyz = getPosition(atom_ind_hashmap[aID.index()]);
            delta_r = pbcDistance(ref_xyz, step_xyz).modulo();
            if(delta_r >= r_tolerance[n]){
              stale_tolerance[n] = true;
              //logMsg("Stale tolerance triggered for atom: " + std::to_string(aID.index()), "Calculate");
              //logMsg("delta_r: " + std::to_string(delta_r) + "; r_tolerance: " + std::to_string(r_tolerance[n]), "Calculate");
            }
          }
        }
      }
      //timer.stop("toleranceCheck");
      const Vector zeroVec(0., 0., 0.);
      //timer.start("setToZero");
      #pragma omp parallel for schedule(static)
      for (unsigned j = 0; j < ann_deriv.size(); j++) {
        for (unsigned i = 0; i < ann_deriv[j].size(); i++) {
          ann_deriv[j][i] = zeroVec;
        }
      }

      #pragma omp parallel for schedule(static)
      for (unsigned n = 0; n < N_Blocks; n++) {
        PIV[n].resize(block_lengths[n]);
        for (unsigned i = 0; i < block_lengths[n]; i++) {
          PIV[n][i] = 0.;
        }
      }


      //timer.stop("setToZero");



      //timer.start("bigCalc");
      int PINES_element = 0;
      if (freeze_selection) {
        for (int n = 0; n < N_Blocks; ++n) {
	  //ensureBlockUpdated(n);
	  logMsg("***************N: "+std::to_string(n)+"***************\n","bigCalc");
          for (int i = 0; i < block_lengths[n]; ++i) {
            const auto& pr = latched_pairs[n][i];
            int a0 = atom_ind_hashmap[pr.first.index()];
            int a1 = atom_ind_hashmap[pr.second.index()];
            Vector r0 = getPosition(a0);
	    Vector r1 = getPosition(a1);
	    double mod_dist = pbcDistance(r0, r1).modulo();
            Vector dr = pbcDistance(r0, r1) / mod_dist;
	    logMsg("P: "+std::to_string(vecMaxHeapVecs[n][i].first)+"\n","bigCalc");

            double dfunc = 0.0;
            PIV[n][i] = sfs[n].calculate(mod_dist, dfunc);
            double ds = dfunc * mod_dist;
            ann_deriv[a0][PINES_element] = -ds * dr;
            ann_deriv[a1][PINES_element] =  ds * dr;
            PINES_element += 1;
          }
        }
      } else{
      for (unsigned n = 0; n < N_Blocks; n++) {

        ensureBlockUpdated(n);
	    bool updated = true;
        // }
	// //timer.stop("updateBPL");
        //logMsg("Switching Function Parameters: " + sfs[n].description(), "Calculate");
        //logMsg("Steps Since Update: " + std::to_string(steps_since_update[n]), "Calculate");
        for (int i=0; i < block_lengths[n]; i++) {
          double ds_element = 0.;
          double dfunc = 0.;
          int local_aid0, local_aid1;
          int piv_ind;
          if (updated){
            piv_ind = i;
          }
          else{
            piv_ind = i + Buffer_Pairs[n];
          }
          //logMsg("Pair Distance: Dist[" + std::to_string(n) + "][" + std::to_string(i) + "]= " + std::to_string(vecMaxHeapVecs[n][i].first), "Calculate");
	  //timer.start("sfs");
          PIV[n][i] = sfs[n].calculate(vecMaxHeapVecs[n][piv_ind].first, dfunc);
          //timer.stop("sfs");
	  //logMsg("PIV Values: PIV[" + std::to_string(n) + "][" + std::to_string(i) + "]= " + std::to_string(PIV[n][i]), "Calculate");
	  //timer.start("distCalc");
          local_aid0 = atom_ind_hashmap[vecMaxHeapVecs[n][piv_ind].second.first.index()];
          local_aid1 = atom_ind_hashmap[vecMaxHeapVecs[n][piv_ind].second.second.index()];
          Vector Pos0 = getPosition(local_aid0);
          Vector Pos1 = getPosition(local_aid1);
          Vector dr_dcoord = pbcDistance(Pos0, Pos1) / vecMaxHeapVecs[n][piv_ind].first;
	  //timer.stop("distCalc");
	  //timer.start("derivCalc");
          ds_element = dfunc * vecMaxHeapVecs[n][piv_ind].first;
          ann_deriv[local_aid0][PINES_element] = -ds_element * dr_dcoord;
          ann_deriv[local_aid1][PINES_element] = ds_element * dr_dcoord;
	  //timer.stop("derivCalc");
          //logMsg("ann_deriv[" + std::to_string(local_aid0) + "][" + std::to_string(PINES_element) + "]: ", "Calculate");
          //logMsg(ann_deriv[local_aid0][PINES_element], "Calculate");
          //logMsg("ann_deriv[" + std::to_string(local_aid1) + "][" + std::to_string(PINES_element) + "]: ", "Calculate");
          //logMsg(ann_deriv[local_aid1][PINES_element], "Calculate");
          PINES_element += 1;
        }
      }
      }

      //timer.stop("bigCalc");
      //timer.start("commStuff");
      //timer.start("commStuff");

      if (comm.initialized()) {
        // Flatten PIV
        std::vector<double> flat_PIV;
        std::vector<int> piv_sizes;
        for (const auto& vec : PIV) {
          piv_sizes.push_back(vec.size());
          flat_PIV.insert(flat_PIV.end(), vec.begin(), vec.end());
        }
      
        // Flatten ann_deriv
        std::vector<double> flat_deriv;
        std::vector<std::pair<int, int>> deriv_shapes;
        for (const auto& atom_deriv : ann_deriv) {
          deriv_shapes.emplace_back(atom_deriv.size(), 3);
          for (const auto& val : atom_deriv) {
            flat_deriv.insert(flat_deriv.end(), val[0]);
            flat_deriv.insert(flat_deriv.end(), val[1]);
            flat_deriv.insert(flat_deriv.end(), val[2]);
          }
        }
      
        // Perform global summation ONCE per array
        comm.Sum(flat_PIV);
        comm.Sum(flat_deriv);
      
        // Unflatten PIV
        {
          size_t idx = 0;
          for (size_t j = 0; j < N_Blocks; ++j) {
            for (int i = 0; i < piv_sizes[j]; ++i) {
              PIV[j][i] = flat_PIV[idx++];
            }
          }
        }
      
        // Unflatten ann_deriv
        {
          size_t idx = 0;
          for (size_t i = 0; i < ann_deriv.size(); ++i) {
            int ncols = deriv_shapes[i].first;
            for (int j = 0; j < ncols; ++j) {
              ann_deriv[i][j][0] = flat_deriv[idx++];
              ann_deriv[i][j][1] = flat_deriv[idx++];
              ann_deriv[i][j][2] = flat_deriv[idx++];
            }
          }
        }
      }


      //timer.stop("commStuff");



      //timer.start("valuePass");
      // Pass values and derivates to next stage
      unsigned total_count = 0;
      for (unsigned j = 0; j < N_Blocks; j++)
      {
        for (unsigned i = 0; i < block_lengths[j]; i++)
        {
          string comp = "ELEMENT-" + to_string(total_count);
          Value *valueNew = getPntrToComponent(comp);
          valueNew->set(PIV[j][i]);
          for (unsigned k = 0; k < ann_deriv.size(); k++)
          {
            setAtomsDerivatives(valueNew, k, ann_deriv[k][total_count]);
          }
          total_count += 1;
        }
      }
      //timer.stop("valuePass");

      const int Natoms = static_cast<int>(ann_deriv.size());

// Precompute current box and fractional (scaled) coords once per step
      Tensor box = getPbc().getBox();
      std::vector<Vector> scaled(Natoms);
      for (int a = 0; a < Natoms; ++a) {
        scaled[a] = getPbc().realToScaled(getPosition(a));
      }

// Build one virial tensor per component (column)
      std::vector<Tensor> virials(total_PIV_length);
      for (int col = 0; col < total_PIV_length; ++col) {
        Tensor dSdB; dSdB.zero();                 // ∂S/∂B at fixed scaled coords

  // dSdB = sum_i  s_i ⊗ g_i  with s_i = scaled (fractional), g_i = ∂S/∂r_i
        for (int a = 0; a < Natoms; ++a) {
          const Vector& gi = ann_deriv[a][col];   // atomic gradient for this component
          if (gi[0]==0.0 && gi[1]==0.0 && gi[2]==0.0) continue;
          const Vector& si = scaled[a];

    // Outer product: (∂S/∂B)_{ik} += s_i[i] * g_i[k]
          dSdB(0,0) += si[0]*gi[0];  dSdB(0,1) += si[0]*gi[1];  dSdB(0,2) += si[0]*gi[2];
          dSdB(1,0) += si[1]*gi[0];  dSdB(1,1) += si[1]*gi[1];  dSdB(1,2) += si[1]*gi[2];
          dSdB(2,0) += si[2]*gi[0];  dSdB(2,1) += si[2]*gi[1];  dSdB(2,2) += si[2]*gi[2];
        }

  // PLUMED numerical uses: virial = - B^T * (∂S/∂B)
        virials[col] = - matmul(box.transpose(), dSdB);
      }

      for (int col = 0; col < total_PIV_length; ++col) {
  // If you already have the Value* via getPntrToComponent("ELEMENT-<col>"), reuse it.
  // Otherwise copyOutput(col) is fine here.
        Value* v = copyOutput(col);

        for (int i = 0; i < 3; ++i) {      // i = column
          for (int k = 0; k < 3; ++k) {    // k = row
            v->addDerivative(3*Natoms + 3*k + i, virials[col](k,i));
          }
        }
      }
      //log << timer;
    }
  }
}
