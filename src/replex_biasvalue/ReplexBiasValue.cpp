/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2026 Olivier Languin-Cattoën
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
#include "bias/Bias.h"
#include "core/ActionRegister.h"
#include "core/ActionSet.h"
#include "core/PlumedMain.h"
#include "tools/Communicator.h"
#include "tools/IFile.h"
#include "tools/Threefry.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <numeric>

namespace PLMD {
namespace replex_biasvalue {

//+PLUMEDOC BIAS REPLEX_BIASVALUE
/*
In a multi-simulation context, takes the values of variables and uses them as
exchangeable biases between replicas.

This is analogous to multiple [BIASVALUE](BIASVALUE.md) where biases are
periodically exchanged using the Metropolis algorithm. The design is consistent
with replica exchange in GROMACS, except that the replicas exchange biases
instead of system coordinates. As a result, output trajectories are
time-continuous and must be sorted ("multiplexed") to recover the correct biased
ensembles (as opposed to the usual "demultiplexing", or "demuxing", which recovers
continuous trajectories from discontinuous ones).

The exchange scheme matches the GROMACS one, where exchanges between neighboring
biases are attempted alternating odd pairs ($0 \leftrightarrow 1$,
$2 \leftrightarrow 3$, …) and even pairs ($1 \leftrightarrow 2$, …). Exchanges
happen at every step that is a multiple of PACE, except at step 0. Multiple
attempts can be made on the same exchange step (NEX keyword), alternating odd
and even pairs, to approximate independence sampling in state space. This incurs
little additional cost since replicas already compute all biases on exchange
steps. Using NEX=2 ensures all neighboring pairs are attempted once per exchange
step.

By default the number of variables $N$ provided with ARG should match the number
of replicas $M$ in the multi-simulation, and will correspond to a different bias
each. Alternatively, the bias ladder can be constructed as any linear
combination of the $N$ input variables. In that case, the COEFFICIENTS keyword
should specify the appropriate $M \times N$ matrix $(c_{ij})$ (see examples below)
and biases $i = 1, ..., M$ will be computed as

$$
V_i = \sum_j^N c_{ij} x_j.
$$

One can also perform multidimensional replica exchange on a grid. In this case
one must provide a SHAPE keyword specifying the size of each of the $D$
dimensions. The product of SHAPE should match the number of replicas $M$. As for
the 1D case, COEFFICIENTS for each bias should be given as an $M \times N$
matrix, with one row per bias. Biases are numbered with C-style, row-major
ordering on the grid (last dimension is the fastest-varying index), see examples
below. For each exchange attempt (odd or even), all dimensions are attempted in
order starting from the last one.

!!! warning
    All replicas are assumed to share the *same* PLUMED biases, forcefield
    topologies and temperatures. REPLEX_BIASVALUE is unaware of differences
    outside its own biases and thus cannot include them in the Metropolis
    acceptance. Keywords ARG, COEFFICIENTS, SHAPE, PACE, NEX, SEED, and TEMP
    must be identical in all replicas and will be checked. We recommend using a
    single common PLUMED input file and avoiding the `@replicas:` syntax.

!!! note
    Only one REPLEX_BIASVALUE action is allowed in each PLUMED input. Combine
    all the bias terms in its ARG list and COEFFICIENTS matrix.

!!! note
    This action cannot be combined with replica exchange provided by the MD
    code. For example, when running with GROMACS, use
    ```bash
    mdrun -multidir … -plumed …
    ```
    *without* `-replex`.


## Examples

One use case is the interpolation between two biases:

```plumed
#SETTINGS NREPLICAS=4
#HIDDEN
cv1: DISTANCE ATOMS=1,2
cv2: DISTANCE ATOMS=3,4
#ENDHIDDEN

# 4 replicas switching from cv1 to cv2
bias: REPLEX_BIASVALUE ...
  ARG=cv1,cv2
  COEFFICIENTS={
    {1.0 0.0}
    {0.7 0.3}
    {0.1 0.9}
    {0.0 1.0}
  }
  PACE=100
  NEX=10
  SEED=42
...
```

If $N=1$, a vector can be given to imply different scalings of the same bias
function:

```plumed
#SETTINGS NREPLICAS=6
#HIDDEN
cv: DISTANCE ATOMS=1,2
#ENDHIDDEN

bias: REPLEX_BIASVALUE ARG=cv COEFFICIENTS=0.0,0.2,0.4,0.6,0.8,1.0 PACE=100
# equivalent to COEFFICIENTS={{0.0} {0.2} {0.4} {0.6} {0.8} {1.0}}
```

If both ARG and COEFFICIENTS are vectors of $M$ elements, a diagonal matrix is
implied:

```plumed
#SETTINGS NREPLICAS=5
#HIDDEN
cv1: DISTANCE ATOMS=1,2
cv2: DISTANCE ATOMS=3,4
cv3: DISTANCE ATOMS=5,6
cv4: DISTANCE ATOMS=7,8
cv5: DISTANCE ATOMS=9,10
#ENDHIDDEN

bias: REPLEX_BIASVALUE ...
  ARG=cv1,cv2,cv3,cv4,cv5
  COEFFICIENTS=0.8,1.2,0.5,0.9,1.0
  # equivalent to COEFFICIENTS={
  #   {0.8    0    0    0    0}
  #   {  0  1.2    0    0    0}
  #   {  0    0  0.5    0    0}
  #   {  0    0    0  0.9    0}
  #   {  0    0    0    0  1.0}
  # }
  PACE=100
  NEX=100
  SEED=1312
...
```


Two-dimensional ladder—the coefficients are given for each bias in row-major
ordering:

```plumed
#SETTINGS NREPLICAS=6
#HIDDEN
cv1: DISTANCE ATOMS=1,2
cv2: DISTANCE ATOMS=3,4
#ENDHIDDEN

bias: REPLEX_BIASVALUE ...
  ARG=cv1,cv2
  SHAPE=2,3
  COEFFICIENTS={
    {0   0}
    {0 0.5}
    {0   1}
    {1   0}
    {1 0.5}
    {1   1}
  }
  # row-major ordering:
  # bias[0,0]
  # bias[0,1]
  # bias[0,2]
  # bias[1,0]
  # bias[1,1]
  # bias[1,2]
  PACE=100
...
```


## Output files

The first replica (index 0) keeps track of the position of the biases in the
output FILE (`bias_index.xvg` by default). Column `replica-i` is the index of
the bias on replica $i$.

The bias indices can be used to "multiplex" the continuous trajectories (with
variable biases) into discontinuous trajectories (with fixed biases). Please
note that exchanges are attempted *after* calculating forces: this means that
the PLUMED outputs at time $t$ correspond to bias indices *strictly before* $t$
in the `bias_index.xvg` file. This is consistent with the output of the
`demux.pl` script provided with GROMACS and can be seamlessly used for
multiplexing:

```bash
gmx trjcat -f rep?/traj.xtc -demux rep0/bias_index.xvg -o rep?/muxed.xtc
```

!!! warning
    For long simulations the time reported in `bias_index.xvg` might be inexact
    due to floating point precision and may be inconsistent with GROMACS
    trajectories.

For convenience, it is also possible to output the replica indices in bias space
by setting REPINDEX_FILE (note that going from replica indices to bias indices
and vice versa corresponds to finding the inverse permutation, that is the
permutation that sorts the indices, e.g. `np.argsort`).

```bash
# This should be the inverse operation to the previous one
gmx trjcat -f rep?/muxed.xtc -demux rep0/replica_index.xvg -o rep?/traj.xtc
```

## Log output

Each exchange step logs the new indices. At the end of the run, per-pair
statistics are given (accepted/attempted, acceptance ratio, average acceptance
probability), counted over this run only. The VERBOSE flag can be used to output
statistics after each exchange step (average over NEX attempts).


## Restart

When restarting a simulation (as indicated internally by the MD engine or
enforced with RESTART=YES) the content of FILE is automatically read to recover
the state at a time consistent with the checkpoint. The random numbers used to
accept or reject the exchanges depend only on SEED and on the step, so the
restarted simulation makes the same decisions as an uninterrupted one would
(assuming deterministic and bit-identical energies). This also holds when
restarting from a PLUMED [CHECKPOINT](CHECKPOINT.md), which needs nothing
specific to this action.

If needed one can run a "fresh" simulation starting from a state provided by
BIASINDEX_RFILE, in which case the time consistency check is ignored and the
very last line is used.
*/
//+ENDPLUMEDOC

class ReplexBiasValue : public bias::Bias {
private:
  // Statistics of the exchange attempts of each pair of neighbors: index
  // ibias + d*nrep is the pair (ibias, ibias + strides[d]) along dimension d
  struct PairStats {
    std::vector<double> sum_proba; // sum of the acceptance probabilities
    std::vector<long long> attempts;
    std::vector<long long> accepted;

    void resize(unsigned n) {
      sum_proba.assign(n, 0.0);
      attempts.assign(n, 0);
      accepted.assign(n, 0);
    }

    void record(unsigned ipair, double pacc, bool acc) {
      plumed_dbg_assert(ipair < sum_proba.size());
      sum_proba[ipair] += pacc;
      attempts[ipair]++;
      if (acc) {
        accepted[ipair]++;
      }
    }

    PairStats operator-(const PairStats &other) const {
      plumed_assert(sum_proba.size() == other.sum_proba.size());
      PairStats out(*this);
      for (size_t i = 0; i < sum_proba.size(); i++) {
        out.sum_proba[i] -= other.sum_proba[i];
        out.attempts[i] -= other.attempts[i];
        out.accepted[i] -= other.accepted[i];
      }
      return out;
    }
  };

  // Settings, identical on all the replicas (checked in the constructor)
  int nrep = 1;
  int narg = 0;
  int repid = 0; // index of this replica
  double kbt = 0;
  int pace = 0;
  int nex = 1;
  int seed = 0;
  bool verbose = false;
  // coefficients[ibias*narg + j]: coefficient of argument j in bias ibias
  std::vector<double> coefficients;
  // grid of the biases (multidimensional replica exchange), product = nrep
  std::vector<int> shape;
  int dim = 1; // shape.size()
  // row-major strides: ibias = sum over d of coordinate(ibias, d) * strides[d]
  std::vector<int> strides;

  // State of the exchanges
  std::vector<int> bias_indices;    // [irep]: bias currently applied on replica irep
  std::vector<int> replica_indices; // [ibias]: replica that has bias ibias (inverse permutation)
  Threefry2x64 rng;
  PairStats stats; // over this run

  // Flags of the current step
  bool first_step = true;
  bool is_exchange_step = false;
  bool update_dependencies = true;

  // Input/output (the files are written by replica 0)
  std::string biasindex_file;
  std::string biasindex_rfile;
  std::string repindex_file;
  bool read_initial_state = false;
  OFile biasindex_ofile;
  OFile repindex_ofile;
  Value *value_index = nullptr;

  // Cartesian coordinate of bias ibias along dimension d
  int coordinate(int ibias, int d) const {
    return (ibias / strides[d]) % shape[d];
  }
  // Broadcast from replica 0 to all the ranks of all the replicas
  template <typename T> void bcastFromReplica0(T &x) {
    if (comm.Get_rank() == 0) {
      multi_sim_comm.Bcast(x, 0);
    }
    comm.Bcast(x, 0);
  }
  void parseCoefficients();
  void readInitialState();
  std::string readBiasIndices();
  void doExchanges();
  void writeIndices();
  void logIndices() const;
  void logStats(const PairStats &s) const;
  void logGrid(const std::vector<std::string> &cells) const;
public:
  static void registerKeywords(Keywords &keys);
  explicit ReplexBiasValue(const ActionOptions &ao);
  ~ReplexBiasValue() override;
  void prepare() override;
  void calculate() override;
  void update() override;
};

PLUMED_REGISTER_ACTION(ReplexBiasValue, "REPLEX_BIASVALUE")

void ReplexBiasValue::registerKeywords(Keywords &keys) {
  Bias::registerKeywords(keys);
  keys.addInputKeyword("compulsory", "ARG", "scalar",
                       "the labels of the scalars combined to define the biases");
  keys.add("compulsory", "COEFFICIENTS", "1.0",
           "the coefficients of the arguments in each bias, one row per bias or "
           "one value per bias (see examples)");
  keys.add("optional", "SHAPE",
           "the shape of the grid of the biases, for multidimensional replica exchange");
  keys.add("compulsory", "PACE", "the frequency of the exchange attempts");
  keys.add("compulsory", "NEX", "1",
           "the number of sweeps over the pairs of neighbors on each exchange step, "
           "alternating odd and even pairs. Many sweeps (e.g. M^3, where M is the "
           "number of replicas) approximate independence sampling");
  keys.add("compulsory", "SEED", "0",
           "the seed of the random numbers used to accept or reject the exchanges, "
           "the same on all replicas");
  keys.add("optional", "TEMP", "the temperature, if not given by the MD code");
  keys.add("compulsory", "FILE", "bias_index.xvg", "the output file of the bias indices");
  keys.add("optional", "BIASINDEX_RFILE",
           "a file of bias indices to start from, if different from FILE: "
           "its last line is used, without checking the time");
  keys.add("optional", "REPINDEX_FILE", "an output file of the replica indices");
  keys.addFlag("VERBOSE", false,
               "write in the log the statistics of every exchange step");
  keys.use("RESTART");
  keys.addOutputComponent("_bias", "default", "scalar",
                          "the contribution of each argument to the bias, "
                          "scaled by the coefficient, "
                          "named after the argument");
  keys.addOutputComponent("index", "default", "scalar",
                          "the index of the bias currently applied on this replica");
  keys.addDOI("10.1063/1.3660669");
  keys.addDOI("10.1063/1.3086832");
}

ReplexBiasValue::ReplexBiasValue(const ActionOptions &ao)
  : PLUMED_BIAS_INIT(ao) {
  if (!plumed.getActionSet().select<ReplexBiasValue*>().empty()) {
    error("only one REPLEX_BIASVALUE action is allowed; combine its arguments in a single action");
  }
  nrep = getNumberOfReplicas();
  narg = getNumberOfArguments();
  repid = getReplicaIndex();

  parseCoefficients();

  parseVector("SHAPE", shape);
  if (shape.empty()) {
    shape.push_back(nrep);
  }
  int grid_size = 1;
  for (int k : shape) {
    if (k < 1) {
      error("shape elements should be strictly positive, not " +
            std::to_string(k));
    }
    grid_size *= k;
  }
  if (grid_size != nrep) {
    error("multidimensional grid shape incompatible with number of replicas: " +
          std::to_string(grid_size) + " != " + std::to_string(nrep));
  }

  parse("PACE", pace);
  if (pace <= 0) {
    error("frequency for exchanges should be strictly positive, not " +
          std::to_string(pace));
  }
  parse("NEX", nex);
  if (nex <= 0) {
    error("number of exchanges attempts should be strictly positive, not " +
          std::to_string(nex));
  }
  parse("SEED", seed);
  rng = Threefry2x64(static_cast<std::uint64_t>(seed));

  kbt = getkBT();
  if (!(kbt > 0)) {
    error("the temperature is not known: set TEMP or use an MD code that "
          "passes it to PLUMED (kBT=" + std::to_string(kbt) + ")");
  }

  parse("FILE", biasindex_file);
  parse("BIASINDEX_RFILE", biasindex_rfile);
  if (biasindex_rfile.empty()) {
    biasindex_rfile = biasindex_file;
  } else if (getRestart() && biasindex_rfile != biasindex_file) {
    error("when restarting a simulation FILE and BIASINDEX_RFILE should be "
          "identical");
  }
  read_initial_state = getRestart() || biasindex_rfile != biasindex_file;
  parse("REPINDEX_FILE", repindex_file);
  parseFlag("VERBOSE", verbose);
  checkRead();

  dim = shape.size();
  strides.assign(dim, 1);
  for (int d = dim - 1; d > 0; d--) {
    strides[d - 1] = strides[d] * shape[d];
  }
  replica_indices.resize(nrep);
  std::iota(replica_indices.begin(), replica_indices.end(), 0);
  bias_indices = replica_indices;
  // One pair for each bias along each dimension (the last ones are unused)
  stats.resize(nrep * dim);

  // Check that all the replicas agree
  std::vector<double> settings{double(narg), double(pace), double(nex), double(seed), kbt};
  settings.insert(settings.end(), shape.begin(), shape.end());
  settings.insert(settings.end(), coefficients.begin(), coefficients.end());
  int mismatch = 0; // int: MPI has no bool type in Communicator
  if (comm.Get_rank() == 0) {
    unsigned size = settings.size();
    multi_sim_comm.Bcast(size, 0);
    mismatch = (size != settings.size());
    multi_sim_comm.Sum(mismatch);
    if (!mismatch) {
      std::vector<double> ref(settings);
      multi_sim_comm.Bcast(ref, 0);
      mismatch = (ref != settings);
      multi_sim_comm.Sum(mismatch);
      if (!mismatch) {
        // Include separators so that different ordered ARG lists cannot alias.
        std::string arguments;
        for (int j = 0; j < narg; j++) {
          arguments += getPntrToArgument(j)->getName();
          arguments += '\0';
        }
        unsigned length = arguments.size();
        multi_sim_comm.Bcast(length, 0);
        std::string reference = arguments;
        reference.resize(length);
        multi_sim_comm.Bcast(reference, 0);
        mismatch = (reference != arguments);
        multi_sim_comm.Sum(mismatch);
      }
    }
  }
  comm.Bcast(mismatch, 0);
  if (mismatch) {
    error("ARG, PACE, NEX, SEED, SHAPE, COEFFICIENTS or temperature differ between replicas");
  }

  for (int j = 0; j < narg; j++) {
    const std::string name = getPntrToArgument(j)->getName() + "_bias";
    addComponent(name);
    componentIsNotPeriodic(name);
  }
  addComponent("index");
  componentIsNotPeriodic("index");
  value_index = getPntrToComponent("index");

  log.printf("  number of replicas: %d, index of this replica: %d\n", nrep, repid);
  log.printf("  coefficients of each argument in the biases:\n");
  for (int j = 0; j < narg; j++) {
    log.printf("    %s:", getPntrToArgument(j)->getName().c_str());
    for (int ibias = 0; ibias < nrep; ibias++) {
      log.printf("  %f", coefficients[ibias * narg + j]);
    }
    log.printf("\n");
  }
  if (dim > 1) {
    log.printf("  %d-dimensional grid of biases with shape:", dim);
    for (int k : shape) {
      log.printf(" %d", k);
    }
    log.printf("\n");
  }
  log.printf("  exchange attempts every %d steps, %d sweep(s) per exchange step",
             pace, nex);
  if (nex > 1) {
    log << " " << cite("Chodera and Shirts, J. Chem. Phys. 135, 194110 (2011)");
  }
  log.printf("\n");
  log.printf("  seed of the random numbers: %d\n", seed);
  log.printf("  kBT: %f\n", kbt);
  log.printf("  bias indices written to file: %s\n", biasindex_file.c_str());
  if (read_initial_state) {
    log.printf("  initial state read from file: %s\n", biasindex_rfile.c_str());
  }
  if (!repindex_file.empty()) {
    log.printf("  replica indices written to file: %s\n", repindex_file.c_str());
  }
  if (verbose) {
    log.printf("  statistics of every exchange step written in the log\n");
  }

  if (repid == 0) {
    biasindex_ofile.link(*this);
    biasindex_ofile.enforceSuffix("");
    biasindex_ofile.open(biasindex_file);
    biasindex_ofile.fmtField(" %f");
    if (!repindex_file.empty()) {
      repindex_ofile.link(*this);
      repindex_ofile.enforceSuffix("");
      repindex_ofile.open(repindex_file);
      repindex_ofile.fmtField(" %f");
    }
  }
}

void ReplexBiasValue::parseCoefficients() {
  // One row per bias (a single value, like the default, applies to all the biases)
  std::vector<std::string> rows(nrep);
  parseVector("COEFFICIENTS", rows);
  auto toDouble = [this](const std::string &word) {
    double c = 0.0;
    if (!Tools::convertNoexcept(word, c)) {
      error("cannot parse coefficient: " + word);
    }
    return c;
  };
  for (int ibias = 0; ibias < nrep; ibias++) {
    const std::vector<std::string> words = Tools::getWords(rows[ibias], " ,");
    const int nwords = words.size();
    if (nwords == narg) {
      for (const auto &word : words) {
        coefficients.push_back(toDouble(word));
      }
    } else if (nwords == 1 && nrep == narg) {
      // diagonal matrix
      const double c = toDouble(words[0]);
      for (int j = 0; j < narg; j++) {
        coefficients.push_back(j == ibias ? c : 0.0);
      }
    } else {
      error("COEFFICIENTS for bias " + std::to_string(ibias) + " has " +
            std::to_string(nwords) + " numbers: expected one per argument (" +
            std::to_string(narg) + ")" +
            (nrep == narg ? " or a single one (diagonal matrix)" : ""));
    }
  }
}

ReplexBiasValue::~ReplexBiasValue() {
  log.printf("%s %s: statistics of the exchanges over this run\n",
             getName().c_str(), getLabel().c_str());
  logStats(stats);
}

void ReplexBiasValue::prepare() {
  if (first_step) {
    if (read_initial_state) {
      readInitialState();
    }
    // A new run starts with a line for its initial state
    if (!getRestart()) {
      writeIndices();
    }
    first_step = false;
  }

  // Abort if external exchanges (e.g. gmx -replex)
  if (getExchangeStep()) {
    error("cannot be used together with the replica exchange of the MD code "
          "(e.g. GROMACS -replex)");
  }

  is_exchange_step = (getStep() != 0 && getStep() % pace == 0);
  if (is_exchange_step) {
    // All the biases are computed on exchange steps
    clearDependencies();
    for (int j = 0; j < narg; j++) {
      addDependency(getPntrToArgument(j)->getPntrToAction());
    }
  } else if (update_dependencies) {
    // Otherwise only the arguments of the current bias are needed
    clearDependencies();
    const int ibias = bias_indices[repid];
    for (int j = 0; j < narg; j++) {
      if (coefficients[ibias * narg + j] != 0.0) {
        addDependency(getPntrToArgument(j)->getPntrToAction());
      }
    }
    update_dependencies = false;
  }
}

void ReplexBiasValue::calculate() {
  double bias = 0.0;
  const int ibias = bias_indices[repid];
  for (int j = 0; j < narg; j++) {
    const double coef = coefficients[ibias * narg + j];
    // Arguments with a zero coefficient are not computed on non-exchange steps
    // (see prepare), so their values must not be used
    if (coef != 0.0) {
      const double val = coef * getArgument(j);
      getPntrToComponent(j + 1)->set(val); // NB: component 0 is total bias
      setOutputForce(j, -coef);
      bias += val;
    } else {
      getPntrToComponent(j + 1)->set(0.0);
    }
  }
  value_index->set(ibias);
  setBias(bias);
}

void ReplexBiasValue::update() {
  if (is_exchange_step) {
    doExchanges();
    update_dependencies = true;
  }
}

void ReplexBiasValue::doExchanges() {
  // all_reduced_biases[irep*nrep + ibias]: bias ibias on replica irep, over kBT
  std::vector<double> reduced_biases(nrep);
  for (int ibias = 0; ibias < nrep; ibias++) {
    for (int j = 0; j < narg; j++) {
      reduced_biases[ibias] += coefficients[ibias * narg + j] * getArgument(j) / kbt;
    }
  }
  std::vector<double> all_reduced_biases(reduced_biases.size() * nrep);
  if (comm.Get_rank() == 0) {
    multi_sim_comm.Allgather(reduced_biases, all_reduced_biases);
  }
  comm.Bcast(all_reduced_biases, 0);

  // Counter-based random numbers, as in GROMACS replica exchange: they depend
  // only on SEED and on the step (exact restarts without storing a state)
  rng.restart(static_cast<std::uint64_t>(getStep()));

  PairStats before;
  if (verbose) {
    before = stats;
  }

  // Sweeps alternate between the pairs (0,1),(2,3),... and (1,2),(3,4),...
  // starting with the former: (step/pace - 1)*nex sweeps were done before
  bool even_pairs = ((getStep() / pace) - 1) * nex % 2 == 0;
  for (int k = 0; k < nex; k++) {
    // Along each dimension, starting from the last one
    for (int d = dim - 1; d >= 0; d--) {
      for (int ibias = 0; ibias < nrep; ibias++) {
        const int pos = coordinate(ibias, d);
        // Is ibias the first bias of a pair attempted in this sweep?
        if (pos < shape[d] - 1 && (pos % 2 == 0) == even_pairs) {
          const int jbias = ibias + strides[d];
          const int irep = replica_indices[ibias];
          const int jrep = replica_indices[jbias];
          // Change of the total reduced energy if irep and jrep swap their biases
          const double dE = all_reduced_biases[irep * nrep + jbias] +
                            all_reduced_biases[jrep * nrep + ibias] -
                            all_reduced_biases[irep * nrep + ibias] -
                            all_reduced_biases[jrep * nrep + jbias];
          const double pacc = std::isnan(dE) ? 0.0 : std::min(1.0, std::exp(-dE));
          const bool acc = rng.U01() < pacc;
          stats.record(ibias + d * nrep, pacc, acc);
          if (acc) {
            replica_indices[ibias] = jrep;
            replica_indices[jbias] = irep;
          }
        }
      }
    }
    even_pairs = !even_pairs;
  }
  for (int ibias = 0; ibias < nrep; ibias++) {
    bias_indices[replica_indices[ibias]] = ibias;
  }

  log.printf("%s %s: exchange at step %lld, time %f\n", getName().c_str(),
             getLabel().c_str(), static_cast<long long>(getStep()), getTime());
  logIndices();
  if (verbose) {
    log.printf("  statistics of this exchange step (NEX=%d):\n", nex);
    logStats(stats - before);
  }
  writeIndices();
}

/// Reads the bias indices on replica 0 and broadcasts them to all the replicas.
void ReplexBiasValue::readInitialState() {
  log.printf("%s %s: reading the initial state\n", getName().c_str(),
             getLabel().c_str());
  std::string message;
  if (repid == 0) {
    message = readBiasIndices();
  }
  // The error is raised on all the replicas, otherwise the others would wait
  int failed = !message.empty();
  bcastFromReplica0(failed);
  if (failed) {
    if (repid == 0) {
      // The other replicas refer to this log: the error message itself might
      // never be printed, if another replica stops the whole run first
      log.printf("  ERROR: %s\n", message.c_str());
      log.flush();
      error(message);
    }
    error("cannot read the initial state from file " + biasindex_rfile +
          " (see the log of replica 0)");
  }
  bcastFromReplica0(bias_indices);
  for (int irep = 0; irep < nrep; irep++) {
    replica_indices[bias_indices[irep]] = irep;
  }
  logIndices();
}

/// Reads the bias indices from BIASINDEX_RFILE (on replica 0 only).
/// Returns an error message, empty on success.
std::string ReplexBiasValue::readBiasIndices() {
  try {
    // Exchanges happen in update(), after the output of the step: the state to
    // restart from is the last line strictly before the current time, and at
    // most one PACE earlier. With a different BIASINDEX_RFILE, the last line.
    const bool check_time = (biasindex_rfile == biasindex_file);
    const double now = getTime();
    const double eps = 0.5 * getTimeStep();
    const double window = pace * getTimeStep();

    IFile ifile;
    ifile.link(*this);
    ifile.enforceSuffix("");
    if (!ifile.FileExist(biasindex_rfile)) {
      return "cannot find file " + biasindex_rfile;
    }
    ifile.open(biasindex_rfile);

    bool found = false;
    double time = 0, rtime = 0;
    std::vector<int> indices(nrep), line(nrep);
    // All the fields of every line are read, so that a file written with a
    // different number of replicas is always detected
    while (ifile.scanField("time", time)) {
      for (int i = 0; i < nrep; i++) {
        ifile.scanField("replica-" + std::to_string(i), line[i]);
      }
      ifile.scanField();
      if (!check_time || (time + eps < now && time + window + eps > now)) {
        indices = line;
        rtime = time;
        found = true;
      }
    }
    ifile.close();

    if (!found) {
      return "no state compatible with time " + std::to_string(now) +
             " found in file " + biasindex_rfile;
    }
    std::vector<int> identity(nrep);
    std::iota(identity.begin(), identity.end(), 0);
    if (!std::is_permutation(indices.begin(), indices.end(), identity.begin())) {
      return "the indices read from file " + biasindex_rfile +
             " are not a permutation of " + std::to_string(nrep) + " replicas";
    }
    bias_indices = indices;
    log.printf("  read from file %s, line at time %f\n", biasindex_rfile.c_str(),
               rtime);
  } catch (const std::exception &e) {
    return "cannot read the state from file " + biasindex_rfile + ": " +
           e.what();
  }
  return "";
}

/// Writes the current indices to the output files (on replica 0 only).
void ReplexBiasValue::writeIndices() {
  if (repid != 0) {
    return;
  }
  auto write = [this](OFile &file, const std::vector<int> &indices,
  const std::string &prefix) {
    file.printField("time", getTime());
    for (size_t i = 0; i < indices.size(); i++) {
      file.printField(prefix + std::to_string(i), indices[i]);
    }
    file.printField();
  };
  write(biasindex_ofile, bias_indices, "replica-");
  if (repindex_ofile.isOpen()) {
    write(repindex_ofile, replica_indices, "bias-");
  }
}

void ReplexBiasValue::logIndices() const {
  log.printf("  bias indices:");
  for (int i : bias_indices) {
    log.printf(" %d", i);
  }
  log.printf("\n  replica indices:");
  for (int i : replica_indices) {
    log.printf(" %d", i);
  }
  log.printf("\n");
}

/// Logs the statistics of each pair of neighbors: as grids for 1D and 2D
/// ladders, as a flat list (attempted pairs only) for higher dimensions.
void ReplexBiasValue::logStats(const PairStats &s) const {
  auto fmt = [](double x) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.3f", x);
    return std::string(buf);
  };
  const unsigned npairs = s.sum_proba.size();
  std::vector<std::string> fraction(npairs), ratio(npairs), average(npairs);
  for (unsigned ipair = 0; ipair < npairs; ipair++) {
    fraction[ipair] = std::to_string(s.accepted[ipair]) + "/" +
                      std::to_string(s.attempts[ipair]);
    if (s.attempts[ipair] > 0) {
      ratio[ipair] = fmt(static_cast<double>(s.accepted[ipair]) / s.attempts[ipair]);
      average[ipair] = fmt(s.sum_proba[ipair] / s.attempts[ipair]);
    } else {
      ratio[ipair] = "-";
      average[ipair] = "-";
    }
  }

  if (dim <= 2) {
    log.printf("  accepted/attempted:\n");
    logGrid(fraction);
    log.printf("  acceptance ratio:\n");
    logGrid(ratio);
    log.printf("  average probability:\n");
    logGrid(average);
    return;
  }

  // Flat list of the pairs (i, j = i + strides[d]), sorted by (i, j)
  if (std::accumulate(s.attempts.begin(), s.attempts.end(), 0LL) == 0) {
    log.printf("  no exchange attempted\n");
    return;
  }
  const int wi = std::to_string(nrep - 1).size();
  int wf = 0;
  for (const auto &f : fraction) {
    wf = std::max(wf, static_cast<int>(f.size()));
  }
  log.printf("  # [i] [j] accepted/attempted ratio average_proba\n");
  for (int ibias = 0; ibias < nrep; ibias++) {
    for (int d = dim - 1; d >= 0; d--) {
      const int ipair = ibias + d * nrep;
      if (coordinate(ibias, d) < shape[d] - 1 && s.attempts[ipair] > 0) {
        log.printf("  [%*d] [%*d] %*s %s %s\n", wi, ibias, wi, ibias + strides[d],
                   wf, fraction[ipair].c_str(), ratio[ipair].c_str(),
                   average[ipair].c_str());
      }
    }
  }
}

/// Logs one value per pair of neighbors (cells[ibias + d*nrep] is the pair
/// (ibias, ibias + strides[d])) on the grid of biases, for 1D and 2D ladders:
///   [0]  c01   [1]
///   c02        c13
///   [2]  c23   [3]
void ReplexBiasValue::logGrid(const std::vector<std::string> &cells) const {
  plumed_assert(dim <= 2);
  const int ncols = shape[dim - 1];
  const int nrows = dim > 1 ? shape[0] : 1;
  const int hpair = (dim - 1) * nrep; // offset of the pairs within a row
  // A single width for all the slots, so that columns line up
  int w = std::to_string(nrep - 1).size() + 2;
  for (const auto &c : cells) {
    w = std::max(w, static_cast<int>(c.size()));
  }
  for (int r = 0; r < nrows; r++) {
    // Biases, with the pairs along the row in between
    log.printf(" ");
    for (int c = 0; c < ncols; c++) {
      const int ibias = r * ncols + c;
      log.printf(" %*s", w, ("[" + std::to_string(ibias) + "]").c_str());
      if (c < ncols - 1) {
        log.printf(" %*s", w, cells[ibias + hpair].c_str());
      }
    }
    log.printf("\n");
    // Pairs between this row and the next one, under each bias
    if (r < nrows - 1) {
      log.printf(" ");
      for (int c = 0; c < ncols; c++) {
        const int ibias = r * ncols + c;
        log.printf(" %*s", w, cells[ibias].c_str());
        if (c < ncols - 1) {
          log.printf(" %*s", w, "");
        }
      }
      log.printf("\n");
    }
  }
}

} // namespace replex_biasvalue
} // namespace PLMD
