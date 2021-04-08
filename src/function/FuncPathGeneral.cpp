/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2020 The plumed team
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

#include "Function.h"
#include "ActionRegister.h"
#include "tools/IFile.h"

namespace PLMD {
namespace function {

//+PLUMEDOC FUNCTION FUNCPATHGENERAL
/*
This function calculates path collective variables (PCVs) using an arbitrary combination of collective variables (see \cite Hovan2019).

This variable computes the progress along a given set of frames that is provided in an input file ("s" component) and the distance from them ("z" component).
The input file could be a colvar file generated with plumed driver on a trajectory containing the frames.

The metric for the path collective variables takes the following form:

\f[
R[X - X_i] = \sum_{j=1}^M c_j^2 (x_j - x_{i,j})^2\,.
\f]

Here, the coefficients \f$c_j\f$ determine the relative weights of the collective variables \f$c_j\f$ in the metric.
A value for the lambda coefficient also needs to be provided, typically chosen in such a way that it ensures a smooth variation of the "s" component.

\par Examples

This command calculates the PCVs using the values from the file COLVAR_TRAJ and the provided values for the lambda and the coefficients.
Since the columns in the file were not specified, the first one will be ignored (assumed to correspond to the time) and the rest used.

\plumedfile
FUNCPATHGENERAL ...
LABEL=path
LAMBDA=12.2
REFERENCE=COLVAR_TRAJ
COEFFICIENTS=0.3536,0.3536,0.3536,0.3536,0.7071
ARG=d1,d2,d,t,drmsd
... FUNCPATHGENERAL
\endplumedfile

The command below is a variation of the previous one, specifying a subset of the collective variables and using a neighbor list.
The columns are zero-indexed.
The neighbor list will include the 10 closest frames and will be recalculated every 20 steps.

\plumedfile
FUNCPATHGENERAL ...
LABEL=path
LAMBDA=5.0
REFERENCE=COLVAR_TRAJ
COLUMNS=2,3,4
COEFFICIENTS=0.3536,0.3536,0.3536
ARG=d2,d,t
NEIGH_SIZE=10
NEIGH_STRIDE=20
... FUNCPATHGENERAL
\endplumedfile

*/
//+ENDPLUMEDOC

class FuncPathGeneral : public Function {
  double lambda;
  int neigh_size;
  double neigh_stride;

  std::vector<double> coefficients;
  std::vector< std::vector<double> > path_cv_values;

  // For faster calculation
  std::vector<double> expdists;

  // For calculating derivatives
  std::vector< std::vector<double> > numerators;
  std::vector<double> s_path_ders;
  std::vector<double> z_path_ders;

  // For handling periodicity
  std::vector<double> domains;

  std::string reference;
  std::vector<int> columns;

  std::vector< std::pair<int,double> > neighpair;
  std::vector <Value*> allArguments;

  // Methods
  void loadReference();

  struct pairordering {
    bool operator ()(std::pair<int, double> const& a, std::pair<int, double> const& b) {
      return (a).second < (b).second;
    }
  };

public:
  explicit FuncPathGeneral(const ActionOptions&);
// Active methods:
  virtual void calculate();
  virtual void prepare();
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(FuncPathGeneral, "FUNCPATHGENERAL")

void FuncPathGeneral::loadReference() {
  IFile input;
  input.open(reference);
  if (!input)
    plumed_merror("Could not open the reference file!");
  while (input)
  {
    std::vector<std::string> strings;
    Tools::getParsedLine(input, strings);
    if (strings.empty())
      continue;
    std::vector<double> colvarLine;
    double value;
    int max = columns.empty() ? strings.size() : columns.size();
    for (int i = 0; i < max; ++i)
    {
      int col = columns.empty() ? i : columns[i];
      // If no columns have been entered, ignore the first (time) and take the rest
      if (columns.empty() && i == 0)
        continue;

      Tools::convert(strings[col], value);
      colvarLine.push_back(value);
    }
    path_cv_values.push_back(colvarLine);
  }
}

void FuncPathGeneral::registerKeywords(Keywords& keys) {
  Function::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory", "LAMBDA", "Lambda parameter required for smoothing");
  keys.add("compulsory", "COEFFICIENTS", "Coefficients to be assigned to the CVs");
  keys.add("compulsory", "REFERENCE", "Colvar file needed to provide the CV milestones");
  keys.add("optional", "COLUMNS", "List of columns in the reference colvar file specifying the CVs");
  keys.add("optional", "NEIGH_SIZE", "Size of the neighbor list");
  keys.add("optional", "NEIGH_STRIDE", "How often the neighbor list needs to be calculated in time units");
  componentsAreNotOptional(keys);
  keys.addOutputComponent("s", "default", "Position on the path");
  keys.addOutputComponent("z", "default", "Distance from the path");
}

FuncPathGeneral::FuncPathGeneral(const ActionOptions&ao):
  Action(ao),
  Function(ao),
  neigh_size(-1),
  neigh_stride(-1.)
{
  parse("LAMBDA", lambda);
  parse("NEIGH_SIZE", neigh_size);
  parse("NEIGH_STRIDE", neigh_stride);
  parse("REFERENCE", reference);
  parseVector("COEFFICIENTS", coefficients);
  parseVector("COLUMNS", columns);
  checkRead();
  log.printf("  lambda is %f\n", lambda);
  if (getNumberOfArguments() != coefficients.size())
    plumed_merror("The numbers of coefficients and CVs are different!");
  if (!columns.empty()) {
    if (columns.size() != coefficients.size())
      plumed_merror("The numbers of coefficients and columns are different!");
  }
  log.printf("  Consistency check completed! Your path cvs look good!\n");

  // Load the reference colvar file
  loadReference();

  // Do some neighbour printout
  if (neigh_stride > 0. || neigh_size > 0) {
    if (static_cast<unsigned>(neigh_size) > path_cv_values.size()) {
      log.printf(" List size required ( %d ) is too large: resizing to the maximum number of arg required: %d  \n", neigh_size, getNumberOfArguments());
      neigh_size = path_cv_values.size();
    }
    log.printf("  Neighbour list enabled: \n");
    log.printf("                 size   :  %d elements\n", neigh_size);
    log.printf("                 stride :  %f time \n", neigh_stride);
  } else {
    log.printf("  Neighbour list NOT enabled \n");
  }

  addComponentWithDerivatives("s"); componentIsNotPeriodic("s");
  addComponentWithDerivatives("z"); componentIsNotPeriodic("z");

  // Initialise vectors
  std::vector<double> temp (coefficients.size());
  for (unsigned i = 0; i < path_cv_values.size(); ++i) {
    numerators.push_back(temp);
    expdists.push_back(0.);
    s_path_ders.push_back(0.);
    z_path_ders.push_back(0.);
  }

  // Store the arguments
  for (unsigned i=0; i<getNumberOfArguments(); i++)
    allArguments.push_back(getPntrToArgument(i));

  // Get periodic domains, negative for not periodic, stores half the domain length (maximum difference)
  for (unsigned i = 0; i < allArguments.size(); ++i) {
    if (allArguments[i]->isPeriodic()) {
      double min_lim, max_lim;
      allArguments[i]->getDomain(min_lim, max_lim);
      domains.push_back((max_lim - min_lim) / 2);
    }
    else
      domains.push_back(-1.);
  }
}

// Calculator
void FuncPathGeneral::calculate() {
  double s_path = 0.;
  double partition = 0.;
  double tmp, value, diff, expdist, s_der, z_der;
  int ii;

  typedef std::vector< std::pair< int,double> >::iterator pairiter;

  for (pairiter it = neighpair.begin(); it != neighpair.end(); ++it) {
    (*it).second = 0.;
  }

  if (neighpair.empty()) {
    // Resize at the first step
    neighpair.resize(path_cv_values.size());
    for (unsigned i = 0; i < path_cv_values.size(); ++i)
      neighpair[i].first = i;
  }

  Value* val_s_path=getPntrToComponent("s");
  Value* val_z_path=getPntrToComponent("z");

  for(unsigned j = 0; j < allArguments.size(); ++j) {
    value = allArguments[j]->get();
    for (pairiter it = neighpair.begin(); it != neighpair.end(); ++it) {
      diff = (value - path_cv_values[(*it).first][j]);
      if (domains[j] > 0) {
        if (diff > domains[j])
          diff -= 2 * domains[j];
        if (diff < -domains[j])
          diff += 2 * domains[j];
      }
      (*it).second += Tools::fastpow(coefficients[j] * diff, 2);
      numerators[(*it).first][j] = 2 * Tools::fastpow(coefficients[j], 2) * diff;
    }
  }

  for (pairiter it = neighpair.begin(); it != neighpair.end(); ++it) {
    expdist = std::exp(-lambda * (*it).second);
    expdists[(*it).first] = expdist;
    s_path += ((*it).first + 1) * expdist;
    partition += expdist;
  }

  s_path /= partition;
  val_s_path->set(s_path);
  val_z_path->set(-(1. / lambda) * std::log(partition));

  // Derivatives
  for (pairiter it = neighpair.begin(); it != neighpair.end(); ++it) {
    ii = (*it).first;
    tmp = lambda * expdists[ii] * (s_path - (ii + 1)) / partition;
    s_path_ders[ii] = tmp;
    z_path_ders[ii] = expdists[ii] / partition;
  }
  for (unsigned i = 0; i < coefficients.size(); ++i) {
    s_der = 0.;
    z_der = 0.;
    for (pairiter it = neighpair.begin(); it != neighpair.end(); ++it) {
      ii = (*it).first;
      s_der += s_path_ders[ii] * numerators[ii][i];
      z_der += z_path_ders[ii] * numerators[ii][i];
    }
    setDerivative(val_s_path, i, s_der);
    setDerivative(val_z_path, i, z_der);
  }
}

// Prepare the required arguments
void FuncPathGeneral::prepare() {
  // Neighbour list: rank and activate the chain for the next step

  // Neighbour list: if neigh_size < 0 never sort and keep the full vector
  // Neighbour list: if neigh_size > 0
  //                 if the size is full -> sort the vector and decide the dependencies for next step
  //                 if the size is not full -> check if next step will need the full dependency otherwise keep these dependencies

  if (neigh_size > 0) {
    if (neighpair.size() == path_cv_values.size()) {
      // The complete round has been done: need to sort, shorten and give it a go
      // Sort the values
      std::sort(neighpair.begin(), neighpair.end(), pairordering());
      // Resize the effective list
      neighpair.resize(neigh_size);
      log.printf("  NEIGHBOUR LIST NOW INCLUDES INDICES: ");
      for (int i = 0; i < neigh_size; ++i)
        log.printf(" %i ",neighpair[i].first);
      log.printf(" \n");
    } else {
      if (int(getStep()) % int(neigh_stride / getTimeStep()) == 0) {
        log.printf(" Time %f : recalculating full neighbour list \n", getStep() * getTimeStep());
        neighpair.resize(path_cv_values.size());
        for (unsigned i = 0; i < path_cv_values.size(); ++i)
          neighpair[i].first = i;
      }
    }
  }

  requestArguments(allArguments);
}

}
}
