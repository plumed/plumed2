/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#include "bias/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include <math.h>
#include <fstream>
#include <sstream>
#include <map>


using namespace std;

namespace PLMD {
namespace isdb {

//+PLUMEDOC ISDB_BIAS ROSETTA_COEVOLUTION
/*


*/
//+ENDPLUMEDOC


class RosettaCoEvolutionRestraint : public bias::Bias
{
  // "cutoff", gamma, and weights parameters;
  vector<double> R0_;
  vector<double> gamma_;
  vector<double> weights_;

  // parallel stuff
  unsigned rank_;
  unsigned nrep_;

  void setup_restraint(double kappa, bool doScaleProb, string res_file);

public:
  RosettaCoEvolutionRestraint(const ActionOptions&);
  void calculate();
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(RosettaCoEvolutionRestraint,"ROSETTA_COEVOLUTION")

void RosettaCoEvolutionRestraint::registerKeywords(Keywords& keys) {
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","RES_FILE","file with residue ids for each argument");
  keys.add("compulsory","KAPPA","overall weight");
  keys.addFlag("SCALE_PROB",false,"scale weights by probability");
  componentsAreNotOptional(keys);
}

RosettaCoEvolutionRestraint::RosettaCoEvolutionRestraint(const ActionOptions&ao):
  PLUMED_BIAS_INIT(ao)
{
  string res_file;
  parse("RES_FILE", res_file);

  // check if scaling by probability
  bool doScaleProb = false;
  parseFlag("SCALE_PROB", doScaleProb);

  // read overall weight
  double kappa;
  parse("KAPPA", kappa);

  checkRead();

  // setup restraint
  setup_restraint(kappa, doScaleProb, res_file);
  if(R0_.size()!=getNumberOfArguments()) error("Check residue pair file");

  log.printf("  overall weight %f\n",kappa);
  if(doScaleProb) log.printf("  scale weights by probability\n");

  // initialize parallel stuff
  rank_ = comm.Get_rank();
  nrep_ = comm.Get_size();

}

void  RosettaCoEvolutionRestraint::setup_restraint(double kappa, bool doScaleProb, string res_file)
{
  // reference distance and slope for each pairs of residues
  map< pair<string,string>, pair<double,double> > DREF_;
  // load from Baker's paper - in Angstrom
  DREF_[make_pair("G","G")]=make_pair(4.467,0.017); DREF_[make_pair("G","A")]=make_pair(5.201,0.269); DREF_[make_pair("G","S")]=make_pair(5.510,0.153);
  DREF_[make_pair("G","V")]=make_pair(5.671,0.107); DREF_[make_pair("G","C")]=make_pair(5.777,0.129); DREF_[make_pair("G","T")]=make_pair(5.619,0.120);
  DREF_[make_pair("G","P")]=make_pair(6.140,0.245); DREF_[make_pair("G","D")]=make_pair(6.135,0.193); DREF_[make_pair("G","N")]=make_pair(6.321,0.169);
  DREF_[make_pair("G","I")]=make_pair(6.413,0.179); DREF_[make_pair("G","L")]=make_pair(6.554,0.125); DREF_[make_pair("G","E")]=make_pair(7.036,0.249);
  DREF_[make_pair("G","Q")]=make_pair(7.297,0.216); DREF_[make_pair("G","M")]=make_pair(7.383,0.255); DREF_[make_pair("G","H")]=make_pair(7.472,0.206);
  DREF_[make_pair("G","K")]=make_pair(8.216,0.358); DREF_[make_pair("G","F")]=make_pair(7.966,0.219); DREF_[make_pair("G","Y")]=make_pair(9.098,0.267);
  DREF_[make_pair("G","R")]=make_pair(9.166,0.334); DREF_[make_pair("G","W")]=make_pair(8.966,0.239); DREF_[make_pair("A","A")]=make_pair(5.381,0.262);
  DREF_[make_pair("A","S")]=make_pair(5.829,0.291); DREF_[make_pair("A","V")]=make_pair(5.854,0.312); DREF_[make_pair("A","C")]=make_pair(6.057,0.394);
  DREF_[make_pair("A","T")]=make_pair(5.982,0.378); DREF_[make_pair("A","P")]=make_pair(6.412,0.399); DREF_[make_pair("A","D")]=make_pair(6.388,0.289);
  DREF_[make_pair("A","N")]=make_pair(6.766,0.349); DREF_[make_pair("A","I")]=make_pair(6.587,0.214); DREF_[make_pair("A","L")]=make_pair(6.707,0.250);
  DREF_[make_pair("A","E")]=make_pair(7.124,0.340); DREF_[make_pair("A","Q")]=make_pair(7.583,0.356); DREF_[make_pair("A","M")]=make_pair(7.605,0.394);
  DREF_[make_pair("A","H")]=make_pair(7.591,0.380); DREF_[make_pair("A","K")]=make_pair(8.327,0.550); DREF_[make_pair("A","F")]=make_pair(8.162,0.260);
  DREF_[make_pair("A","Y")]=make_pair(9.121,0.443); DREF_[make_pair("A","R")]=make_pair(9.365,0.485); DREF_[make_pair("A","W")]=make_pair(9.252,0.290);
  DREF_[make_pair("S","S")]=make_pair(6.190,0.292); DREF_[make_pair("S","V")]=make_pair(6.567,0.205); DREF_[make_pair("S","C")]=make_pair(6.590,0.240);
  DREF_[make_pair("S","T")]=make_pair(6.450,0.214); DREF_[make_pair("S","P")]=make_pair(6.937,0.321); DREF_[make_pair("S","D")]=make_pair(6.760,0.323);
  DREF_[make_pair("S","N")]=make_pair(7.081,0.305); DREF_[make_pair("S","I")]=make_pair(7.142,0.342); DREF_[make_pair("S","L")]=make_pair(7.394,0.287);
  DREF_[make_pair("S","E")]=make_pair(7.483,0.446); DREF_[make_pair("S","Q")]=make_pair(7.807,0.408); DREF_[make_pair("S","M")]=make_pair(8.010,0.369);
  DREF_[make_pair("S","H")]=make_pair(8.051,0.435); DREF_[make_pair("S","K")]=make_pair(8.792,0.445); DREF_[make_pair("S","F")]=make_pair(8.694,0.394);
  DREF_[make_pair("S","Y")]=make_pair(9.594,0.467); DREF_[make_pair("S","R")]=make_pair(9.753,0.483); DREF_[make_pair("S","W")]=make_pair(9.770,0.497);
  DREF_[make_pair("V","V")]=make_pair(6.759,0.145); DREF_[make_pair("V","C")]=make_pair(6.941,0.173); DREF_[make_pair("V","T")]=make_pair(6.791,0.138);
  DREF_[make_pair("V","P")]=make_pair(7.063,0.298); DREF_[make_pair("V","D")]=make_pair(6.972,0.287); DREF_[make_pair("V","N")]=make_pair(7.219,0.232);
  DREF_[make_pair("V","I")]=make_pair(7.441,0.242); DREF_[make_pair("V","L")]=make_pair(7.633,0.179); DREF_[make_pair("V","E")]=make_pair(7.404,0.510);
  DREF_[make_pair("V","Q")]=make_pair(8.008,0.359); DREF_[make_pair("V","M")]=make_pair(8.335,0.295); DREF_[make_pair("V","H")]=make_pair(8.179,0.383);
  DREF_[make_pair("V","K")]=make_pair(8.077,0.634); DREF_[make_pair("V","F")]=make_pair(9.057,0.246); DREF_[make_pair("V","Y")]=make_pair(9.442,0.535);
  DREF_[make_pair("V","R")]=make_pair(9.513,0.514); DREF_[make_pair("V","W")]=make_pair(10.021,0.271); DREF_[make_pair("C","C")]=make_pair(6.426,0.178);
  DREF_[make_pair("C","T")]=make_pair(6.801,0.181); DREF_[make_pair("C","P")]=make_pair(7.157,0.259); DREF_[make_pair("C","D")]=make_pair(6.985,0.299);
  DREF_[make_pair("C","N")]=make_pair(7.205,0.240); DREF_[make_pair("C","I")]=make_pair(7.476,0.295); DREF_[make_pair("C","L")]=make_pair(7.685,0.206);
  DREF_[make_pair("C","E")]=make_pair(7.449,0.538); DREF_[make_pair("C","Q")]=make_pair(7.962,0.347); DREF_[make_pair("C","M")]=make_pair(8.265,0.439);
  DREF_[make_pair("C","H")]=make_pair(8.422,0.203); DREF_[make_pair("C","K")]=make_pair(8.494,0.521); DREF_[make_pair("C","F")]=make_pair(9.026,0.286);
  DREF_[make_pair("C","Y")]=make_pair(9.362,0.585); DREF_[make_pair("C","R")]=make_pair(9.460,0.491); DREF_[make_pair("C","W")]=make_pair(9.752,0.417);
  DREF_[make_pair("T","T")]=make_pair(6.676,0.188); DREF_[make_pair("T","P")]=make_pair(7.062,0.320); DREF_[make_pair("T","D")]=make_pair(6.971,0.307);
  DREF_[make_pair("T","N")]=make_pair(7.159,0.262); DREF_[make_pair("T","I")]=make_pair(7.442,0.259); DREF_[make_pair("T","L")]=make_pair(7.642,0.190);
  DREF_[make_pair("T","E")]=make_pair(7.628,0.409); DREF_[make_pair("T","Q")]=make_pair(8.055,0.378); DREF_[make_pair("T","M")]=make_pair(8.397,0.292);
  DREF_[make_pair("T","H")]=make_pair(8.221,0.417); DREF_[make_pair("T","K")]=make_pair(8.715,0.464); DREF_[make_pair("T","F")]=make_pair(9.030,0.264);
  DREF_[make_pair("T","Y")]=make_pair(9.813,0.430); DREF_[make_pair("T","R")]=make_pair(9.764,0.477); DREF_[make_pair("T","W")]=make_pair(9.980,0.315);
  DREF_[make_pair("P","P")]=make_pair(7.288,0.339); DREF_[make_pair("P","D")]=make_pair(7.321,0.416); DREF_[make_pair("P","N")]=make_pair(7.497,0.334);
  DREF_[make_pair("P","I")]=make_pair(7.554,0.336); DREF_[make_pair("P","L")]=make_pair(7.751,0.317); DREF_[make_pair("P","E")]=make_pair(7.938,0.475);
  DREF_[make_pair("P","Q")]=make_pair(8.308,0.410); DREF_[make_pair("P","M")]=make_pair(8.247,0.388); DREF_[make_pair("P","H")]=make_pair(8.537,0.457);
  DREF_[make_pair("P","K")]=make_pair(9.198,0.550); DREF_[make_pair("P","F")]=make_pair(8.895,0.425); DREF_[make_pair("P","Y")]=make_pair(9.965,0.506);
  DREF_[make_pair("P","R")]=make_pair(10.266,0.506); DREF_[make_pair("P","W")]=make_pair(9.719,0.462); DREF_[make_pair("D","D")]=make_pair(8.001,0.392);
  DREF_[make_pair("D","N")]=make_pair(7.672,0.337); DREF_[make_pair("D","I")]=make_pair(7.472,0.341); DREF_[make_pair("D","L")]=make_pair(7.696,0.348);
  DREF_[make_pair("D","E")]=make_pair(8.945,0.354); DREF_[make_pair("D","Q")]=make_pair(8.601,0.357); DREF_[make_pair("D","M")]=make_pair(8.401,0.361);
  DREF_[make_pair("D","H")]=make_pair(8.634,0.325); DREF_[make_pair("D","K")]=make_pair(9.306,0.343); DREF_[make_pair("D","F")]=make_pair(9.111,0.351);
  DREF_[make_pair("D","Y")]=make_pair(9.979,0.676); DREF_[make_pair("D","R")]=make_pair(10.123,0.327); DREF_[make_pair("D","W")]=make_pair(9.867,0.475);
  DREF_[make_pair("N","N")]=make_pair(7.682,0.249); DREF_[make_pair("N","I")]=make_pair(7.631,0.341); DREF_[make_pair("N","L")]=make_pair(7.889,0.279);
  DREF_[make_pair("N","E")]=make_pair(8.485,0.423); DREF_[make_pair("N","Q")]=make_pair(8.502,0.373); DREF_[make_pair("N","M")]=make_pair(8.550,0.310);
  DREF_[make_pair("N","H")]=make_pair(8.672,0.289); DREF_[make_pair("N","K")]=make_pair(9.319,0.398); DREF_[make_pair("N","F")]=make_pair(9.168,0.393);
  DREF_[make_pair("N","Y")]=make_pair(10.039,0.586); DREF_[make_pair("N","R")]=make_pair(10.135,0.372); DREF_[make_pair("N","W")]=make_pair(9.976,0.458);
  DREF_[make_pair("I","I")]=make_pair(8.096,0.321); DREF_[make_pair("I","L")]=make_pair(8.342,0.261); DREF_[make_pair("I","E")]=make_pair(7.949,0.453);
  DREF_[make_pair("I","Q")]=make_pair(8.302,0.406); DREF_[make_pair("I","M")]=make_pair(8.874,0.327); DREF_[make_pair("I","H")]=make_pair(8.523,0.379);
  DREF_[make_pair("I","K")]=make_pair(8.329,0.582); DREF_[make_pair("I","F")]=make_pair(9.602,0.347); DREF_[make_pair("I","Y")]=make_pair(9.719,0.589);
  DREF_[make_pair("I","R")]=make_pair(9.746,0.557); DREF_[make_pair("I","W")]=make_pair(10.470,0.397); DREF_[make_pair("L","L")]=make_pair(8.522,0.198);
  DREF_[make_pair("L","E")]=make_pair(8.077,0.475); DREF_[make_pair("L","Q")]=make_pair(8.480,0.411); DREF_[make_pair("L","M")]=make_pair(9.122,0.318);
  DREF_[make_pair("L","H")]=make_pair(8.676,0.401); DREF_[make_pair("L","K")]=make_pair(8.479,0.591); DREF_[make_pair("L","F")]=make_pair(9.900,0.260);
  DREF_[make_pair("L","Y")]=make_pair(9.889,0.611); DREF_[make_pair("L","R")]=make_pair(9.852,0.578); DREF_[make_pair("L","W")]=make_pair(10.707,0.331);
  DREF_[make_pair("E","E")]=make_pair(9.863,0.389); DREF_[make_pair("E","Q")]=make_pair(9.328,0.450); DREF_[make_pair("E","M")]=make_pair(8.870,0.511);
  DREF_[make_pair("E","H")]=make_pair(9.454,0.443); DREF_[make_pair("E","K")]=make_pair(9.842,0.434); DREF_[make_pair("E","F")]=make_pair(9.403,0.512);
  DREF_[make_pair("E","Y")]=make_pair(10.544,0.469); DREF_[make_pair("E","R")]=make_pair(10.713,0.363); DREF_[make_pair("E","W")]=make_pair(10.303,0.493);
  DREF_[make_pair("Q","Q")]=make_pair(9.074,0.436); DREF_[make_pair("Q","M")]=make_pair(9.102,0.498); DREF_[make_pair("Q","H")]=make_pair(9.391,0.401);
  DREF_[make_pair("Q","K")]=make_pair(9.667,0.521); DREF_[make_pair("Q","F")]=make_pair(9.506,0.451); DREF_[make_pair("Q","Y")]=make_pair(10.534,0.547);
  DREF_[make_pair("Q","R")]=make_pair(10.610,0.535); DREF_[make_pair("Q","W")]=make_pair(10.429,0.490); DREF_[make_pair("M","M")]=make_pair(9.530,0.457);
  DREF_[make_pair("M","H")]=make_pair(9.396,0.342); DREF_[make_pair("M","K")]=make_pair(9.096,0.611); DREF_[make_pair("M","F")]=make_pair(10.253,0.377);
  DREF_[make_pair("M","Y")]=make_pair(10.400,0.661); DREF_[make_pair("M","R")]=make_pair(10.250,0.641); DREF_[make_pair("M","W")]=make_pair(11.110,0.397);
  DREF_[make_pair("H","H")]=make_pair(10.606,0.333); DREF_[make_pair("H","K")]=make_pair(9.582,0.714); DREF_[make_pair("H","F")]=make_pair(9.602,0.542);
  DREF_[make_pair("H","Y")]=make_pair(10.843,0.554); DREF_[make_pair("H","R")]=make_pair(10.879,0.595); DREF_[make_pair("H","W")]=make_pair(10.661,0.458);
  DREF_[make_pair("K","K")]=make_pair(10.662,0.738); DREF_[make_pair("K","F")]=make_pair(9.344,0.441); DREF_[make_pair("K","Y")]=make_pair(10.627,0.704);
  DREF_[make_pair("K","R")]=make_pair(11.322,0.648); DREF_[make_pair("K","W")]=make_pair(10.136,0.470); DREF_[make_pair("F","F")]=make_pair(10.903,0.460);
  DREF_[make_pair("F","Y")]=make_pair(10.999,0.767); DREF_[make_pair("F","R")]=make_pair(10.577,0.738); DREF_[make_pair("F","W")]=make_pair(11.758,0.447);
  DREF_[make_pair("Y","Y")]=make_pair(11.536,0.855); DREF_[make_pair("Y","R")]=make_pair(11.615,0.822); DREF_[make_pair("Y","W")]=make_pair(11.807,0.684);
  DREF_[make_pair("R","R")]=make_pair(12.050,0.704); DREF_[make_pair("R","W")]=make_pair(11.355,0.889); DREF_[make_pair("W","W")]=make_pair(12.806,0.473);

// open residue pair file
  ifstream rfile;
  rfile.open(res_file);
  string line, res0, res1;
  string score, prob;
// iterator for residue map
  map< pair<string,string>, pair<double,double> >::iterator it;
// read file
  if (rfile.is_open()) {
    // read line by line
    while ( getline (rfile, line) )
    {
      // split line into strings separated by a space
      stringstream ss(line);
      // read residue pair, score and prob
      ss >> res0; ss >> res1; ss >> score; ss >> prob;
      // find entry in DREF_
      pair<string,string> key1 = make_pair(res0, res1);
      pair<string,string> key2 = make_pair(res1, res0);
      // look for key1 in DREF_
      it = DREF_.find(key1);
      // and add to reference distances and slopes - convert to nm
      if (it != DREF_.end()) {
        R0_.push_back( 0.1 * DREF_[key1].first );
        gamma_.push_back( 0.1 * DREF_[key1].second );
      } else {
        R0_.push_back( 0.1 * DREF_[key2].first );
        gamma_.push_back( 0.1 * DREF_[key2].second );
      }
      // calculate weight
      double w = kappa * atof(score.c_str());
      // scale by probability, in case
      if(doScaleProb) w *= atof(prob.c_str());
      weights_.push_back(w);
    }
    rfile.close();
  }
  else error("Unable to open residue file");
}


void RosettaCoEvolutionRestraint::calculate()
{
  // allocate force vector
  vector<double> force(getNumberOfArguments(), 0.0);

  // calculate energy
  double ene = 0.0;
  // cycle on arguments
  for(unsigned i=rank_; i<getNumberOfArguments(); i=i+nrep_) {
    // get distance
    double dist = getArgument(i);
    // calculate temporary factor
    double tmp0 = exp( - ( dist - R0_[i] ) / gamma_[i] );
    // increment energy
    ene += weights_[i] / ( 1.0 + tmp0 );
    // calculate force
    force[i] = -weights_[i] / ( 1.0 + tmp0 ) / ( 1.0 + tmp0 ) * tmp0 / gamma_[i];
  }

  // sum energy and forces
  comm.Sum(&force[0], force.size());
  comm.Sum(&ene, 1);

  // apply forces
  for(unsigned i=0; i<force.size(); ++i) setOutputForce(i, force[i]);

  // set value of the bias
  setBias(ene);

}


}
}


