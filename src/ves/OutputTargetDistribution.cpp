/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2018 The VES code team
   (see the PEOPLE-VES file at the root of this folder for a list of names)

   See http://www.ves-code.org for more information.

   This file is part of VES code module.

   The VES code module is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   The VES code module is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with the VES code module.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#include "TargetDistribution.h"

#include "GridIntegrationWeights.h"
#include "VesTools.h"

#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Value.h"
#include "tools/File.h"
#include "tools/Grid.h"


namespace PLMD {
namespace ves {

//+PLUMEDOC VES_UTILS VES_OUTPUT_TARGET_DISTRIBUTION
/*
Output target distribution to file.

This action can be used to output target distributions to a grid file,
for example to see how they look like before using them in a VES bias.
This action only support static target distributions.

This action is normally used through the \ref driver.


\par Examples

In the following input we define a target distribution that is uniform for
argument 1 and a Gaussian for argument 2 and then output it to a file
called targetdist-1.data.
\plumedfile
t1_1: TD_UNIFORM  MINIMA=-4.0  MAXIMA=+4.0
t1_2: TD_GAUSSIAN  CENTER1=-2.0  SIGMA1=0.5
t1: TD_PRODUCT_DISTRIBUTION  DISTRIBUTIONS=t1_1,t1_2

VES_OUTPUT_TARGET_DISTRIBUTION ...
 GRID_MIN=-4.0,-4.0
 GRID_MAX=+4.0,+4.0
 GRID_BINS=100,100
 TARGET_DISTRIBUTION=t1
 TARGETDIST_FILE=targetdist-1.data
 LOG_TARGETDIST_FILE=targetdist-1.log.data
 FMT_GRIDS=%11.6f
... VES_OUTPUT_TARGET_DISTRIBUTION
\endplumedfile

This input should be run through the driver by using a command similar to the
following one where the trajectory/configuration file configuration.gro is needed to
trick the code to exit correctly.
\verbatim
plumed driver --plumed plumed.dat --igro configuration.gro
\endverbatim

*/
//+ENDPLUMEDOC


class OutputTargetDistribution :
  public Action
{
public:
  explicit OutputTargetDistribution(const ActionOptions&);
  void calculate() override {}
  void apply() override {}
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(OutputTargetDistribution,"VES_OUTPUT_TARGET_DISTRIBUTION")

void OutputTargetDistribution::registerKeywords(Keywords& keys) {
  Action::registerKeywords(keys);
  keys.add("compulsory","GRID_MIN","the lower bounds for the grid");
  keys.add("compulsory","GRID_MAX","the upper bounds for the grid");
  keys.add("compulsory","GRID_BINS","the number of bins used for the grid.");
  keys.add("optional","GRID_PERIODICITY","specify if the individual arguments should be made periodic (YES) or not (NO). By default all arguments are taken as not periodic.");
  keys.add("compulsory","TARGETDIST_FILE","filename of the file for writing the target distribution");
  keys.add("optional","LOG_TARGETDIST_FILE","filename of the file for writing the log of the target distribution");
  keys.add("compulsory","TARGET_DISTRIBUTION","the target distribution to be used.");
  keys.add("optional","FMT_GRIDS","the numerical format of the target distribution grids written to file. By default it is %14.9f");
  keys.addFlag("DO_1D_PROJECTIONS",false,"Also output the one-dimensional marginal distributions for multi-dimensional target distribution.");
}

OutputTargetDistribution::OutputTargetDistribution(const ActionOptions&ao):
  Action(ao)
{

  std::string targetdist_fname;
  parse("TARGETDIST_FILE",targetdist_fname);
  std::string log_targetdist_fname;
  parse("LOG_TARGETDIST_FILE",log_targetdist_fname);
  if(targetdist_fname==log_targetdist_fname) {
    plumed_merror("error in " + getName() + ":TARGETDIST_FILE and LOG_TARGETDIST_FILE cannot be the same");
  }

  std::vector<unsigned int> grid_bins;
  parseVector("GRID_BINS",grid_bins);
  unsigned int nargs = grid_bins.size();

  std::vector<std::string> grid_min(nargs);
  parseVector("GRID_MIN",grid_min);
  std::vector<std::string> grid_max(nargs);
  parseVector("GRID_MAX",grid_max);

  std::vector<std::string> grid_periodicity(nargs);
  parseVector("GRID_PERIODICITY",grid_periodicity);
  if(grid_periodicity.size()==0) {grid_periodicity.assign(nargs,"NO");}

  std::string fmt_grids="%14.9f";
  parse("FMT_GRIDS",fmt_grids);

  bool do_1d_proj = false;
  parseFlag("DO_1D_PROJECTIONS",do_1d_proj);
  if(do_1d_proj && nargs==1) {
    plumed_merror("doesn't make sense to use the DO_1D_PROJECTIONS keyword for a one-dimensional distribution");
  }

  plumed_massert(grid_min.size()==nargs,"mismatch between number of values given for grid parameters");
  plumed_massert(grid_max.size()==nargs,"mismatch between number of values given for grid parameters");
  plumed_massert(grid_periodicity.size()==nargs,"mismatch between number of values given for grid parameters");

  std::string targetdist_label;
  parse("TARGET_DISTRIBUTION",targetdist_label);
  checkRead();
  //
  std::vector<Value*> arguments(nargs);
  for(unsigned int i=0; i < nargs; i++) {
    std::string is; Tools::convert(i+1,is);
    if(nargs==1) {is="";}
    arguments[i]= new Value(NULL,"arg"+is,false);
    if(grid_periodicity[i]=="YES") {
      arguments[i]->setDomain(grid_min[i],grid_max[i]);
    }
    else if(grid_periodicity[i]=="NO") {
      arguments[i]->setNotPeriodic();
    }
    else {
      plumed_merror("wrong value given in GRID_PERIODICITY, either specify YES or NO");
    }
  }

  std::string error_msg = "";
  TargetDistribution* targetdist_pntr = VesTools::getPointerFromLabel<TargetDistribution*>(targetdist_label,plumed.getActionSet(),error_msg);
  if(error_msg.size()>0) {plumed_merror("Error in keyword TARGET_DISTRIBUTION of "+getName()+": "+error_msg);}
  //
  if(targetdist_pntr->isDynamic()) {
    plumed_merror(getName() + " only works for static target distributions");
  }
  targetdist_pntr->setupGrids(arguments,grid_min,grid_max,grid_bins);
  targetdist_pntr->updateTargetDist();
  Grid* targetdist_grid_pntr = targetdist_pntr->getTargetDistGridPntr();
  Grid* log_targetdist_grid_pntr = targetdist_pntr->getLogTargetDistGridPntr();


  double sum_grid = TargetDistribution::integrateGrid(targetdist_grid_pntr);
  log.printf("  target distribution integrated over the grid: %16.12f\n",sum_grid);
  log.printf("                                                (%30.16e)\n",sum_grid);
  //
  OFile ofile;
  ofile.link(*this);
  ofile.enforceBackup();
  ofile.open(targetdist_fname);
  targetdist_grid_pntr->setOutputFmt(fmt_grids);
  targetdist_grid_pntr->writeToFile(ofile);
  ofile.close();
  if(log_targetdist_fname.size()>0) {
    OFile ofile2;
    ofile2.link(*this);
    ofile2.enforceBackup();
    ofile2.open(log_targetdist_fname);
    log_targetdist_grid_pntr->setOutputFmt(fmt_grids);
    log_targetdist_grid_pntr->writeToFile(ofile2);
    ofile2.close();
  }

  if(do_1d_proj) {
    for(unsigned int i=0; i<nargs; i++) {
      std::vector<std::string> arg1d(1);
      arg1d[0] = arguments[i]->getName();
      Grid marginal_grid = targetdist_pntr->getMarginal(arg1d);
      //
      std::string suffix;
      Tools::convert(i+1,suffix);
      suffix = "proj-" + suffix;
      std::string marginal_fname = FileBase::appendSuffix(targetdist_fname,"."+suffix);
      //
      OFile ofile3;
      ofile3.link(*this);
      ofile3.enforceBackup();
      ofile3.open(marginal_fname);
      marginal_grid.setOutputFmt(fmt_grids);
      marginal_grid.writeToFile(ofile3);
    }
  }

  //
  for(unsigned int i=0; i < nargs; i++) {
    delete arguments[i];
  }
  arguments.clear();


}





}
}
