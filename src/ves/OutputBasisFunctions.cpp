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

#include "BasisFunctions.h"
#include "TargetDistribution.h"

#include "CoeffsVector.h"
#include "VesTools.h"

#include "core/ActionRegister.h"
#include "core/ActionSet.h"
#include "core/PlumedMain.h"
#include "core/Value.h"
#include "tools/File.h"
#include "tools/Grid.h"


namespace PLMD {
namespace ves {

//+PLUMEDOC VES_UTILS VES_OUTPUT_BASISFUNCTIONS
/*
Output basis functions to file.

This action can be used to write out to a grid file the values and derivatives of
given basis functions. This is normally used for debugging when programming new
types of basis functions. For example, it is possible to calculate the
derivatives numerically and compare to the analytically calculated
derivatives.

This action is normally used through the \ref driver.

\par Examples

In the following input we define a Legendre polynomials basis functions
of order 14 over the interval -4.0 to 4.0 and output their values
and derivatives to files called bfL.values.data and bfL.derivs.data.
\plumedfile
BF_LEGENDRE ...
 ORDER=14
 MINIMUM=-4.0
 MAXIMUM=4.0
 LABEL=bfL
... BF_LEGENDRE

VES_OUTPUT_BASISFUNCTIONS ...
 BASIS_FUNCTIONS=bfL
 GRID_BINS=200
 FORMAT_VALUES_DERIVS=%13.6f
... VES_OUTPUT_BASISFUNCTIONS
\endplumedfile

This input should be run through the driver by using a command similar to the
following one where the trajectory/configuration file configuration.gro is needed to
trick the code to exit correctly.
\verbatim
plumed driver --plumed plumed.dat --igro configuration.gro
\endverbatim

*/
//+ENDPLUMEDOC


class OutputBasisFunctions :
  public Action
{
  std::vector<BasisFunctions*> bf_pntrs;
public:
  explicit OutputBasisFunctions(const ActionOptions&);
  TargetDistribution* setupTargetDistPntr(std::string keyword) const;
  void calculate() override {}
  void apply() override {}
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(OutputBasisFunctions,"VES_OUTPUT_BASISFUNCTIONS")

void OutputBasisFunctions::registerKeywords(Keywords& keys) {
  Action::registerKeywords(keys);
  keys.add("compulsory","BASIS_FUNCTIONS","the label of the basis functions that you want to use");
  keys.add("optional","GRID_BINS","the number of bins used for the grid for writing the basis function values and derivatives. The default value is 1000.");
  keys.add("optional","GRID_MIN","the minimum of the grid for writing the basis function values and derivatives. By default it is the minimum of the interval on which the basis functions are defined.");
  keys.add("optional","GRID_MAX","the maximum of the grid for writing the basis function values and derivatives. By default it is the maximum of the interval on which the basis functions are defined.");
  keys.add("optional","FILE_VALUES","filename of the file on which the basis function values are written. By default it is BF_LABEL.values.data.");
  keys.add("optional","FILE_DERIVS","filename of the file on which the basis function derivatives are written. By default it is BF_LABEL.derivs.data.");
  keys.add("optional","FORMAT_VALUES_DERIVS","the numerical format of the basis function values and derivatives written to file. By default it is %15.8f.\n");
  keys.add("optional","FILE_TARGETDIST_AVERAGES","filename of the file on which the averages over the target distributions are written. By default it is BF_LABEL.targetdist-averages.data.");
  keys.add("optional","FORMAT_TARGETDIST_AVERAGES","the numerical format of the target distribution averages written to file. By default it is %15.8f.\n");
  keys.add("optional","FILE_TARGETDIST","filename of the files on which the target distributions are written. By default it is BF_LABEL.targetdist-#.data.");
  keys.add("numbered","TARGET_DISTRIBUTION","the target distribution to be used.");
  keys.addFlag("IGNORE_PERIODICITY",false,"if the periodicity of the basis functions should be ignored.");
  keys.addFlag("NUMERICAL_DERIVATIVES",false,"if the derivatives of the basis functions should be calculated numerically.");
}

OutputBasisFunctions::OutputBasisFunctions(const ActionOptions&ao):
  Action(ao),
  bf_pntrs(0)
{
  std::vector<std::string> basisset_labels(0);
  parseVector("BASIS_FUNCTIONS",basisset_labels);
  if(basisset_labels.size()>1) {plumed_merror("Only one basis set label allowed in keyword BASIS_FUNCTIONS of "+getName());}

  std::string error_msg = "";
  bf_pntrs = VesTools::getPointersFromLabels<BasisFunctions*>(basisset_labels,plumed.getActionSet(),error_msg);
  if(error_msg.size()>0) {plumed_merror("Error in keyword BASIS_FUNCTIONS of "+getName()+": "+error_msg);}

  unsigned int nbins = 1000;
  parse("GRID_BINS",nbins);

  std::string min_str = bf_pntrs[0]->intervalMinStr();
  std::string max_str = bf_pntrs[0]->intervalMaxStr();
  parse("GRID_MIN",min_str);
  parse("GRID_MAX",max_str);

  std::string fname_values = bf_pntrs[0]->getLabel()+".values.data";
  parse("FILE_VALUES",fname_values);
  std::string fname_derives = bf_pntrs[0]->getLabel()+".derivs.data";
  parse("FILE_DERIVS",fname_derives);
  std::string fname_targetdist_aver = bf_pntrs[0]->getLabel()+".targetdist-averages.data";
  parse("FILE_TARGETDIST_AVERAGES",fname_targetdist_aver);
  std::string fname_targetdist = bf_pntrs[0]->getLabel()+".targetdist-.data";
  parse("FILE_TARGETDIST",fname_targetdist);

  std::string fmt_values_derivs = "%15.8f";
  parse("FORMAT_VALUES_DERIVS",fmt_values_derivs);
  std::string fmt_targetdist_aver = "%15.8f";
  parse("FORMAT_TARGETDIST_AVERAGES",fmt_targetdist_aver);

  bool ignore_periodicity = false;
  parseFlag("IGNORE_PERIODICITY",ignore_periodicity);

  bool numerical_deriv = false;
  parseFlag("NUMERICAL_DERIVATIVES",numerical_deriv);

  std::vector<TargetDistribution*> targetdist_pntrs;
  targetdist_pntrs.push_back(NULL);
  std::string targetdist_label="";
  for(int i=1;; i++) {
    if(!parseNumbered("TARGET_DISTRIBUTION",i,targetdist_label)) {break;}
    std::string error_msg = "";
    TargetDistribution* pntr_tmp = VesTools::getPointerFromLabel<TargetDistribution*>(targetdist_label,plumed.getActionSet(),error_msg);
    if(error_msg.size()>0) {plumed_merror("Error in keyword TARGET_DISTRIBUTION of "+getName()+": "+error_msg);}
    targetdist_pntrs.push_back(pntr_tmp);
  }
  checkRead();
  //
  OFile ofile_values;
  ofile_values.link(*this);
  ofile_values.enforceBackup();
  ofile_values.open(fname_values);
  OFile ofile_derivs;
  ofile_derivs.link(*this);
  ofile_derivs.enforceBackup();
  ofile_derivs.open(fname_derives);
  bf_pntrs[0]->writeBasisFunctionsToFile(ofile_values,ofile_derivs,min_str,max_str,nbins,ignore_periodicity,fmt_values_derivs,numerical_deriv);
  ofile_values.close();
  ofile_derivs.close();
  //
  std::vector<std::string> grid_min(1); grid_min[0]=bf_pntrs[0]->intervalMinStr();
  std::vector<std::string> grid_max(1); grid_max[0]=bf_pntrs[0]->intervalMaxStr();
  std::vector<unsigned int> grid_bins(1); grid_bins[0]=nbins;
  std::vector<Value*> arguments(1);
  arguments[0]= new Value(NULL,"arg",false);
  if(bf_pntrs[0]->arePeriodic() && !ignore_periodicity) {
    arguments[0]->setDomain(bf_pntrs[0]->intervalMinStr(),bf_pntrs[0]->intervalMaxStr());
  }
  else {
    arguments[0]->setNotPeriodic();
  }

  OFile ofile_targetdist_aver;
  ofile_targetdist_aver.link(*this);
  ofile_targetdist_aver.enforceBackup();
  ofile_targetdist_aver.open(fname_targetdist_aver);

  for(unsigned int i=0; i<targetdist_pntrs.size(); i++) {
    std::string is; Tools::convert(i,is);
    //
    if(targetdist_pntrs[i]!=NULL) {
      targetdist_pntrs[i]->setupGrids(arguments,grid_min,grid_max,grid_bins);
      plumed_massert(targetdist_pntrs[i]->getDimension()==1,"the target distribution must be one dimensional");
      targetdist_pntrs[i]->updateTargetDist();
    }
    //
    std::vector<double> bf_integrals = bf_pntrs[0]->getTargetDistributionIntegrals(targetdist_pntrs[i]);
    CoeffsVector targetdist_averages = CoeffsVector("aver.targetdist-"+is,arguments,bf_pntrs,comm,false);
    targetdist_averages.setValues(bf_integrals);
    if(fmt_targetdist_aver.size()>0) {targetdist_averages.setOutputFmt(fmt_targetdist_aver);}
    targetdist_averages.writeToFile(ofile_targetdist_aver,true);
    if(targetdist_pntrs[i]!=NULL) {
      Grid* targetdist_grid_pntr = targetdist_pntrs[i]->getTargetDistGridPntr();
      std::string fname = FileBase::appendSuffix(fname_targetdist,is);
      OFile ofile;
      ofile.link(*this);
      ofile.enforceBackup();
      ofile.open(fname);
      targetdist_grid_pntr->writeToFile(ofile);
      ofile.close();
    }
  }
  ofile_targetdist_aver.close();
  delete arguments[0]; arguments.clear();



}



}
}
