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

#include "CoeffsVector.h"
#include "CoeffsMatrix.h"
#include "BasisFunctions.h"

#include "tools/Tools.h"
#include "core/Value.h"
#include "tools/File.h"
#include "tools/Exception.h"
#include "tools/Random.h"
#include "tools/Communicator.h"

#include <vector>
#include <cmath>
#include <iostream>
#include <sstream>
#include <cstdio>
#include <cfloat>


namespace PLMD {
namespace ves {

CoeffsVector::CoeffsVector(
  const std::string& label,
  const std::vector<std::string>& dimension_labels,
  const std::vector<unsigned int>& indices_shape,
  Communicator& cc,
  const bool use_iteration_counter):
  CoeffsBase(label,dimension_labels,indices_shape,use_iteration_counter),
  data(0),
  averaging_counter(0),
  averaging_exp_decay_(0),
  mycomm(cc)
{
  clear();
}


CoeffsVector::CoeffsVector(
  const std::string& label,
  std::vector<Value*>& args,
  std::vector<BasisFunctions*>& basisf,
  Communicator& cc,
  const bool use_iteration_counter):
  CoeffsBase(label,args,basisf,use_iteration_counter),
  data(0),
  averaging_counter(0),
  averaging_exp_decay_(0),
  mycomm(cc)
{
  clear();
}


CoeffsVector::CoeffsVector(
  const std::string& label,
  std::vector<std::vector<Value*> >& argsv,
  std::vector<std::vector<BasisFunctions*> >& basisfv,
  Communicator& cc,
  const bool use_iteration_counter,
  const std::string& multicoeffs_label):
  CoeffsBase(label,argsv,basisfv,use_iteration_counter,multicoeffs_label),
  data(0),
  averaging_counter(0),
  averaging_exp_decay_(0),
  mycomm(cc)
{
  clear();
}


CoeffsVector::CoeffsVector(
  const std::string& label,
  CoeffsMatrix* coeffsMat,
  Communicator& cc):
  CoeffsBase( *(static_cast<CoeffsBase*>(coeffsMat)) ),
  data(0),
  averaging_counter(0),
  averaging_exp_decay_(0),
  mycomm(cc)
{
  clear();
}


CoeffsVector::~CoeffsVector() {}


void CoeffsVector::clear() {
  data.resize(getSize());
  for(size_t i=0; i<data.size(); i++) {
    data[i]=0.0;
  }
}


void CoeffsVector::setAllValuesToZero() {
  for(size_t i=0; i<data.size(); i++) {
    data[i]=0.0;
  }
}


bool CoeffsVector::sameShape(CoeffsVector& coeffsvector_in) const {
  return CoeffsBase::sameShape( *(static_cast<CoeffsBase*>(&coeffsvector_in)) );
}


bool CoeffsVector::sameShape(CoeffsMatrix& coeffsmat_in) const {
  return CoeffsBase::sameShape( *(static_cast<CoeffsBase*>(&coeffsmat_in)) );
}


bool CoeffsVector::sameShape(CoeffsVector& coeffsvec0, CoeffsVector& coeffsvec1) {
  return coeffsvec0.sameShape(coeffsvec1);
}


void CoeffsVector::resizeCoeffs(const std::vector<unsigned int>& indices_shape_new) {
  CoeffsVector coeffsVecOld(*this);
  resizeIndices(indices_shape_new);
  clear();
  setValuesFromDifferentShape(coeffsVecOld);
}


void CoeffsVector::resizeCoeffs(std::vector<BasisFunctions*>& basisf_new) {
  CoeffsVector coeffsVecOld(*this);
  resizeIndices(basisf_new);
  clear();
  setValuesFromDifferentShape(coeffsVecOld);
}


void CoeffsVector::sumCommMPI() {
  mycomm.Sum(data);
}


void CoeffsVector::sumCommMPI(Communicator& cc) {
  cc.Sum(data);
}


void CoeffsVector::sumMultiSimCommMPI(Communicator& multi_sim_cc) {
  if(mycomm.Get_rank()==0) {
    double nwalkers = static_cast<double>(multi_sim_cc.Get_size());
    multi_sim_cc.Sum(data);
    scaleAllValues(1.0/nwalkers);
  }
  mycomm.Bcast(data,0);
}


double& CoeffsVector::operator[](const size_t index) {
  plumed_dbg_assert(index<data.size());
  return data[index];
}


const double& CoeffsVector::operator[](const size_t index) const {
  plumed_dbg_assert(index<data.size());
  return data[index];
}


double& CoeffsVector::operator[](const std::vector<unsigned int>& indices) {
  return data[getIndex(indices)];
}


const double& CoeffsVector::operator[](const std::vector<unsigned int>& indices) const {
  return data[getIndex(indices)];
}


double& CoeffsVector::operator()(const size_t index) {
  plumed_dbg_assert(index<data.size());
  return data[index];
}


const double& CoeffsVector::operator()(const size_t index) const {
  plumed_dbg_assert(index<data.size());
  return data[index];
}


double& CoeffsVector::operator()(const std::vector<unsigned int>& indices) {
  return data[getIndex(indices)];
}


const double& CoeffsVector::operator()(const std::vector<unsigned int>& indices) const {
  return data[getIndex(indices)];
}


void CoeffsVector::setValue(const size_t index, const double value) {
  plumed_dbg_assert(index<data.size());
  data[index]=value;
}


void CoeffsVector::setValue(const std::vector<unsigned int>& indices, const double value) {
  setValue(getIndex(indices),value);
}


void CoeffsVector::addToValue(const size_t index, const double value) {
  plumed_dbg_assert(index<data.size());
  data[index]+=value;
}


void CoeffsVector::addToValue(const std::vector<unsigned int>& indices, const double value) {
  addToValue(getIndex(indices),value);
}


void CoeffsVector::scaleAllValues(const double scalef) {
  for(size_t i=0; i<data.size(); i++) {
    data[i]*=scalef;
  }
}


CoeffsVector& CoeffsVector::operator*=(const double scalef) {
  scaleAllValues(scalef);
  return *this;
}


CoeffsVector operator*(const double scalef, const CoeffsVector& coeffsvector) {
  return CoeffsVector(coeffsvector)*=scalef;
}


CoeffsVector operator*(const CoeffsVector& coeffsvector, const double scalef) {
  return scalef*coeffsvector;
}


CoeffsVector& CoeffsVector::operator*=(const CoeffsVector& other_coeffsvector) {
  plumed_massert(data.size()==other_coeffsvector.getSize(),"Coeffs vectors do not have the same size");
  for(size_t i=0; i<data.size(); i++) {
    data[i]*=other_coeffsvector.data[i];
  }
  return *this;
}


CoeffsVector CoeffsVector::operator*(const CoeffsVector& other_coeffsvector) const {
  return CoeffsVector(*this)*=other_coeffsvector;
}


void CoeffsVector::setValues(const double value) {
  for(size_t i=0; i<data.size(); i++) {
    data[i]=value;
  }
}


void CoeffsVector::setValues(const std::vector<double>& values) {
  plumed_massert( data.size()==values.size(), "Incorrect size");
  for(size_t i=0; i<data.size(); i++) {
    data[i]=values[i];
  }
}


void CoeffsVector::setValues(const CoeffsVector& other_coeffsvector) {
  plumed_massert( data.size()==other_coeffsvector.getSize(), "Incorrect size");
  for(size_t i=0; i<data.size(); i++) {
    data[i]=other_coeffsvector.data[i];
  }
}


CoeffsVector& CoeffsVector::operator=(const double value) {
  setValues(value);
  return *this;
}


CoeffsVector& CoeffsVector::operator=(const std::vector<double>& values) {
  setValues(values);
  return *this;
}


// CoeffsVector& CoeffsVector::operator=(const CoeffsVector& other_coeffsvector) {
//   setValues(other_coeffsvector);
//   return *this;
// }


CoeffsVector CoeffsVector::operator+() const {
  return *this;
}


CoeffsVector CoeffsVector::operator-() const {
  return CoeffsVector(*this)*=-1.0;
}


void CoeffsVector::addToValues(const double value) {
  for(size_t i=0; i<data.size(); i++) {
    data[i]+=value;
  }
}


void CoeffsVector::addToValues(const std::vector<double>& values) {
  plumed_massert( data.size()==values.size(), "Incorrect size");
  for(size_t i=0; i<data.size(); i++) {
    data[i]+=values[i];
  }
}


void CoeffsVector::addToValues(const CoeffsVector& other_coeffsvector) {
  plumed_massert( data.size()==other_coeffsvector.getSize(), "Incorrect size");
  for(size_t i=0; i<data.size(); i++) {
    data[i]+=other_coeffsvector.data[i];
  }
}


void CoeffsVector::subtractFromValues(const double value) {
  for(size_t i=0; i<data.size(); i++) {
    data[i]-=value;
  }
}


void CoeffsVector::subtractFromValues(const std::vector<double>& values) {
  plumed_massert( data.size()==values.size(), "Incorrect size");
  for(size_t i=0; i<data.size(); i++) {
    data[i]-=values[i];
  }
}


void CoeffsVector::subtractFromValues(const CoeffsVector& other_coeffsvector) {
  plumed_massert( data.size()==other_coeffsvector.getSize(), "Incorrect size");
  for(size_t i=0; i<data.size(); i++) {
    data[i]-=other_coeffsvector.data[i];
  }
}


CoeffsVector& CoeffsVector::operator+=(const double value) {
  addToValues(value);
  return *this;
}


CoeffsVector operator+(const double value, const CoeffsVector& coeffsvector) {
  return coeffsvector+value;
}


CoeffsVector operator+(const CoeffsVector& coeffsvector, const double value) {
  return CoeffsVector(coeffsvector)+=value;
}


CoeffsVector& CoeffsVector::operator+=(const std::vector<double>& values) {
  addToValues(values);
  return *this;
}


CoeffsVector operator+(const std::vector<double>& values, const CoeffsVector& coeffsvector) {
  return coeffsvector+values;
}


CoeffsVector operator+(const CoeffsVector& coeffsvector, const std::vector<double>& values) {
  return CoeffsVector(coeffsvector)+=values;
}


CoeffsVector& CoeffsVector::operator-=(const double value) {
  subtractFromValues(value);
  return *this;
}


CoeffsVector operator-(const double value, const CoeffsVector& coeffsvector) {
  return -1.0*coeffsvector+value;
}


CoeffsVector operator-(const CoeffsVector& coeffsvector, const double value) {
  return CoeffsVector(coeffsvector)-=value;
}


CoeffsVector& CoeffsVector::operator-=(const std::vector<double>& values) {
  subtractFromValues(values);
  return *this;
}


CoeffsVector operator-(const std::vector<double>& values, const CoeffsVector& coeffsvector) {
  return -1.0*coeffsvector+values;
}


CoeffsVector operator-(const CoeffsVector& coeffsvector, const std::vector<double>& values) {
  return CoeffsVector(coeffsvector)-=values;
}


CoeffsVector& CoeffsVector::operator+=(const CoeffsVector& other_coeffsvector) {
  addToValues(other_coeffsvector);
  return *this;
}


CoeffsVector CoeffsVector::operator+(const CoeffsVector& other_coeffsvector) const {
  return CoeffsVector(*this)+=other_coeffsvector;
}


CoeffsVector& CoeffsVector::operator-=(const CoeffsVector& other_coeffsvector) {
  subtractFromValues(other_coeffsvector);
  return *this;
}


CoeffsVector CoeffsVector::operator-(const CoeffsVector& other_coeffsvector) const {
  return CoeffsVector(*this)-=other_coeffsvector;
}


void CoeffsVector::setValuesFromDifferentShape(const CoeffsVector& other_coeffsvector) {
  plumed_massert(numberOfDimensions()==other_coeffsvector.numberOfDimensions(),"both coeffs vector need to have the same dimension");
  for(size_t i=0; i<data.size(); i++) {
    std::vector<unsigned int> indices=getIndices(i);
    if(other_coeffsvector.indicesExist(indices)) {
      size_t oidx = other_coeffsvector.getIndex(indices);
      data[i] = other_coeffsvector.data[oidx];
    }
  }
}


void CoeffsVector::averageVectors(CoeffsVector& coeffsvec0, CoeffsVector& coeffsvec1) {
  plumed_massert(sameShape(coeffsvec0,coeffsvec1),"both CoeffsVector objects need to have the same shape");
  for(size_t i=0; i<coeffsvec0.getSize(); i++) {
    coeffsvec0.data[i] = coeffsvec1.data[i] = 0.5 * (coeffsvec0.data[i]+coeffsvec1.data[i]);
  }
}


void CoeffsVector::averageVectors(const std::vector<CoeffsVector*>& coeffsvecSet) {
  const double norm_factor = 1.0/static_cast<double>(coeffsvecSet.size());
  for(unsigned int k=1; k<coeffsvecSet.size(); k++) {
    plumed_massert(coeffsvecSet[0]->sameShape(*coeffsvecSet[k]),"All CoeffsVector objects need to have the same shape");
  }
  for(size_t i=0; i<coeffsvecSet[0]->getSize(); i++) {
    double value = 0.0;
    for(unsigned int k=0; k<coeffsvecSet.size(); k++) {
      value += coeffsvecSet[k]->data[i];
    }
    value *= norm_factor;
    for(unsigned int k=0; k<coeffsvecSet.size(); k++) {
      coeffsvecSet[k]->data[i] = value;
    }
  }
}


double CoeffsVector::getMinValue() const {
  size_t min_index=0;
  return getMinValue(min_index);
}


double CoeffsVector::getMinValue(size_t& min_index) const {
  min_index=0;
  double min_value=DBL_MAX;
  for(size_t i=0; i<data.size(); i++) {
    if(data[i]<min_value) {
      min_value=data[i];
      min_index=i;
    }
  }
  return min_value;
}


double CoeffsVector::getMinAbsValue() const {
  size_t min_index=0;
  return getMinAbsValue(min_index);
}


double CoeffsVector::getMinAbsValue(size_t& min_index) const {
  min_index=0;
  double min_value=DBL_MAX;
  for(size_t i=0; i<data.size(); i++) {
    if(std::abs(data[i])<min_value) {
      min_value=std::abs(data[i]);
      min_index=i;
    }
  }
  return min_value;
}


double CoeffsVector::getMaxValue() const {
  size_t max_index=0;
  return getMaxValue(max_index);
}


double CoeffsVector::getMaxValue(size_t& max_index) const {
  max_index=0;
  double max_value=DBL_MIN;
  for(size_t i=0; i<data.size(); i++) {
    if(data[i]>max_value) {
      max_value=data[i];
      max_index=i;
    }
  }
  return max_value;
}


double CoeffsVector::getMaxAbsValue() const {
  size_t max_index=0;
  return getMaxAbsValue(max_index);
}


double CoeffsVector::getMaxAbsValue(size_t& max_index) const {
  max_index=0;
  double max_value=0.0;
  for(size_t i=0; i<data.size(); i++) {
    if(std::abs(data[i])>max_value) {
      max_value=std::abs(data[i]);
      max_index=i;
    }
  }
  return max_value;
}


double CoeffsVector::getNorm() const {
  return getL2Norm();
}


double CoeffsVector::getL1Norm() const {
  double norm=0.0;
  for(size_t i=0; i<data.size(); i++) {
    norm+=std::abs(data[i]);
  }
  return norm;
}


double CoeffsVector::getL2Norm() const {
  double norm=0.0;
  for(size_t i=0; i<data.size(); i++) {
    norm+=data[i]*data[i];
  }
  norm=sqrt(norm);
  return norm;
}


double CoeffsVector::getLpNorm(const double p) const {
  double norm=0.0;
  for(size_t i=0; i<data.size(); i++) {
    norm+=pow(data[i],p);
  }
  norm=pow(norm,(1.0/p));
  return norm;
}


double CoeffsVector::getRMS() const {
  return getNorm()/sqrt(numberOfCoeffs());
}


void CoeffsVector::normalizeCoeffs() {
  double norm=getNorm();
  scaleAllValues(norm);
}


void CoeffsVector::randomizeValuesGaussian(int randomSeed) {
  Random rnd;
  if (randomSeed<0) {randomSeed = -randomSeed;}
  rnd.setSeed(-randomSeed);
  for(size_t i=0; i<data.size(); i++) {
    data[i]=rnd.Gaussian();
  }
}


void CoeffsVector::resetAveraging() {
  clear();
  resetAveragingCounter();
}


void CoeffsVector::addToAverage(const CoeffsVector& coeffsvec) {
  plumed_massert( data.size()==coeffsvec.getSize(), "Incorrect size");
  //
  double aver_decay = 1.0 / ( static_cast<double>(averaging_counter) + 1.0 );
  if(averaging_exp_decay_>0 &&  (averaging_counter+1 > averaging_exp_decay_) ) {
    aver_decay = 1.0 / static_cast<double>(averaging_exp_decay_);
  }
  //
  for(size_t i=0; i<data.size(); i++) {
    data[i]+=(coeffsvec.data[i]-data[i])*aver_decay;
  }
  averaging_counter++;
}


size_t CoeffsVector::countValues(const double value) const {
  size_t numvalues=0;
  for(size_t i=0; i<data.size(); i++) {
    if(data[i]==value) {
      numvalues++;
    }
  }
  return numvalues;
}


void CoeffsVector::writeToFile(const std::string& filepath, const bool print_coeffs_descriptions, const bool append_file, Action* action_pntr) {
  OFile file;
  if(action_pntr!=NULL) {
    file.link(*action_pntr);
  }
  else {
    file.link(mycomm);
  }
  if(append_file) { file.enforceRestart(); }
  file.open(filepath);
  writeToFile(file,print_coeffs_descriptions);
  file.close();
}


void CoeffsVector::writeToFile(OFile& ofile, const bool print_coeffs_descriptions) {
  std::vector<CoeffsVector*> CoeffsSetTmp;
  CoeffsSetTmp.push_back(this);
  writeHeaderToFile(ofile);
  writeDataToFile(ofile,CoeffsSetTmp,print_coeffs_descriptions);
}


void CoeffsVector::writeToFile(OFile& ofile, CoeffsVector* aux_coeffsvector, const bool print_coeffs_descriptions) {
  std::vector<CoeffsVector*> CoeffsSetTmp;
  CoeffsSetTmp.push_back(this);
  CoeffsSetTmp.push_back(aux_coeffsvector);
  writeHeaderToFile(ofile);
  writeDataToFile(ofile,CoeffsSetTmp,print_coeffs_descriptions);
}


void CoeffsVector::writeToFile(const std::string& filepath, const std::vector<CoeffsVector*>& coeffsvecSet, const bool print_coeffs_descriptions, const bool append_file, Action* action_pntr) {
  OFile file;
  if(action_pntr!=NULL) {
    file.link(*action_pntr);
  }
  else {
    file.link(coeffsvecSet[0]->getCommunicator());
  }
  if(append_file) { file.enforceRestart(); }
  file.open(filepath);
  writeToFile(file,coeffsvecSet,print_coeffs_descriptions);
  file.close();
}


void CoeffsVector::writeToFile(OFile& ofile, const std::vector<CoeffsVector*>& coeffsvecSet, const bool print_coeffs_descriptions) {
  for(unsigned int k=1; k<coeffsvecSet.size(); k++) {
    plumed_massert(coeffsvecSet[k]->sameShape(*coeffsvecSet[0]),"Error in writing a set of coeffs to file: The coeffs do not have the same shape and size");
  }
  coeffsvecSet[0]->writeHeaderToFile(ofile);
  writeDataToFile(ofile,coeffsvecSet, print_coeffs_descriptions);
}


void CoeffsVector::writeHeaderToFile(OFile& ofile) const {
  ofile.clearFields();
  if(isIterationCounterActive()) {
    writeIterationCounterAndTimeToFile(ofile);
  }
  writeCoeffsInfoToFile(ofile);
}


void CoeffsVector::writeDataToFile(OFile& ofile, const std::vector<CoeffsVector*>& coeffsvecSet, const bool print_coeffs_descriptions) {
  //
  std::string field_indices_prefix = "idx_";
  std::string field_index = "index";
  std::string field_description = "description";
  //
  std::string int_fmt = "%8d";
  std::string str_separate = "#!-------------------";
  //
  unsigned int numvec = coeffsvecSet.size();
  unsigned int numdim = coeffsvecSet[0]->numberOfDimensions();
  unsigned int numcoeffs = coeffsvecSet[0]->getSize();
  std::vector<std::string> coeffs_descriptions = coeffsvecSet[0]->getAllCoeffsDescriptions();
  std::string output_fmt = coeffsvecSet[0]->getOutputFmt();
  std::vector<std::string> coeffs_datalabels(numvec);
  for(unsigned int k=0; k<numvec; k++) {
    coeffs_datalabels[k] = coeffsvecSet[k]->getDataLabel();
  }
  //
  char* s1 = new char[20];
  std::vector<unsigned int> indices(numdim);
  std::vector<std::string> ilabels(numdim);
  for(unsigned int k=0; k<numdim; k++) {
    ilabels[k]=field_indices_prefix+coeffsvecSet[0]->getDimensionLabel(k);
  }
  //
  for(size_t i=0; i<numcoeffs; i++) {
    indices=coeffsvecSet[0]->getIndices(i);
    for(unsigned int k=0; k<numdim; k++) {
      sprintf(s1,int_fmt.c_str(),indices[k]);
      ofile.printField(ilabels[k],s1);
    }
    for(unsigned int l=0; l<numvec; l++) {
      ofile.fmtField(" "+output_fmt).printField(coeffs_datalabels[l],coeffsvecSet[l]->getValue(i));
    }
    sprintf(s1,int_fmt.c_str(),i); ofile.printField(field_index,s1);
    if(print_coeffs_descriptions) { ofile.printField(field_description,"  "+coeffs_descriptions[i]);}
    ofile.printField();
  }
  ofile.fmtField();
  // blank line between iterations to allow proper plotting with gnuplot
  ofile.printf("%s\n",str_separate.c_str());
  ofile.printf("\n");
  ofile.printf("\n");
  delete [] s1;
}


size_t CoeffsVector::readFromFile(IFile& ifile, const bool ignore_missing_coeffs, const bool ignore_header) {
  ifile.allowIgnoredFields();
  size_t ncoeffs_read=0;
  while(ifile) {
    if(!ignore_header) {readHeaderFromFile(ifile);}
    if(ifile) {
      ncoeffs_read=readDataFromFile(ifile,ignore_missing_coeffs);
    }
  }
  return ncoeffs_read;
}


size_t CoeffsVector::readOneSetFromFile(IFile& ifile, const bool ignore_header) {
  ifile.allowIgnoredFields();
  size_t ncoeffs_read=0;
  if(ifile) {
    if(!ignore_header) {readHeaderFromFile(ifile);}
    if(ifile) {ncoeffs_read=readDataFromFile(ifile,false);}
  }
  return ncoeffs_read;
}


size_t CoeffsVector::readFromFile(const std::string& filepath, const bool ignore_missing_coeffs, const bool ignore_header) {
  IFile file;
  file.link(mycomm);
  file.open(filepath);
  size_t ncoeffs_read=readFromFile(file,ignore_missing_coeffs, ignore_header);
  file.close();
  return ncoeffs_read;
}


void CoeffsVector::readHeaderFromFile(IFile& ifile, const bool ignore_coeffs_info) {
  if(ifile && isIterationCounterActive()) {
    getIterationCounterAndTimeFromFile(ifile);
  }
  if(ifile) {
    getCoeffsInfoFromFile(ifile,ignore_coeffs_info);
  }
}


size_t CoeffsVector::readDataFromFile(IFile& ifile, const bool ignore_missing_coeffs) {
  //
  std::string field_indices_prefix = "idx_";
  std::string field_coeffs = getDataLabel();
  std::string field_index = "index";
  std::string field_description = "description";
  //
  std::vector<std::string> ilabels(numberOfDimensions());
  for(unsigned int k=0; k<numberOfDimensions(); k++) {
    ilabels[k]=field_indices_prefix+getDimensionLabel(k);
  }
  //
  std::vector<unsigned int> indices(numberOfDimensions());
  double coeff_tmp=0.0;
  std::string str_tmp;
  size_t ncoeffs_read=0;
  //
  while(ifile.scanField(field_coeffs,coeff_tmp)) {
    int idx_tmp;
    for(unsigned int k=0; k<numberOfDimensions(); k++) {
      ifile.scanField(ilabels[k],idx_tmp);
      indices[k] = static_cast<unsigned int>(idx_tmp);
    }
    data[getIndex(indices)] = coeff_tmp;
    ifile.scanField(field_index,idx_tmp);
    if(getIndex(indices)!=static_cast<unsigned int>(idx_tmp)) {
      std::string is1; Tools::convert(idx_tmp,is1);
      std::string msg="ERROR: problem with indices at index " + is1 + " when reading coefficients from file";
      plumed_merror(msg);
    }
    if(ifile.FieldExist(field_description)) { ifile.scanField(field_description,str_tmp); }
    //
    ifile.scanField();
    ncoeffs_read++;
    if(ncoeffs_read==numberOfCoeffs()) {
      if((static_cast<unsigned int>(idx_tmp)+1)!=numberOfCoeffs()) {
        plumed_merror("something strange about the coefficient file that is being read in, perhaps multiple entries or missing values");
      }
      break;
    }
  }
  // checks on the coeffs read
  if(ncoeffs_read>0 &&!ignore_missing_coeffs && ncoeffs_read < numberOfCoeffs()) {
    plumed_merror("ERROR: missing coefficients when reading from file");
  }
  //
  return ncoeffs_read;
}


}
}
