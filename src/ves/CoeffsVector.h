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
#ifndef __PLUMED_ves_CoeffsVector_h
#define __PLUMED_ves_CoeffsVector_h

#include "CoeffsBase.h"

#include <vector>
#include <string>
#include <cmath>


namespace PLMD {

class Action;
class Value;
class IFile;
class OFile;
class Communicator;

namespace ves {

class BasisFunctions;
class CoeffsMatrix;


class CoeffsVector:
  public CoeffsBase
{
public:
private:
  std::vector<double> data;
  //
  unsigned int averaging_counter;
  unsigned int averaging_exp_decay_;
  //
  Communicator& mycomm;
  //
  CoeffsVector& operator=(const CoeffsVector&);
public:
  explicit CoeffsVector(
    const std::string&,
    const std::vector<std::string>&,
    const std::vector<unsigned int>&,
    Communicator&,
    const bool use_counter=false);
  //
  explicit CoeffsVector(
    const std::string&,
    std::vector<Value*>&,
    std::vector<BasisFunctions*>&,
    Communicator&,
    const bool use_counter=false);
  //
  explicit CoeffsVector(
    const std::string&,
    std::vector<std::vector<Value*> >&,
    std::vector<std::vector<BasisFunctions*> >&,
    Communicator&,
    const bool use_counter=false,
    const std::string& multicoeffs_label="bias");
  //
  explicit CoeffsVector(
    const std::string&,
    CoeffsMatrix*,
    Communicator&);
  //
  ~CoeffsVector();
  //
  size_t getSize() const {return numberOfCoeffs();}
  // clear coeffs
  void clear();
  void setAllValuesToZero();
  //
  std::vector<double> getDataAsVector() const {return data;}
  //
  bool sameShape(CoeffsVector&) const;
  bool sameShape(CoeffsMatrix&) const;
  static bool sameShape(CoeffsVector&, CoeffsVector&);
  //
  void resizeCoeffs(const std::vector<unsigned int>&);
  void resizeCoeffs(std::vector<BasisFunctions*>&);
  //
  void sumCommMPI();
  void sumCommMPI(Communicator& cc);
  //
  void sumMultiSimCommMPI(Communicator&);
  // get value
  double getValue(const size_t) const;
  double getValue(const std::vector<unsigned int>&) const;
  double& operator[](const size_t index);
  const double& operator[](const size_t index) const;
  double& operator[](const std::vector<unsigned int>&);
  const double& operator[](const std::vector<unsigned int>&) const;
  double& operator()(const size_t index);
  const double& operator()(const size_t index) const;
  double& operator()(const std::vector<unsigned int>&);
  const double& operator()(const std::vector<unsigned int>&) const;
  // set value
  void setValue(const size_t, const double);
  void setValue(const std::vector<unsigned int>&, const double);
  // add to value
  void addToValue(const size_t, const double);
  void addToValue(const std::vector<unsigned int>&, const double);
  // scale all values
  void scaleAllValues(const double);
  CoeffsVector& operator*=(const double);
  friend CoeffsVector operator*(const double, const CoeffsVector&);
  friend CoeffsVector operator*(const CoeffsVector&, const double);
  CoeffsVector& operator*=(const CoeffsVector&);
  CoeffsVector operator*(const CoeffsVector&) const;
  // set all values
  void setValues(const double);
  void setValues(const std::vector<double>&);
  void setValues(const CoeffsVector&);
  CoeffsVector& operator=(const double);
  CoeffsVector& operator=(const std::vector<double>&);
  // CoeffsVector& operator=(const CoeffsVector&);
  // add to all values
  CoeffsVector operator+() const;
  CoeffsVector operator-() const;
  void addToValues(const double);
  void addToValues(const std::vector<double>&);
  void addToValues(const CoeffsVector&);
  void subtractFromValues(const double);
  void subtractFromValues(const std::vector<double>&);
  void subtractFromValues(const CoeffsVector&);
  CoeffsVector& operator+=(const double);
  friend CoeffsVector operator+(const double, const CoeffsVector&);
  friend CoeffsVector operator+(const CoeffsVector&, const double);
  CoeffsVector& operator+=(const std::vector<double>&);
  friend CoeffsVector operator+(const std::vector<double>&, const CoeffsVector&);
  friend CoeffsVector operator+(const CoeffsVector&, const std::vector<double>&);
  CoeffsVector& operator-=(const double);
  friend CoeffsVector operator-(const double, const CoeffsVector&);
  friend CoeffsVector operator-(const CoeffsVector&, const double);
  CoeffsVector& operator-=(const std::vector<double>&);
  friend CoeffsVector operator-(const std::vector<double>&, const CoeffsVector&);
  friend CoeffsVector operator-(const CoeffsVector&, const std::vector<double>&);
  CoeffsVector& operator+=(const CoeffsVector&);
  CoeffsVector operator+(const CoeffsVector&) const;
  CoeffsVector& operator-=(const CoeffsVector&);
  CoeffsVector operator-(const CoeffsVector&) const;
  //
  void setValuesFromDifferentShape(const CoeffsVector&);
  //
  static void averageVectors(CoeffsVector&, CoeffsVector&);
  static void averageVectors(const std::vector<CoeffsVector*>&);
  //
  double getMinValue() const;
  double getMinValue(size_t&) const;
  double getMinAbsValue() const;
  double getMinAbsValue(size_t&) const;
  //
  double getMaxValue() const;
  double getMaxValue(size_t&) const;
  double getMaxAbsValue() const;
  double getMaxAbsValue(size_t&) const;
  //
  double getNorm() const;
  double getL1Norm() const;
  double getL2Norm() const;
  double getLpNorm(const double) const;
  double getRMS() const;
  //
  void normalizeCoeffs();
  // Random values
  void randomizeValuesGaussian(int);
  //
  void resetAveragingCounter() {averaging_counter=0;}
  void setupExponentiallyDecayingAveraging(const unsigned int averaging_exp_decay_in) {averaging_exp_decay_=averaging_exp_decay_in;}
  void turnOffExponentiallyDecayingAveraging() {averaging_exp_decay_=0;}
  void resetAveraging();
  void addToAverage(const CoeffsVector&);
  //
  size_t countValues(const double) const;

  // file input/output stuff
  void writeToFile(const std::string&, const bool print_description=false, const bool append_file=false, Action* action_pntr=NULL);
  void writeToFile(OFile&, const bool print_description=false);
  void writeToFile(OFile& ofile, CoeffsVector*, const bool print_coeffs_descriptions=false);
  static void writeToFile(const std::string&, const std::vector<CoeffsVector*>&, const bool print_description=false, const bool append_file=false, Action* action_pntr=NULL);
  static void writeToFile(OFile&, const std::vector<CoeffsVector*>&, const bool print_description=false);
private:
  void writeHeaderToFile(OFile&) const;
  static void writeDataToFile(OFile&, const std::vector<CoeffsVector*>&, const bool print_description=false);
public:
  size_t readFromFile(IFile&, const bool ignore_missing_coeffs=false, const bool ignore_header=false);
  size_t readFromFile(const std::string&, const bool ignore_missing_coeffs=false, const bool ignore_header=false);
  size_t readOneSetFromFile(IFile& ifile, const bool ignore_header=false);
private:
  void readHeaderFromFile(IFile&, const bool ignore_coeffs_info=false);
  size_t readDataFromFile(IFile&, const bool ignore_missing_coeffs=false);
public:
  Communicator& getCommunicator() const {return mycomm;}
};


inline
double CoeffsVector::getValue(const size_t index) const {
  return data[index];
}


inline
double CoeffsVector::getValue(const std::vector<unsigned int>& indices) const {
  return data[getIndex(indices)];
}


}

}


#endif
