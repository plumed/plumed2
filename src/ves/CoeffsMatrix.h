/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2017 The VES code team
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
#ifndef __PLUMED_ves_CoeffsMatrix_h
#define __PLUMED_ves_CoeffsMatrix_h

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
class CoeffsVector;


class CoeffsMatrix:
  public CoeffsBase
{
public:
private:
  std::vector<double> data;
  //
  size_t size_;
  size_t nrows_;
  size_t ncolumns_;
  //
  bool diagonal_;
  //
  unsigned int averaging_counter;
  unsigned int averaging_exp_decay_;
  //
  Communicator& mycomm;
  //
  void setupMatrix();
  //
  CoeffsMatrix& operator=(const CoeffsMatrix&);
public:
  explicit CoeffsMatrix(
    const std::string&,
    const std::vector<std::string>&,
    const std::vector<unsigned int>&,
    Communicator& cc,
    const bool diagonal=true,
    const bool use_iteration_counter=false);
  //
  explicit CoeffsMatrix(
    const std::string&,
    std::vector<Value*>&,
    std::vector<BasisFunctions*>&,
    Communicator& cc,
    const bool diagonal=true,
    const bool use_iteration_counter=false);
  //
  explicit CoeffsMatrix(
    const std::string&,
    std::vector<std::vector<Value*> >& argsv,
    std::vector<std::vector<BasisFunctions*> >& basisfv,
    Communicator& cc,
    const bool diagonal=true,
    const bool use_iteration_counter=false,
    const std::string& multicoeffs_label="bias");
  //
  explicit CoeffsMatrix(
    const std::string&,
    CoeffsVector*,
    Communicator& cc,
    const bool diagonal=true);
  //
  ~CoeffsMatrix();
  //
  size_t getSize() const;
  //
  bool isSymmetric() const;
  bool isDiagonal() const;
  //
  bool sameShape(CoeffsVector&) const;
  bool sameShape(CoeffsMatrix&) const;
  static bool sameShape(CoeffsMatrix&, CoeffsMatrix&);
  static bool sameShape(CoeffsVector&, CoeffsMatrix&);
  static bool sameShape(CoeffsMatrix&, CoeffsVector&);
  //
  void sumCommMPI();
  void sumCommMPI(Communicator&);
  //
  void sumMultiSimCommMPI(Communicator&);
  //
  size_t getMatrixIndex(const size_t, const size_t) const;
  //
  // clear coeffs
  void clear();
  void setAllValuesToZero();
  //
  std::vector<double> getDataAsVector() const {return data;}
  // get value
  double getValue(const size_t, const size_t) const;
  double getValue(const std::vector<unsigned int>&, const std::vector<unsigned int>&) const;
  // set value
  void setValue(const size_t, const size_t, const double);
  void setValue(const std::vector<unsigned int>&, const std::vector<unsigned int>&, const double);
  double& operator()(const size_t, const size_t);
  const double& operator()(const size_t, const size_t) const;
  double& operator()(const std::vector<unsigned int>&, const std::vector<unsigned int>&);
  const double& operator()(const std::vector<unsigned int>&, const std::vector<unsigned int>&) const;
  //
  friend CoeffsVector operator*(const CoeffsMatrix&, const CoeffsVector&);
  // add to value
  void addToValue(const size_t, const size_t, const double);
  void addToValue(const std::vector<unsigned int>&, const std::vector<unsigned int>&, const double);
  // scale all values
  void scaleAllValues(const double);
  CoeffsMatrix& operator*=(const double);
  friend CoeffsMatrix operator*(const double, const CoeffsMatrix&);
  friend CoeffsMatrix operator*(const CoeffsMatrix&, const double);
  CoeffsMatrix& operator*=(const CoeffsMatrix&);
  CoeffsMatrix operator*(const CoeffsMatrix&) const;
  // set all values
  void setValues(const double);
  void setValues(const std::vector<double>&);
  void setValues(const CoeffsMatrix&);
  CoeffsMatrix& operator=(const double);
  CoeffsMatrix& operator=(const std::vector<double>&);
  // CoeffsMatrix& operator=(const CoeffsMatrix&);
  // add to all values
  CoeffsMatrix operator+() const;
  CoeffsMatrix operator-() const;
  void addToValues(const double);
  void addToValues(const std::vector<double>&);
  void addToValues(const CoeffsMatrix&);
  void subtractFromValues(const double);
  void subtractFromValues(const std::vector<double>&);
  void subtractFromValues(const CoeffsMatrix&);
  CoeffsMatrix& operator+=(const double);
  friend CoeffsMatrix operator+(const double, const CoeffsMatrix&);
  friend CoeffsMatrix operator+(const CoeffsMatrix&, const double);
  CoeffsMatrix& operator+=(const std::vector<double>&);
  friend CoeffsMatrix operator+(const std::vector<double>&, const CoeffsMatrix&);
  friend CoeffsMatrix operator+(const CoeffsMatrix&, const std::vector<double>&);
  CoeffsMatrix& operator-=(const double);
  friend CoeffsMatrix operator-(const double, const CoeffsMatrix&);
  friend CoeffsMatrix operator-(const CoeffsMatrix&, const double);
  CoeffsMatrix& operator-=(const std::vector<double>&);
  friend CoeffsMatrix operator-(const std::vector<double>&, const CoeffsMatrix&);
  friend CoeffsMatrix operator-(const CoeffsMatrix&, const std::vector<double>&);
  CoeffsMatrix& operator+=(const CoeffsMatrix&);
  CoeffsMatrix operator+(const CoeffsMatrix&) const;
  CoeffsMatrix& operator-=(const CoeffsMatrix&);
  CoeffsMatrix operator-(const CoeffsMatrix&) const;
  //
  static void averageMatrices(CoeffsMatrix&, CoeffsMatrix&);
  static void averageMatrices(const std::vector<CoeffsMatrix*>&);
  //
  double getMinValue() const;
  double getMaxValue() const;
  //
  void randomizeValuesGaussian(int);
  //
  void resetAveragingCounter() {averaging_counter=0;}
  void setupExponentiallyDecayingAveraging(const unsigned int averaging_exp_decay_in) {averaging_exp_decay_=averaging_exp_decay_in;}
  void turnOffExponentiallyDecayingAveraging() { averaging_exp_decay_=0;}
  void resetAveraging();
  void addToAverage(const CoeffsMatrix&);
  void addToAverage(const CoeffsMatrix&, const unsigned int);
  //
  // file input/output stuff
  void writeToFile(OFile&);
  void writeToFile(const std::string&, const bool append_file=false, Action* action_pntr=NULL);
private:
  void writeDataToFile(OFile&);
  void writeMatrixInfoToFile(OFile&);
  void writeHeaderToFile(OFile&);
  void writeDataDiagonalToFile(OFile&);
  void writeDataFullToFile(OFile&);
public:
  Communicator& getCommunicator() const {return mycomm;}

};
}
}


#endif
