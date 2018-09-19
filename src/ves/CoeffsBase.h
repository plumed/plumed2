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
#ifndef __PLUMED_ves_CoeffsBase_h
#define __PLUMED_ves_CoeffsBase_h

#include <vector>
#include <string>


namespace PLMD {

class Action;
class Value;
class IFile;
class OFile;

namespace ves {

class BasisFunctions;
class VesBias;

/// \ingroup TOOLBOX
class CoeffsBase
{
public:
  // the type of 1D index
  // typedef size_t index_t;
  // typedef unsigned int index_t;
private:
  std::string label_;
  std::string data_label_;
  enum CoeffsType {
    Generic,
    LinearBasisSet,
    MultiCoeffs_LinearBasisSet
  } coeffs_type_;
  //
  bool iteration_and_time_active_;
  unsigned int iteration_opt;
  double time_md;
  //
  bool active;
  //
  Action* action_pntr_;
  VesBias* vesbias_pntr_;
  //
  unsigned int ndimensions_;
  std::vector<unsigned int> indices_shape_;
  size_t ncoeffs_;
  std::vector<std::string> coeffs_descriptions_;
  std::vector<std::string> dimension_labels_;
  //
  std::vector<Value*> args_;
  std::vector<BasisFunctions*> basisf_;
  //
  bool multicoeffs_;
  std::vector<std::vector<Value*> > multicoeffs_args_;
  std::vector<std::vector<BasisFunctions*> >multicoeffs_basisf_;
  // Labels for fields in output/input files
  const std::string field_type_;
  const std::string field_ndimensions_;
  const std::string field_ncoeffs_total_;
  const std::string field_shape_prefix_;
  const std::string field_time_;
  const std::string field_iteration_;
  //
  std::string output_fmt_;
  //
  void initializeIndices(const std::vector<unsigned int>&, const std::vector<std::string>&);
  void reinitializeIndices(const std::vector<unsigned int>&);
public:
  explicit CoeffsBase();
  //
  explicit CoeffsBase(
    const std::string&,
    const std::vector<std::string>&,
    const std::vector<unsigned int>&,
    const bool use_iteration_counter=false);
  //
  explicit CoeffsBase(
    const std::string&,
    std::vector<Value*>&,
    std::vector<BasisFunctions*>&,
    const bool use_iteration_counter=false);
  //
  explicit CoeffsBase(
    const std::string&,
    std::vector<std::vector<Value*> >&,
    std::vector<std::vector<BasisFunctions*> >&,
    const bool use_iteration_counter=false,
    const std::string& multicoeffs_label="bias"
  );
  //
  ~CoeffsBase();
  //
  std::string getLabel() const {return label_;}
  void setLabel(const std::string&);
  std::string getDataLabel() const {return data_label_;};
  void setDataLabel(const std::string&);
  void setLabels(const std::string&);
  void setLabels(const std::string&, const std::string&);
  //
  CoeffsType getType() const {return coeffs_type_;}
  std::string getTypeStr() const;
  void setType(const CoeffsType coeffs_type);
  void linkVesBias(VesBias*);
  void linkAction(Action*);
  VesBias* getPntrToVesBias() const {return vesbias_pntr_;}
  Action* getPntrToAction() const {return action_pntr_;}
  bool isGenericCoeffs() const {return coeffs_type_==Generic;}
  bool isLinearBasisSetCoeffs() const {return coeffs_type_==LinearBasisSet;}
  bool isMultiLinearBasisSetCoeffs() const {return coeffs_type_==MultiCoeffs_LinearBasisSet;}
  //
  std::vector<unsigned int> shapeOfIndices() const {return indices_shape_;}
  unsigned int shapeOfIndices(const unsigned int dim_index) const {return indices_shape_[dim_index];}
  size_t numberOfCoeffs() const {return ncoeffs_;}
  unsigned int numberOfDimensions() const {return ndimensions_;}
  //
  bool isActive() const {return active;}
  void activate() {active=true;}
  void deactivate() {active=false;}
  //
  size_t getIndex(const std::vector<unsigned int>&) const;
  std::vector<unsigned int> getIndices(const size_t) const;
  bool indicesExist(const std::vector<unsigned int>&) const;
  //
  std::string getCoeffDescription(const size_t index) const {return coeffs_descriptions_[index];}
  std::string getCoeffDescription(const std::vector<unsigned int>&) const;
  std::vector<std::string> getAllCoeffsDescriptions() const {return coeffs_descriptions_;}
  void setCoeffDescription(const size_t, const std::string&);
  void setCoeffDescription(const std::vector<unsigned int>&, const std::string&);
  void setAllCoeffsDescriptions(const std::string& description_prefix="C");
  void setAllCoeffsDescriptions(const std::vector<std::string>&);
  //
  std::string getDimensionLabel(const unsigned int) const;
  std::vector<std::string> getAllDimensionLabels() const {return dimension_labels_;}
  void setDimensionLabel(const unsigned int, const std::string&);
  void setAllDimensionLabels(const std::string&);
  void setAllDimensionLabels(const std::vector<std::string>&);
  void writeCoeffsInfoToFile(OFile&) const;
  void writeTimeInfoToFile(OFile&, const double) const;
  void getCoeffsInfoFromFile(IFile&, const bool ignore_coeffs_info=false);
  void checkCoeffsInfo(const std::string&, const std::string&, const unsigned int, const size_t, const std::vector<unsigned int>&);
  //
  void turnOnIterationCounter() {iteration_and_time_active_=true;}
  void turnOffIterationCounter() {iteration_and_time_active_=false;}
  bool isIterationCounterActive() const {return iteration_and_time_active_;}
  void setIterationCounter(const unsigned int);
  void setTime(const double);
  void setIterationCounterAndTime(const unsigned int, const double);
  unsigned int getIterationCounter() const {return iteration_opt;}
  double getTimeValue() const {return time_md;}
  //
  void setOutputFmt(const std::string& ss) { output_fmt_=ss; }
  void resetOutputFmt() {output_fmt_="%30.16e";}
  std::string getOutputFmt() const {return output_fmt_;}
  //
protected:
  void setupBasisFunctionsInfo();
  void resizeIndices(const std::vector<unsigned int>&);
  void resizeIndices(std::vector<BasisFunctions*>&);
  bool sameShape(const CoeffsBase&) const;
  //
  void writeIterationCounterAndTimeToFile(OFile&) const;
  bool getIterationCounterAndTimeFromFile(IFile&);
  //

};

inline
void CoeffsBase::setIterationCounter(const unsigned int iteration_opt_in) {
  iteration_opt=iteration_opt_in;
}

inline
void CoeffsBase::setTime(const double time_md_in) {
  time_md=time_md_in;
}

inline
void CoeffsBase::setIterationCounterAndTime(const unsigned int iteration_opt_in, const double time_md_in) {
  iteration_opt=iteration_opt_in;
  time_md=time_md_in;
}

inline
std::string CoeffsBase::getCoeffDescription(const std::vector<unsigned int>& indices) const {
  return getCoeffDescription(getIndex(indices));
}

inline
std::string CoeffsBase::getDimensionLabel(const unsigned int dim_index) const {
  // plumed_massert(dim_index<numberOfDimensions(),"Trying to get the label of a dimension outside the number of dimensions");
  return dimension_labels_[dim_index];
}


// we are flattening arrays using a column-major order
inline
size_t CoeffsBase::getIndex(const std::vector<unsigned int>& indices) const {
  // plumed_dbg_assert(indices.size()==ndimensions_);
  // for(unsigned int i=0; i<ndimensions_; i++){
  //   if(indices[i]>=indices_shape_[i]){
  //     std::string is;
  //     Tools::convert(i,is);
  //     std::string msg="ERROR: the system is looking for a value outside the indices along the " + is + "dimension!";
  //     plumed_merror(msg);
  //   }
  // }
  size_t index=indices[ndimensions_-1];
  for(unsigned int i=ndimensions_-1; i>0; --i) {
    index=index*indices_shape_[i-1]+indices[i-1];
  }
  return index;
}

// we are flattening arrays using a column-major order
inline
std::vector<unsigned int> CoeffsBase::getIndices(const size_t index) const {
  std::vector<unsigned int> indices(ndimensions_);
  size_t kk=index;
  indices[0]=(index%indices_shape_[0]);
  for(unsigned int i=1; i<ndimensions_-1; ++i) {
    kk=(kk-indices[i-1])/indices_shape_[i-1];
    indices[i]=(kk%indices_shape_[i]);
  }
  if(ndimensions_>=2) {
    indices[ndimensions_-1]=((kk-indices[ndimensions_-2])/indices_shape_[ndimensions_-2]);
  }
  return indices;
}




}
}

#endif
