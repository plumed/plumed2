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

#include "CoeffsBase.h"
#include "BasisFunctions.h"
#include "VesBias.h"

#include "tools/Tools.h"
#include "tools/File.h"
#include "tools/Exception.h"
#include "core/Value.h"

#include <vector>
#include <string>


namespace PLMD {
namespace ves {

CoeffsBase::CoeffsBase(
  const std::string& label,
  const std::vector<std::string>& dimension_labels,
  const std::vector<unsigned int>& indices_shape,
  const bool use_iteration_counter):
  label_(label),
  data_label_(label),
  coeffs_type_(Generic),
  iteration_and_time_active_(use_iteration_counter),
  iteration_opt(0),
  time_md(-1.0),
  active(true),
  action_pntr_(NULL),
  vesbias_pntr_(NULL),
  ndimensions_(0),
  indices_shape_(0),
  ncoeffs_(0),
  coeffs_descriptions_(0),
  dimension_labels_(0),
  args_(0),
  basisf_(0),
  multicoeffs_(false),
  multicoeffs_args_(0),
  multicoeffs_basisf_(0),
  field_type_("type"),
  field_ndimensions_("ndimensions"),
  field_ncoeffs_total_("ncoeffs_total"),
  field_shape_prefix_("shape_"),
  field_time_("time"),
  field_iteration_("iteration"),
  output_fmt_("%30.16e")
{
  initializeIndices(indices_shape,dimension_labels);
  setAllCoeffsDescriptions();
}


CoeffsBase::CoeffsBase(
  const std::string& label,
  std::vector<Value*>& args,
  std::vector<BasisFunctions*>& basisf,
  const bool use_iteration_counter):
  label_(label),
  data_label_(label),
  coeffs_type_(LinearBasisSet),
  iteration_and_time_active_(use_iteration_counter),
  iteration_opt(0),
  time_md(-1.0),
  active(true),
  action_pntr_(NULL),
  vesbias_pntr_(NULL),
  ndimensions_(0),
  indices_shape_(0),
  ncoeffs_(0),
  coeffs_descriptions_(0),
  dimension_labels_(0),
  args_(args),
  basisf_(basisf),
  multicoeffs_(false),
  multicoeffs_args_(0),
  multicoeffs_basisf_(0),
  field_type_("type"),
  field_ndimensions_("ndimensions"),
  field_ncoeffs_total_("ncoeffs_total"),
  field_shape_prefix_("shape_"),
  field_time_("time"),
  field_iteration_("iteration"),
  output_fmt_("%30.16e")
{
  plumed_massert(args_.size()==basisf_.size(),"CoeffsBase: number of arguments do not match number of basis functions");
  std::vector<std::string> dimension_labels(args_.size());
  std::vector<unsigned int> indices_shape(args_.size());
  for(unsigned int i=0; i<args_.size(); i++) {
    dimension_labels[i]=args_[i]->getName();
    indices_shape[i]=basisf_[i]->getNumberOfBasisFunctions();
  }
  initializeIndices(indices_shape,dimension_labels);
  setupBasisFunctionsInfo();
}


CoeffsBase::CoeffsBase(
  const std::string& label,
  std::vector<std::vector<Value*> >& multicoeffs_args,
  std::vector<std::vector<BasisFunctions*> >& multicoeffs_basisf,
  const bool use_iteration_counter,
  const std::string& multicoeffs_label):
  label_(label),
  data_label_(label),
  coeffs_type_(MultiCoeffs_LinearBasisSet),
  iteration_and_time_active_(use_iteration_counter),
  iteration_opt(0),
  time_md(-1.0),
  active(true),
  action_pntr_(NULL),
  vesbias_pntr_(NULL),
  ndimensions_(0),
  indices_shape_(0),
  ncoeffs_(0),
  coeffs_descriptions_(0),
  dimension_labels_(0),
  args_(0),
  basisf_(0),
  multicoeffs_(true),
  multicoeffs_args_(multicoeffs_args),
  multicoeffs_basisf_(multicoeffs_basisf),
  field_type_("type"),
  field_ndimensions_("ndimensions"),
  field_ncoeffs_total_("ncoeffs_total"),
  field_shape_prefix_("shape_"),
  field_time_("time"),
  field_iteration_("iteration"),
  output_fmt_("%30.16e")
{
  plumed_massert(multicoeffs_args.size()==multicoeffs_basisf.size(),"Multi Coeffs: number of arguments vectors does not match number of basis functions vectors");
  unsigned int num_args = multicoeffs_args[0].size();
  unsigned int dim = num_args+1;
  std::vector<std::string> dimension_labels(dim);
  std::vector<unsigned int> indices_shape(dim);
  for(unsigned int i=0; i<num_args; i++) {
    std::string ip;
    Tools::convert(i+1,ip);
    dimension_labels[i] = "bf" + ip;
    indices_shape[i] = multicoeffs_basisf[0][i]->getNumberOfBasisFunctions();
  }
  indices_shape[dim-1] = multicoeffs_args.size();
  dimension_labels[dim-1] = multicoeffs_label;
  for(unsigned int k=0; k<multicoeffs_args.size(); k++) {
    plumed_massert(multicoeffs_args[k].size()==num_args && multicoeffs_basisf[k].size()==num_args,"Multi Coeffs: arguments and basis functions vectors for each bias should be of the same size");
    for(unsigned int i=0; i<num_args; i++) {
      plumed_massert(indices_shape[i]==multicoeffs_basisf[k][i]->getNumberOfBasisFunctions(),"Multi Coeffs: the coeffs shape for each bias should be identical");
    }
  }
  initializeIndices(indices_shape,dimension_labels);
  setupBasisFunctionsInfo();
}


CoeffsBase::~CoeffsBase() {}


void CoeffsBase::initializeIndices(const std::vector<unsigned int>& indices_shape, const std::vector<std::string>& dimension_labels) {
  plumed_massert(indices_shape.size()==dimension_labels.size(),"indices shape and dimension labels must be of the same size");
  ndimensions_=indices_shape.size();
  indices_shape_=indices_shape;
  dimension_labels_=dimension_labels;
  ncoeffs_=1;
  for(unsigned int i=0; i<ndimensions_; i++) {
    ncoeffs_*=indices_shape_[i];
  }
  coeffs_descriptions_.resize(ncoeffs_);
}


void CoeffsBase::reinitializeIndices(const std::vector<unsigned int>& indices_shape_new) {
  plumed_massert(indices_shape_.size()>0,"indices must have been previously initialized before using this function");
  plumed_massert(dimension_labels_.size()>0,"indices must have been previously initialized before using this function");
  plumed_massert(indices_shape_new.size()==numberOfDimensions(),"when resizeing Coeffs the dimension must be constant");
  indices_shape_=indices_shape_new;
  ncoeffs_=1;
  for(unsigned int i=0; i<ndimensions_; i++) {
    ncoeffs_*=indices_shape_[i];
  }
  coeffs_descriptions_.clear();
  coeffs_descriptions_.resize(ncoeffs_);
}


void CoeffsBase::setupBasisFunctionsInfo() {
  plumed_massert(indices_shape_.size()>0,"indices must be initialized before running this function");
  if(coeffs_type_==LinearBasisSet) {
    for(unsigned int i=0; i<numberOfCoeffs(); i++) {
      std::vector<unsigned int> indices=getIndices(i);
      std::string desc;
      desc=basisf_[0]->getBasisFunctionLabel(indices[0]);
      for(unsigned int k=1; k<numberOfDimensions(); k++) {
        desc+="*"+basisf_[k]->getBasisFunctionLabel(indices[k]);
      }
      setCoeffDescription(i,desc);
    }
  }
  else if(coeffs_type_==MultiCoeffs_LinearBasisSet) {
    for(unsigned int i=0; i<numberOfCoeffs(); i++) {
      std::vector<unsigned int> indices=getIndices(i);
      unsigned int mc_id = indices[ndimensions_-1];
      std::string mc_idstr;
      Tools::convert(mc_id,mc_idstr);
      // std::string mc_label = getDimensionLabel(ndimensions_-1);
      std::string postfix = ":" + mc_idstr;
      std::string desc ="";
      desc+=multicoeffs_basisf_[mc_id][0]->getBasisFunctionLabel(indices[0]);
      for(unsigned int k=1; k<(numberOfDimensions()-1); k++) {
        desc+="*"+multicoeffs_basisf_[mc_id][k]->getBasisFunctionLabel(indices[k]);
      }
      desc+=postfix;
      setCoeffDescription(i,desc);
    }
  }
}


void CoeffsBase::resizeIndices(const std::vector<unsigned int>& indices_shape_new) {
  plumed_massert(coeffs_type_==Generic,"Coeffs type must be Generic when resizeing based on a new indices shape vector");
  reinitializeIndices(indices_shape_new);
  setAllCoeffsDescriptions();
}


void CoeffsBase::resizeIndices(std::vector<BasisFunctions*>& basisf_new) {
  plumed_massert(coeffs_type_==LinearBasisSet,"Coeffs type must be LinearBasisSet when resizeing based on a new basis function set");
  basisf_=basisf_new;
  std::vector<unsigned int> indices_shape_new(basisf_new.size());
  for(unsigned int i=0; i<basisf_new.size(); i++) {
    indices_shape_new[i]=basisf_new[i]->getNumberOfBasisFunctions();
  }
  reinitializeIndices(indices_shape_new);
  setupBasisFunctionsInfo();
}


bool CoeffsBase::sameShape(const CoeffsBase& coeffsbase_in) const {
  if(numberOfDimensions()!=coeffsbase_in.numberOfDimensions()) {
    return false;
  }
  if(numberOfCoeffs()!=coeffsbase_in.numberOfCoeffs()) {
    return false;
  }
  for(unsigned int k=0; k<numberOfDimensions(); k++) {
    if(shapeOfIndices(k)!=coeffsbase_in.shapeOfIndices(k)) {
      return false;
    }
  }
  return true;
}


void CoeffsBase::setLabel(const std::string& label) {
  label_=label;
}


void CoeffsBase::setDataLabel(const std::string& data_label) {
  data_label_=data_label;
}


void CoeffsBase::setLabels(const std::string& label) {
  label_=label;
  data_label_=label;
}


void CoeffsBase::setLabels(const std::string& label, const std::string& data_label) {
  label_=label;
  data_label_=data_label;
}


std::string CoeffsBase::getTypeStr() const {
  std::string type_str="";
  if(coeffs_type_==Generic) {
    type_str = "Generic";
  }
  else if(coeffs_type_==LinearBasisSet) {
    type_str = "LinearBasisSet";
  }
  else if(coeffs_type_==MultiCoeffs_LinearBasisSet) {
    type_str = "MultiCoeffs_LinearBasisSet";
  }
  return type_str;
}


void CoeffsBase::setType(const CoeffsType coeffs_type) {
  coeffs_type_=coeffs_type;
}


void CoeffsBase::linkVesBias(VesBias* vesbias_pntr_in) {
  vesbias_pntr_ = vesbias_pntr_in;
  action_pntr_ = static_cast<Action*>(vesbias_pntr_in);
}


void CoeffsBase::linkAction(Action* action_pntr_in) {
  action_pntr_ = action_pntr_in;
}


bool CoeffsBase::indicesExist(const std::vector<unsigned int>& indices) const {
  plumed_dbg_assert(indices.size()==ndimensions_);
  for(unsigned int k=0; k<ndimensions_; k++) {
    if(indices[k]>=indices_shape_[k]) {
      return false;
    }
  }
  return true;
}


void CoeffsBase::setCoeffDescription(const size_t index, const std::string& description) {
  coeffs_descriptions_[index]=description;
}


void CoeffsBase::setCoeffDescription(const std::vector<unsigned int>& indices, const std::string& description) {
  setCoeffDescription(getIndex(indices), description);
}


void CoeffsBase::setAllCoeffsDescriptions(const std::string& description_prefix) {
  for(size_t i=0; i<numberOfCoeffs(); i++) {
    std::vector<unsigned int> indices=getIndices(i);
    std::string is; Tools::convert(indices[0],is);
    std::string desc=description_prefix+"("+is;
    for(unsigned int k=1; k<numberOfDimensions(); k++) {
      Tools::convert(indices[k],is); desc+=","+is;
    }
    desc+=")";
    coeffs_descriptions_[i]=desc;
  }
}


void CoeffsBase::setAllCoeffsDescriptions(const std::vector<std::string>& coeffs_descriptions) {
  plumed_massert(coeffs_descriptions.size()==numberOfCoeffs(),"The coeffs description vector doesn't match the number of coeffs");
  for(size_t i=0; i<numberOfCoeffs(); i++) {
    coeffs_descriptions_[i]=coeffs_descriptions[i];
  }
}


void CoeffsBase::setDimensionLabel(const unsigned int dim_index, const std::string& label) {
  plumed_massert(dim_index<numberOfDimensions(),"Trying to set the label of a dimension outside the number of dimensions");
  dimension_labels_[dim_index]=label;
}


void CoeffsBase::setAllDimensionLabels(const std::string& label_prefix) {
  for(unsigned int i=0; i<numberOfDimensions(); i++) {
    std::string is; Tools::convert(i,is);
    dimension_labels_[i]=label_prefix + is;
  }
}


void CoeffsBase::setAllDimensionLabels(const std::vector<std::string>& labels) {
  for(unsigned int i=0; i<numberOfDimensions(); i++) {
    dimension_labels_[i]=labels[i];
  }
}


void CoeffsBase::writeCoeffsInfoToFile(OFile& ofile) const {
  ofile.addConstantField(field_type_).printField(field_type_,getTypeStr());
  ofile.addConstantField(field_ndimensions_).printField(field_ndimensions_,(int) numberOfDimensions());
  ofile.addConstantField(field_ncoeffs_total_).printField(field_ncoeffs_total_,(int) numberOfCoeffs());
  for(unsigned int k=0; k<numberOfDimensions(); k++) {
    ofile.addConstantField(field_shape_prefix_+getDimensionLabel(k));
    ofile.printField(field_shape_prefix_+getDimensionLabel(k),(int) shapeOfIndices(k));
  }
}


void CoeffsBase::getCoeffsInfoFromFile(IFile& ifile, const bool ignore_coeffs_info) {
  int int_tmp;
  // label
  std::string coeffs_type_f;
  if(ifile.scanField(field_type_,coeffs_type_f)) {
    // empty for now
  }
  else {
    return;
  }
  // number of dimensions
  unsigned int ndimensions_f = 0;
  if(ifile.scanField(field_ndimensions_,int_tmp)) {
    ndimensions_f=(unsigned int) int_tmp;
  }
  else {
    return;
  }
  // total number of coeffs
  size_t ncoeffs_total_f = 0;
  if(ifile.scanField(field_ncoeffs_total_,int_tmp)) {
    ncoeffs_total_f=(size_t) int_tmp;
  }
  else {
    return;
  }
  // shape of indices
  std::vector<unsigned int> indices_shape_f(numberOfDimensions());
  for(unsigned int k=0; k<numberOfDimensions(); k++) {
    if(ifile.scanField(field_shape_prefix_+getDimensionLabel(k),int_tmp)) {
      indices_shape_f[k]=(unsigned int) int_tmp;
    }
    else {
      return;
    }
  }
  if(!ignore_coeffs_info) {
    std::string msg_header="Error when reading in coeffs from file " + ifile.getPath() + ": ";
    checkCoeffsInfo(msg_header, coeffs_type_f, ndimensions_f, ncoeffs_total_f, indices_shape_f);
  }
}


void CoeffsBase::checkCoeffsInfo(const std::string& msg_header, const std::string& coeffs_type_f, const unsigned int ndimensions_f, const size_t ncoeffs_total_f, const std::vector<unsigned int>& indices_shape_f) {

  if(coeffs_type_f != getTypeStr()) {
    std::string msg = msg_header + " coeffs type " + coeffs_type_f + " from file doesn't match the defined value " + getTypeStr();
    plumed_merror(msg);
  }
  if(ndimensions_f != numberOfDimensions() ) {
    std::string s1; Tools::convert(ndimensions_f,s1);
    std::string s2; Tools::convert(numberOfDimensions(),s2);
    std::string msg = msg_header + " the number of dimensions " + s1 + " in file doesn't match the defined value " + s2;
    plumed_merror(msg);
  }
  if(ncoeffs_total_f != numberOfCoeffs() ) {
    std::string s1; Tools::convert(ncoeffs_total_f,s1);
    std::string s2; Tools::convert(numberOfCoeffs(),s2);
    std::string msg = msg_header + " the number of coeffs " + s1 + " in file doesn't match the defined value " + s2;
    plumed_merror(msg);
  }
  for(unsigned int k=0; k<numberOfDimensions(); k++) {
    if(indices_shape_f[k] != shapeOfIndices(k) ) {
      std::string s1; Tools::convert(indices_shape_f[k],s1);
      std::string s2; Tools::convert(shapeOfIndices(k),s2);
      std::string msg = msg_header + " for dimension labeled " + getDimensionLabel(k) + " the shape of indices " + s1 + " in file doesn't match defined value " + s2;
      plumed_merror(msg);
    }
  }
}


void CoeffsBase::writeIterationCounterAndTimeToFile(OFile& ofile) const {
  if(time_md>=0.0) {
    ofile.fmtField("%f");
    ofile.addConstantField(field_time_).printField(field_time_,time_md);
    ofile.fmtField();
  }
  ofile.addConstantField(field_iteration_).printField(field_iteration_,(int) iteration_opt);
}


bool CoeffsBase::getIterationCounterAndTimeFromFile(IFile& ifile) {
  bool field_found=false;
  if(ifile.FieldExist(field_time_)) {
    field_found=true;
    double time_tmp;
    ifile.scanField(field_time_,time_tmp);
    time_md=time_tmp;
  }
  if(ifile.FieldExist(field_iteration_)) {
    field_found=true;
    int iter_tmp;
    ifile.scanField(field_iteration_,iter_tmp);
    iteration_opt=(unsigned int) iter_tmp;
  }
  return field_found;
}


}
}
