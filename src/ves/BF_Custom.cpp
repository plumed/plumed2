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

#include "core/ActionRegister.h"
#include "lepton/Lepton.h"


namespace PLMD {
namespace ves {

//+PLUMEDOC VES_BASISF BF_CUSTOM
/*
Basis functions given by arbitrary mathematical expressions.

This allows you to define basis functions using arbitrary mathematical expressions
that are parsed using the lepton library.
The basis functions
\f$f_{i}(x)\f$ are given in mathematical expressions with _x_ as a variable using
the numbered FUNC keywords that start from
FUNC1. Consistent with other basis functions is \f$f_{0}(x)=1\f$ defined as
the constant. The interval on which the basis functions are defined is
given using the MINIMUM and MAXIMUM keywords.

Using the TRANSFORM keyword it is possible to define a function \f$x(t)\f$ that
is used to transform the argument before calculating the basis functions
values. The variables _min_ and _max_ can be used to indicate the minimum
and the maximum of the interval. By default the arguments are not transformed,
i.e. \f$x(t)=t\f$.

For periodic basis functions you should use the PERIODIC flag to indicate
that they are periodic.

The basis functions \f$f_{i}(x)\f$ and the transform function \f$x(t)\f$ need
to be well behaved in the interval on which the basis functions are defined,
e.g. not result in a not a number (nan) or infinity (inf).
The code will not perform checks to make sure that this is the case unless the
flag CHECK_NAN_INF is enabled.

\par Examples

Defining Legendre polynomial basis functions of order 6 using BF_CUSTOM
where the appropriate transform function is given by the TRANSFORM keyword.
This is just an example of what can be done, in practice you should use
\ref BF_LEGENDRE for Legendre polynomial basis functions.
\plumedfile
BF_CUSTOM ...
 TRANSFORM=(t-(min+max)/2)/((max-min)/2)
 FUNC1=x
 FUNC2=(1/2)*(3*x^2-1)
 FUNC3=(1/2)*(5*x^3-3*x)
 FUNC4=(1/8)*(35*x^4-30*x^2+3)
 FUNC5=(1/8)*(63*x^5-70*x^3+15*x)
 FUNC6=(1/16)*(231*x^6-315*x^4+105*x^2-5)
 MINIMUM=-4.0
 MAXIMUM=4.0
 LABEL=bf1
... BF_CUSTOM
\endplumedfile


Defining Fourier basis functions of order 3 using BF_CUSTOM where the
periodicity is indicated using the PERIODIC flag. This is just an example
of what can be done, in practice you should use \ref BF_FOURIER
for Fourier basis functions.
\plumedfile
BF_CUSTOM ...
 FUNC1=cos(x)
 FUNC2=sin(x)
 FUNC3=cos(2*x)
 FUNC4=sin(2*x)
 FUNC5=cos(3*x)
 FUNC6=sin(3*x)
 MINIMUM=-pi
 MAXIMUM=+pi
 LABEL=bf1
 PERIODIC
... BF_CUSTOM
\endplumedfile


*/
//+ENDPLUMEDOC

class BF_Custom : public BasisFunctions {
private:
  lepton::CompiledExpression transf_value_expression_;
  lepton::CompiledExpression transf_deriv_expression_;
  double* transf_value_lepton_ref_;
  double* transf_deriv_lepton_ref_;
  std::vector<lepton::CompiledExpression> bf_values_expressions_;
  std::vector<lepton::CompiledExpression> bf_derivs_expressions_;
  std::vector<double*> bf_values_lepton_ref_;
  std::vector<double*> bf_derivs_lepton_ref_;
  std::string variable_str_;
  std::string transf_variable_str_;
  bool do_transf_;
  bool check_nan_inf_;
public:
  static void registerKeywords( Keywords&);
  explicit BF_Custom(const ActionOptions&);
  ~BF_Custom() {};
  void getAllValues(const double, double&, bool&, std::vector<double>&, std::vector<double>&) const override;
};

PLUMED_REGISTER_ACTION(BF_Custom,"BF_CUSTOM")

void BF_Custom::registerKeywords(Keywords& keys) {
  BasisFunctions::registerKeywords(keys);
  keys.remove("ORDER");
  keys.add("numbered","FUNC","The basis functions f_i(x) given in mathematical expressions using _x_ as a variable.");
  keys.add("optional","TRANSFORM","An optional function that can be used to transform the argument before calculating the basis function values. You should use _t_ as a variable. You can use the variables _min_ and _max_ to give the minimum and the maximum of the interval.");
  keys.addFlag("PERIODIC",false,"Indicate that the basis functions are periodic.");
  keys.addFlag("CHECK_NAN_INF",false,"Check that the basis functions do not result in a not a number (nan) or infinity (inf).");
  keys.remove("NUMERICAL_INTEGRALS");
}


BF_Custom::BF_Custom(const ActionOptions&ao):
  PLUMED_VES_BASISFUNCTIONS_INIT(ao),
  transf_value_lepton_ref_(nullptr),
  transf_deriv_lepton_ref_(nullptr),
  bf_values_expressions_(0),
  bf_derivs_expressions_(0),
  bf_values_lepton_ref_(0,nullptr),
  bf_derivs_lepton_ref_(0,nullptr),
  variable_str_("x"),
  transf_variable_str_("t"),
  do_transf_(false),
  check_nan_inf_(false)
{
  std::vector<std::string> bf_str;
  std::string str_t1="1";
  bf_str.push_back(str_t1);
  for(int i=1;; i++) {
    std::string str_t2;
    if(!parseNumbered("FUNC",i,str_t2)) {break;}
    std::string is; Tools::convert(i,is);
    addKeywordToList("FUNC"+is,str_t2);
    bf_str.push_back(str_t2);
  }
  //
  if(bf_str.size()==1) {plumed_merror(getName()+" with label "+getLabel()+": No FUNC keywords given");}

  setOrder(bf_str.size()-1);
  setNumberOfBasisFunctions(getOrder()+1);
  setIntrinsicInterval(intervalMin(),intervalMax());
  bool periodic = false;
  parseFlag("PERIODIC",periodic); addKeywordToList("PERIODIC",periodic);
  if(periodic) {setPeriodic();}
  else {setNonPeriodic();}
  setIntervalBounded();
  setType("custom_functions");
  setDescription("Custom Functions");
  //
  std::vector<std::string> bf_values_parsed(getNumberOfBasisFunctions());
  std::vector<std::string> bf_derivs_parsed(getNumberOfBasisFunctions());
  bf_values_parsed[0] = "1";
  bf_derivs_parsed[0] = "0";
  //
  bf_values_expressions_.resize(getNumberOfBasisFunctions());
  bf_derivs_expressions_.resize(getNumberOfBasisFunctions());
  bf_values_lepton_ref_.resize(getNumberOfBasisFunctions());
  bf_derivs_lepton_ref_.resize(getNumberOfBasisFunctions());
  //
  for(unsigned int i=1; i<getNumberOfBasisFunctions(); i++) {
    std::string is; Tools::convert(i,is);
    try {
      lepton::ParsedExpression pe_value = lepton::Parser::parse(bf_str[i]).optimize(lepton::Constants());
      std::ostringstream tmp_stream; tmp_stream << pe_value;
      bf_values_parsed[i] = tmp_stream.str();
      bf_values_expressions_[i] = pe_value.createCompiledExpression();
    }
    catch(PLMD::lepton::Exception& exc) {
      plumed_merror("There was some problem in parsing the function "+bf_str[i]+" given in FUNC"+is + " with lepton");
    }

    std::vector<std::string> var_str;
    for(auto &p: bf_values_expressions_[i].getVariables()) {
      var_str.push_back(p);
    }
    if(var_str.size()!=1) {
      plumed_merror("Problem with function "+bf_str[i]+" given in FUNC"+is+": there should only be one variable");
    }
    if(var_str[0]!=variable_str_) {
      plumed_merror("Problem with function "+bf_str[i]+" given in FUNC"+is+": you should use "+variable_str_+" as a variable");
    }

    try {
      lepton::ParsedExpression pe_deriv = lepton::Parser::parse(bf_str[i]).differentiate(variable_str_).optimize(lepton::Constants());
      std::ostringstream tmp_stream2; tmp_stream2 << pe_deriv;
      bf_derivs_parsed[i] = tmp_stream2.str();
      bf_derivs_expressions_[i] = pe_deriv.createCompiledExpression();
    }
    catch(PLMD::lepton::Exception& exc) {
      plumed_merror("There was some problem in parsing the derivative of the function "+bf_str[i]+" given in FUNC"+is + " with lepton");
    }

    try {
      bf_values_lepton_ref_[i] = &bf_values_expressions_[i].getVariableReference(variable_str_);
    } catch(PLMD::lepton::Exception& exc) {}

    try {
      bf_derivs_lepton_ref_[i] = &bf_derivs_expressions_[i].getVariableReference(variable_str_);
    } catch(PLMD::lepton::Exception& exc) {}

  }

  std::string transf_value_parsed;
  std::string transf_deriv_parsed;
  std::string transf_str;
  parse("TRANSFORM",transf_str);
  if(transf_str.size()>0) {
    do_transf_ = true;
    addKeywordToList("TRANSFORM",transf_str);
    for(unsigned int k=0;; k++) {
      if(transf_str.find("min")!=std::string::npos) {transf_str.replace(transf_str.find("min"), std::string("min").length(),intervalMinStr());}
      else {break;}
    }
    for(unsigned int k=0;; k++) {
      if(transf_str.find("max")!=std::string::npos) {transf_str.replace(transf_str.find("max"), std::string("max").length(),intervalMaxStr());}
      else {break;}
    }

    try {
      lepton::ParsedExpression pe_value = lepton::Parser::parse(transf_str).optimize(lepton::Constants());;
      std::ostringstream tmp_stream; tmp_stream << pe_value;
      transf_value_parsed = tmp_stream.str();
      transf_value_expression_ = pe_value.createCompiledExpression();
    }
    catch(PLMD::lepton::Exception& exc) {
      plumed_merror("There was some problem in parsing the function "+transf_str+" given in TRANSFORM with lepton");
    }

    std::vector<std::string> var_str;
    for(auto &p: transf_value_expression_.getVariables()) {
      var_str.push_back(p);
    }
    if(var_str.size()!=1) {
      plumed_merror("Problem with function "+transf_str+" given in TRANSFORM: there should only be one variable");
    }
    if(var_str[0]!=transf_variable_str_) {
      plumed_merror("Problem with function "+transf_str+" given in TRANSFORM: you should use "+transf_variable_str_+" as a variable");
    }

    try {
      lepton::ParsedExpression pe_deriv = lepton::Parser::parse(transf_str).differentiate(transf_variable_str_).optimize(lepton::Constants());;
      std::ostringstream tmp_stream2; tmp_stream2 << pe_deriv;
      transf_deriv_parsed = tmp_stream2.str();
      transf_deriv_expression_ = pe_deriv.createCompiledExpression();
    }
    catch(PLMD::lepton::Exception& exc) {
      plumed_merror("There was some problem in parsing the derivative of the function "+transf_str+" given in TRANSFORM with lepton");
    }

    try {
      transf_value_lepton_ref_ = &transf_value_expression_.getVariableReference(transf_variable_str_);
    } catch(PLMD::lepton::Exception& exc) {}

    try {
      transf_deriv_lepton_ref_ = &transf_deriv_expression_.getVariableReference(transf_variable_str_);
    } catch(PLMD::lepton::Exception& exc) {}

  }
  //
  log.printf("  Using the following functions [lepton parsed function and derivative]:\n");
  for(unsigned int i=0; i<getNumberOfBasisFunctions(); i++) {
    log.printf("   %u:  %s   [   %s   |   %s   ] \n",i,bf_str[i].c_str(),bf_values_parsed[i].c_str(),bf_derivs_parsed[i].c_str());

  }
  //
  if(do_transf_) {
    log.printf("  Arguments are transformed using the following function [lepton parsed function and derivative]:\n");
    log.printf("   %s   [   %s   |   %s   ] \n",transf_str.c_str(),transf_value_parsed.c_str(),transf_deriv_parsed.c_str());
  }
  else {
    log.printf("  Arguments are not transformed\n");
  }
  //

  parseFlag("CHECK_NAN_INF",check_nan_inf_); addKeywordToList("CHECK_NAN_INF",check_nan_inf_);
  if(check_nan_inf_) {
    log.printf("  The code will check that values given are numercially stable, e.g. do not result in a not a number (nan) or infinity (inf).\n");
  }
  else {
    log.printf("  The code will NOT check that values given are numercially stable, e.g. do not result in a not a number (nan) or infinity (inf).\n");
  }

  setupBF();
  checkRead();
}


void BF_Custom::getAllValues(const double arg, double& argT, bool& inside_range, std::vector<double>& values, std::vector<double>& derivs) const {
  inside_range=true;
  argT=checkIfArgumentInsideInterval(arg,inside_range);
  double transf_derivf=1.0;
  //
  if(do_transf_) {

    if(transf_value_lepton_ref_) {*transf_value_lepton_ref_ = argT;}
    if(transf_deriv_lepton_ref_) {*transf_deriv_lepton_ref_ = argT;}

    argT = transf_value_expression_.evaluate();
    transf_derivf = transf_deriv_expression_.evaluate();

    if(check_nan_inf_ && (std::isnan(argT) || std::isinf(argT)) ) {
      std::string vs; Tools::convert(argT,vs);
      plumed_merror(getName()+" with label "+getLabel()+": problem with the transform function, it gives " + vs);
    }

    if(check_nan_inf_ && (std::isnan(transf_derivf) || std::isinf(transf_derivf)) ) {
      std::string vs; Tools::convert(transf_derivf,vs);
      plumed_merror(getName()+" with label "+getLabel()+": problem with the transform function, its derivative gives " + vs);
    }
  }
  //
  values[0]=1.0;
  derivs[0]=0.0;
  for(unsigned int i=1; i < getNumberOfBasisFunctions(); i++) {

    if(bf_values_lepton_ref_[i]) {*bf_values_lepton_ref_[i] = argT;}
    if(bf_derivs_lepton_ref_[i]) {*bf_derivs_lepton_ref_[i] = argT;}

    values[i] = bf_values_expressions_[i].evaluate();
    derivs[i] = bf_derivs_expressions_[i].evaluate();

    if(do_transf_) {derivs[i]*=transf_derivf;}
    // NaN checks
    if(check_nan_inf_ && (std::isnan(values[i]) || std::isinf(values[i])) ) {
      std::string vs; Tools::convert(values[i],vs);
      std::string is; Tools::convert(i,is);
      plumed_merror(getName()+" with label "+getLabel()+": problem with the basis function given in FUNC"+is+", it gives "+vs);
    }
    //
    if(check_nan_inf_ && (std::isnan(derivs[i])|| std::isinf(derivs[i])) ) {
      std::string vs; Tools::convert(derivs[i],vs);
      std::string is; Tools::convert(i,is);
      plumed_merror(getName()+" with label "+getLabel()+": problem with derivative of the basis function given in FUNC"+is+", it gives "+vs);
    }
  }
  if(!inside_range) {for(unsigned int i=0; i<derivs.size(); i++) {derivs[i]=0.0;}}
}




}
}
