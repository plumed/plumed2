/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2017 The plumed team
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
#include "ActionWithArguments.h"
#include "ActionWithValue.h"
#include "ActionAtomistic.h"
#include "tools/PDB.h"
#include "PlumedMain.h"
#include "ActionSet.h"
#include <iostream>
#ifdef __PLUMED_HAS_CREGEX
#include <cstring>
#include <regex.h>
#endif

using namespace std;
namespace PLMD {

void ActionWithArguments::registerKeywords(Keywords& keys) {
  keys.reserve("numbered","ARG","the input for this action is the scalar output from one or more other actions. The particular scalars that you will use "
               "are referenced using the label of the action. If the label appears on its own then it is assumed that the Action calculates "
               "a single scalar value.  The value of this scalar is thus used as the input to this new action.  If * or *.* appears the "
               "scalars calculated by all the proceding actions in the input file are taken.  Some actions have multi-component outputs and "
               "each component of the output has a specific label.  For example a \\ref DISTANCE action labelled dist may have three componets "
               "x, y and z.  To take just the x component you should use dist.x, if you wish to take all three components then use dist.*."
               "More information on the referencing of Actions can be found in the section of the manual on the PLUMED \\ref Syntax.  "
               "Scalar values can also be "
               "referenced using POSIX regular expressions as detailed in the section on \\ref Regex. To use this feature you you must compile "
               "PLUMED with the appropriate flag.");
}

void ActionWithArguments::parseArgumentList(const std::string&key,std::vector<Value*>&arg) {
  vector<string> c; arg.clear(); parseVector(key,c);
  if( c.size()==0 && (keywords.style(key,"compulsory") || keywords.style(key,"hidden")) ) {
    std::string def; if( keywords.getDefaultValue(key,def) ) c.push_back( def );
  }
  interpretArgumentList(c,arg);
}

bool ActionWithArguments::parseArgumentList(const std::string&key,int i,std::vector<Value*>&arg) {
  vector<string> c;
  arg.clear();
  if(parseNumberedVector(key,i,c)) {
    interpretArgumentList(c,arg);
    return true;
  } else return false;
}

void ActionWithArguments::interpretArgumentList(const std::vector<std::string>& c, std::vector<Value*>&arg) {
  for(unsigned i=0; i<c.size(); i++) {
    // is a regex? then just interpret it. The signal is ()
    if(!c[i].compare(0,1,"(")) {
      unsigned l=c[i].length();
      if(!c[i].compare(l-1,1,")")) {
        // start regex parsing
#ifdef __PLUMED_HAS_CREGEX
        // take the string enclosed in quotes and put in round brackets
        std::string myregex=c[i];
        log.printf("  Evaluating regexp for this action: %s \n",myregex.c_str());
        int errcode;
        regex_t *preg = (regex_t*)malloc(sizeof(regex_t)); // pointer to the regular expression
        regmatch_t *pmatch;
        if ((errcode=regcomp(preg, myregex.c_str(),REG_EXTENDED|REG_NEWLINE))) {  // compile the regular expression
          char* errbuf;
          size_t errbuf_size;
          // one can check the errors asking to regerror
          errbuf_size = regerror(errcode, preg, NULL, 0);
          if (!(errbuf=(char*)malloc(errbuf_size))) {
            plumed_merror("cannot allocate the buffer for error detection in regexp!");
          };
          regerror(errcode, preg, errbuf, errbuf_size);
          error(errbuf);
        }
        plumed_massert(preg->re_nsub==1,"I can parse with only one subexpression");
        pmatch = (regmatch_t*)malloc(sizeof(regmatch_t)*preg->re_nsub);
        // select all the actions that have a value
        std::vector<ActionWithValue*> all=plumed.getActionSet().select<ActionWithValue*>();
        if( all.empty() ) error("your input file is not telling plumed to calculate anything");
        for(unsigned j=0; j<all.size(); j++) {
          std::vector<std::string> ss=all[j]->getComponentsVector();
          for(unsigned  k=0; k<ss.size(); ++k) {
            unsigned ll=strlen(ss[k].c_str())+1;
            std::vector<char> str(ll);
            strcpy(&str[0],ss[k].c_str());
            const char *ppstr=&str[0];
            if(!regexec(preg, ppstr, preg->re_nsub, pmatch, 0)) {
              log.printf("  Something matched with \"%s\" : ",ss[k].c_str());
              do {
                if (pmatch[0].rm_so != -1) {	/* The regex is matching part of a string */
                  char *submatch;
                  size_t matchlen = pmatch[0].rm_eo - pmatch[0].rm_so;
                  submatch = (char*)malloc(matchlen+1);
                  strncpy(submatch, ppstr+pmatch[0].rm_so, matchlen+1);
                  submatch[matchlen]='\0';
                  log.printf("  subpattern %s\n", submatch);
                  // this is the match: try to see if it is a valid action
                  std::string putativeVal(submatch);
                  if( all[j]->exists(putativeVal) ) {
                     all[j]->interpretDataLabel( putativeVal, this, arg );
                     log.printf("  Action %s added! \n",putativeVal.c_str());
                  }
                  free(submatch);
                };
                ppstr += pmatch[0].rm_eo;	/* Restart from last match */
              } while(!regexec(preg,ppstr,preg->re_nsub,pmatch,0));
            }
          }
        };
        regfree(preg);
        free(preg);
        free(pmatch);
#else
        plumed_merror("Regexp support not compiled!");
#endif
      } else {
        plumed_merror("did you want to use regexp to input arguments? enclose it between two round braces (...) with no spaces!");
      }
    } else {
      std::size_t dot=c[i].find_first_of('.');
      string a=c[i].substr(0,dot);
      string name=c[i].substr(dot+1);
      if(c[i].find(".")!=string::npos) {   // if it contains a dot:
        if(a=="*") {
          // Take all values from all actions
          std::vector<ActionWithValue*> all=plumed.getActionSet().select<ActionWithValue*>();
          if( all.empty() ) error("your input file is not telling plumed to calculate anything");
          unsigned carg = arg.size();
          for(unsigned j=0; j<all.size(); j++){
              if( name=="*" || all[j]->exists(all[j]->getLabel() + "." + name) ) all[j]->interpretDataLabel( all[j]->getLabel() + "." + name, this, arg );
          }
          if( arg.size()==carg ) error("found no actions with a component called " + name );
        } else {   
          // Take all the values from an action with a specific name
          ActionWithValue* action=plumed.getActionSet().selectWithLabel<ActionWithValue*>(a);
          if( !action ){ 
            std::string str=" (hint! the actions in this ActionSet are: ";
            str+=plumed.getActionSet().getLabelList()+")";
            error("cannot find action named " + a +str);
          }
          unsigned carg = arg.size(); action->interpretDataLabel( c[i], this, arg );
          // if( arg.size()==carg && name=="*" ) error("found " + a +".* indicating use all components calculated by action with label " + a + " but this action has no components");
          if( arg.size()==carg ) {
             std::string str=" (hint! the components in this actions are: ";
             str+=action->getComponentsList()+")";
             error("action " + a + " has no component named " + name + str);
          }
        }
      } else {    // if it doesn't contain a dot
        if(c[i]=="*") {
          // Take all values from all actions
          std::vector<ActionWithValue*> all=plumed.getActionSet().select<ActionWithValue*>();
          if( all.empty() ) error("your input file is not telling plumed to calculate anything");
          unsigned carg = arg.size();
          for(unsigned j=0; j<all.size(); j++) all[j]->interpretDataLabel( all[j]->getLabel() + ".*", this, arg );
          if( arg.size()==carg ) error("found no actions with a component called " + name ); 
        } else {
          ActionWithValue* action=plumed.getActionSet().selectWithLabel<ActionWithValue*>(c[i]);
          if(!action) {
            std::string str=" (hint! the actions in this ActionSet are: ";
            str+=plumed.getActionSet().getLabelList()+")";
            error("cannot find action named " + c[i] + str );
          }
          unsigned carg=arg.size(); action->interpretDataLabel( "", this, arg );
          if( arg.size()!=carg+1 ){
             std::string str=" (hint! the components in this actions are: ";
             str+=action->getComponentsList()+")";
             error("action " + c[i] + " has no component named " + c[i] +str);
          };
        }
      }
    }
  }
}

void ActionWithArguments::expandArgKeywordInPDB( PDB& pdb ) {
  std::vector<std::string> pdb_remark=pdb.getRemark();
  std::vector<std::string> arg_names;
  bool found=Tools::parseVector(pdb_remark,"ARG",arg_names);
  if( found ) {
    std::vector<Value*> arg_vals;
    interpretArgumentList( arg_names, arg_vals );
    std::string new_args="ARG=" + arg_vals[0]->getName();
    for(unsigned i=1; i<arg_vals.size(); ++i) new_args = new_args + "," + arg_vals[i]->getName();
    pdb.setArgKeyword( new_args );
  }
}

void ActionWithArguments::requestArguments(const vector<Value*> &arg, const bool& allow_streams ) {
  plumed_massert(!lockRequestArguments,"requested argument list can only be changed in the prepare() method");
  bool firstcall=(arguments.size()==0);
  arguments=arg;
  clearDependencies();
  std::string fullname,name;
  std::vector<ActionWithValue*> f_actions;
  for(unsigned i=0; i<arguments.size(); i++) {
    fullname=arguments[i]->getName();
    if(fullname.find(".")!=string::npos) {
      std::size_t dot=fullname.find_first_of('.');
      name=fullname.substr(0,dot);
    } else {
      name=fullname;
    }
    ActionWithValue* action=plumed.getActionSet().selectWithLabel<ActionWithValue*>(name);
    plumed_massert(action,"cannot find action named (in requestArguments - this is weird)" + name);
    addDependency(action);
    // Check if we already have this argument in the stream
    if( arguments[i]->getRank()>0 ){
        bool found=false; ActionWithValue* myact = (arguments[i]->getPntrToAction())->getActionThatCalculates();
        for(unsigned k=0;k<f_actions.size();++k){
            if( f_actions[k]==myact ){ found=true; break; }
        }   
        if( !found ) f_actions.push_back( myact );
    }   
  }
  // This is a way of checking if we are in an ActionWithValue by looking at the keywords -- is there better fix?
  if( firstcall ) {
      if( !keywords.exists("SERIAL") ){
          for(unsigned i=0;i<arg.size();++i){ if( arg[i]->getRank()>0 ) arg[i]->buildDataStore(); }
          return;
      }
  } else {
      ActionWithValue* av = dynamic_cast<ActionWithValue*>(this);
      if(!av) return;
  }

  if( !allow_streams ){
      done_over_stream=false;
  } else if( f_actions.size()>1 ){
      done_over_stream=true;
      for(unsigned i=1;i<f_actions.size();++i){
          if( f_actions[0]->getFullNumberOfTasks()!=f_actions[i]->getFullNumberOfTasks() ){ done_over_stream=false; break; }
      }
      if( done_over_stream ){
          std::vector<std::string> empty(1); empty[0] = f_actions[0]->getLabel();
          for(unsigned i=1;i<f_actions.size();++i){ f_actions[0]->addActionToChain( empty, f_actions[i] ); f_actions[i]->addDependency( f_actions[0] ); }
      } else {
          for(unsigned i=0;i<arg.size();++i){ if( arg[i]->getRank()>0 ) arg[i]->buildDataStore(); }
      }
  } else if( f_actions.size()==1 ) done_over_stream=true; 

  if( done_over_stream ) {
      // Get the action where this argument should be applied 
      ActionWithArguments* aa=dynamic_cast<ActionWithArguments*>( arguments[0]->getPntrToAction() );
      if( !aa || aa->mustBeTreatedAsDistinctArguments() ) distinct_arguments.push_back( arguments[0]->getPntrToAction() );
      else distinct_arguments.push_back( aa->getFirstNonStream() );
      // Build vector with locations to keep derivatives of arguments 
      arg_deriv_starts.clear(); arg_deriv_starts.resize(0);
      arg_deriv_starts.push_back(0); unsigned nder = distinct_arguments[0]->getNumberOfDerivatives();
      
      for(unsigned i=1;i<getNumberOfArguments();++i){
          // Work out what action applies forces
          ActionWithValue* myval; ActionWithArguments* aa=dynamic_cast<ActionWithArguments*>( arguments[i]->getPntrToAction() );
          if( !aa || aa->mustBeTreatedAsDistinctArguments() ) myval = arguments[i]->getPntrToAction(); 
          else myval = aa->getFirstNonStream();
          
          // Check we haven't already dealt with this argument
          int argno=-1;;
          for(unsigned j=0;j<distinct_arguments.size();++j){
             if( myval==distinct_arguments[j] ){ argno=j; break; }
          }  
          if( argno>=0 ){
              arg_deriv_starts.push_back( arg_deriv_starts[argno] );
          } else {
              arg_deriv_starts.push_back( nder ); distinct_arguments.push_back( myval ); nder += myval->getNumberOfDerivatives();
          }   
      }
  } else {
      for(unsigned i=0;i<getNumberOfArguments();++i){ if( arg[i]->getRank()>0 ) arg[i]->buildDataStore(); } 
  }

}

ActionWithArguments::ActionWithArguments(const ActionOptions&ao):
  Action(ao),
  lockRequestArguments(false),
  done_over_stream(false),
  numberedkeys(false)
{
  if( keywords.exists("ARG") && !keywords.exists("DATA") ) {
    vector<Value*> arg;
    parseArgumentList("ARG",arg);

    if(!arg.empty()) {
      log.printf("  with arguments"); arg_ends.push_back(0); numberedkeys=false;
      for(unsigned i=0; i<arg.size(); i++){ log.printf(" %s",arg[i]->getName().c_str()); arg_ends.push_back(1+i); }
      log.printf("\n");
    } else {
      unsigned narg=0; arg_ends.push_back(0); numberedkeys=true;
      for(unsigned i=1;;++i){
         vector<Value*> argn; parseArgumentList("ARG",i,argn);
         if( argn.size()==0 ) break;
         unsigned nargt=0; log.printf("  %dth set of arguments",i);
         for(unsigned j=0;j<argn.size();++j){
             log.printf(" %s",argn[j]->getName().c_str());
             arg.push_back( argn[j] ); nargt += argn[j]->getNumberOfValues();
         }
         arg_ends.push_back( arg.size() ); log.printf("\n"); 
         if( i==1 ) narg = nargt;
         else if( narg!=nargt ) error("mismatch between number of arguments specified for different numbered ARG values");
      }
    }
    requestArguments(arg,true);
  }
}

ActionWithValue* ActionWithArguments::getFirstNonStream() {
  plumed_massert( getNumberOfArguments()<2, "cannot use functions with multiple arguments in this way " + getLabel() );
  ActionWithArguments* aa=dynamic_cast<ActionWithArguments*>( getPntrToArgument(0)->getPntrToAction() );
  if( !aa || aa->mustBeTreatedAsDistinctArguments() ) return getPntrToArgument(0)->getPntrToAction();
  else return aa->getFirstNonStream();  
}

void ActionWithArguments::createTasksFromArguments(){
  ActionWithValue* av = dynamic_cast<ActionWithValue*>(this); plumed_assert( av );
  unsigned ntasks=0; for(unsigned j=arg_ends[0];j<arg_ends[1];++j) ntasks += getPntrToArgument(j)->getNumberOfValues();
  for(unsigned i=1;i<arg_ends.size()-1;++i){
      unsigned nt = 0; for(unsigned j=arg_ends[i];j<arg_ends[i+1];++j) nt += getPntrToArgument(j)->getNumberOfValues();
      plumed_assert( nt==ntasks );
  } 
  for(unsigned i=0;i<ntasks;++i) av->addTaskToList( i );
}

void ActionWithArguments::calculateNumericalDerivatives( ActionWithValue* a ) {
  if( done_over_stream ) error("cannot use numerical derivatives if calculation is done over stream");
  if(!a) {
    a=dynamic_cast<ActionWithValue*>(this);
    plumed_massert(a,"cannot compute numerical derivatives for an action without values");
  }

  std::vector<Value*> myvals; a->retrieveAllScalarValuesInLoop( myvals );
  const int npar=arguments.size();
  std::vector<double> value (myvals.size()*npar);
  for(int i=0; i<npar; i++) {
    double arg0=arguments[i]->get();
    arguments[i]->set(arg0+sqrt(epsilon));
    a->calculate();
    arguments[i]->set(arg0);
    for(int j=0; j<myvals.size(); j++) {
      value[i*myvals.size()+j]=myvals[j]->get(); 
    }
  }
  a->calculate();
  a->clearDerivatives();
  for(int j=0; j<myvals.size(); j++) {
    if( myvals[j]->hasDerivatives() ) for(int i=0; i<npar; i++) myvals[j]->addDerivative(i,(value[i*myvals.size()+j]-a->getOutputQuantity(j))/sqrt(epsilon));
  }
}

double ActionWithArguments::getProjection(unsigned i,unsigned j)const {
  plumed_massert(i<arguments.size()," making projections with an index which  is too large");
  plumed_massert(j<arguments.size()," making projections with an index which  is too large");
  plumed_massert(arguments[i]->getRank()==0 && arguments[j]->getRank()==0,"cannot calculate projection for data stream input"); 
  const Value* v1=arguments[i];
  const Value* v2=arguments[j];
  return Value::projection(*v1,*v2);
}

void ActionWithArguments::retrieveArguments( const MultiValue& myvals, std::vector<double>& args ) const {
  if( done_over_stream ){
      for(unsigned i=0;i<args.size();++i) args[i]=myvals.get( arguments[i]->streampos );
      return;
  }
  if( arguments[0]->getRank()>0 ){
      for(unsigned i=0;i<arg_ends.size()-1;++i){
          unsigned nt=0, nn=0; 
          for(unsigned j=arg_ends[i];j<arg_ends[i+1];++j){
              nt += arguments[j]->getNumberOfValues();
              if( myvals.getTaskIndex()<nt ) break;
              nn += arguments[j]->getNumberOfValues();
          }
          args[i]=arguments[i]->get( myvals.getTaskIndex() - nn );
      }
  } else {
      for(unsigned i=0;i<arg_ends.size()-1;++i) args[i]=arguments[arg_ends[i] + myvals.getTaskIndex()]->get(); 
  }
}

void ActionWithArguments::setForcesOnArguments( const std::vector<double>& forces, unsigned& start ) {
  if( done_over_stream ){
      for(unsigned i=0;i<distinct_arguments.size();++i){
          ActionWithArguments* aarg = dynamic_cast<ActionWithArguments*>( distinct_arguments[i] );
          if( aarg ) aarg->setForcesOnArguments( forces, start ); 
          ActionAtomistic* aat = dynamic_cast<ActionAtomistic*>( distinct_arguments[i] );
          if( aat ){ aat->setForcesOnAtoms( forces, start ); start += distinct_arguments[i]->getNumberOfDerivatives(); } 
      } 
  } else {
      for(unsigned k=0;k<getNumberOfScalarArguments();++k) {
          unsigned nt = 0, nn = 0, j=0;
          for(unsigned i=0;i<arguments.size();++i){
             nt += arguments[i]->getNumberOfValues();
             if( k<nt ){ j=i; break ; }
             nn += arguments[i]->getNumberOfValues();
          }
          arguments[j]->addForce( k-nn, forces[start + k] );
      }
      start += getNumberOfScalarArguments();
  }
}

}
