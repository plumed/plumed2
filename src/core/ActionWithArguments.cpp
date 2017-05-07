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
#include "tools/PDB.h"
#include "PlumedMain.h"
#include "ActionSet.h"
#include <iostream>
#ifdef __PLUMED_HAS_CREGEX
#include <cstring>
#include "regex.h"
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
    std::size_t found1 = c[i].find("(");
    if(found1!=std::string::npos) {
      std::size_t found2=c[i].find(")",found1+1,1); // find it again
      if(found2!=std::string::npos) {
        // start regex parsing
#ifdef __PLUMED_HAS_CREGEX
        // take the string enclosed in quotes and put in round brackets
        std::string myregex=c[i].substr(found1,found2-found1+1);
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
                    arg.push_back(all[j]->copyOutput(putativeVal));
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
        if(a=="*" && name=="*") {
          // Take all values from all actions
          std::vector<ActionWithValue*> all=plumed.getActionSet().select<ActionWithValue*>();
          if( all.empty() ) error("your input file is not telling plumed to calculate anything");
          for(unsigned j=0; j<all.size(); j++) {
            for(int k=0; k<all[j]->getNumberOfComponents(); ++k) arg.push_back(all[j]->copyOutput(k));
          }
        } else if ( name=="*") {
          // Take all the values from an action with a specific name
          ActionWithValue* action=plumed.getActionSet().selectWithLabel<ActionWithValue*>(a);
          if(!action) {
            std::string str=" (hint! the actions in this ActionSet are: ";
            str+=plumed.getActionSet().getLabelList()+")";
            error("cannot find action named " + a + str);
          }
          if( action->getNumberOfComponents()==0 ) error("found " + a +".* indicating use all components calculated by action with label " + a + " but this action has no components");
          for(int k=0; k<action->getNumberOfComponents(); ++k) arg.push_back(action->copyOutput(k));
        } else if ( a=="*" ) {
          // Take components from all actions with a specific name
          std::vector<ActionWithValue*> all=plumed.getActionSet().select<ActionWithValue*>();
          if( all.empty() ) error("your input file is not telling plumed to calculate anything");
          unsigned nval=0;
          for(unsigned j=0; j<all.size(); j++) {
            std::string flab; flab=all[j]->getLabel() + "." + name;
            if( all[j]->exists(flab) ) { arg.push_back(all[j]->copyOutput(flab)); nval++; }
          }
          if(nval==0) error("found no actions with a component called " + name );
        } else {
          // Take values with a specific name
          ActionWithValue* action=plumed.getActionSet().selectWithLabel<ActionWithValue*>(a);
          if(!action) {
            std::string str=" (hint! the actions in this ActionSet are: ";
            str+=plumed.getActionSet().getLabelList()+")";
            error("cannot find action named " + a +str);
          }
          if( !(action->exists(c[i])) ) {
            std::string str=" (hint! the components in this actions are: ";
            str+=action->getComponentsList()+")";
            error("action " + a + " has no component named " + name + str);
          } ;
          arg.push_back(action->copyOutput(c[i]));
        }
      } else {    // if it doesn't contain a dot
        if(c[i]=="*") {
          // Take all values from all actions
          std::vector<ActionWithValue*> all=plumed.getActionSet().select<ActionWithValue*>();
          if( all.empty() ) error("your input file is not telling plumed to calculate anything");
          for(unsigned j=0; j<all.size(); j++) {
            for(int k=0; k<all[j]->getNumberOfComponents(); ++k) arg.push_back(all[j]->copyOutput(k));
          }
        } else {
          ActionWithValue* action=plumed.getActionSet().selectWithLabel<ActionWithValue*>(c[i]);
          if(!action) {
            std::string str=" (hint! the actions in this ActionSet are: ";
            str+=plumed.getActionSet().getLabelList()+")";
            error("cannot find action named " + c[i] + str );
          }
          if( !(action->exists(c[i])) ) {
            std::string str=" (hint! the components in this actions are: ";
            str+=action->getComponentsList()+")";
            error("action " + c[i] + " has no component named " + c[i] +str);
          };
          arg.push_back(action->copyOutput(c[i]));
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

void ActionWithArguments::requestArguments(const vector<Value*> &arg) {
  plumed_massert(!lockRequestArguments,"requested argument list can only be changed in the prepare() method");
  arguments=arg;
  clearDependencies();
  std::string fullname,name;
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
  }
}

ActionWithArguments::ActionWithArguments(const ActionOptions&ao):
  Action(ao),
  lockRequestArguments(false)
{
  if( keywords.exists("ARG") && !keywords.exists("DATA") ) {
    vector<Value*> arg;
    parseArgumentList("ARG",arg);

    if(!arg.empty()) {
      log.printf("  with arguments");
      for(unsigned i=0; i<arg.size(); i++) log.printf(" %s",arg[i]->getName().c_str());
      log.printf("\n");
    }
    requestArguments(arg);
  }
}

void ActionWithArguments::calculateNumericalDerivatives( ActionWithValue* a ) {
  if(!a) {
    a=dynamic_cast<ActionWithValue*>(this);
    plumed_massert(a,"cannot compute numerical derivatives for an action without values");
  }

  const int nval=a->getNumberOfComponents();
  const int npar=arguments.size();
  std::vector<double> value (nval*npar);
  for(int i=0; i<npar; i++) {
    double arg0=arguments[i]->get();
    arguments[i]->set(arg0+sqrt(epsilon));
    a->calculate();
    arguments[i]->set(arg0);
    for(int j=0; j<nval; j++) {
      value[i*nval+j]=a->getOutputQuantity(j);
    }
  }
  a->calculate();
  a->clearDerivatives();
  for(int j=0; j<nval; j++) {
    Value* v=a->copyOutput(j);
    if( v->hasDerivatives() ) for(int i=0; i<npar; i++) v->addDerivative(i,(value[i*nval+j]-a->getOutputQuantity(j))/sqrt(epsilon));
  }
}

double ActionWithArguments::getProjection(unsigned i,unsigned j)const {
  plumed_massert(i<arguments.size()," making projections with an index which  is too large");
  plumed_massert(j<arguments.size()," making projections with an index which  is too large");
  const Value* v1=arguments[i];
  const Value* v2=arguments[j];
  return Value::projection(*v1,*v2);
}

void ActionWithArguments::addForcesOnArguments( const std::vector<double>& forces ) {
  for(unsigned i=0; i<arguments.size(); ++i) arguments[i]->addForce( forces[i] );
}

}
