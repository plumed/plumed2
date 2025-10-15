/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#include "ActionForInterface.h"
#include "ActionWithVector.h"
#include "ActionWithVirtualAtom.h"
#include "ActionShortcut.h"
#include "tools/PDB.h"
#include "PlumedMain.h"
#include "ActionSet.h"
#include <iostream>
#include <regex>

namespace PLMD {

void ActionWithArguments::registerKeywords(Keywords& keys) {
//  keys.reserve("numbered","ARG","the input for this action is the scalar output from one or more other actions. The particular scalars that you will use "
//               "are referenced using the label of the action. If the label appears on its own then it is assumed that the Action calculates "
//               "a single scalar value.  The value of this scalar is thus used as the input to this new action.  If * or *.* appears the "
//               "scalars calculated by all the proceeding actions in the input file are taken.  Some actions have multi-component outputs and "
//               "each component of the output has a specific label.  For example a \\ref DISTANCE action labelled dist may have three components "
//               "x, y and z.  To take just the x component you should use dist.x, if you wish to take all three components then use dist.*."
//               "More information on the referencing of Actions can be found in the section of the manual on the PLUMED \\ref Syntax.  "
//               "Scalar values can also be "
//               "referenced using POSIX regular expressions as detailed in the section on \\ref Regex. To use this feature you you must compile "
//               "PLUMED with the appropriate flag.");
}

void ActionWithArguments::parseArgumentList(const std::string&key,std::vector<Value*>&arg) {
  if( keywords.getArgumentType(key).length()==0 ) {
    warning("keyword " + key + " for reading arguments is registered using Keyword::add rather than Keyword::addInputKeyword.  The keyword will thus not appear in the correct place in the manual");
  }
  std::string def;
  std::vector<std::string> c;
  arg.clear();
  parseVector(key,c);
  if( c.size()==0 && (keywords.style(key,"compulsory") || keywords.style(key,"hidden")) ) {
    if( keywords.getDefaultValue(key,def) ) {
      c.push_back( def );
    } else {
      return;
    }
  }
  interpretArgumentList(c,plumed.getActionSet(),this,arg);
}

bool ActionWithArguments::parseArgumentList(const std::string&key,int i,std::vector<Value*>&arg) {
  if( keywords.getArgumentType(key).length()==0 ) {
    warning("keyword " + key + " for reading argument is registered using Keyword::add rather than Keyword::addInputKeyword.  The keyword will thus not appear in the correct place in the manual");
  }
  std::vector<std::string> c;
  arg.clear();
  if(parseNumberedVector(key,i,c)) {
    interpretArgumentList(c,plumed.getActionSet(),this,arg);
    return true;
  } else {
    return false;
  }
}

void ActionWithArguments::interpretArgumentList(const std::vector<std::string>& c, const ActionSet& as, Action* readact, std::vector<Value*>&arg) {
  for(unsigned i=0; i<c.size(); i++) {
    // is a regex? then just interpret it. The signal is ()
    if(!c[i].compare(0,1,"(")) {
      unsigned l=c[i].length();
      if(!c[i].compare(l-1,1,")")) {
        // start regex parsing
        bool found_something=false;
        // take the string enclosed in quotes and put in round brackets
        std::string myregex=c[i];
        std::vector<ActionWithValue*> all=as.select<ActionWithValue*>();
        if( all.empty() ) {
          readact->error("your input file is not telling plumed to calculate anything");
        }

        try {
          std::regex txt_regex(myregex,std::regex::extended);
          plumed_massert(txt_regex.mark_count()==1,"I can parse with only one subexpression");
          for(unsigned j=0; j<all.size(); j++) {
            std::vector<std::string> ss=all[j]->getComponentsVector();
            for(unsigned  k=0; k<ss.size(); ++k) {
              if(std::regex_match(ss[k],txt_regex)) {
                arg.push_back(all[j]->copyOutput(ss[k]));
                found_something=true;
              }
            }
          }
        } catch(std::regex_error & e) {
          plumed_error()<<"Error parsing regular expression: "<<e.what();
        }
        if(!found_something) {
          plumed_error()<<"There isn't any action matching your regex " << myregex;
        }
      } else {
        plumed_merror("did you want to use regexp to input arguments? enclose it between two round braces (...) with no spaces!");
      }
    } else {
      std::size_t dot=c[i].find_first_of('.');
      std::string a=c[i].substr(0,dot);
      std::string name=c[i].substr(dot+1);
      if(c[i].find(".")!=std::string::npos) {   // if it contains a dot:
        if(a=="*" && name=="*") {
          // Take all values from all actions
          std::vector<ActionWithValue*> all=as.select<ActionWithValue*>();
          if( all.empty() ) {
            readact->error("your input file is not telling plumed to calculate anything");
          }
          for(unsigned j=0; j<all.size(); j++) {
            plumed_assert(all[j]); // needed for following calls, see #1046
            ActionForInterface* ap=all[j]->castToActionForInterface();
            if( ap ) {
              continue;
            }
            for(unsigned k=0; k<all[j]->getNumberOfComponents(); ++k) {
              arg.push_back(all[j]->copyOutput(k));
            }
          }
        } else if ( name=="*") {
          unsigned carg=arg.size();
          // Take all the values from an action with a specific name
          ActionShortcut* shortcut=as.getShortcutActionWithLabel(a);
          if( shortcut ) {
            shortcut->interpretDataLabel( a + "." + name, readact, arg );
          }
          if( arg.size()==carg ) {
            // Take all the values from an action with a specific name
            ActionWithValue* action=as.selectWithLabel<ActionWithValue*>(a);
            if(!action) {
              std::string str=" (hint! the actions with value in this ActionSet are: ";
              str+=as.getLabelList<ActionWithValue*>()+")";
              readact->error("cannot find action named " + a + str);
            }
            if( action->getNumberOfComponents()==0 ) {
              readact->error("found " + a +".* indicating use all components calculated by action with label " + a + " but this action has no components");
            }
            for(unsigned k=0; k<action->getNumberOfComponents(); ++k) {
              arg.push_back(action->copyOutput(k));
            }
          }
        } else if ( a=="*" ) {
          std::vector<ActionShortcut*> shortcuts=as.select<ActionShortcut*>();
          // Take components from all actions with a specific name
          std::vector<ActionWithValue*> all=as.select<ActionWithValue*>();
          if( all.empty() ) {
            readact->error("your input file is not telling plumed to calculate anything");
          }
          unsigned carg=arg.size();
          for(unsigned j=0; j<shortcuts.size(); ++j) {
            shortcuts[j]->interpretDataLabel( shortcuts[j]->getShortcutLabel() + "." + name, readact, arg );
          }
          unsigned nval=0;
          for(unsigned j=0; j<all.size(); j++) {
            std::string flab;
            flab=all[j]->getLabel() + "." + name;
            if( all[j]->exists(flab) ) {
              arg.push_back(all[j]->copyOutput(flab));
              nval++;
            }
          }
          if(nval==0 && arg.size()==carg) {
            readact->error("found no actions with a component called " + name );
          }
        } else {
          // Take values with a specific name
          ActionWithValue* action=as.selectWithLabel<ActionWithValue*>(a);
          ActionShortcut* shortcut=as.getShortcutActionWithLabel(a);
          if( !shortcut && !action ) {
            std::string str=" (hint! the actions with value in this ActionSet are: ";
            str+=as.getLabelList<ActionWithValue*>()+")";
            readact->error("cannot find action named " + a +str);
          } else if( action && action->exists(c[i]) ) {
            arg.push_back(action->copyOutput(c[i]));
          } else if( shortcut ) {
            unsigned narg=arg.size();
            shortcut->interpretDataLabel( a + "." + name, readact, arg );
            if( arg.size()==narg ) {
              readact->error("found no element in " + a + " with label " + name );
            }
          } else {
            std::string str=" (hint! the components in this actions are: ";
            str+=action->getComponentsList()+")";
            readact->error("action " + a + " has no component named " + name + str);
          }
        }
      } else {    // if it doesn't contain a dot
        if(c[i]=="*") {
          // Take all values from all actions
          std::vector<ActionWithValue*> all=as.select<ActionWithValue*>();
          if( all.empty() ) {
            readact->error("your input file is not telling plumed to calculate anything");
          }
          for(unsigned j=0; j<all.size(); j++) {
            plumed_assert(all[j]); // needed for following calls, see #1046
            ActionWithVirtualAtom* av=all[j]->castToActionWithVirtualAtom();
            if( av ) {
              continue;
            }
            ActionForInterface* ap=all[j]->castToActionForInterface();
            if( ap && all[j]->getName()!="ENERGY" ) {
              continue;
            }
            for(unsigned k=0; k<all[j]->getNumberOfComponents(); ++k) {
              arg.push_back(all[j]->copyOutput(k));
            }
          }
        } else {
          ActionWithValue* action=as.selectWithLabel<ActionWithValue*>(c[i]);
          if(!action) {
            std::string str=" (hint! the actions with value in this ActionSet are: ";
            str+=as.getLabelList<ActionWithValue*>()+")";
            readact->error("cannot find action named " + c[i] + str );
          }
          if( !(action->exists(c[i])) ) {
            std::string str=" (hint! the components in this actions are: ";
            str+=action->getComponentsList()+")";
            readact->error("action " + c[i] + " has no component named " + c[i] +str);
          };
          arg.push_back(action->copyOutput(c[i]));
        }
      }
    }
  }
  for(unsigned i=0; i<arg.size(); ++i) {
    if( !readact->keywords.checkArgumentType( arg[i]->getRank(), arg[i]->hasDerivatives() ) ) {
      readact->error("documentation for input type is not provided in " + readact->getName() );
    }
  }
  if( readact->keywords.exists("MASKED_INPUT_ALLOWED") || readact->keywords.exists("IS_SHORTCUT") || readact->keywords.exists("MASK") ) {
    return;
  }
  for(unsigned i=0; i<arg.size(); ++i) {
    if( arg[i]->getRank()==0 ) {
      continue;
    }
    ActionWithVector* av=dynamic_cast<ActionWithVector*>( arg[i]->getPntrToAction() );
    if( av && av->getNumberOfMasks()>=0 ) {
      readact->error("cannot use argument " + arg[i]->getName() + " in input as not all elements are computed");
    }
  }
}

void ActionWithArguments::expandArgKeywordInPDB( const PDB& pdb ) {
  std::vector<std::string> arg_names = pdb.getArgumentNames();
  if( arg_names.size()>0 ) {
    std::vector<Value*> arg_vals;
    interpretArgumentList( arg_names, plumed.getActionSet(), this, arg_vals );
  }
}

void ActionWithArguments::requestArguments(const std::vector<Value*> &arg) {
  plumed_massert(!lockRequestArguments,"requested argument list can only be changed in the prepare() method");
  arguments=arg;
  clearDependencies();
  std::string fullname;
  std::string argName;
  for(unsigned i=0; i<arguments.size(); i++) {
    fullname=arguments[i]->getName();
    if(fullname.find(".")!=std::string::npos) {
      std::size_t dot=fullname.find_first_of('.');
      argName=fullname.substr(0,dot);
    } else {
      argName=fullname;
    }
    ActionWithValue* action=plumed.getActionSet().selectWithLabel<ActionWithValue*>(argName);
    plumed_massert(action,"cannot find action named (in requestArguments - this is weird)" + argName);
    addDependency(action);
  }
  ActionWithValue* av=dynamic_cast<ActionWithValue*>(this);
  if(av) {
    av->firststep=true;
  }
}

void ActionWithArguments::requestExtraDependencies(const std::vector<Value*> &extra) {
  plumed_massert(!lockRequestArguments,"requested argument list can only be changed in the prepare() method");
  std::string fullname;
  std::string argName;
  for(unsigned i=0; i<extra.size(); i++) {
    fullname=extra[i]->getName();
    if(fullname.find(".")!=std::string::npos) {
      std::size_t dot=fullname.find_first_of('.');
      argName=fullname.substr(0,dot);
    } else {
      argName=fullname;
    }
    ActionWithValue* action=plumed.getActionSet().selectWithLabel<ActionWithValue*>(argName);
    plumed_massert(action,"cannot find action named (in requestArguments - this is weird)" + argName);
    addDependency(action);
  }
}

ActionWithArguments::ActionWithArguments(const ActionOptions&ao):
  Action(ao),
  lockRequestArguments(false) {
  if( keywords.exists("ARG") ) {
    std::vector<Value*> arg;
    parseArgumentList("ARG",arg);

    if(!arg.empty()) {
      log.printf("  with arguments : \n");
      for(unsigned i=0; i<arg.size(); i++) {
        if( arg[i]->hasDerivatives() && arg[i]->getRank()>0 ) {
          log.printf(" function on grid with label %s \n",arg[i]->getName().c_str());
        } else if( arg[i]->getRank()==2 ) {
          log.printf("   matrix with label %s \n",arg[i]->getName().c_str());
        } else if( arg[i]->getRank()==1 ) {
          log.printf("   vector with label %s \n",arg[i]->getName().c_str());
        } else if( arg[i]->getRank()==0 ) {
          log.printf("   scalar with label %s \n",arg[i]->getName().c_str());
        } else {
          error("type of argument does not make sense");
        }
      }
    }
    requestArguments(arg);
  }
}

void ActionWithArguments::calculateNumericalDerivatives( ActionWithValue* a ) {
  if(!a) {
    a=castToActionWithValue();
    plumed_massert(a,"cannot compute numerical derivatives for an action without values");
  }

  const size_t nval=a->getNumberOfComponents();
  const size_t npar=arguments.size();
  std::vector<double> value (nval*npar);
  for(unsigned i=0; i<npar; i++) {
    double arg0=arguments[i]->get();
    arguments[i]->set(arg0+std::sqrt(epsilon));
    a->calculate();
    arguments[i]->set(arg0);
    for(unsigned j=0; j<nval; j++) {
      value[i*nval+j]=a->getOutputQuantity(j);
    }
  }
  a->calculate();
  a->clearDerivatives();
  for(unsigned j=0; j<nval; j++) {
    Value* v=a->copyOutput(j);
    if( v->hasDerivatives() )
      for(unsigned i=0; i<npar; i++) {
        v->addDerivative(i,(value[i*nval+j]-a->getOutputQuantity(j))/std::sqrt(epsilon));
      }
  }
}

double ActionWithArguments::getProjection(unsigned i,unsigned j)const {
  plumed_massert(i<arguments.size()," making projections with an index which  is too large");
  plumed_massert(j<arguments.size()," making projections with an index which  is too large");
  const Value* v1=arguments[i];
  const Value* v2=arguments[j];
  return Value::projection(*v1,*v2);
}

void ActionWithArguments::addForcesOnArguments( const unsigned& argstart, const std::vector<double>& forces, unsigned& ind ) {
  unsigned nargs=arguments.size();
  const ActionWithVector* av=dynamic_cast<const ActionWithVector*>( this );
  if( av && av->getNumberOfMasks()>0 ) {
    nargs=nargs-av->getNumberOfMasks();
  }
  for(unsigned i=0; i<nargs; ++i) {
    ind += arguments[i]->addForces(View(&forces[ind],forces.size()-ind));
  }
}

void ActionWithArguments::setGradients( Value* myval, unsigned& start ) const {
  if( !myval->hasDeriv ) {
    return;
  }
  plumed_assert( myval->getRank()==0 );

  bool scalar=true;
  for(unsigned i=0; i<arguments.size(); ++i ) {
    if( arguments[i]->getRank()!=0 ) {
      scalar=false;
      break;
    }
  }
  if( !scalar ) {
    bool constant=true;
    for(unsigned i=0; i<arguments.size(); ++i ) {
      if( !arguments[i]->isConstant() ) {
        constant=false;
        break;
      } else {
        start += arguments[i]->getNumberOfValues();
      }
    }
    if( !constant ) {
      error("cannot set gradient as unable to handle non-constant actions that take vectors/matrices/grids in input");
    }
  }
  // Now pass the gradients
  for(unsigned i=0; i<arguments.size(); ++i ) {
    arguments[i]->passGradients( myval->getDerivative(i), myval->gradients );
  }
}

bool ActionWithArguments::calculateConstantValues( const bool& haveatoms ) {
  ActionWithValue* awval = castToActionWithValue();
  if( !awval || arguments.size()==0 ) {
    return false;
  }
  bool constant = true, atoms=false;
  for(unsigned i=0; i<arguments.size(); ++i) {
    auto * ptr=arguments[i]->getPntrToAction();
    plumed_assert(ptr); // needed for following calls, see #1046
    ActionAtomistic* aa=ptr->castToActionAtomistic();
    if( aa ) {
      ActionWithVector* awvec=dynamic_cast<ActionWithVector*>( arguments[i]->getPntrToAction() );
      if( !awvec || aa->getNumberOfAtoms()>0 ) {
        atoms=true;
      }
    }
    if( !arguments[i]->isConstant() ) {
      constant=false;
      break;
    }
  }
  if( constant ) {
    // Set everything constant first as we need to set the shape
    for(unsigned i=0; i<awval->getNumberOfComponents(); ++i) {
      (awval->copyOutput(i))->setConstant();
    }
    if( !haveatoms ) {
      log.printf("  values stored by this action are computed during startup and stay fixed during the simulation\n");
    }
    if( atoms ) {
      return haveatoms;
    }
  }
  // Now do the calculation and store the values if we don't need anything from the atoms
  if( constant && !haveatoms ) {
    plumed_assert( !atoms );
    activate();
    calculate();
    deactivate();
    for(unsigned i=0; i<awval->getNumberOfComponents(); ++i) {
      unsigned nv = awval->copyOutput(i)->getNumberOfValues();
      log.printf("  %d values stored in component labelled %s are : ", nv, (awval->copyOutput(i))->getName().c_str() );
      for(unsigned j=0; j<nv; ++j) {
        log.printf(" %f", (awval->copyOutput(i))->get(j) );
      }
      log.printf("\n");
    }
  }
  return constant;
}

}
