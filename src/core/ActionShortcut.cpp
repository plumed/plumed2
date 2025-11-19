/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2018-2023 The plumed team
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
#include "ActionShortcut.h"
#include "PlumedMain.h"
#include "ActionWithValue.h"
#include "ActionRegister.h"
#include "ActionSet.h"

namespace PLMD {

void ActionShortcut::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  keys.add("hidden","IS_SHORTCUT","hidden keyword to tell if actions are shortcuts so that example generator can provide expansions of shortcuts");
  keys.add("hidden","HAS_VALUES","this is used in json output to determine those actions that have values");
}

void ActionShortcut::readShortcutKeywords( const Keywords& keys, std::map<std::string,std::string>& keymap ) {
  for (auto& keyname:keys.getKeys()) {
    std::string t;
    bool def;
    if( keys.getLogicalDefault(keyname,def) ) {
      bool found=false;
      parseFlag(keyname,found);
      if( found ) {
        keymap.insert(std::pair<std::string,std::string>(keyname,""));
      }
    } else if( keys.style( keyname, "optional") || keys.style( keyname, "compulsory") || keys.style( keyname, "deprecated") ) {
      parse(keyname,t);
      if( t.length()>0 ) {
        keymap.insert(std::pair<std::string,std::string>(keyname,t));
      } else if( keys.numbered( keyname ) ) {
        for(unsigned i=1;; ++i) {
          std::string istr;
          Tools::convert( i, istr );
          if( !parseNumbered(keyname,i,t) ) {
            break ;
          }
          keymap.insert(std::pair<std::string,std::string>(keyname + istr,t));
        }
      }
    } else {
      plumed_merror("shortcut keywords should be optional, compulsory or flags");
    }
  }
}

ActionShortcut::ActionShortcut(const ActionOptions&ao):
  Action(ao),
  shortcutlabel(actionLabel) {
  std::string s;
  Tools::convert(plumed.getActionSet().size(),s);
  if( shortcutlabel==("@" + s) ) {
    std::string t;
    Tools::convert(plumed.getActionSet().size(),t);
    shortcutlabel="@" + t;
  } else {
    actionLabel = ("@s" + s);
  }
}

void ActionShortcut::readInputLine( const std::string& input, bool saveline ) {
  std::vector<std::string> words=Tools::getWords(input);
  Tools::interpretLabel(words);
  // Check if this action name has been registered
  bool founds=false;
  bool found = keywords.isActionNeeded(words[0]);
  // Check if we are just calling something like SUM_VECTOR using just SUM.
  if( !found && words[0].find(getName())!=std::string::npos ) {
    found =keywords.isActionSuffixed(words[0],getName());
    founds=true;
  }
  if( found ) {
    if( !founds && saveline ) {
      addToSavedInputLines( input );
    }
    if( keywords.exists("RESTART") ) {
      if( restart ) {
        words.push_back("RESTART=YES");
      }
      if( !restart ) {
        words.push_back("RESTART=NO");
      }
    }
    plumed.readInputWords( words, false);
    if( !founds ) {
      ActionWithValue* av=NULL;
      for(auto pp=plumed.getActionSet().rbegin(); pp!=plumed.getActionSet().rend(); ++pp) {
        av = pp->get()->castToActionWithValue();
        if( !av ) {
          continue ;
        }
        if( std::find(savedOutputs.begin(), savedOutputs.end(), av->getLabel() )!=savedOutputs.end() ) {
          av=NULL;
        }
        break;
      }
      if( av ) {
        std::string av_label = av->getLabel();
        if( av_label == getShortcutLabel() && av->getNumberOfComponents()==1 ) {
          savedOutputs.push_back( av_label );
          plumed_massert( keywords.componentHasCorrectType(".#!value",
                          (av->copyOutput(0))->getRank(), (av->copyOutput(0))->hasDerivatives() ),
                          "documentation for type of value is incorrect");
        } else {
          for(const auto& cname: keywords.getOutputComponents()) {
            if( av_label == getShortcutLabel() + "_" + cname ) {
              savedOutputs.push_back( av_label );
              plumed_massert( keywords.componentHasCorrectType(cname,
                              (av->copyOutput(0))->getRank(),
                              (av->copyOutput(0))->hasDerivatives() ),
                              "documentation for type of component " + cname + " is incorrect");
            } else if( keywords.getOutputComponentFlag(cname)!="default" ) {
              std::string thisflag = keywords.getOutputComponentFlag(cname);
              if( keywords.numbered(thisflag) && av_label.find(getShortcutLabel() + "_" + cname)!=std::string::npos ) {
                savedOutputs.push_back( av_label );
                plumed_massert( keywords.componentHasCorrectType(cname,
                                (av->copyOutput(0))->getRank(),
                                (av->copyOutput(0))->hasDerivatives() ),
                                "documentation for type of component " + cname + " is incorrect");
              }
            }
          }
        }
      }
    } else {
      ActionWithValue* av = plumed.getActionSet()[plumed.getActionSet().size()-1]->castToActionWithValue();
      if( !av ) {
        error("shortcut is using suffix but action created is not ActionWithValue");
      }
      Keywords thiskeys;
      actionRegister().getKeywords( av->getName(), thiskeys );
      if( thiskeys.getDisplayName()!=getName() ) {
        error("mismatch between display name of hidden action " + thiskeys.getDisplayName() + " and shortcut that creates it " + getName() );
      }
    }
  } else {
    error("requirement for action " + words[0] + " should be registered in registerKeywords function for shortcut action using keys.useAction");
  }
}

void ActionShortcut::addCommentToShortcutOutput( const std::string& input ) {
  savedInputLines.push_back( input );
}

std::string ActionShortcut::getUpdateLimits() const {
  std::string f_input="";
  if( update_from!=std::numeric_limits<double>::max() ) {
    std::string ufrom;
    Tools::convert( update_from, ufrom );
    f_input += " UPDATE_FROM=" + ufrom;
  }
  if( update_until!=std::numeric_limits<double>::max() ) {
    std::string util;
    Tools::convert( update_until, util );
    f_input += " UPDATE_UNTIL=" + util;
  }
  return f_input;
}

void ActionShortcut::addToSavedInputLines( const std::string& line ) {
  std::vector<std::string> words = Tools::getWords(line);
  std::string actname;
  if( words[0].find_first_of(":")!=std::string::npos) {
    actname = words[1];
  } else {
    actname = words[0];
  }
  if( !actionRegister().check(actname) ) {
    error("found no action with name " + actname + " to create shortcut");
  }
  Keywords thiskeys;
  actionRegister().getKeywords( actname, thiskeys );
  std::vector<std::string> numberedkeys;
  for (auto& keyname: thiskeys.getKeys()) {
    if( thiskeys.numbered( keyname ) ) {
      numberedkeys.push_back( keyname );
    }
  }
  if( numberedkeys.size()>0 && actname!="CONCATENATE" ) {
    std::string reducedline;
    for(unsigned i=0; i<words.size(); ++i) {
      bool notnumbered=true;
      for(unsigned j=0; j<numberedkeys.size(); ++j) {
        if( words[i].find(numberedkeys[j])!=std::string::npos && words[i].substr(0,numberedkeys[j].length()+1)!=numberedkeys[j]+"=" ) {
          notnumbered=false;
          break;
        }
      }
      if( notnumbered || words[i]==actname ) {
        if( words[i].find(" ")!=std::string::npos) {
          std::size_t eq=words[i].find_first_of("=");
          reducedline += words[i].substr(0,eq) + "={" + words[i].substr(eq+1) + "} ";
        } else {
          reducedline += words[i] + " ";
        }
      }
    }
    std::vector<unsigned> ninstances( numberedkeys.size(), 0 );
    for(unsigned j=0; j<numberedkeys.size(); ++j) {
      for(unsigned i=1;; ++i) {
        std::string num, val;
        Tools::convert(i, num);
        bool found = Tools::parse(words, numberedkeys[j] + num, val );
        if( !found) {
          break ;
        }
        if( i<6 ) {
          reducedline += numberedkeys[j] + num + "=" + val + " ";
        } else {
          ninstances[j]++;
        }
      }
    }
    bool outputcomment=false;
    for(unsigned j=0; j<numberedkeys.size(); ++j) {
      if( ninstances[j]>0 ) {
        outputcomment=true;
        break;
      }
    }
    if( outputcomment ) {
      reducedline += "    # Action input conctinues with ";
      for(unsigned  j=0; j<numberedkeys.size(); ++j) {
        std::string num;
        Tools::convert( ninstances[j], num );
        if( ninstances[j]>0 ) {
          reducedline += num + " further " + numberedkeys[j] + "n keywords, ";
        }
      }
      savedInputLines.push_back( reducedline );
    } else {
      savedInputLines.push_back( line );
    }
  } else {
    savedInputLines.push_back( line );
  }
}

const std::string & ActionShortcut::getShortcutLabel() const {
  return shortcutlabel;
}

const std::vector<std::string>& ActionShortcut::getSavedInputLines() const {
  return savedInputLines;
}

const std::vector<std::string>& ActionShortcut::getSavedOutputs() const {
  return savedOutputs;
}

std::string ActionShortcut::convertInputLineToString() {
  return linemap.convertToString(true);
}

void ActionShortcut::interpretDataLabel( const std::string& mystr, Action* myuser, std::vector<Value*>& arg ) const {
  std::size_t dot=mystr.find_first_of('.');
  std::string a=mystr.substr(0,dot);
  std::string dataName=mystr.substr(dot+1);
  // Retrieve the keywords for the shortcut
  Keywords skeys;
  actionRegister().getKeywords( getName(), skeys );
  std::vector<std::string> out_comps( skeys.getOutputComponents() );
  // Now get the output components
  if( dataName=="*" ) {
    for(unsigned k=0; k<out_comps.size(); ++k) {
      if( out_comps[k]=="" ) {
        ActionWithValue* action=plumed.getActionSet().selectWithLabel<ActionWithValue*>( a );
        if( action ) {
          if( action->getNumberOfComponents()!=1 ) {
            myuser->error("action named " + a + " has more than one component");
          }
          arg.push_back(action->copyOutput(0));
        }
      } else {
        ActionWithValue* action=plumed.getActionSet().selectWithLabel<ActionWithValue*>( a + "_" + out_comps[k] );
        if( action ) {
          if( action->getNumberOfComponents()!=1 ) {
            myuser->error("action named " + a + "_" + out_comps[k] + " has more than one component");
          }
          arg.push_back(action->copyOutput(0));
        } else {
          for(unsigned j=1;; ++j) {
            std::string numstr;
            Tools::convert( j, numstr );
            ActionWithValue* act=plumed.getActionSet().selectWithLabel<ActionWithValue*>( a + "_" + out_comps[k] + "-" + numstr );
            if( act ) {
              for(unsigned n=0; n<act->getNumberOfComponents(); ++n ) {
                arg.push_back(act->copyOutput(n));
              }
            } else if( j>1 ) {
              break;  // This ensures that * syntax works with moments, which normally start from 2
            }
          }
        }
      }
    }
  } else {
    // Check for an action that has action.component
    ActionWithValue* act=plumed.getActionSet().selectWithLabel<ActionWithValue*>( a );
    if( act && act->exists(mystr) ) {
      return;
    }
    // Get components that are actually actions
    for(unsigned k=0; k<out_comps.size(); ++k) {
      if(dataName.find_first_of(out_comps[k])!=std::string::npos ) {
        ActionWithValue* action=plumed.getActionSet().selectWithLabel<ActionWithValue*>( a + "_" + dataName );
        if( action ) {
          arg.push_back(action->copyOutput(a+"_"+dataName));
        }
        break;
      }
    }
  }
}

}
