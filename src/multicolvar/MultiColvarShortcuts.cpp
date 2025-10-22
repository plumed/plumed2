/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2023 The plumed team
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
#include "MultiColvarShortcuts.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/Group.h"

namespace PLMD {
namespace multicolvar {

void MultiColvarShortcuts::shortcutKeywords( Keywords& keys ) {
  keys.add("numbered","LESS_THAN","calculate the number of variables that are less than a certain target value. "
           "This quantity is calculated using \\f$\\sum_i \\sigma(s_i)\\f$, where \\f$\\sigma(s)\\f$ "
           "is a \\ref switchingfunction.");
  keys.linkActionInDocs("LESS_THAN","LESS_THAN");
  keys.reset_style("LESS_THAN","deprecated");
  keys.addOutputComponent("lessthan","LESS_THAN","scalar","the number of colvars that have a value less than a threshold");
  keys.add("numbered","MORE_THAN","calculate the number of variables that are more than a certain target value. "
           "This quantity is calculated using \\f$\\sum_i 1 - \\sigma(s_i)\\f$, where \\f$\\sigma(s)\\f$ "
           "is a \\ref switchingfunction.");
  keys.linkActionInDocs("MORE_THAN","MORE_THAN");
  keys.reset_style("MORE_THAN","deprecated");
  keys.addOutputComponent("morethan","MORE_THAN","scalar","the number of colvars that have a value more than a threshold");
  keys.add("optional","ALT_MIN","calculate the minimum value. "
           "To make this quantity continuous the minimum is calculated using "
           "\\f$ \\textrm{min} = -\\frac{1}{\\beta} \\log \\sum_i \\exp\\left( -\\beta s_i \\right)  \\f$ "
           "The value of \\f$\\beta\\f$ in this function is specified using (BETA=\\f$\\beta\\f$).");
  keys.reset_style("ALT_MIN","deprecated");
  keys.addOutputComponent("altmin","ALT_MIN","scalar","the minimum value of the cv");
  keys.add("optional","MIN","calculate the minimum value. "
           "To make this quantity continuous the minimum is calculated using "
           "\\f$ \\textrm{min} = \\frac{\\beta}{ \\log \\sum_i \\exp\\left( \\frac{\\beta}{s_i} \\right) } \\f$ "
           "The value of \\f$\\beta\\f$ in this function is specified using (BETA=\\f$\\beta\\f$)");
  keys.reset_style("MIN","deprecated");
  keys.addOutputComponent("min","MIN","scalar","the minimum colvar");
  keys.add("optional","MAX","calculate the maximum value. "
           "To make this quantity continuous the maximum is calculated using "
           "\\f$ \\textrm{max} = \\beta \\log \\sum_i \\exp\\left( \\frac{s_i}{\\beta}\\right) \\f$ "
           "The value of \\f$\\beta\\f$ in this function is specified using (BETA=\\f$\\beta\\f$)");
  keys.reset_style("MAX","deprecated");
  keys.addOutputComponent("max","MAX","scalar","the maximum colvar");
  keys.add("numbered","BETWEEN","calculate the number of values that are within a certain range. "
           "These quantities are calculated using kernel density estimation as described on "
           "\\ref histogrambead.");
  keys.linkActionInDocs("BETWEEN","BETWEEN");
  keys.reset_style("BETWEEN","deprecated");
  keys.addOutputComponent("between","BETWEEN","scalar","the number of colvars that have a value that lies in a particular interval");
  keys.addFlag("HIGHEST",false,"this flag allows you to recover the highest of these variables.");
  keys.reset_style("HIGHEST","deprecated");
  keys.addOutputComponent("highest","HIGHEST","scalar","the largest of the colvars");
  keys.add("optional","HISTOGRAM","calculate a discretized histogram of the distribution of values. "
           "This shortcut allows you to calculates NBIN quantites like BETWEEN.");
  keys.reset_style("HISTOGRAM","deprecated");
  keys.addFlag("LOWEST",false,"this flag allows you to recover the lowest of these variables.");
  keys.reset_style("LOWEST","deprecated");
  keys.addOutputComponent("lowest","LOWEST","scalar","the smallest of the colvars");
  keys.addFlag("SUM",false,"calculate the sum of all the quantities.");
  keys.reset_style("SUM","deprecated");
  keys.addOutputComponent("sum","SUM","scalar","the sum of the colvars");
  keys.addFlag("MEAN",false,"calculate the mean of all the quantities.");
  keys.reset_style("MEAN","deprecated");
  keys.addOutputComponent("mean","MEAN","scalar","the mean of the colvars");
  keys.needsAction("SUM");
  keys.needsAction("MEAN");
  keys.needsAction("CUSTOM");
  keys.needsAction("HIGHEST");
  keys.needsAction("LOWEST");
  keys.needsAction("LESS_THAN");
  keys.needsAction("MORE_THAN");
  keys.needsAction("BETWEEN");
}

void MultiColvarShortcuts::expandFunctions( const std::string& labout, const std::string& argin, const std::string& weights, ActionShortcut* action ) {
  std::map<std::string,std::string> keymap;
  readShortcutKeywords( keymap, action );
  expandFunctions( labout, argin, weights, keymap, action );
}

void MultiColvarShortcuts::readShortcutKeywords( std::map<std::string,std::string>& keymap, ActionShortcut* action ) {
  Keywords keys;
  shortcutKeywords( keys );
  action->readShortcutKeywords( keys, keymap );
}

void MultiColvarShortcuts::parseAtomList( const std::string& key, std::vector<std::string>& atoms, ActionShortcut* action ) {
  std::vector<std::string> astr;
  action->parseVector(key,astr);
  if( astr.size()==0 ) {
    return ;
  }
  Tools::interpretRanges( astr );
  for(unsigned i=0; i<astr.size(); ++i) {
    Group* mygr=action->plumed.getActionSet().selectWithLabel<Group*>(astr[i]);
    if( mygr ) {
      std::vector<std::string> grstr( mygr->getGroupAtoms() );
      for(unsigned j=0; j<grstr.size(); ++j) {
        atoms.push_back(grstr[j]);
      }
    } else {
      Group* mygr2=action->plumed.getActionSet().selectWithLabel<Group*>(astr[i] + "_grp");
      if( mygr2 ) {
        std::vector<std::string> grstr( mygr2->getGroupAtoms() );
        for(unsigned j=0; j<grstr.size(); ++j) {
          atoms.push_back(grstr[j]);
        }
      } else {
        atoms.push_back(astr[i]);
      }
    }
  }
}

void MultiColvarShortcuts::expandFunctions( const std::string& labout, const std::string& argin, const std::string& weights,
    const std::map<std::string,std::string>& keymap, ActionShortcut* action ) {
  if( keymap.empty() ) {
    return;
  }
  // Parse LESS_THAN
  if( keymap.count("LESS_THAN") ) {
    std::string sum_arg = labout + "_lt", lt_string = keymap.find("LESS_THAN")->second;
    action->readInputLine( labout + "_lt: LESS_THAN ARG=" + argin + " SWITCH={" + lt_string + "}");
    if( weights.length()>0 ) {
      sum_arg = labout + "_wlt";
      action->readInputLine( labout + "_wlt: CUSTOM ARG=" + weights + "," + labout + "_lt FUNC=x*y PERIODIC=NO");
    }
    action->readInputLine( labout + "_lessthan: SUM ARG=" + sum_arg + " PERIODIC=NO");
  }
  if( keymap.count("LESS_THAN1") ) {
    for(unsigned i=1;; ++i) {
      std::string istr;
      Tools::convert( i, istr );
      if( !keymap.count("LESS_THAN" + istr ) ) {
        break;
      }
      std::string sum_arg = labout + "_lt" + istr, lt_string1 = keymap.find("LESS_THAN" + istr)->second;
      action->readInputLine( labout + "_lt" + istr + ": LESS_THAN ARG=" + argin + " SWITCH={" + lt_string1 + "}");
      if( weights.length()>0 ) {
        sum_arg = labout + "_wlt" + istr;
        action->readInputLine( labout + "_wlt" + istr + ": CUSTOM ARG=" + weights + "," + labout + "_lt" + istr + " FUNC=x*y PERIODIC=NO");
      }
      action->readInputLine( labout + "_lessthan-" + istr + ": SUM ARG=" + sum_arg + " PERIODIC=NO");
    }
  }
  // Parse MORE_THAN
  if( keymap.count("MORE_THAN") ) {
    std::string sum_arg=labout + "_mt", mt_string = keymap.find("MORE_THAN")->second;
    action->readInputLine( labout + "_mt: MORE_THAN ARG=" + argin + " SWITCH={" + mt_string + "}");
    if( weights.length()>0 ) {
      sum_arg = labout + "_wmt";
      action->readInputLine( labout + "_wmt: CUSTOM ARG=" + weights + "," + labout + "_mt FUNC=x*y PERIODIC=NO" );
    }
    action->readInputLine( labout + "_morethan: SUM ARG=" + sum_arg + " PERIODIC=NO");
  }
  if(  keymap.count("MORE_THAN1") ) {
    for(unsigned i=1;; ++i) {
      std::string istr;
      Tools::convert( i, istr );
      if( !keymap.count("MORE_THAN" + istr ) ) {
        break;
      }
      std::string sum_arg = labout + "_mt" + istr, mt_string1 = keymap.find("MORE_THAN" + istr)->second;
      action->readInputLine( labout + "_mt" + istr + ": MORE_THAN ARG=" + argin + " SWITCH={" + mt_string1 + "}");
      if( weights.length()>0 ) {
        sum_arg = labout + "_wmt" + istr;
        action->readInputLine( labout + "_wmt" + istr + ": CUSTOM ARG=" + weights + "," + labout + "_lt" + istr + " FUNC=x*y PERIODIC=NO");
      }
      action->readInputLine( labout + "_morethan-" + istr + ": SUM ARG=" + sum_arg + " PERIODIC=NO");
    }
  }
  // Parse ALT_MIN
  if( keymap.count("ALT_MIN") ) {
    if( weights.length()>0 ) {
      plumed_merror("cannot use ALT_MIN with this shortcut");
    }
    std::string amin_string = keymap.find("ALT_MIN")->second;
    std::size_t dd = amin_string.find("BETA");
    std::string beta_str = amin_string.substr(dd+5);
    beta_str.erase(std::remove_if(beta_str.begin(), beta_str.end(), ::isspace), beta_str.end());
    action->readInputLine( labout + "_me_altmin: CUSTOM ARG=" + argin + " FUNC=exp(-x*" + beta_str + ") PERIODIC=NO");
    action->readInputLine( labout + "_mec_altmin: SUM ARG=" + labout + "_me_altmin PERIODIC=NO");
    action->readInputLine( labout + "_altmin: CUSTOM ARG=" + labout + "_mec_altmin FUNC=-log(x)/" + beta_str + " PERIODIC=NO");
  }
  // Parse MIN
  if( keymap.count("MIN") ) {
    if( weights.length()>0 ) {
      plumed_merror("cannot use MIN with this shortcut");
    }
    std::string min_string = keymap.find("MIN")->second;
    std::size_t dd = min_string.find("BETA");
    std::string beta_str = min_string.substr(dd+5);
    beta_str.erase(std::remove_if(beta_str.begin(), beta_str.end(), ::isspace), beta_str.end());
    action->readInputLine( labout + "_me_min: CUSTOM ARG=" + argin + " FUNC=exp(" + beta_str + "/x) PERIODIC=NO");
    action->readInputLine( labout + "_mec_min: SUM ARG=" + labout + "_me_min PERIODIC=NO");
    action->readInputLine( labout + "_min: CUSTOM ARG=" + labout + "_mec_min FUNC=" + beta_str + "/log(x) PERIODIC=NO");
  }
  // Parse MAX
  if( keymap.count("MAX") ) {
    if( weights.length()>0 ) {
      plumed_merror("cannot use MAX with this shortcut");
    }
    std::string max_string = keymap.find("MAX")->second;
    std::size_t dd = max_string.find("BETA");
    std::string beta_str = max_string.substr(dd+5);
    beta_str.erase(std::remove_if(beta_str.begin(), beta_str.end(), ::isspace), beta_str.end());
    action->readInputLine( labout + "_me_max: CUSTOM ARG=" + argin + " FUNC=exp(x/" + beta_str + ") PERIODIC=NO");
    action->readInputLine( labout + "_mec_max: SUM ARG=" + labout + "_me_max PERIODIC=NO");
    action->readInputLine( labout + "_max: CUSTOM ARG=" + labout + "_mec_max FUNC=" + beta_str  + "*log(x) PERIODIC=NO");
  }
  // Parse HIGHEST
  if( keymap.count("HIGHEST") ) {
    if( weights.length()>0 ) {
      plumed_merror("cannot use HIGHEST with this shortcut");
    }
    action->readInputLine( labout + "_highest: HIGHEST ARG=" + argin );
  }
  // Parse LOWEST
  if( keymap.count("LOWEST") ) {
    if( weights.length()>0 ) {
      plumed_merror("cannot use LOWEST with this shortcut");
    }
    action->readInputLine( labout + "_lowest: LOWEST ARG=" + argin );
  }
  // Parse SUM
  if( keymap.count("SUM") ) {
    std::string sum_arg=argin;
    if( weights.length()>0 ) {
      sum_arg = labout + "_wsum";
      action->readInputLine( labout + "_wsum: CUSTOM ARG=" + weights + "," + argin + " FUNC=x*y PERIODIC=NO");
    }
    action->readInputLine( labout + "_sum: SUM ARG=" + sum_arg + " PERIODIC=NO");
  }
  // Parse MEAN
  if( keymap.count("MEAN") ) {
    if( weights.length()>0 ) {
      plumed_merror("cannot use MEAN with this shortcut");
    }
    action->readInputLine( labout + "_mean: MEAN ARG=" + argin + " PERIODIC=NO");
  }
  // Parse BETWEEN
  if( keymap.count("BETWEEN") ) {
    std::string sum_arg=labout + "_bt";
    std::string bt_string = keymap.find("BETWEEN")->second;
    action->readInputLine( labout + "_bt: BETWEEN ARG=" + argin + " SWITCH={" + bt_string + "}" );
    if( weights.length()>0 ) {
      sum_arg = labout + "_wbt";
      action->readInputLine( labout + "_wbt: CUSTOM ARG=" + weights + "," + labout + "_bt FUNC=x*y PERIODIC=NO");
    }
    action->readInputLine( labout + "_between: SUM ARG=" + sum_arg + " PERIODIC=NO");
  }

  if( keymap.count("BETWEEN1") ) {
    for(unsigned i=1;; ++i) {
      std::string istr;
      Tools::convert( i, istr );
      if( !keymap.count("BETWEEN" + istr) ) {
        break;
      }
      std::string sum_arg=labout + "_bt" + istr;
      std::string bt_string1 = keymap.find("BETWEEN" + istr)->second;
      action->readInputLine( labout + "_bt" + istr + ": BETWEEN ARG=" + argin + " SWITCH={" + bt_string1 + "}" );
      if( weights.length()>0 ) {
        sum_arg = labout + "_wbt" + istr;
        action->readInputLine( labout + "_wbt" + istr + ": CUSTOM ARG=" + weights + "," + labout + "_bt" + istr + " FUNC=x*y PERIODIC=NO");
      }
      action->readInputLine( labout + "_between-" + istr + ": SUM ARG=" + sum_arg + " PERIODIC=NO");
    }
  }
  // Parse HISTOGRAM
  if( keymap.count("HISTOGRAM") ) {
    std::vector<std::string> words=Tools::getWords( keymap.find("HISTOGRAM")->second );
    unsigned nbins;
    bool found=Tools::parse(words,"NBINS",nbins,0); // Need replica index
    if( !found ) {
      plumed_merror("did not find NBINS in specification for HISTOGRAM");
    }
    double lower;
    found=Tools::parse(words,"LOWER",lower,0);
    if( !found ) {
      plumed_merror("did not find LOWER in specification for HISTOGRAM");
    }
    double upper;
    found=Tools::parse(words,"UPPER",upper,0);
    if( !found ) {
      plumed_merror("did not find UPPER in specification for HISTOGRAM");
    }
    double delr = ( upper - lower ) / static_cast<double>( nbins );
    double smear=0.5;
    found=Tools::parse(words,"SMEAR",smear,0);
    if( !found ) {
      smear = 0.5;
    }
    for(unsigned i=0; i<nbins; ++i) {
      std::string smstr, istr;
      Tools::convert( i+1, istr );
      Tools::convert( smear, smstr );
      std::string sum_arg=labout + "_bt" + istr;
      std::string low_str, high_str;
      Tools::convert( lower + i*delr, low_str );
      Tools::convert( lower + (i+1)*delr, high_str );
      action->readInputLine( labout + "_bt" + istr + ": BETWEEN ARG=" + argin + " SWITCH={" + words[0] + " LOWER=" + low_str + " UPPER=" + high_str + " SMEAR=" + smstr + "}");
      if( weights.length()>0 ) {
        sum_arg = labout + "_wbt" + istr;
        action->readInputLine( labout + "_wbt" + istr + ": CUSTOM ARG=" + weights + "," + labout + "_bt" + istr + " FUNC=x*y PERIODIC=NO");
      }
      action->readInputLine( labout + "_between-" + istr + ": SUM ARG=" + sum_arg + " PERIODIC=NO");
    }
  }
}

}
}
