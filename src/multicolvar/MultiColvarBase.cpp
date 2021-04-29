/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
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
#include "MultiColvarBase.h"
#include "core/ActionShortcut.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include <vector>
#include <string>
#include <limits>

namespace PLMD {
namespace multicolvar {

void MultiColvarBase::shortcutKeywords( Keywords& keys ) {
  keys.add("numbered","LESS_THAN","calculate the number of variables that are less than a certain target value. "
           "This quantity is calculated using \\f$\\sum_i \\sigma(s_i)\\f$, where \\f$\\sigma(s)\\f$ "
           "is a \\ref switchingfunction.");
  keys.addOutputComponent("lessthan","LESS_THAN","the number of colvars that have a value less than a threshold");
  keys.add("numbered","MORE_THAN","calculate the number of variables that are more than a certain target value. "
           "This quantity is calculated using \\f$\\sum_i 1 - \\sigma(s_i)\\f$, where \\f$\\sigma(s)\\f$ "
           "is a \\ref switchingfunction.");
  keys.addOutputComponent("morethan","MORE_THAN","the number of colvars that have a value more than a threshold");
  keys.add("optional","ALT_MIN","calculate the minimum value. "
           "To make this quantity continuous the minimum is calculated using "
           "\\f$ \\textrm{min} = -\\frac{1}{\\beta} \\log \\sum_i \\exp\\left( -\\beta s_i \\right)  \\f$ "
           "The value of \\f$\\beta\\f$ in this function is specified using (BETA=\\f$\\beta\\f$).");
  keys.addOutputComponent("altmin","ALT_MIN","the minimum value of the cv");
  keys.add("optional","MIN","calculate the minimum value. "
           "To make this quantity continuous the minimum is calculated using "
           "\\f$ \\textrm{min} = \\frac{\\beta}{ \\log \\sum_i \\exp\\left( \\frac{\\beta}{s_i} \\right) } \\f$ "
           "The value of \\f$\\beta\\f$ in this function is specified using (BETA=\\f$\\beta\\f$)");
  keys.addOutputComponent("min","MIN","the minimum colvar");
  keys.add("optional","MAX","calculate the maximum value. "
           "To make this quantity continuous the maximum is calculated using "
           "\\f$ \\textrm{max} = \\beta \\log \\sum_i \\exp\\left( \\frac{s_i}{\\beta}\\right) \\f$ "
           "The value of \\f$\\beta\\f$ in this function is specified using (BETA=\\f$\\beta\\f$)");
  keys.addOutputComponent("max","MAX","the maximum colvar");
  keys.add("numbered","BETWEEN","calculate the number of values that are within a certain range. "
           "These quantities are calculated using kernel density estimation as described on "
           "\\ref histogrambead.");
  keys.addOutputComponent("between","BETWEEN","the number of colvars that have a value that lies in a particular interval");
  keys.addFlag("HIGHEST",false,"this flag allows you to recover the highest of these variables.");
  keys.addOutputComponent("highest","HIGHEST","the largest of the colvars");
  keys.add("optional","HISTOGRAM","calculate a discretized histogram of the distribution of values. "
           "This shortcut allows you to calculates NBIN quantites like BETWEEN.");
  keys.addFlag("LOWEST",false,"this flag allows you to recover the lowest of these variables.");
  keys.addOutputComponent("lowest","LOWEST","the smallest of the colvars");
  keys.addFlag("SUM",false,"calculate the sum of all the quantities.");
  keys.addOutputComponent("sum","SUM","the sum of the colvars");
  keys.addFlag("MEAN",false,"calculate the mean of all the quantities.");
  keys.addOutputComponent("mean","MEAN","the mean of the colvars");
}

void MultiColvarBase::expandFunctions( const std::string& labout, const std::string& argin, const std::string& weights, ActionShortcut* action ) {
  std::map<std::string,std::string> keymap; readShortcutKeywords( keymap, action ); expandFunctions( labout, argin, weights, keymap, action );
}

void MultiColvarBase::readShortcutKeywords( std::map<std::string,std::string>& keymap, ActionShortcut* action ) {
  Keywords keys; MultiColvarBase::shortcutKeywords( keys ); action->readShortcutKeywords( keys, keymap );
}

void MultiColvarBase::expandFunctions( const std::string& labout, const std::string& argin, const std::string& weights, 
                                       const std::map<std::string,std::string>& keymap, ActionShortcut* action ) {
  // Parse LESS_THAN
  if( keymap.count("LESS_THAN") ) {
    std::string sum_arg = labout + "_lt", lt_string = keymap.find("LESS_THAN")->second;
    action->readInputLine( labout + "_lt: LESS_THAN ARG1=" + argin + " SWITCH={" + lt_string + "}"); 
    if( weights.length()>0 ) {
        sum_arg = labout + "_wlt";
        action->readInputLine( labout + "_wlt: MATHEVAL ARG1=" + weights + " ARG2=" + labout + "_lt FUNC=x*y PERIODIC=NO");   
    }
    action->readInputLine( labout + "_lessthan: COMBINE ARG=" + sum_arg + " PERIODIC=NO"); 
  }
  if( keymap.count("LESS_THAN1") ) { 
    for(unsigned i=1;; ++i) {
      std::string istr; Tools::convert( i, istr ); 
      if( !keymap.count("LESS_THAN" + istr ) ) { break; }
      std::string sum_arg = labout + "_lt" + istr, lt_string1 = keymap.find("LESS_THAN" + istr)->second;
      action->readInputLine( labout + "_lt" + istr + ": LESS_THAN ARG1=" + argin + " SWITCH={" + lt_string1 + "}");  
      if( weights.length()>0 ) {
          sum_arg = labout + "_wlt" + istr;
          action->readInputLine( labout + "_wlt" + istr + ": MATHEVAL ARG1=" + weights + "ARG2=" + labout + "_lt" + istr + " FUNC=x*y PERIODIC=NO");
      }
      action->readInputLine( labout + "_lessthan" + istr + ": COMBINE ARG=" + sum_arg + " PERIODIC=NO");
    }
  }
  // Parse MORE_THAN
  if( keymap.count("MORE_THAN") ) {
    std::string sum_arg=labout + "_mt", mt_string = keymap.find("MORE_THAN")->second;
    action->readInputLine( labout + "_mt: MORE_THAN ARG1=" + argin + " SWITCH={" + mt_string + "}");
    if( weights.length()>0 ) {
        sum_arg = labout + "_wmt";
        action->readInputLine( labout + "_wmt: MATHEVAL ARG1=" + weights + " ARG2=" + labout + "_mt FUNC=x*y PERIODIC=NO" );
    }
    action->readInputLine( labout + "_morethan: COMBINE ARG=" + sum_arg + " PERIODIC=NO");
  }
  if(  keymap.count("MORE_THAN1") ) {
    for(unsigned i=1;; ++i) {
      std::string istr; Tools::convert( i, istr ); 
      if( !keymap.count("MORE_THAN" + istr ) ) { break; }
      std::string sum_arg = labout + "_mt" + istr, mt_string1 = keymap.find("MORE_THAN" + istr)->second;
      action->readInputLine( labout + "_mt" + istr + ": MORE_THAN ARG1=" + argin + " SWITCH={" + mt_string1 + "}");  
      if( weights.length()>0 ) {
          sum_arg = labout + "_wmt" + istr;
          action->readInputLine( labout + "_wmt" + istr + ": MATHEVAL ARG1=" + weights + "ARG2=" + labout + "_lt" + istr + " FUNC=x*y PERIODIC=NO");
      }
      action->readInputLine( labout + "_morethan" + istr + ": COMBINE ARG=" + sum_arg + " PERIODIC=NO");
      if( !action->parseNumbered("MORE_THAN",i+1,mt_string1) ) { break; }
    }
  }
  // Parse ALT_MIN
  if( keymap.count("ALT_MIN") ) {
    if( weights.length()>0 ) plumed_merror("cannot use ALT_MIN with this shortcut");
    std::string amin_string = keymap.find("ALT_MIN")->second;
    std::size_t dd = amin_string.find("BETA"); std::string beta_str = amin_string.substr(dd+5);
    beta_str.erase(std::remove_if(beta_str.begin(), beta_str.end(), ::isspace), beta_str.end());
    action->readInputLine( labout + "_me_altmin: MATHEVAL ARG1=" + argin + " FUNC=exp(-x*" + beta_str + ") PERIODIC=NO");
    action->readInputLine( labout + "_mec_altmin: COMBINE ARG=" + labout + "_me_altmin PERIODIC=NO");
    action->readInputLine( labout + "_altmin: MATHEVAL ARG=" + labout + "_mec_altmin FUNC=-log(x)/" + beta_str + " PERIODIC=NO");
  }
  // Parse MIN
  if( keymap.count("MIN") ) {
    if( weights.length()>0 ) plumed_merror("cannot use MIN with this shortcut");
    std::string min_string = keymap.find("MIN")->second;
    std::size_t dd = min_string.find("BETA"); std::string beta_str = min_string.substr(dd+5);
    beta_str.erase(std::remove_if(beta_str.begin(), beta_str.end(), ::isspace), beta_str.end());
    action->readInputLine( labout + "_me_min: MATHEVAL ARG1=" + argin + " FUNC=exp(" + beta_str + "/x) PERIODIC=NO");  
    action->readInputLine( labout + "_mec_min: COMBINE ARG=" + labout + "_me_min PERIODIC=NO"); 
    action->readInputLine( labout + "_min: MATHEVAL ARG=" + labout + "_mec_min FUNC=" + beta_str + "/log(x) PERIODIC=NO");
  }
  // Parse MAX
  if( keymap.count("MAX") ) {
    if( weights.length()>0 ) plumed_merror("cannot use MAX with this shortcut");
    std::string max_string = keymap.find("MAX")->second; 
    std::size_t dd = max_string.find("BETA"); std::string beta_str = max_string.substr(dd+5);
    beta_str.erase(std::remove_if(beta_str.begin(), beta_str.end(), ::isspace), beta_str.end());
    action->readInputLine( labout + "_me_max: MATHEVAL ARG1=" + argin + " FUNC=exp(x/" + beta_str + ") PERIODIC=NO");
    action->readInputLine( labout + "_mec_max: COMBINE ARG=" + labout + "_me_max PERIODIC=NO"); 
    action->readInputLine( labout + "_max: MATHEVAL ARG=" + labout + "_mec_max FUNC=" + beta_str  + "*log(x) PERIODIC=NO");
  }
  // Parse HIGHEST
  if( keymap.count("HIGHEST") ) {
    if( weights.length()>0 ) plumed_merror("cannot use HIGHEST with this shortcut");
    action->readInputLine( labout + "_highest: HIGHEST ARG=" + argin );
  }
  // Parse LOWEST
  if( keymap.count("LOWEST") ) {
    if( weights.length()>0 ) plumed_merror("cannot use LOWEST with this shortcut");
    action->readInputLine( labout + "_lowest: LOWEST ARG=" + argin ); 
  }
  // Parse SUM
  if( keymap.count("SUM") ) {
    std::string sum_arg=argin;
    if( weights.length()>0 ) {
      sum_arg = labout + "_wsum";
      action->readInputLine( labout + "_wsum: MATHEVAL ARG1=" + weights + " ARG2=" + argin + " FUNC=x*y PERIODIC=NO");
    }
    action->readInputLine( labout + "_sum: COMBINE ARG=" + sum_arg + " PERIODIC=NO");
  }
  // Parse MEAN
  if( keymap.count("MEAN") ) {
    if( weights.length()>0 ) plumed_merror("cannot use LOWEST with this shortcut");
    action->readInputLine( labout + "_mean: COMBINE ARG=" + argin + " NORMALIZE PERIODIC=NO");
  }
  // Parse BETWEEN
  if( keymap.count("BETWEEN") ) {
    std::string sum_arg=labout + "_bt", bt_string = keymap.find("BETWEEN")->second;
    action->readInputLine( labout + "_bt: BETWEEN ARG1=" + argin + " SWITCH={" + bt_string + "}" );
    if( weights.length()>0 ) {
      sum_arg = labout + "_wbt";
      action->readInputLine( labout + "_wbt: MATHEVAL ARG1=" + weights + " ARG2=" + labout + "_bt FUNC=x*y PERIODIC=NO");
    }
    action->readInputLine( labout + "_between: COMBINE ARG=" + sum_arg + " PERIODIC=NO");
  }
  std::string bt_string1;
  if( keymap.count("BETWEEN1") ) {
    for(unsigned i=1;; ++i) {
      std::string istr; Tools::convert( i, istr ); 
      if( !keymap.count("BETWEEN" + istr) ) break;
      std::string sum_arg=labout + "_bt" + istr, bt_string1 = keymap.find("BETWEEN" + istr)->second;
      action->readInputLine( labout + "_bt" + istr + ": BETWEEN ARG1=" + argin + " SWITCH={" + bt_string1 + "}" );
      if( weights.length()>0 ) {
        sum_arg = labout + "_wbt" + istr;
        action->readInputLine( labout + "_wbt" + istr + ": MATHEVAL ARG1=" + weights + " ARG2=" + labout + "_bt" + istr + " FUNC=x*y PERIODIC=NO");
      }
      action->readInputLine( labout + "_between" + istr + ": COMBINE ARG=" + sum_arg + " PERIODIC=NO");
    }
  }
  // Parse HISTOGRAM
  if( keymap.count("HISTOGRAM") ) {
    std::vector<std::string> words=Tools::getWords( keymap.find("HISTOGRAM")->second );
    unsigned nbins; bool found=Tools::parse(words,"NBINS",nbins,0); // Need replica index
    if( !found ) plumed_merror("did not find NBINS in specification for HISTOGRAM");
    double lower; found=Tools::parse(words,"LOWER",lower,0);
    if( !found ) plumed_merror("did not find LOWER in specification for HISTOGRAM");
    double upper; found=Tools::parse(words,"UPPER",upper,0);
    if( !found ) plumed_merror("did not find UPPER in specification for HISTOGRAM");
    double delr = ( upper - lower ) / static_cast<double>( nbins );
    double smear=0.5; found=Tools::parse(words,"SMEAR",smear,0);
    if( !found ) smear = 0.5;
    for(unsigned i=0; i<nbins; ++i) {
      std::string smstr, istr; Tools::convert( i+1, istr ); Tools::convert( smear, smstr ); std::string sum_arg=labout + "_bt" + istr;
      std::string low_str, high_str; Tools::convert( lower + i*delr, low_str ); Tools::convert( lower + (i+1)*delr, high_str );
      action->readInputLine( labout + "_bt" + istr + ": BETWEEN ARG1=" + argin + " SWITCH={" + words[0] + " LOWER=" + low_str + " UPPER=" + high_str + " SMEAR=" + smstr + "}");
      if( weights.length()>0 ) {
        sum_arg = labout + "_wbt" + istr;
        action->readInputLine( labout + "_wbt" + istr + ": MATHEVAL ARG1=" + weights + " ARG2=" + labout + "_bt" + istr + " FUNC=x*y PERIODIC=NO");
      }
      action->readInputLine( labout + "_between" + istr + ": COMBINE ARG=" + sum_arg + " PERIODIC=NO"); 
    }
  }
}

void MultiColvarBase::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  keys.add("numbered","ATOMS","the atoms involved in each of the colvars you wish to calculate. "
           "Keywords like ATOMS1, ATOMS2, ATOMS3,... should be listed and one or more scalars will be "
           "calculated for each ATOM keyword you specify");
  keys.add("numbered","LOCATION","the location at which the CV is assumed to be in space");
  keys.addFlag("NOLOCATION",false,"do not create a virtual atom for this location");
  keys.reset_style("ATOMS","atoms"); keys.reset_style("LOCATION","atoms");
}

MultiColvarBase::MultiColvarBase(const ActionOptions& ao):
  Action(ao),
  ActionAtomistic(ao),
  ActionWithValue(ao),
  usepbc(true),
  nolocation(false)
{
  if( keywords.exists("NOPBC") ) {
    bool nopbc=!usepbc; parseFlag("NOPBC",nopbc);
    usepbc=!nopbc;
  }
  if( usepbc ) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");
  parseFlag("NOLOCATION",nolocation);

  std::vector<AtomNumber> catoms, all_atoms; parseAtomList( "ATOMS", all_atoms );
  if( getName()=="TORSION" ) {
      std::vector<AtomNumber> v1, v2, axis; parseAtomList("VECTORA", v1 ); parseAtomList("VECTORB", v2 ); parseAtomList("AXIS", axis );
      if( v1.size()>0 ) {
          if( all_atoms.size()>0 ) error("cannot mix ATOMS with VECTORA/VECTORB/AXIS"); 
          if( v1.size()!=2 || v2.size()!=2 || axis.size()!=2 ) error("wrong number of atoms specified to VECTORA, VECTORB or AXIS keyword");
          all_atoms.resize(6); all_atoms[0]=v1[1]; all_atoms[1]=v1[0]; all_atoms[2]=axis[0]; all_atoms[3]=axis[1]; all_atoms[4]=v2[0]; all_atoms[5]=v2[1];
      } else if( all_atoms.size()==4 ) {
          all_atoms.resize(6); all_atoms[5]=all_atoms[3]; all_atoms[4]=all_atoms[2]; all_atoms[3]=all_atoms[2]; all_atoms[2]=all_atoms[1];
      } else if( all_atoms.size()!=0 ) error("wrong number of atoms specified to torsion");
  }
  if( all_atoms.size()>0 ) {
    ablocks.resize(all_atoms.size());
    log.printf("  Colvar is calculated from atoms : ");
    for(unsigned j=0; j<ablocks.size(); ++j) { ablocks[j].push_back(j); log.printf("%d ",all_atoms[j].serial() ); }
    log.printf("\n"); parseAtomList("LOCATION",catoms);
    if( catoms.size()>0 ) {
      if( catoms.size()!=1 ) error("should provide position of one atom only for location");
      log.printf("  CV is located on position of atom : %d \n", catoms[0].serial() );
      catom_indices.push_back( all_atoms.size() ); mygroup.push_back( catoms[0] );
    } else if( !nolocation ) {
      log.printf("  CV is located at center of mass for atoms : ");
      for(unsigned j=0; j<ablocks.size(); ++j) log.printf("%d ", all_atoms[j].serial() );
      log.printf("\n"); mygroup.push_back( atoms.addVirtualAtom( this ) );
    }
  } else {
    std::vector<AtomNumber> t;
    for(int i=1;; ++i ) {
      parseAtomList("ATOMS", i, t );
      if( getName()=="TORSION" ) {
          if( t.empty() ) {
              std::vector<AtomNumber> v1; parseAtomList("VECTORA", i, v1 );
              if( v1.empty() ) break; 
              std::vector<AtomNumber> v2; parseAtomList("VECTORB", i, v2 );               
              std::vector<AtomNumber> axis; parseAtomList("AXIS", i, axis );
              if( v1.size()!=2 || v2.size()!=2 || axis.size()!=2 ) error("wrong number of atoms specified to VECTORA, VECTORB or AXIS keyword");
              t.resize(6); t[0]=v1[1]; t[1]=v1[0]; t[2]=axis[0]; t[3]=axis[1]; t[4]=v2[0]; t[5]=v2[1];
          } else if( t.size()==4 ) {
              std::vector<AtomNumber> v1, v2, axis; parseAtomList("VECTORA", v1 ); parseAtomList("VECTORB", v2 ); parseAtomList("AXIS", axis );
              if( v1.size()>0 || v2.size()>0 || axis.size()>0 ) error("cannot mix ATOMS with VECTORA/VECTORB/AXIS");
              t.resize(6); t[5]=t[3]; t[4]=t[2]; t[3]=t[2]; t[2]=t[1];
          } else plumed_error();
      }
      if( t.empty() ) break;

      log.printf("  Colvar %d is calculated from atoms : ", i);
      for(unsigned j=0; j<t.size(); ++j) log.printf("%d ",t[j].serial() );
      log.printf("\n");

      if( i==1 ) { ablocks.resize(t.size()); }
      if( t.size()!=ablocks.size() ) {
        std::string ss; Tools::convert(i,ss);
        error("ATOMS" + ss + " keyword has the wrong number of atoms");
      }
      for(unsigned j=0; j<ablocks.size(); ++j) {
        ablocks[j].push_back( ablocks.size()*(i-1)+j ); all_atoms.push_back( t[j] );
      }
      t.resize(0);
    }
    parseAtomList("LOCATION", 1, catoms );
    if( catoms.size()>0 ) {
      if( catoms.size()!=1 ) error("should provide position of one atom only for location");
      log.printf("  CV 1 is located on position of atom : %d \n", catoms[0].serial() );
      catom_indices.push_back( all_atoms.size() ); mygroup.push_back( catoms[0] );

      for(int i=2; i<=ablocks[0].size(); ++i) {
        std::vector<AtomNumber> cc; parseAtomList("LOCATION", i, cc );
        if( cc.empty() ) error("LOCATION should be specified for all or none of the atoms in your CV");

        log.printf("  CV %d is located on position of atom : %d \n", i, cc[0].serial() );
        catom_indices.push_back( all_atoms.size() + i ); catoms.push_back( cc[0] ); mygroup.push_back( cc[0] );
      }
    } else if( !nolocation) {
      for(int i=0; i<ablocks[0].size(); ++i) mygroup.push_back( atoms.addVirtualAtom( this ) );
    }
  }
  std::vector<AtomNumber> atoms_for_request(all_atoms); // atoms.insertGroup( getLabel(), mygroup );
  if( catoms.size()>0 ) atoms_for_request.insert( atoms_for_request.end(),catoms.begin(),catoms.end() );
  requestAtoms(atoms_for_request); forcesToApply.resize( getNumberOfDerivatives() );
  if( all_atoms.size()>0 ) {
    for(unsigned i=0; i<ablocks[0].size(); ++i) addTaskToList( i );
  }
  if( catom_indices.size()==0 && !nolocation ) vatom_forces.resize( getNumberOfAtoms() );
  // Create a useful group to hold these atoms
  if( !nolocation ) {
      std::string num; Tools::convert( mygroup[0].serial(), num ); std::string grp_str = getLabel() + "_grp: GROUP ATOMS=" + num;
      for(unsigned i=1;i<mygroup.size();++i) { Tools::convert( mygroup[i].serial(), num ); grp_str += "," + num; }
      plumed.readInputLine( grp_str );
  }
}

MultiColvarBase::~MultiColvarBase() {
  if(catom_indices.size()==0 && !nolocation) atoms.removeVirtualAtom( this );
  // atoms.removeGroup( getLabel() );
}

void MultiColvarBase::addValueWithDerivatives() {
  if( getFullNumberOfTasks()==1 ) { ActionWithValue::addValueWithDerivatives(); }
  else addValue();
}

void MultiColvarBase::addValue() {
  std::vector<unsigned> shape;
  if( getFullNumberOfTasks()>1 ) { shape.resize(1); shape[0]=getFullNumberOfTasks(); }
  ActionWithValue::addValue( shape );
}

void MultiColvarBase::addComponentWithDerivatives( const std::string& name ) {
  if( getFullNumberOfTasks()==1 ) { ActionWithValue::addComponentWithDerivatives(name); }
  else addComponent( name );
}

void MultiColvarBase::addComponent( const std::string& name ) {
  std::vector<unsigned> shape;
  if( getFullNumberOfTasks()>1 ) { shape.resize(1); shape[0]=getFullNumberOfTasks(); }
  ActionWithValue::addComponent( name, shape );
}

void MultiColvarBase::useFourAtomsForEachCV() {
  std::vector<std::vector<unsigned> > tblocks( 4 );
  for(unsigned i=0; i<getFullNumberOfTasks(); ++i) {
    tblocks[0].push_back(ablocks[0][i]); tblocks[1].push_back(ablocks[1][i]);
    tblocks[2].push_back(ablocks[1][i]); tblocks[3].push_back(ablocks[2][i]);
  }
  ablocks.resize(0); ablocks.resize(4);
  for(unsigned i=0; i<getFullNumberOfTasks(); ++i) {
    for(unsigned j=0; j<4; ++j) ablocks[j].push_back(tblocks[j][i]);
  }
}

Vector MultiColvarBase::getSeparation( const Vector& vec1, const Vector& vec2 ) const {
  if(usepbc) { return pbcDistance( vec1, vec2 ); }
  else { return delta( vec1, vec2 ); }
}

void MultiColvarBase::calculate() {
  // Set positions of all virtual atoms
  if( catom_indices.size()==0 && !nolocation ) {
    unsigned stride=comm.Get_size();
    unsigned rank=comm.Get_rank();
    if( runInSerial() ) { stride=1; rank=0; }
    std::vector<Vector> catomp( getFullNumberOfTasks() );
    for(unsigned i=rank; i<getFullNumberOfTasks(); i+=stride) {
      catomp[i].zero(); double wsum = 0.0;
      for(unsigned j=0; j<ablocks.size(); ++j) {
        bool newi=true;
        for(unsigned k=0; k<j; ++k) {
          if( ablocks[j][i]==ablocks[k][i] ) { newi=false; break; }
        }
        if( !newi ) continue;
        wsum += 1.0; catomp[i] += getPosition( ablocks[j][i] );
      }
      double normaliz = 1 / wsum; catomp[i] *= normaliz;
    }
    if( !runInSerial() ) comm.Sum( catomp );
    for(unsigned i=0; i<getFullNumberOfTasks(); ++i) atoms.setVatomPosition( mygroup[i], catomp[i] );
  }
  runAllTasks();
}

void MultiColvarBase::performTask( const unsigned& task_index, MultiValue& myvals ) const {
  // Set the positions for this particular CV
  std::vector<Vector> & fpositions( myvals.getFirstAtomVector() );
  if( fpositions.size()!=ablocks.size() ) fpositions.resize( ablocks.size() );
  for(unsigned i=0; i<ablocks.size(); ++i) fpositions[i] = getPosition( ablocks[i][task_index] );
  // If we are using pbc make whole
  if( usepbc ) {
    for(unsigned j=0; j<fpositions.size()-1; ++j) {
      const Vector & first (fpositions[j]); Vector & second (fpositions[j+1]);
      second=first+pbcDistance(first,second);
    }
  }
  // And compute
  compute( fpositions, myvals );
  // Now update the active derivatives
  if( !doNotCalculateDerivatives() ) {
    for(unsigned i=0; i<ablocks.size(); ++i) {
      // Check for duplicated indices during update to avoid double counting
      bool newi=true;
      for(unsigned j=0; j<i; ++j) {
        if( ablocks[j][task_index]==ablocks[i][task_index] ) { newi=false; break; }
      }
      if( !newi ) continue;
      unsigned base=3*ablocks[i][task_index];
      for(unsigned j=0; j<getNumberOfComponents(); ++j) {
        myvals.updateIndex( getPntrToOutput(j)->getPositionInStream(), base );
        myvals.updateIndex( getPntrToOutput(j)->getPositionInStream(), base + 1 );
        myvals.updateIndex( getPntrToOutput(j)->getPositionInStream(), base + 2 );
      }
    }
    unsigned nvir=3*getNumberOfAtoms();
    for(unsigned j=0; j<getNumberOfComponents(); ++j) {
      for(unsigned i=0; i<9; ++i) myvals.updateIndex( getPntrToOutput(j)->getPositionInStream(), nvir + i );
    }
  }
}

void MultiColvarBase::setBoxDerivativesNoPbc( const unsigned& ival, const std::vector<Vector>& fpositions, MultiValue& myvals ) const {
  if( doNotCalculateDerivatives() ) return;
  Tensor virial; unsigned itask = myvals.getTaskIndex();
  for(unsigned i=0; i<ablocks.size(); i++) {
    // Check for duplicated indices during update to avoid double counting
    bool newi=true;
    for(unsigned j=0; j<i; ++j) {
      if( ablocks[j][itask]==ablocks[i][itask] ) { newi=false; break; }
    }
    if( !newi ) continue;
    virial-=Tensor( fpositions[i], Vector(myvals.getDerivative(ival,3*ablocks[i][itask]+0),
                                          myvals.getDerivative(ival,3*ablocks[i][itask]+1),
                                          myvals.getDerivative(ival,3*ablocks[i][itask]+2)));
  }
  addBoxDerivatives(ival,virial,myvals);
}

void MultiColvarBase::apply() {
  if( doNotCalculateDerivatives() ) return;
  std::fill(forcesToApply.begin(),forcesToApply.end(),0); unsigned mm=0;
  if( getForcesFromValues( forcesToApply ) ) setForcesOnAtoms( forcesToApply, mm );

  // Virtual atom forces
  if( catom_indices.size()==0 && !nolocation ) {
    unsigned stride=comm.Get_size();
    unsigned rank=comm.Get_rank();
    if( runInSerial() ) { stride=1; rank=0; }
    // Clear the forces
    for(unsigned i=0; i<getNumberOfAtoms(); ++i) vatom_forces[i].zero();
    // Accumulate the force on each virtual atom
    for(unsigned i=rank; i<getFullNumberOfTasks(); i+=stride) {
      Vector & f(atoms.getVatomForces(mygroup[i]));
      //printf("FORCES %s %d %f %f %f \n",getLabel().c_str(), i,f[0],f[1],f[2]);
      double wsum=0.0;
      for(unsigned j=0; j<ablocks.size(); ++j) {
        bool newi=true;
        for(unsigned k=0; k<j; ++k) {
          if( ablocks[j][i]==ablocks[k][i] ) { newi=false; break; }
        }
        if( !newi ) continue;
        wsum += 1.0;
      }
      for(unsigned j=0; j<ablocks.size(); ++j) {
        bool newi=true;
        for(unsigned k=0; k<j; ++k) {
          if( ablocks[j][i]==ablocks[k][i] ) { newi=false; break; }
        }
        if( !newi ) continue;
        vatom_forces[ablocks[j][i]] += matmul( (1./wsum)*Tensor::identity(), f );
      }
    }
    if( !runInSerial() ) comm.Sum( vatom_forces );
    // Add the final forces to the atoms
    std::vector<Vector>& final_forces(modifyForces());
    for(unsigned i=0; i<final_forces.size(); ++i) final_forces[i] += vatom_forces[i];
    // Clear the forces on the virtual atoms
    for(unsigned i=0; i<getFullNumberOfTasks(); ++i) atoms.getVatomForces(mygroup[i]).zero();
  }
}

}
}
