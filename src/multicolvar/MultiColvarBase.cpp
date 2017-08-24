/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2017 The plumed team
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
#include "AtomValuePack.h"

namespace PLMD {
namespace multicolvar {

void MultiColvarBase::shortcutKeywords( Keywords& keys ) {
  keys.add("numbered","LESS_THAN","calculate the number of variables that are less than a certain target value. "
                                  "This quantity is calculated using \\f$\\sum_i \\sigma(s_i)\\f$, where \\f$\\sigma(s)\\f$ "
                                  "is a \\ref switchingfunction.");
  keys.addOutputComponent("_lessthan","LESS_THAN","the number of colvars that have a value less than a threshold");
  keys.add("numbered","MORE_THAN","calculate the number of variables that are more than a certain target value. "
                                  "This quantity is calculated using \\f$\\sum_i 1 - \\sigma(s_i)\\f$, where \\f$\\sigma(s)\\f$ "
                                  "is a \\ref switchingfunction.");
  keys.addOutputComponent("_morethan","MORE_THAN","the number of colvars that have a value more than a threshold");
  keys.add("optional","ALT_MIN","calculate the minimum value. "
                                "To make this quantity continuous the minimum is calculated using "
                                "\\f$ \\textrm{min} = -\\frac{1}{\\beta} \\log \\sum_i \\exp\\left( -\\beta s_i \\right)  \\f$ "
                                "The value of \\f$\\beta\\f$ in this function is specified using (BETA=\\f$\\beta\\f$).");
  keys.addOutputComponent("_altmin","ALT_MIN","the minimum value of the cv");
  keys.add("optional","MIN","calculate the minimum value. "
                            "To make this quantity continuous the minimum is calculated using "
                            "\\f$ \\textrm{min} = \\frac{\\beta}{ \\log \\sum_i \\exp\\left( \\frac{\\beta}{s_i} \\right) } \\f$ "
                            "The value of \\f$\\beta\\f$ in this function is specified using (BETA=\\f$\\beta\\f$)");
  keys.addOutputComponent("_min","MIN","the minimum colvar");
  keys.add("optional","MAX","calculate the maximum value. "
                            "To make this quantity continuous the maximum is calculated using "
                            "\\f$ \\textrm{max} = \\beta \\log \\sum_i \\exp\\left( \\frac{s_i}{\\beta}\\right) \\f$ "
                            "The value of \\f$\\beta\\f$ in this function is specified using (BETA=\\f$\\beta\\f$)");
  keys.addOutputComponent("_max","MAX","the maximum colvar");
  keys.add("numbered","BETWEEN","calculate the number of values that are within a certain range. "
                                "These quantities are calculated using kernel density estimation as described on "
                                "\\ref histogrambead."); 
  keys.addOutputComponent("_between","BETWEEN","the number of colvars that have a value that lies in a particular interval");
  keys.addFlag("HIGHEST",false,"this flag allows you to recover the highest of these variables.");
  keys.addOutputComponent("_highest","HIGHEST","the largest of the colvars");
  keys.add("optional","HISTOGRAM","calculate a discretized histogram of the distribution of values. "
                                  "This shortcut allows you to calculates NBIN quantites like BETWEEN.");
  keys.addFlag("LOWEST",false,"this flag allows you to recover the lowest of these variables.");
  keys.addOutputComponent("_lowest","LOWEST","the smallest of the colvars");
  keys.addFlag("SUM",false,"calculate the sum of all the quantities.");
  keys.addOutputComponent("_sum","SUM","the sum of the colvars");
  keys.addFlag("MEAN",false,"calculate the mean of all the quantities.");
  keys.addOutputComponent("_mean","MEAN","the mean of the colvars");
}

void MultiColvarBase::expandFunctions( const std::string& labout, const std::string& argin, 
                                       const std::vector<std::string>& words,
                                       const std::map<std::string,std::string>& keys,
                                       std::vector<std::vector<std::string> >& actions ){
  // Parse LESS_THAN
  if( keys.count("LESS_THAN") ){
      std::vector<std::string> input; input.push_back( labout + "_lt:" ); input.push_back("LESS_THAN"); 
      input.push_back("ARG1=" + argin ); 
      input.push_back("SWITCH=" + keys.find("LESS_THAN")->second  );
      actions.push_back( input );
      std::vector<std::string> sum_inp; sum_inp.push_back( labout + "_lessthan:" ); 
      sum_inp.push_back("COMBINE"); sum_inp.push_back("ARG=" + labout + "_lt");
      sum_inp.push_back("PERIODIC=NO"); actions.push_back( sum_inp ); 
  }
  if( keys.count("LESS_THAN1") ){
      for(unsigned i=1;; ++i){
          std::string istr; Tools::convert( i, istr );  
          if( !keys.count("LESS_THAN" + istr ) ){ break; } 
          std::vector<std::string> input; input.push_back( labout + "_lt" + istr + ":" ); input.push_back("LESS_THAN");
          input.push_back("ARG1=" + argin );
          input.push_back("SWITCH=" + keys.find("LESS_THAN" + istr)->second  );
          actions.push_back( input );
          std::vector<std::string> sum_inp; sum_inp.push_back( labout + "_lessthan" + istr + ":" ); 
          sum_inp.push_back("COMBINE"); sum_inp.push_back("ARG=" + labout + "_lt" + istr );
          sum_inp.push_back("PERIODIC=NO"); actions.push_back( sum_inp );
      }
  }
  // Parse MORE_THAN
  if( keys.count("MORE_THAN") ){
      std::vector<std::string> input; input.push_back( labout + "_mt:" ); input.push_back("MORE_THAN");
      input.push_back("ARG1=" + argin );
      input.push_back("SWITCH=" + keys.find("MORE_THAN")->second  );
      actions.push_back( input );
      std::vector<std::string> sum_inp; sum_inp.push_back( labout + "_morethan:" ); 
      sum_inp.push_back("COMBINE"); sum_inp.push_back("ARG=" + labout + "_mt");
      sum_inp.push_back("PERIODIC=NO"); actions.push_back( sum_inp );
  }
  if( keys.count("MORE_THAN1") ){
      for(unsigned i=1;; ++i){ 
          std::string istr; Tools::convert( i, istr );    
          if( !keys.count("MORE_THAN" + istr ) ){ break; }
          std::vector<std::string> input; input.push_back( labout + "_mt" + istr + ":" ); input.push_back("MORE_THAN");
          input.push_back("ARG1=" + argin );
          input.push_back("SWITCH=" + keys.find("MORE_THAN" + istr)->second  );
          actions.push_back( input );
          std::vector<std::string> sum_inp; sum_inp.push_back( labout + "_morethan" + istr + ":" ); 
          sum_inp.push_back("COMBINE"); sum_inp.push_back("ARG=" + labout + "_mt" + istr );
          sum_inp.push_back("PERIODIC=NO"); actions.push_back( sum_inp );
      }
  }
  // Parse ALT_MIN
  if( keys.count("ALT_MIN") ){
      std::vector<std::string> input; input.push_back( labout + "_altmin:" ); input.push_back("ALT_MIN");
      input.push_back("ARG=" + argin ); std::size_t dd = keys.find("ALT_MIN")->second.find("BETA");
      input.push_back( keys.find("ALT_MIN")->second.substr(dd) ); 
      actions.push_back( input );
  }
  // Parse MIN
  if( keys.count("MIN") ){
      std::vector<std::string> input; input.push_back( labout + "_min:" ); input.push_back("MIN");
      input.push_back("ARG=" + argin ); std::size_t dd = keys.find("MIN")->second.find("BETA"); 
      input.push_back( keys.find("MIN")->second.substr(dd) ); actions.push_back( input );
  }
  // Parse MAX
  if( keys.count("MAX") ){
      std::vector<std::string> input; input.push_back( labout + "_max:" ); input.push_back("MAX");
      input.push_back("ARG=" + argin ); std::size_t dd = keys.find("MAX")->second.find("BETA"); 
      input.push_back( keys.find("MAX")->second.substr(dd) ); actions.push_back( input );
  }
  // Parse HIGHEST
  if( keys.count("HIGHEST") ){
      std::vector<std::string> input; input.push_back( labout + "_highest:" ); input.push_back("HIGHEST");
      input.push_back("ARG=" + argin ); actions.push_back( input );
  }
  // Parse LOWEST
  if( keys.count("LOWEST") ){
      std::vector<std::string> input; input.push_back( labout + "_lowest:" ); input.push_back("LOWEST");
      input.push_back("ARG=" + argin ); actions.push_back( input );
  }
  // Parse SUM
  if( keys.count("SUM") ){
      std::vector<std::string> input; input.push_back( labout + "_sum:" ); 
      input.push_back("COMBINE"); input.push_back("ARG=" + argin ); 
      input.push_back("PERIODIC=NO"); actions.push_back( input );
  }
  // Parse MEAN
  if( keys.count("MEAN") ){
      std::vector<std::string> input; input.push_back( labout + "_mean:" ); input.push_back("COMBINE");
      input.push_back("ARG=" + argin ); input.push_back("NORMALIZE"); 
      input.push_back("PERIODIC=NO"); actions.push_back( input );
  }
  // Parse BETWEEN
  if( keys.count("BETWEEN") ){
      std::vector<std::string> input; input.push_back( labout + "_bt:" ); input.push_back("BETWEEN");
      input.push_back("ARG1=" + argin );
      input.push_back("SWITCH=" + keys.find("BETWEEN")->second  );
      actions.push_back( input );
      std::vector<std::string> sum_inp; sum_inp.push_back( labout + "_between:" ); 
      sum_inp.push_back("COMBINE"); sum_inp.push_back("ARG=" + labout + "_bt");
      sum_inp.push_back("PERIODIC=NO"); actions.push_back( sum_inp );
  }
  if( keys.count("BETWEEN1") ){
      for(unsigned i=1;; ++i){ 
          std::string istr; Tools::convert( i, istr );    
          if( !keys.count("BETWEEN" + istr ) ){ break; }
          std::vector<std::string> input; input.push_back( labout + "_bt" + istr + ":" ); input.push_back("BETWEEN");
          input.push_back("ARG1=" + argin );
          input.push_back("SWITCH=" + keys.find("BETWEEN" + istr)->second  );
          actions.push_back( input );
          std::vector<std::string> sum_inp; sum_inp.push_back( labout + "_between" + istr + ":" ); 
          sum_inp.push_back("COMBINE"); sum_inp.push_back("ARG=" + labout + "_bt" + istr );
          sum_inp.push_back("PERIODIC=NO"); actions.push_back( sum_inp );
      }
  }
  // Parse HISTOGRAM  
  if( keys.count("HISTOGRAM") ){
      std::vector<std::string> words=Tools::getWords( keys.find("HISTOGRAM")->second );
      unsigned nbins; bool found=Tools::parse(words,"NBINS",nbins,0); // Need replica index
      if( !found ) plumed_merror("did not find NBINS in specification for HISTOGRAM");
      double lower; found=Tools::parse(words,"LOWER",lower,0);
      if( !found ) plumed_merror("did not find LOWER in specification for HISTOGRAM"); 
      double upper; found=Tools::parse(words,"UPPER",upper,0);
      if( !found ) plumed_merror("did not find UPPER in specification for HISTOGRAM");
      double delr = ( upper - lower ) / static_cast<double>( nbins );
      double smear=0.5; found=Tools::parse(words,"SMEAR",smear,0);
      if( !found ) smear = 0.5; 
      for(unsigned i=0;i<nbins;++i){
          std::string smstr, istr; Tools::convert( i+1, istr ); Tools::convert( smear, smstr );
          std::vector<std::string> input; input.push_back( labout + "_bt" + istr + ":" ); input.push_back("BETWEEN");
          input.push_back("ARG1=" + argin ); std::string low_str, high_str; 
          Tools::convert( lower + i*delr, low_str ); Tools::convert( lower + (i+1)*delr, high_str );
          input.push_back("SWITCH= " + words[0] + " LOWER=" + low_str + " UPPER=" + high_str + " SMEAR=" + smstr );  actions.push_back( input );   
          std::vector<std::string> sum_inp; sum_inp.push_back( labout + "_between" + istr + ":" ); 
          sum_inp.push_back("COMBINE"); sum_inp.push_back("ARG=" + labout + "_bt" + istr );
          sum_inp.push_back("PERIODIC=NO"); actions.push_back( sum_inp );
      }
  }
}

void MultiColvarBase::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  keys.add("numbered","ATOMS","the atoms involved in each of the colvars you wish to calculate. "
           "Keywords like ATOMS1, ATOMS2, ATOMS3,... should be listed and one or more scalars will be "
           "calculated for each ATOM keyword you specify");
  keys.add("numbered","LOCATION","the location at which the CV is assumed to be in space");
  keys.reset_style("ATOMS","atoms"); keys.reset_style("LOCATION","atoms");
}

MultiColvarBase::MultiColvarBase(const ActionOptions& ao):
Action(ao),
ActionAtomistic(ao),
ActionWithValue(ao),
usepbc(true)
{
  if( keywords.exists("NOPBC") ) {
    bool nopbc=!usepbc; parseFlag("NOPBC",nopbc);
    usepbc=!nopbc;
  }
  if( usepbc ) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");

  std::vector<AtomNumber> catoms, all_atoms; parseAtomList( "ATOMS", all_atoms );
  if( all_atoms.size()>0 ){
      ablocks.resize(all_atoms.size()); 
      log.printf("  Colvar is calculated from atoms : ");
      for(unsigned j=0; j<ablocks.size(); ++j){ ablocks[j].push_back(j); log.printf("%d ",all_atoms[j].serial() ); }
      log.printf("\n"); parseAtomList("LOCATION",catoms);
      if( catoms.size()>0 ){
           if( catoms.size()!=1 ) error("should provide position of one atom only for location");
           log.printf("  CV is located on position of atom : %d \n", catoms[0].serial() ); 
           catom_indices.push_back( all_atoms.size() ); mygroup.push_back( catoms[0] );
      } else {
           log.printf("  CV is located at center of mass for atoms : ");
           for(unsigned j=0; j<ablocks.size(); ++j) log.printf("%d ", all_atoms[j].serial() );
           log.printf("\n"); mygroup.push_back( atoms.addVirtualAtom( this ) );
      }
  } else {
      std::vector<AtomNumber> t;
      for(int i=1;; ++i ) {
        parseAtomList("ATOMS", i, t );
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
      if( catoms.size()>0 ){
          if( catoms.size()!=1 ) error("should provide position of one atom only for location");
          log.printf("  CV 1 is located on position of atom : %d \n", catoms[0].serial() );
          catom_indices.push_back( all_atoms.size() ); mygroup.push_back( catoms[0] );

          for(int i=2;i<ablocks[0].size();++i){
              std::vector<AtomNumber> cc; parseAtomList("LOCATION", 1, cc );
              if( cc.empty() ) error("LOCATION should be specified for all or none of the atoms in your CV");

              log.printf("  CV %d is located on position of atom : %d \n", i, cc[0].serial() );
              catom_indices.push_back( all_atoms.size() + i ); catoms.push_back( cc[0] ); mygroup.push_back( cc[0] );
          }
      } else {
          for(int i=0;i<ablocks[0].size();++i) mygroup.push_back( atoms.addVirtualAtom( this ) );
      }
  }
  std::vector<AtomNumber> atoms_for_request(all_atoms); atoms.insertGroup( getLabel(), mygroup );
  if( catoms.size()>0 ) atoms_for_request.insert( atoms_for_request.end(),catoms.begin(),catoms.end() );
  requestAtoms(atoms_for_request); forcesToApply.resize( 3*all_atoms.size()+9 );
  if( all_atoms.size()>0 ) {
      for(unsigned i=0; i<ablocks[0].size(); ++i) addTaskToList( i );
  }
  if( catom_indices.size()==0 ) vatom_forces.resize( getNumberOfAtoms() );
}

MultiColvarBase::~MultiColvarBase(){
  atoms.removeVirtualAtom( this ); atoms.removeGroup( getLabel() );
}

void MultiColvarBase::addValueWithDerivatives(){
  if( getFullNumberOfTasks()==1 ){ ActionWithValue::addValueWithDerivatives(); }
  else addValue();
}

void MultiColvarBase::addValue(){
  std::vector<unsigned> shape;
  if( getFullNumberOfTasks()>1 ){ shape.resize(1); shape[0]=getFullNumberOfTasks(); }
  ActionWithValue::addValue( shape );
}

void MultiColvarBase::addComponentWithDerivatives( const std::string& name ){
  if( getFullNumberOfTasks()==1 ){ ActionWithValue::addComponentWithDerivatives(name); }
  else addComponent( name );
}

void MultiColvarBase::addComponent( const std::string& name ){
  std::vector<unsigned> shape;
  if( getFullNumberOfTasks()>1 ){ shape.resize(1); shape[0]=getFullNumberOfTasks(); }
  ActionWithValue::addComponent( name, shape );
}

Vector MultiColvarBase::getSeparation( const Vector& vec1, const Vector& vec2 ) const {
  if(usepbc) { return pbcDistance( vec1, vec2 ); }
  else { return delta( vec1, vec2 ); }
}

void MultiColvarBase::buildCurrentTaskList( std::vector<unsigned>& tflags ) {
  tflags.assign(tflags.size(),1);
}

void MultiColvarBase::calculate(){
  // Set positions of all virtual atoms
  if( catom_indices.size()==0 ) {
      unsigned stride=comm.Get_size();
      unsigned rank=comm.Get_rank();
      if( runInSerial() ) { stride=1; rank=0; }
      std::vector<Vector> catomp( getFullNumberOfTasks() ); 
      double normaliz = 1. / static_cast<double>( ablocks.size() );
      for(unsigned i=rank;i<getFullNumberOfTasks();i+=stride){
          catomp[i].zero(); 
          for(unsigned j=0;j<ablocks.size();++j) catomp[i] += getPosition( ablocks[j][i] );
          catomp[i] *= normaliz;
      }
      if( !runInSerial() ) comm.Sum( catomp );
      for(unsigned i=0;i<getFullNumberOfTasks();++i) atoms.setVatomPosition( mygroup[i], catomp[i] );
  }
  runAllTasks();
}

void MultiColvarBase::performTask( const unsigned& task_index, MultiValue& myvals ) const {
  // Set the atoms pack up for the calculation 
  AtomValuePack myatoms( myvals, this ); myatoms.setNumberOfAtoms( ablocks.size() );
  for(unsigned i=0; i<ablocks.size(); ++i) myatoms.setAtom( i, ablocks[i][task_index] );
  // If we are using pbc make whole
  if(usepbc) myatoms.makeWhole();
  // And compute
  compute( task_index, myatoms ); 
  // Now update the active derivatives
  if( !doNotCalculateDerivatives() ) myatoms.updateUsingIndices();
}

void MultiColvarBase::apply(){
  if( doNotCalculateDerivatives() ) return;
  std::fill(forcesToApply.begin(),forcesToApply.end(),0);
  if( getForcesFromValues( forcesToApply ) ) setForcesOnAtoms( forcesToApply );
  
  // Virtual atom forces
  if( catom_indices.size()==0 ) {
      Tensor deriv; deriv = (1./static_cast<double>( ablocks.size() ))*Tensor::identity();
      unsigned stride=comm.Get_size();
      unsigned rank=comm.Get_rank();
      if( runInSerial() ) { stride=1; rank=0; }
      // Clear the forces
      for(unsigned i=0;i<getNumberOfAtoms();++i) vatom_forces[i].zero();
      // Accumulate the force on each virtual atom
      for(unsigned i=rank;i<getFullNumberOfTasks();i+=stride){
          Vector & f(atoms.getVatomForces(mygroup[i]));  
          //printf("FORCES %s %d %f %f %f \n",getLabel().c_str(), i,f[0],f[1],f[2]);
          for(unsigned j=0;j<ablocks.size();++j) vatom_forces[ablocks[j][i]] += matmul( deriv, f );
      }       
      if( !runInSerial() ) comm.Sum( vatom_forces ); 
      // Add the final forces to the atoms
      std::vector<Vector>& final_forces(modifyForces());
      for(unsigned i=0;i<final_forces.size();++i) final_forces[i] += vatom_forces[i]; 
      // Clear the forces on the virtual atoms
      for(unsigned i=0;i<getFullNumberOfTasks();++i) atoms.getVatomForces(mygroup[i]).zero(); 
  }
}

}
}
