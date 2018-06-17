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
#include "Function.h"
#include "ActionRegister.h"
#include "tools/PDB.h"

#include <cmath>

using namespace std;

namespace PLMD {
namespace function {

//+PLUMEDOC FUNCTION DIFFERENCE
/*
Use a switching function to determine how many of the input variables are less than a certain cutoff.

\par Examples

*/
//+ENDPLUMEDOC


class Difference :
  public Function
{
public:
  static void shortcutKeywords( Keywords& keys );
  static void expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                              const std::map<std::string,std::string>& keys,
                              std::vector<std::vector<std::string> >& actions );
  explicit Difference(const ActionOptions&);
  void calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const ;
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(Difference,"DIFFERENCE")
PLUMED_REGISTER_SHORTCUT(Difference,"DRMSD")

void Difference::shortcutKeywords( Keywords& keys ) {
  keys.add("compulsory","REFERENCE","a file in pdb format containing the reference structure and the atoms involved in the CV.");
  keys.add("compulsory","LOWER_CUTOFF","only pairs of atoms further than LOWER_CUTOFF are considered in the calculation.");
  keys.add("compulsory","UPPER_CUTOFF","only pairs of atoms closer than UPPER_CUTOFF are considered in the calculation.");
  keys.add("compulsory","TYPE","DRMSD","what kind of DRMSD would you like to calculate.  You can use either the normal DRMSD involving all the distances between "
           "the atoms in your molecule.  Alternatively, if you have multiple molecules you can use the type INTER-DRMSD "
           "to compute DRMSD values involving only those distances between the atoms at least two molecules or the type INTRA-DRMSD "
           "to compute DRMSD values involving only those distances between atoms in the same molecule");
  keys.addFlag("SQUARED",false,"This should be setted if you want MSD instead of RMSD ");
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
}

void Difference::expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                                 const std::map<std::string,std::string>& keys,
                                 std::vector<std::vector<std::string> >& actions ) {
  if( words[0]=="DRMSD" ) {
      // Read in the reference configuration
      std::vector<std::string> ref_str; ref_str.push_back(lab + "_atoms:"); ref_str.push_back("READ_ATOMS");
      ref_str.push_back("REFERENCE=" + keys.find("REFERENCE")->second ); actions.push_back( ref_str );
      // First bit of input for reference values
      std::vector<std::string> ref_vals; ref_vals.push_back(lab + "_ref:"); ref_vals.push_back("CALCULATE_REFERENCE");
      ref_vals.push_back("ATOMS=" + lab + "_atoms"); std::string din = "INPUT=DISTANCE NOPBC";
      // First bit of input for the instantaneous distances
      std::vector<std::string> mat_inp; mat_inp.push_back( lab + "_mat:"); mat_inp.push_back("DISTANCE");
      if( keys.count("NOPBC") ) mat_inp.push_back("NOPBC");
      // Get cutoff information
      double lcut, ucut=std::numeric_limits<double>::max();
      if( keys.count("LOWER_CUTOFF") ) Tools::convert( keys.find("LOWER_CUTOFF")->second, lcut );
      if( keys.count("UPPER_CUTOFF") ) Tools::convert( keys.find("UPPER_CUTOFF")->second, ucut ); 
      PDB pdb;  std::string reference = keys.find("REFERENCE")->second;
      // This is a potential problem if the user is using units
      if( !pdb.read(reference,false,0.1) ) plumed_merror("missing input file " + reference ); 
      std::vector<AtomNumber> atoms( pdb.getAtomNumbers() ); std::vector<Vector> pos( pdb.getPositions() );
      // Work out what distances we need to calculate from the reference configuration
      std::string drmsd_type = keys.find("TYPE")->second; unsigned nn=1; std::string istr, jstr, num;
      if( drmsd_type=="DRMSD" ) {
          for(unsigned i=0;i<atoms.size()-1;++i) {
              Tools::convert( atoms[i].serial(), istr );
              for(unsigned j=i+1; j<atoms.size(); ++j) {
                  Tools::convert( atoms[j].serial(), jstr );
                  double distance = delta( pos[i], pos[j] ).modulo();
                  if( distance < ucut && distance > lcut ) { 
                      Tools::convert( nn, num );
                      std::string dinstr = "ATOMS" + num + "=" + istr + "," + jstr;
                      mat_inp.push_back( dinstr ); din += " " + dinstr; nn++;
                  }
              }
          }
      } else { 
          unsigned nblocks = pdb.getNumberOfAtomBlocks(); std::vector<unsigned> blocks( nblocks+1 );
          if( nblocks==1 ) plumed_merror("Trying to compute intermolecular rmsd but found no TERs in input PDB");
          blocks[0]=0; for(unsigned i=0; i<nblocks; ++i) blocks[i+1]=pdb.getAtomBlockEnds()[i];
          if( drmsd_type=="INTRA-DRMSD" ) { 
              for(unsigned i=0; i<nblocks; ++i) {
                for(unsigned iatom=blocks[i]+1; iatom<blocks[i+1]; ++iatom) {
                  Tools::convert( atoms[iatom].serial(), istr );
                  for(unsigned jatom=blocks[i]; jatom<iatom; ++jatom) {
                    Tools::convert( atoms[jatom].serial(), jstr );
                    double distance = delta( pos[iatom], pos[jatom] ).modulo();
                    if(distance < ucut && distance > lcut ) {
                       Tools::convert( nn, num );
                       std::string dinstr = "ATOMS" + num + "=" + istr + "," + jstr;
                       mat_inp.push_back( dinstr ); din += " " + dinstr; nn++;
                    }
                  }
                }
              }
          } else if( drmsd_type=="INTER-DRMSD" ) {
              for(unsigned i=1; i<nblocks; ++i) {
                for(unsigned j=0; j<i; ++j) {
                  for(unsigned iatom=blocks[i]; iatom<blocks[i+1]; ++iatom) {
                    Tools::convert( atoms[iatom].serial(), istr );
                    for(unsigned jatom=blocks[j]; jatom<blocks[j+1]; ++jatom) {
                      Tools::convert( atoms[jatom].serial(), jstr );
                      double distance = delta( pos[iatom], pos[jatom] ).modulo();
                      if(distance < ucut && distance > lcut ) {
                         Tools::convert( nn, num );
                         std::string dinstr = "ATOMS" + num + "=" + istr + "," + jstr;
                         mat_inp.push_back( dinstr ); din += " " + dinstr; nn++; 
                      }
                    }
                  }
                } 
              }
          } else plumed_error();
      }
      // Put this information into the reference matrix
      ref_vals.push_back( din ); actions.push_back( ref_vals ); actions.push_back( mat_inp );
      // And the difference between these two sets of matrices
      std::vector<std::string> diff_mat; diff_mat.push_back( lab + "_diffm:"); diff_mat.push_back("DIFFERENCE");
      diff_mat.push_back("ARG1=" + lab + "_mat"); diff_mat.push_back("ARG2=" + lab + "_ref");
      actions.push_back( diff_mat );
      // And the total difference
      std::vector<std::string> comb_inp; 
      if( keys.count("SQUARED") ) comb_inp.push_back( lab + ":"); else comb_inp.push_back( lab + "_2:"); 
      comb_inp.push_back("COMBINE"); comb_inp.push_back("ARG=" + lab + "_diffm"); std::string powstr ="POWERS=2";
      for(unsigned i=1;i<nn-1;++i) powstr += ",2"; comb_inp.push_back( powstr );
      comb_inp.push_back("NORMALIZE"); comb_inp.push_back("PERIODIC=NO"); actions.push_back( comb_inp );
      // And the square root of the distance if required
      if( !keys.count("SQUARED") ) {
          std::vector<std::string> fval; fval.push_back( lab + ":" ); fval.push_back("MATHEVAL");
          fval.push_back("ARG=" + lab + "_2"); fval.push_back("FUNC=sqrt(x)"); fval.push_back("PERIODIC=NO");
          actions.push_back( fval );
      }
  }
}

void Difference::registerKeywords(Keywords& keys) {
  Function::registerKeywords(keys); keys.use("ARG");
}

Difference::Difference(const ActionOptions&ao):
  Action(ao),
  Function(ao)
{
  if( getNumberOfArguments()!=2 ) error("difference can only take two arguments as input");
  if( getPntrToArgument(0)->isPeriodic() ) {
      if( !getPntrToArgument(1)->isPeriodic() ) error("period for input variables should be the same"); 
      std::string min0, max0; getPntrToArgument(0)->getDomain( min0, max0 );
      std::string min1, max1; getPntrToArgument(1)->getDomain( min1, max1 );
      if( min0!=max0 || max0!=max1 ) error("period for input variables should be the same");
  } else if( getPntrToArgument(1)->isPeriodic() ) error("period for input variables should be the same");

  getPeriodFromArg=true; addValueWithDerivatives();
  checkRead();
}

void Difference::calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const {
  plumed_dbg_assert( args.size()==2 ); 
  addValue( 0, getPntrToArgument(0)->difference( args[1], args[0] ), myvals ); 
  addDerivative( 0, 0, 1.0, myvals ); addDerivative( 0, 1, -1.0, myvals );
}

}
}


