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
#include "RDF.h"
#include "KDE.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace gridtools {

PLUMED_REGISTER_ACTION(RDF,"RDF_CALC")

void RDF::createX2ReferenceObject( const std::string& lab, const std::string& grid_setup, const bool& calc_dens, ActionShortcut* action ) {
  // Create grid with normalizing function
  action->readInputLine( lab  + "_x2: REFERENCE_FUNCTION PERIODIC=NO FUNC=x*x " + grid_setup );
  // Compute density if required
  if( calc_dens ) action->readInputLine( lab + "_vol: VOLUME" );
}

void RDF::registerKeywords( Keywords& keys ) {
  KDE::registerKeywords( keys ); keys.remove("GRID_MIN"); keys.remove("GRID_MAX");
  keys.add("atoms","GROUP","");
  keys.add("atoms-2","GROUPA","");
  keys.add("atoms-2","GROUPB","");
  keys.add("compulsory","MAXR","the maximum distance to use for the rdf"); 
  keys.add("compulsory","BANDWIDTH","the bandwidths for kernel density esimtation");
  keys.add("optional","DENSITY","the reference density to use when normalizing the RDF");
  keys.add("hidden","REFERENCE","this is the label of the reference objects");
}

RDF::RDF(const ActionOptions&ao):
Action(ao),
ActionShortcut(ao)
{
  // Read in grid extent and number of bins
  std::string maxr, nbins, dens; parse("MAXR",maxr); parse("GRID_BIN",nbins); parse("DENSITY",dens);
  std::string grid_setup = "GRID_MIN=0 GRID_MAX=" + maxr + " GRID_BIN=" + nbins;
  // Create grid with normalizing function on it 
  std::string refstr; parse("REFERENCE",refstr); 
  if( refstr.length()==0 ){
      createX2ReferenceObject( getShortcutLabel(), grid_setup, dens.length()==0, this ); refstr = getShortcutLabel();
  }
  // Read input to histogram
  std::string cutoff, kernel, bandwidth, kernel_data; parse("KERNEL",kernel); 
  if( kernel=="DISCRETE" ) { 
    cutoff = maxr; kernel_data="KERNEL=DISCRETE"; 
  } else { 
    parse("BANDWIDTH",bandwidth); double rcut; parse("CUTOFF",rcut); kernel_data="KERNEL=" + kernel + " BANDWIDTH=" + bandwidth;
    double bw; Tools::convert( bandwidth, bw ); double fcut; Tools::convert( maxr, fcut ); Tools::convert( fcut + sqrt(2.0*rcut)*bw, cutoff );
  }

  // Create contact matrix
  std::string natoms, str_norm_atoms, atom_str, group_str, groupa_str, groupb_str; parse("GROUP",group_str); 
  if( group_str.length()>0 ) { 
    atom_str="GROUP=" + group_str; std::vector<std::string> awords=Tools::getWords(group_str,"\t\n ,"); 
    Tools::interpretRanges( awords ); Tools::convert( awords.size(), natoms ); str_norm_atoms = natoms;
  } else { 
    parse("GROUPA",groupa_str); parse("GROUPB",groupb_str); 
    std::vector<std::string> awords=Tools::getWords(groupb_str,"\t\n ,"); 
    Tools::interpretRanges( awords ); Tools::convert( awords.size(), natoms ); 
    atom_str="GROUPA=" + groupa_str + " GROUPB=" + groupb_str; 
    std::vector<std::string> bwords=Tools::getWords(groupa_str,"\t\n ,"); Tools::interpretRanges( bwords );
    Tools::convert( bwords.size()+1, str_norm_atoms );
  }
  // Retrieve the number of atoms
  readInputLine( getShortcutLabel() + "_mat: DISTANCE_MATRIX CUTOFF=" + cutoff + " " + atom_str);

  // Calculate weights of distances 
  readInputLine( getShortcutLabel() + "_wmat: MATHEVAL ARG1=" + getShortcutLabel() + "_mat.w FUNC=step(" + cutoff + "-x) PERIODIC=NO");
  // Now create a histogram from the contact matrix
  readInputLine( getShortcutLabel() + "_kde: KDE ARG1=" + getShortcutLabel() + "_mat.w HEIGHTS=" + getShortcutLabel() + "_wmat UNORMALIZED " + grid_setup + " " + kernel_data); 
  // Transform the histogram by normalizing factor for rdf
  readInputLine( getShortcutLabel() + "_vrdf: MATHEVAL ARG1=" + getShortcutLabel() + "_kde ARG2=" + refstr + "_x2 FUNC=x/(2*pi*y) PERIODIC=NO");
  // And normalize by density and number of atoms (separated from above to avoid nans)
  std::string func_str; 
  if( dens.length()>0 ) func_str = "FUNC=x/(" + dens + "*" + str_norm_atoms + ")";
  else func_str = "ARG2=" + refstr + "_vol FUNC=x*y/(" + natoms + "*" + str_norm_atoms + ")";
  readInputLine( getShortcutLabel() + "_rdf: MATHEVAL ARG1=" + getShortcutLabel() + "_vrdf PERIODIC=NO " + func_str);
}

class RDFShortcut : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit RDFShortcut(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(RDFShortcut,"RDF")

void RDFShortcut::registerKeywords( Keywords& keys ) {
  RDF::registerKeywords( keys ); HistogramBase::histogramKeywords( keys );
}

RDFShortcut::RDFShortcut(const ActionOptions&ao):
Action(ao),
ActionShortcut(ao)
{
    // Parse the keymap
    std::map<std::string,std::string> keymap; HistogramBase::readHistogramKeywords( keymap, this );
    // Make the kde object
    readInputLine( getShortcutLabel() + "_kde: RDF_CALC " + convertInputLineToString() );
    // And the averaging object
    HistogramBase::createAveragingObject( getShortcutLabel() + "_rdf", getShortcutLabel(), keymap, this );
}

}
}
