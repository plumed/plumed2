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
#include "KDE.h"
#include "RDF.h"
#include "core/ActionRegister.h"
#include "core/ActionShortcut.h"

namespace PLMD {
namespace gridtools {

class PairEntropies : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit PairEntropies(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(PairEntropies,"PAIRENTROPIES")

void PairEntropies::registerKeywords( Keywords& keys ) {
  KDE::registerKeywords( keys ); 
  keys.add("atoms","ATOMS","");
  keys.add("compulsory","MAXR","the maximum distance to use for the rdf");
  keys.add("optional","DENSITY","the reference density to use for the pair entropy");
}

PairEntropies::PairEntropies(const ActionOptions&ao):
Action(ao),
ActionShortcut(ao)
{
   // Read in the atoms involved
   std::string atoms_str; parse("ATOMS",atoms_str); 
   // Create the x2 value
   std::string maxr, nbins, dens; parse("MAXR",maxr); parse("GRID_BIN",nbins); parse("DENSITY",dens);
   std::string grid_setup = "GRID_MIN=0 GRID_MAX=" + maxr + " GRID_BIN=" + nbins;
   RDF::createX2ReferenceObject( getShortcutLabel(), grid_setup, dens.length()==0, this );
   // Create the number of input atoms
   std::string pair_str = "GRID_BIN=" + nbins + " MAXR=" + maxr + " " + convertInputLineToString(); 
   std::vector<std::string> awords=Tools::getWords(atoms_str,"\t\n ,"); Tools::interpretRanges( awords );
   // Now create all the pairentropy objects
   for(unsigned i=0;i<awords.size();++i) {
       std::string ilab, jlab, atoms_str2; Tools::convert( awords[i], ilab ); 
       for(unsigned j=0;j<awords.size();++j) {
           if( awords[i]==awords[j] ) continue; Tools::convert( awords[j], jlab );
           if( atoms_str2.length()==0 ) atoms_str2 = jlab; else atoms_str2 += "," + jlab;
       }
       readInputLine( getShortcutLabel() + ilab + ": PAIRENTROPY GROUPA=" + ilab + " GROUPB=" + atoms_str2 + " " + pair_str + " REFERENCE=" + getShortcutLabel() );
   }
   // And compose a vector containing the values of all the pair entropies
   std::string num, argstr = "ARG=" + getShortcutLabel() + "1"; 
   for(unsigned i=1;i<awords.size();++i){ Tools::convert( i+1, num ); argstr += "," + getShortcutLabel() + num; }
   readInputLine( getShortcutLabel() + ": CONCATENATE " + argstr ); 
}

}
}
