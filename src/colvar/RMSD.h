/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2020 The plumed team
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
#ifndef __PLUMED_colvar_RMSD_h
#define __PLUMED_colvar_RMSD_h

#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"
#include "tools/RMSD.h"

namespace PLMD {

class ActionShortcut;

namespace colvar {

class RMSD : 
public ActionWithArguments,
public ActionWithValue {
private:
  bool firsttime, fixed_reference;
  bool squared;
  bool displacement;
  bool norm_weights;
  bool multiple;
  std::string type;
  std::vector<double> align,displace,sqrtdisplace;
  std::vector<PLMD::RMSD> myrmsd;
  std::vector<double> forcesToApply;
  void setReferenceConfiguration( const unsigned& jconf );
public:
  explicit RMSD(const ActionOptions&);
  unsigned getNumberOfDerivatives() const override;
  unsigned getNumberOfColumns() const override;
  void calculate() override;
  void performTask( const unsigned& task_index, MultiValue& myvals ) const override;
  bool performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const ;
  void apply() override;
  void update() override;
  static void registerRMSD(Keywords& keys );
  static void registerKeywords(Keywords& keys);
  static void createReferenceConfiguration( const std::string& lab, const std::string& input, PlumedMain& plumed, const unsigned number=0 ); 
  static void createPosVector( const std::string& lab, const PDB& pdb, ActionShortcut* action );
};

inline
unsigned RMSD::getNumberOfDerivatives() const {
  return 3*align.size();
}

inline 
unsigned RMSD::getNumberOfColumns() const {
  return 3*align.size();  
}

}
}

#endif

