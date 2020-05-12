/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2018 The VES code team
   (see the PEOPLE-VES file at the root of this folder for a list of names)

   See http://www.ves-code.org for more information.

   This file is part of VES code module.

   The VES code module is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   The VES code module is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with the VES code module.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#include "VesTools.h"

#include "tools/Grid.h"
#include "tools/IFile.h"
#include "tools/Exception.h"


namespace PLMD {
namespace ves {


void VesTools::copyGridValues(GridBase* grid_pntr_orig, GridBase* grid_pntr_copy) {
  // plumed_massert(grid_pntr_orig!=NULL,"grid not defined");
  // plumed_massert(grid_pntr_copy!=NULL,"grid not defined");
  plumed_massert(grid_pntr_orig->getSize()==grid_pntr_copy->getSize(),"the two grids are not of the same size");
  plumed_massert(grid_pntr_orig->getDimension()==grid_pntr_copy->getDimension(),"the two grids are not of the same dimension");
  //
  for(Grid::index_t i=0; i<grid_pntr_orig->getSize(); i++) {
    double value = grid_pntr_orig->getValue(i);
    grid_pntr_copy->setValue(i,value);
  }
}


unsigned int VesTools::getGridFileInfo(const std::string& filepath, std::string& grid_label, std::vector<std::string>& arg_labels, std::vector<std::string>& arg_min, std::vector<std::string>& arg_max, std::vector<bool>& arg_periodic, std::vector<unsigned int>& arg_nbins, bool& derivatives) {

  IFile ifile; ifile.open(filepath);
  std::vector<std::string> fields;
  ifile.scanFieldList(fields);
  ifile.allowIgnoredFields();
  ifile.scanField();

  unsigned int nargs=0;
  for(unsigned int i=0; i<fields.size(); i++) {
    if(fields[i]=="min_"+fields[0]) {
      derivatives = false;
      nargs = i-1;
      break;
    }
    else if(fields[i]=="der_"+fields[0]) {
      derivatives = true;
      nargs = i-1;
      break;
    }
  }

  grid_label = fields[nargs];

  arg_labels.assign(nargs,"");
  arg_min.assign(nargs,"");
  arg_max.assign(nargs,"");
  arg_periodic.assign(nargs,false);
  arg_nbins.assign(nargs,0);
  for(unsigned int i=0; i<nargs; i++) {
    arg_labels[i] = fields[i];
    ifile.scanField("min_"+arg_labels[i],arg_min[i]);
    ifile.scanField("max_"+arg_labels[i],arg_max[i]);
    std::string str_periodic;
    ifile.scanField("periodic_"+arg_labels[i],str_periodic);
    if(str_periodic=="true") {arg_periodic[i]=true;}
    int nbins;
    ifile.scanField("nbins_"+arg_labels[i],nbins);
    arg_nbins[i] = static_cast<unsigned int>(nbins);
  }
  ifile.scanField();
  ifile.close();
  return nargs;
}


}
}
