/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2022,2023 The plumed team
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
#include "CLTool.h"
#include "CLToolRegister.h"
#include "tools/Tools.h"
#include "config/Config.h"
#include "core/ActionRegister.h"
#include "core/GenericMolInfo.h"
#include <cstdio>
#include <string>
#include <iostream>

namespace PLMD {
namespace cltools {

//+PLUMEDOC TOOLS gen_json
/*
gen_json constructs a json file that includes a dictionary of actions, the keywords for those actions and the components and outputs this to standard output

\par Examples

The following command generates the json file
\verbatim
plumed gen_json
\endverbatim


*/
//+ENDPLUMEDOC

class GenJson : public CLTool {
private:
  std::string version;
public:
  static void registerKeywords( Keywords& keys );
  explicit GenJson(const CLToolOptions& co );
  int main(FILE* in, FILE*out,Communicator& pc) override;
  std::string description()const override {
    return "print out a json file that contains the pluemd syntax";
  }
};

PLUMED_REGISTER_CLTOOL(GenJson,"gen_json")

void GenJson::registerKeywords( Keywords& keys ) {
  CLTool::registerKeywords( keys );
  keys.add("compulsory","--actions","a file containing one line descriptions of the various actions");
}

GenJson::GenJson(const CLToolOptions& co ):
  CLTool(co),
  version("master")
{
  inputdata=commandline;
  if( config::getVersionLong().find("dev")==std::string::npos ) version="v"+config::getVersion();
}

int GenJson::main(FILE* in, FILE*out,Communicator& pc) {
  std::string line(""), actionfile; parse("--actions",actionfile);
  IFile myfile; myfile.open(actionfile); bool stat;
  std::map<std::string,std::string> action_map;
  while((stat=myfile.getline(line))) {
    std::size_t col = line.find_first_of(":");
    action_map.insert(std::pair<std::string,std::string>( line.substr(0,col), line.substr(col+1) ) );
  }
  myfile.close();

  // Cycle over all the action names
  std::cout<<"{"<<std::endl;
  // Get the vimlink
  std::cout<<"  \"vimlink\" : \"https://www.plumed.org/doc-"<<version<<"/user-doc/html/_vim_syntax.html\","<<std::endl;
  // And the replicas link
  std::cout<<"  \"replicalink\" : \"https://www.plumed.org/doc-"<<version<<"/user-doc/html/special-replica-syntax.html\","<<std::endl;
  // Get the names of all the actions
  std::vector<std::string> action_names( actionRegister().getActionNames() );
  for(unsigned i=0; i<action_names.size(); ++i) {
    std::cout<<"  \""<<action_names[i]<<'"'<<": {"<<std::endl; std::string action=action_names[i];
    // Handle conversion of action names to links
    std::cout<<"    \"hyperlink\" : \"https://www.plumed.org/doc-"<<version<<"/user-doc/html/";
    std::transform(action.begin(), action.end(), action.begin(), [](unsigned char c) { return std::tolower(c); });
    while(true) {
      std::size_t und=action.find_first_of("_");
      if( und==std::string::npos ) break;
      std::string first=action.substr(0,und);
      for(auto c : first ) { if( isdigit(c) ) std::cout<<c; else std::cout<<"_"<<c; }
      std::cout<<"_"; action=action.substr(und+1);
    }
    for(auto c : action ) { if( isdigit(c) ) std::cout<<c; else std::cout<<"_"<<c; }
    std::cout<<".html\","<<std::endl;
    std::cout<<"    \"description\" : \""<<action_map[action_names[i]]<<"\",\n";
    // Now output keyword information
    Keywords keys; actionRegister().getKeywords( action_names[i], keys );
    std::cout<<"    \"syntax\" : {"<<std::endl;
    for(unsigned j=0; j<keys.size(); ++j) {
      std::string desc = keys.getKeywordDescription( keys.getKeyword(j) );
      if( desc.find("default=")!=std::string::npos ) {
        std::size_t brac=desc.find_first_of(")"); desc = desc.substr(brac+1);
      }
      std::size_t dot=desc.find_first_of(".");
      std::cout<<"       \""<<keys.getKeyword(j)<<"\" : { \"type\": \""<<keys.getStyle(keys.getKeyword(j))<<"\", \"description\": \""<<desc.substr(0,dot)<<"\", \"multiple\": "<<keys.numbered( keys.getKeyword(j) )<<"}";
      if( j==keys.size()-1 && !keys.exists("HAS_VALUES") ) std::cout<<std::endl; else std::cout<<","<<std::endl;
    }
    if( keys.exists("HAS_VALUES") ) {
      std::cout<<"       \"output\" : {"<<std::endl;
      std::vector<std::string> components( keys.getOutputComponents() );
      // Check if we have a value
      bool hasvalue=true;
      for(unsigned k=0; k<components.size(); ++k) {
        if( keys.getOutputComponentFlag( components[k] )=="default" ) { hasvalue=false; break; }
      }
      if( hasvalue ) {
        std::cout<<"         \"value\": {"<<std::endl;
        std::cout<<"           \"flag\": \"value\","<<std::endl;
        std::cout<<"           \"description\": \"a scalar quantity\""<<std::endl;
        if( components.size()==0 ) std::cout<<"         }"<<std::endl; else std::cout<<"         },"<<std::endl;
      }
      for(unsigned k=0; k<components.size(); ++k) {
        std::cout<<"         \""<<components[k]<<"\" : {"<<std::endl;
        std::cout<<"           \"flag\": \""<<keys.getOutputComponentFlag( components[k] )<<"\","<<std::endl;
        std::string desc=keys.getOutputComponentDescription( components[k] ); std::size_t dot=desc.find_first_of(".");
        std::cout<<"           \"description\": \""<<desc.substr(0,dot)<<"\""<<std::endl;
        if( k==components.size()-1 ) std::cout<<"         }"<<std::endl; else std::cout<<"         },"<<std::endl;
      }
      std::cout<<"       }"<<std::endl;

    }
    std::cout<<"    },"<<std::endl;
    // This ensures that \n is replaced by \\n
    std::string unsafen="\n", safen="\\n", helpstr = keys.getHelpString();
    for( std::size_t pos = helpstr.find("\n");
         pos != std::string::npos;
         pos = helpstr.find("\n", pos)
       ) { helpstr.replace(pos, unsafen.size(), safen); pos += safen.size(); }
    std::cout<<"    \"help\" : \""<<helpstr<<"\"\n";
    std::cout<<"  },"<<std::endl;
  }
  // Get all the special groups
  std::cout<<"  \"groups\" : {"<<std::endl;
  std::cout<<"    \"@allatoms\" : { \n"<<std::endl;
  std::cout<<"        \"description\" : \"refers to all the MD codes atoms and PLUMEDs vatoms\","<<std::endl;
  std::cout<<"        \"link\" : \"https://www.plumed.org/doc-"<<version<<"/user-doc/html/_group.html\""<<std::endl;
  std::cout<<"    },"<<std::endl;
  std::cout<<"    \"@mdatoms\" : { \n"<<std::endl;
  std::cout<<"        \"description\" : \"refers to all the MD codes atoms but not PLUMEDs vatoms\","<<std::endl;
  std::cout<<"        \"link\" : \"https://www.plumed.org/doc-"<<version<<"/user-doc/html/_group.html\""<<std::endl;
  // Now print all the special keywords in molinfo
  std::map<std::string,std::string> specials( GenericMolInfo::getSpecialKeywords() );
  for(auto const& s : specials ) {
    std::cout<<"    },"<<std::endl;
    std::cout<<"    \""<<s.first<<"\" : { \n"<<std::endl;
    std::cout<<"        \"description\" : \""<<s.second<<"\","<<std::endl;
    std::cout<<"        \"link\" : \"https://www.plumed.org/doc-"<<version<<"/user-doc/html/_m_o_l_i_n_f_o.html\""<<std::endl;
  }
  std::cout<<"        }"<<std::endl;
  std::cout<<"  }"<<std::endl;
  std::cout<<"}"<<std::endl;
  return 0;
}

} // End of namespace
}
