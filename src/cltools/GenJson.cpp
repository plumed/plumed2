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
#include "core/CLToolRegister.h"
#include "tools/Tools.h"
#include "config/Config.h"
#include "core/ActionRegister.h"
#include "core/CLToolRegister.h"
#include "core/GenericMolInfo.h"
#include "core/ModuleMap.h"
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
  void printHyperlink( const std::string& action );
  void printKeywordDocs( const std::string& k, const std::string& mydescrip, const Keywords& keys );
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
  version("master") {
  inputdata=commandline;
  if( config::getVersionLong().find("dev")==std::string::npos ) {
    version="v"+config::getVersion();
  }
}

void GenJson::printHyperlink( const std::string& action ) {
  std::cout<<"    \"hyperlink\" : \"https://www.plumed.org/doc-"<<version<<"/user-doc/html/"<<action<<"\","<<std::endl;
}

void GenJson::printKeywordDocs( const std::string& k, const std::string& mydescrip, const Keywords& keys ) {
  std::cout<<"       \""<<k<<"\" : { \"type\": \""<<keys.getStyle(k)<<"\", \"description\": \""<<mydescrip<<"\", \"multiple\": "<<keys.numbered(k)<<", \"actionlink\": \""<<keys.getLinkedActions(k)<<"\"";
}

int GenJson::main(FILE* in, FILE*out,Communicator& pc) {
  std::string line(""), actionfile;
  parse("--actions",actionfile);
  IFile myfile;
  myfile.open(actionfile);
  bool stat;
  std::map<std::string,std::string> action_map;
  while((stat=myfile.getline(line))) {
    std::size_t col = line.find_first_of(":");
    std::string docs = line.substr(col+1);
    if( docs.find("\\")!=std::string::npos ) {
      error("found invalid backslash character in first line of documentation for action " + line.substr(0,col) );
    }
    action_map.insert(std::pair<std::string,std::string>( line.substr(0,col), docs ) );
  }
  myfile.close();

  // Cycle over all the action names
  std::cout<<"{"<<std::endl;
  // Get the vimlink
  std::cout<<"  \"vimlink\" : \"https://www.plumed.org/doc-"<<version<<"/user-doc/html/vim\","<<std::endl;
  // And the replicas link
  std::cout<<"  \"replicalink\" : \"https://www.plumed.org/doc-"<<version<<"/user-doc/html/parsing.html\","<<std::endl;
  // Get the names of all the actions
  std::vector<std::string> action_names( actionRegister().getActionNames() );
  std::vector<std::string> allmodules;
  for(unsigned i=0; i<action_names.size(); ++i) {
    std::cout<<"  \""<<action_names[i]<<'"'<<": {"<<std::endl;
    std::string action=action_names[i];
    // Handle conversion of action names to links
    printHyperlink( action );
    std::cout<<"    \"description\" : \""<<action_map[action_names[i]]<<"\",\n";
    bool found=false;
    std::string thismodule = getModuleMap().find(action_names[i])->second;
    for(unsigned i=0; i<allmodules.size(); ++i) {
      if( allmodules[i]==thismodule ) {
        found=true;
        break;
      }
    }
    if( !found ) {
      allmodules.push_back( thismodule );
    }
    std::cout<<"    \"module\" : \""<<thismodule<<"\",\n";
    // Now output keyword information
    Keywords keys;
    actionRegister().getKeywords( action_names[i], keys );
    std::cout<<"    \"displayname\" : \""<<keys.getDisplayName()<<"\",\n";
// This is used for noting actions that have been deprecated
    std::string replacement = keys.getReplacementAction();
    if( replacement!="none" ) {
      bool found_replacement=false;
      for(unsigned j=0; j<action_names.size(); ++j) {
        if( action_names[j]==replacement ) {
          found_replacement=true;
          break;
        }
      }
      if( !found_replacement ) {
        error("could not find action named " + replacement + " that is supposed to be used to replace " + action_names[i] );
      }
      std::cout<<"    \"replacement\" : \""<<replacement<<"\",\n";
    }
    std::cout<<"     \"dois\" : [";
    unsigned ndoi = keys.getDOIList().size();
    if( ndoi>0 ) {
      std::cout<<"\"" + keys.getDOIList()[0] + "\"";
      for(unsigned j=1; j<ndoi; ++j) {
        std::cout<<", \"" + keys.getDOIList()[j] + "\"";
      }
    }
    std::cout<<"],\n";
    std::cout<<"    \"syntax\" : {"<<std::endl;
    for(unsigned j=0; j<keys.size(); ++j) {
      std::string defa = "", desc = keys.getKeywordDescription( keys.getKeyword(j) );
      if( desc.find("default=")!=std::string::npos ) {
        std::size_t defstart = desc.find_first_of("="), brac=desc.find_first_of(")");
        defa = desc.substr(defstart+1,brac-defstart-2);
        desc = desc.substr(brac+1);
      }
      std::size_t dot=desc.find_first_of(".");
      std::string mydescrip = desc.substr(0,dot);
      if( mydescrip.find("\\")!=std::string::npos ) {
        error("found invalid backslash character documentation for keyword " + keys.getKeyword(j) + " in action " + action_names[i] );
      }
      std::string argtype = keys.getArgumentType( keys.getKeyword(j) );
      if( argtype.length()>0 ) {
        printKeywordDocs( keys.getKeyword(j), mydescrip, keys );
        std::cout<<", \"argtype\": \""<<argtype<<"\"}";
      } else if( defa.length()>0 ) {
        printKeywordDocs( keys.getKeyword(j), mydescrip, keys );
        std::cout<<", \"default\": \""<<defa<<"\"}";
      } else {
        printKeywordDocs( keys.getKeyword(j), mydescrip, keys );
        std::cout<<"}";
      }
      if( j==keys.size()-1 && !keys.exists("HAS_VALUES") ) {
        std::cout<<std::endl;
      } else {
        std::cout<<","<<std::endl;
      }
    }
    if( keys.exists("HAS_VALUES") ) {
      std::cout<<"       \"output\" : {"<<std::endl;
      std::vector<std::string> components( keys.getOutputComponents() );
      // Check if we have a value
      bool hasvalue=true;
      for(unsigned k=0; k<components.size(); ++k) {
        if( keys.getOutputComponentFlag( components[k] )=="default" ) {
          hasvalue=false;
          break;
        }
      }
      for(unsigned k=0; k<components.size(); ++k) {
        std::string compname=components[k];
        if( components[k]==".#!value" ) {
          hasvalue=false;
          compname="value";
        }
        std::cout<<"         \""<<compname<<"\" : {"<<std::endl;
        std::cout<<"           \"flag\": \""<<keys.getOutputComponentFlag( components[k] )<<"\","<<std::endl;
        std::cout<<"           \"type\": \""<<keys.getOutputComponentType( components[k] )<<"\","<<std::endl;
        std::string desc=keys.getOutputComponentDescription( components[k] );
        std::size_t dot=desc.find_first_of(".");
        std::string mydescrip = desc.substr(0,dot);
        if( mydescrip.find("\\")!=std::string::npos ) {
          error("found invalid backslash character documentation for output component " + compname + " in action " + action_names[i] );
        }
        std::cout<<"           \"description\": \""<<mydescrip<<"\""<<std::endl;
        if( k==components.size()-1 ) {
          std::cout<<"         }"<<std::endl;
        } else {
          std::cout<<"         },"<<std::endl;
        }
      }
      if( hasvalue && components.size()==0 ) {
        printf("WARNING: no components have been registered for action %s \n", action_names[i].c_str() );
      }
      std::cout<<"       }"<<std::endl;

    }
    std::cout<<"    },"<<std::endl;
    if( keys.getNeededKeywords().size()>0 ) {
      std::vector<std::string> neededActions( keys.getNeededKeywords() );
      std::cout<<"    \"needs\" : ["<<"\""<<neededActions[0]<<"\"";
      for(unsigned j=1; j<neededActions.size(); ++j) {
        std::cout<<", \""<<neededActions[j]<<"\"";
      }
      std::cout<<"],"<<std::endl;
    }
    // This ensures that \n is replaced by \\n
    std::string unsafen="\n", safen="\\n", helpstr = keys.getHelpString();
    for( std::size_t pos = helpstr.find("\n");
         pos != std::string::npos;
         pos = helpstr.find("\n", pos)
       ) {
      helpstr.replace(pos, unsafen.size(), safen);
      pos += safen.size();
    }
    std::cout<<"    \"help\" : \""<<helpstr<<"\"\n";
    std::cout<<"  },"<<std::endl;
  }
  // Get all the special groups
  std::cout<<"  \"groups\" : {"<<std::endl;
  std::cout<<"    \"@allatoms\" : { \n"<<std::endl;
  std::cout<<"        \"description\" : \"refers to all the MD codes atoms and PLUMEDs vatoms\","<<std::endl;
  std::cout<<"        \"link\" : \"https://www.plumed.org/doc-"<<version<<"/user-doc/html/specifying_atoms\""<<std::endl;
  std::cout<<"    },"<<std::endl;
  std::cout<<"    \"@mdatoms\" : { \n"<<std::endl;
  std::cout<<"        \"description\" : \"refers to all the MD codes atoms but not PLUMEDs vatoms\","<<std::endl;
  std::cout<<"        \"link\" : \"https://www.plumed.org/doc-"<<version<<"/user-doc/html/specifying_atoms\""<<std::endl;
  std::cout<<"    },"<<std::endl;
  std::cout<<"    \"@ndx:\" : { \n"<<std::endl;
  std::cout<<"        \"description\" : \"load a group from a GROMACS index file\","<<std::endl;
  std::cout<<"        \"link\" : \"https://www.plumed.org/doc-"<<version<<"/user-doc/html/specifying_atoms.html\""<<std::endl;
  // Now print all the special keywords in molinfo
  std::map<std::string,std::string> specials( GenericMolInfo::getSpecialKeywords() );
  for(auto const& s : specials ) {
    std::cout<<"    },"<<std::endl;
    std::cout<<"    \""<<s.first<<"\" : { \n"<<std::endl;
    std::cout<<"        \"description\" : \""<<s.second<<"\","<<std::endl;
    std::cout<<"        \"link\" : \"https://www.plumed.org/doc-"<<version<<"/user-doc/html/MOLINFO\""<<std::endl;
  }
  std::cout<<"        }"<<std::endl;
  std::cout<<"  },"<<std::endl;
  std::cout<<"  \"modules\" : {"<<std::endl;
  std::cout<<"    \""<<allmodules[0]<<"\" : { "<<std::endl;
  printHyperlink( allmodules[0] );
  std::cout<<"        \"description\" : \"A module that will be used for something\""<<std::endl;
  for(unsigned i=1; i<allmodules.size(); ++i) {
    std::cout<<"    },"<<std::endl;
    std::cout<<"    \""<<allmodules[i]<<"\" : { "<<std::endl;
    printHyperlink( allmodules[i] );
    std::cout<<"        \"description\" : \"A module that will be used for something\""<<std::endl;
  }
  std::cout<<"        }"<<std::endl;
  std::cout<<"  },"<<std::endl;
  std::cout<<"  \"cltools\" : {"<<std::endl;
  std::vector<std::string> cltool_names( cltoolRegister().list() );
  for(unsigned i=0; i<cltool_names.size(); ++i) {
    std::cout<<"  \""<<cltool_names[i]<<'"'<<": {"<<std::endl;
    std::string cltool=cltool_names[i];
    // Handle converstion to link
    printHyperlink( cltool );
    std::string thismodule = getModuleMap().find(cltool_names[i])->second;
    std::cout<<"    \"module:\" : \""<<thismodule<<"\",\n";
    auto mytool = cltoolRegister().create( CLToolOptions(cltool) );
    std::cout<<"    \"description:\" : \""<<mytool->description()<<"\",\n";
    if( mytool->inputdata==commandline ) {
        std::cout<<"    \"inputtype:\" : \"command line args\",\n";
    } else if( mytool->inputdata==ifile ) {
        std::cout<<"    \"inputtype:\" : \"file\",\n";
    } else {
        error("input type for cltool was not specified");
    }
    std::cout<<"    \"syntax\" : {"<<std::endl;
    Keywords keys( cltoolRegister().get(cltool).keys );
    for(unsigned j=0; j<keys.size(); ++j) {
        std::string k = keys.getKeyword(j);
        std::cout<<"     \"" + k + "\": \"" + keys.getKeywordDescription( k ) <<"\"";
        if( j==keys.size()-1 ) {
            std::cout<<std::endl;
        } else {
            std::cout<<","<<std::endl;
        }
    }
    std::cout<<"       }"<<std::endl;
    std::cout<<"     }";
    if( i==cltool_names.size()-1 ) {
       std::cout<<std::endl;
    } else {
       std::cout<<","<<std::endl;
    }
  }
  std::cout<<"  }"<<std::endl;
  std::cout<<"}"<<std::endl;
  return 0;
}

} // End of namespace
}
