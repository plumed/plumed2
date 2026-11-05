/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2026 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "InputBuilderMetadata.h"

#include "ActionDescriptions.h"
#include "ActionRegister.h"
#include "CLTool.h"
#include "CLToolRegister.h"
#include "GenericMolInfo.h"
#include "ModuleMap.h"
#include "config/Config.h"
#include "tools/Exception.h"

#include <algorithm>
#include <iostream>
#include <ostream>
#include <set>
#include <string>
#include <vector>

namespace PLMD {
namespace {

/// Select the manual version matching this kernel.
std::string documentationVersion() {
  if(config::getVersionLong().find("dev")==std::string::npos) {
    return "v"+config::getVersion();
  }
  return "master";
}

/// Quote UTF-8 text as JSON, escaping control characters without changing bytes.
void writeString(std::ostream& output, const std::string& value) {
  static const char hex[]="0123456789abcdef";
  output.put('"');
  for(unsigned char ch : value) {
    switch(ch) {
    case '"':
      output << "\\\"";
      break;
    case '\\':
      output << "\\\\";
      break;
    case '\b':
      output << "\\b";
      break;
    case '\f':
      output << "\\f";
      break;
    case '\n':
      output << "\\n";
      break;
    case '\r':
      output << "\\r";
      break;
    case '\t':
      output << "\\t";
      break;
    default:
      if(ch<0x20) {
        output << "\\u00" << hex[ch>>4] << hex[ch&0x0f];
      } else {
        output.put(static_cast<char>(ch));
      }
    }
  }
  output.put('"');
}

/// Write an absolute link to a page in the matching manual.
void writeLink(std::ostream& output, const std::string& version, const std::string& page) {
  writeString(output,"https://www.plumed.org/doc-"+version+"/user-doc/html/"+page);
}

/// Resolve a built-in module or identify a dynamically loaded registration.
std::string moduleFor(const std::string& name, const std::set<std::string>& dynamicNames) {
  if(dynamicNames.count(name)) {
    return "external";
  }
  const auto& modules=getModuleMap();
  const auto found=modules.find(name);
  // Dynamically loaded extensions may have no entry in the generated module map.
  return found==modules.end() ? "external" : found->second;
}

/// Remove the image-handle prefix from a dynamically registered name.
std::string unqualifiedName(const std::string& registered) {
  const std::size_t colon=registered.find(':');
  return colon==std::string::npos ? registered : registered.substr(colon+1);
}

/// Serialize a keyword while preserving legacy fields and exact numbering.
void writeKeyword(std::ostream& output, const std::string& key,
                  const Keywords& keys, const std::string& indent) {
  std::string description=keys.getKeywordDescription(key);
  std::string defaultValue;
  if(description.find("default=")!=std::string::npos) {
    const std::size_t defaultStart=description.find_first_of('=');
    const std::size_t close=description.find_first_of(')');
    defaultValue=description.substr(defaultStart+1,close-defaultStart-2);
    description=description.substr(close+1);
  }
  const std::size_t dot=description.find('.');
  if(dot!=std::string::npos) {
    description.resize(dot);
  }

  output << indent;
  writeString(output,key);
  output << " : { \"type\": ";
  writeString(output,keys.getStyle(key));
  output << ", \"description\": ";
  writeString(output,description);
  output << ", \"multiple\": " << (keys.numbered(key) ? 1 : 0)
         << ", \"numbered\": " << (keys.explicitlyNumbered(key) ? "true" : "false");
  const std::string actionLink=keys.getLinkedActions(key);
  if(!actionLink.empty()) {
    output << ", \"actionlink\": ";
    writeString(output,actionLink);
  }
  const std::string pageLink=keys.getLinkedPages(key);
  if(!pageLink.empty()) {
    output << ", \"pagelink\": ";
    writeString(output,pageLink);
  }
  const std::string argumentType=keys.getArgumentType(key);
  if(!argumentType.empty()) {
    output << ", \"argtype\": ";
    writeString(output,argumentType);
  } else if(!defaultValue.empty()) {
    output << ", \"default\": ";
    writeString(output,defaultValue);
  }
  output << " }";
}

/// Serialize one action using the registrations visible to the caller.
void writeAction(std::ostream& output, const std::string& action,
                 const std::string& version, const std::string& description,
                 const std::vector<std::string>& actionNames,
                 const std::vector<void*>& images,
                 const std::set<std::string>& dynamicNames,
                 std::vector<std::string>& modules,
                 std::set<std::string>& knownModules) {
  output << "  ";
  writeString(output,action);
  output << ": {\n    \"hyperlink\" : ";
  writeLink(output,version,action);
  output << ",\n    \"description\" : ";
  writeString(output,description);

  const std::string module=moduleFor(action,dynamicNames);
  if(knownModules.insert(module).second) {
    modules.push_back(module);
  }
  output << ",\n    \"module\" : ";
  writeString(output,module);

  Keywords keys;
  actionRegister().getKeywords(images,action,keys);
  output << ",\n    \"displayname\" : ";
  writeString(output,keys.getDisplayName());
  const std::string replacement=keys.getReplacementAction();
  if(replacement!="none") {
    if(std::find(actionNames.begin(),actionNames.end(),replacement)==actionNames.end()) {
      plumed_error() << "could not find action named " << replacement
                     << " that is supposed to be used to replace " << action;
    }
    output << ",\n    \"replacement\" : ";
    writeString(output,replacement);
  }
  output << ",\n     \"dois\" : [";
  const auto& dois=keys.getDOIList();
  for(std::size_t i=0; i<dois.size(); ++i) {
    if(i) {
      output << ", ";
    }
    writeString(output,dois[i]);
  }
  output << "],\n    \"syntax\" : {\n";
  for(unsigned i=0; i<keys.size(); ++i) {
    writeKeyword(output,keys.getKeyword(i),keys,"       ");
    output << (i+1==keys.size() && !keys.exists("HAS_VALUES") ? "\n" : ",\n");
  }
  if(keys.exists("HAS_VALUES")) {
    output << "       \"output\" : {\n";
    const auto& components=keys.getOutputComponents();
    bool hasValue=true;
    for(const auto& component : components) {
      if(keys.getOutputComponentFlag(component)=="default") {
        hasValue=false;
        break;
      }
    }
    for(std::size_t i=0; i<components.size(); ++i) {
      std::string componentName=components[i];
      if(componentName==".#!value") {
        hasValue=false;
        componentName="value";
      }
      output << "         ";
      writeString(output,componentName);
      output << " : {\n           \"flag\": ";
      writeString(output,keys.getOutputComponentFlag(components[i]));
      output << ",\n           \"type\": ";
      writeString(output,keys.getOutputComponentType(components[i]));
      std::string componentDescription=keys.getOutputComponentDescription(components[i]);
      const std::size_t dot=componentDescription.find('.');
      if(dot!=std::string::npos) {
        componentDescription.resize(dot);
      }
      output << ",\n           \"description\": ";
      writeString(output,componentDescription);
      output << "\n         }" << (i+1==components.size() ? "\n" : ",\n");
    }
    if(hasValue && components.empty()) {
      std::cerr << "WARNING: no components have been registered for action " << action << '\n';
    }
    output << "       }\n";
  }
  output << "    },\n";
  const auto& needed=keys.getNeededKeywords();
  if(!needed.empty()) {
    output << "    \"needs\" : [";
    for(std::size_t i=0; i<needed.size(); ++i) {
      if(i) {
        output << ", ";
      }
      writeString(output,needed[i]);
    }
    output << "],\n";
  }
  output << "    \"help\" : ";
  writeString(output,keys.getHelpString());
  output << "\n  }";
}

} // namespace

void writeInputBuilderJSON(std::ostream& output,
                           const std::map<std::string,std::string>* descriptions,
                           const std::vector<void*>& images) {
  const std::string version=documentationVersion();
  const auto& embeddedDescriptions=getActionDescriptions();
  const auto& actionDescriptions=descriptions ? *descriptions : embeddedDescriptions;
  std::set<std::string> dynamicActions;
  std::set<std::string> actionNameSet;
  for(const auto& action : actionRegister().getActionNames()) {
    if(action.find(':')==std::string::npos) {
      actionNameSet.insert(action);
    }
  }
  for(void* image : images) {
    for(const auto& registered : actionRegister().getKeysWithDLHandle(image)) {
      const std::string action=unqualifiedName(registered);
      actionNameSet.insert(action);
      dynamicActions.insert(action);
    }
  }
  std::vector<std::string> actionNames(actionNameSet.begin(),actionNameSet.end());
  std::vector<std::string> modules;
  std::set<std::string> knownModules;

  output << "{\n  \"vimlink\" : ";
  writeLink(output,version,"vim");
  output << ",\n  \"replicalink\" : ";
  writeLink(output,version,"parsing.html");
  output << ",\n";
  for(std::size_t i=0; i<actionNames.size(); ++i) {
    // A loaded implementation can override a built-in name; its summary is
    // unknown unless the caller supplies one explicitly.
    const auto description=(!descriptions && dynamicActions.count(actionNames[i]))
                           ? actionDescriptions.end() : actionDescriptions.find(actionNames[i]);
    writeAction(output,actionNames[i],version,
                description==actionDescriptions.end() ? "" : description->second,
                actionNames,images,dynamicActions,modules,knownModules);
    output << ",\n";
  }

  // Collect CLTools before writing the module table so the external module is
  // defined whenever a dynamically loaded tool refers to it.
  std::set<std::string> cltoolNameSet;
  std::set<std::string> dynamicCLTools;
  for(const auto& tool : cltoolRegister().list()) {
    if(tool.find(':')==std::string::npos) {
      cltoolNameSet.insert(tool);
    }
  }
  for(void* image : images) {
    for(const auto& registered : cltoolRegister().getKeysWithDLHandle(image)) {
      const std::string tool=unqualifiedName(registered);
      cltoolNameSet.insert(tool);
      dynamicCLTools.insert(tool);
    }
  }
  const std::vector<std::string> cltoolNames(cltoolNameSet.begin(),cltoolNameSet.end());
  if(!dynamicCLTools.empty() && knownModules.insert("external").second) {
    modules.push_back("external");
  }

  output << "  \"groups\" : {\n"
         << "    \"@allatoms\" : {\n"
         << "        \"description\" : \"refers to all the MD codes atoms and PLUMEDs vatoms\",\n"
         << "        \"link\" : ";
  writeLink(output,version,"specifying_atoms");
  output << "\n    },\n"
         << "    \"@mdatoms\" : {\n"
         << "        \"description\" : \"refers to all the MD codes atoms but not PLUMEDs vatoms\",\n"
         << "        \"link\" : ";
  writeLink(output,version,"specifying_atoms");
  output << "\n    },\n"
         << "    \"@ndx:\" : {\n"
         << "        \"description\" : \"load a group from a GROMACS index file\",\n"
         << "        \"link\" : ";
  writeLink(output,version,"specifying_atoms.html");
  const auto specials=GenericMolInfo::getSpecialKeywords();
  for(const auto& special : specials) {
    output << "\n    },\n    ";
    writeString(output,special.first);
    output << " : {\n        \"description\" : ";
    writeString(output,special.second);
    output << ",\n        \"link\" : ";
    writeLink(output,version,"MOLINFO");
  }
  output << "\n    }\n  },\n  \"modules\" : {\n";
  for(std::size_t moduleIndex=0; moduleIndex<modules.size(); ++moduleIndex) {
    const auto& module=modules[moduleIndex];
    if(moduleIndex>0) {
      output << ",\n";
    }
    output << "    ";
    writeString(output,module);
    output << " : {\n        \"hyperlink\" : ";
    writeLink(output,version,"module_"+module);
    output << ",\n        \"description\" : \"A module that will be used for something\"\n    }";
  }
  output << "\n  },\n  \"cltools\" : {\n";

  for(std::size_t i=0; i<cltoolNames.size(); ++i) {
    const auto& toolName=cltoolNames[i];
    const auto module=moduleFor(toolName,dynamicCLTools);
    output << "    ";
    writeString(output,toolName);
    output << " : {\n      \"hyperlink\" : ";
    writeLink(output,version,toolName);
    output << ",\n      \"module\" : ";
    writeString(output,module);
    auto tool=cltoolRegister().create(images,CLToolOptions(toolName));
    output << ",\n      \"description\" : ";
    writeString(output,tool->description());
    if(tool->inputdata==CLTool::inputType::commandline) {
      output << ",\n      \"inputtype\" : \"command line args\",\n";
    } else if(tool->inputdata==CLTool::inputType::ifile) {
      output << ",\n      \"inputtype\" : \"file\",\n";
    } else {
      plumed_error() << "input type for cltool was not specified";
    }
    output << "      \"syntax\" : {\n";
    const Keywords& keys=cltoolRegister().get(images,toolName).keys;
    for(unsigned j=0; j<keys.size(); ++j) {
      writeKeyword(output,keys.getKeyword(j),keys,"        ");
      output << (j+1==keys.size() ? "\n" : ",\n");
    }
    output << "      }\n    }" << (i+1==cltoolNames.size() ? "\n" : ",\n");
  }
  output << "  }\n}\n";
  if(!output) {
    plumed_error() << "failed to write input-builder JSON metadata";
  }
}

} // namespace PLMD
