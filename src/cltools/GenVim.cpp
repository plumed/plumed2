/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2026 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "CLTool.h"
#include "core/ActionRegister.h"
#include "core/CLToolRegister.h"
#include "config/Config.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <set>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

namespace PLMD {
namespace cltools {
namespace {

namespace fs = std::filesystem;

/// Return the registered name without the image-handle prefix used by DLLoader.
std::string unqualifiedName(const std::string& registered) {
  const std::size_t colon=registered.find(':');
  return colon==std::string::npos ? registered : registered.substr(colon+1);
}

/// Enumerate base actions and actions registered by the loaded library handles.
std::vector<std::string> actionNames(const std::vector<void*>& images) {
  std::set<std::string> names;
  for(const auto& registered : actionRegister().getActionNames()) {
    if(registered.find(':')==std::string::npos) {
      names.insert(registered);
    }
  }
  for(void* image : images) {
    for(const auto& registered : actionRegister().getKeysWithDLHandle(image)) {
      names.insert(unqualifiedName(registered));
    }
  }
  return {names.begin(),names.end()};
}

/// Escape a string for a JSON-compatible Vim dictionary literal.
std::string quoteVimString(const std::string& value) {
  static const char hex[]="0123456789abcdef";
  std::string result="\"";
  for(unsigned char ch : value) {
    switch(ch) {
    case '"':
      result+="\\\"";
      break;
    case '\\':
      result+="\\\\";
      break;
    case '\b':
      result+="\\b";
      break;
    case '\f':
      result+="\\f";
      break;
    case '\n':
      result+="\\n";
      break;
    case '\r':
      result+="\\r";
      break;
    case '\t':
      result+="\\t";
      break;
    default:
      if(ch<0x20) {
        result+="\\u00";
        result+=hex[ch>>4];
        result+=hex[ch&0x0f];
      } else {
        result+=static_cast<char>(ch);
      }
    }
  }
  result+='"';
  return result;
}

/// Make one dictionary entry in the format consumed by plumed.vim.
std::string dictionaryEntry(const std::string& word, const std::string& menu) {
  return "{\"word\":"+quoteVimString(word)+",\"menu\":"+quoteVimString(menu)+"}";
}

/// Read a template as bytes and fail if it cannot be read completely.
std::string readFile(const fs::path& path) {
  std::ifstream input(path,std::ios::binary);
  if(!input) {
    throw std::runtime_error("cannot open Vim template: "+path.string());
  }
  std::string contents((std::istreambuf_iterator<char>(input)),std::istreambuf_iterator<char>());
  if(input.bad()) {
    throw std::runtime_error("cannot read Vim template: "+path.string());
  }
  return contents;
}

/// Write all bytes to a file and verify the completed write.
void writeFile(const fs::path& path, const std::string& contents) {
  std::ofstream output(path,std::ios::binary|std::ios::trunc);
  if(!output) {
    throw std::runtime_error("cannot create generated Vim file: "+path.string());
  }
  output.write(contents.data(),static_cast<std::streamsize>(contents.size()));
  output.flush();
  output.close();
  if(!output) {
    throw std::runtime_error("cannot write generated Vim file: "+path.string());
  }
}

/// Replace the dictionary marker with all action completion dictionaries.
std::string fillDictionaryTemplate(std::string source,
                                   const std::vector<std::string>& actions,
                                   const std::vector<void*>& images) {
  const std::string marker="@PLUMED_DICTIONARY@";
  const std::size_t markerPos=source.find(marker);
  if(markerPos==std::string::npos || source.find(marker,markerPos+marker.size())!=std::string::npos) {
    throw std::runtime_error("plumed.vim.in must contain exactly one @PLUMED_DICTIONARY@ marker");
  }
  if(markerPos>0 && source[markerPos-1]!='\n') {
    throw std::runtime_error("@PLUMED_DICTIONARY@ must be on its own line");
  }
  std::size_t markerEnd=markerPos+marker.size();
  if(markerEnd<source.size() && source[markerEnd]=='\r') {
    ++markerEnd;
  }
  if(markerEnd<source.size() && source[markerEnd]=='\n') {
    ++markerEnd;
  } else if(markerEnd<source.size()) {
    throw std::runtime_error("@PLUMED_DICTIONARY@ must be on its own line");
  }

  std::string dictionaries;
  for(const auto& action : actions) {
    Keywords keys;
    actionRegister().getKeywords(images,action,keys);
    std::vector<std::string> entries;
    entries.push_back(dictionaryEntry("LABEL=","(label)"));
    for(const auto& key : keys.getKeys()) {
      if(key=="LABEL") {
        continue;
      }
      if(keys.style(key,"flag")) {
        entries.push_back(dictionaryEntry(key,"(flag)"));
      } else if(keys.explicitlyNumbered(key)) {
        entries.push_back(dictionaryEntry(key,"(numbered)"));
      } else {
        entries.push_back(dictionaryEntry(key+"=","(option)"));
      }
    }
    std::sort(entries.begin(),entries.end());
    dictionaries+="let b:plumedDictionary["+quoteVimString(action)+"]=[";
    for(std::size_t i=0; i<entries.size(); ++i) {
      if(i) {
        dictionaries+=',';
      }
      dictionaries+=entries[i];
    }
    dictionaries+="]\n";
  }
  source.replace(markerPos,markerEnd-markerPos,dictionaries);
  return source;
}

/// Remove the temporary directory unless a failed rollback needs its backups.
class StagingDirectory {
  fs::path path_;
  bool preserve_=false;
public:
  /// Create a unique staging directory under the requested output directory.
  explicit StagingDirectory(const fs::path& parent) {
    static std::atomic<unsigned long long> counter{0};
    const auto ticks=std::chrono::steady_clock::now().time_since_epoch().count();
    for(unsigned attempt=0; attempt<100; ++attempt) {
      path_=parent/(".gen_vim_stage_"+std::to_string(ticks)+"_"+
                    std::to_string(counter.fetch_add(1)));
      std::error_code ec;
      if(fs::create_directory(path_,ec)) {
        return;
      }
      if(ec && ec!=std::errc::file_exists) {
        throw fs::filesystem_error("cannot create Vim staging directory",path_,ec);
      }
    }
    throw std::runtime_error("cannot create a unique Vim staging directory in "+parent.string());
  }

  /// Clean staged output on success or after a successful rollback.
  ~StagingDirectory() {
    if(!preserve_) {
      std::error_code ignored;
      fs::remove_all(path_,ignored);
    }
  }

  /// Get the staging path.
  const fs::path& path() const {
    return path_;
  }

  /// Preserve backups if rollback could not restore the original directories.
  void preserve() {
    preserve_=true;
  }
};

/// Replace syntax and help directories, restoring backups on failure.
void installStagedDirectories(StagingDirectory& stage, const fs::path& output) {
  const fs::path& staging=stage.path();
  const fs::path newSyntax=staging/"syntax";
  const fs::path newHelp=staging/"help";
  const fs::path oldSyntax=staging/"old_syntax";
  const fs::path oldHelp=staging/"old_help";
  const fs::path outputSyntax=output/"syntax";
  const fs::path outputHelp=output/"help";
  const bool hadSyntax=fs::exists(outputSyntax);
  const bool hadHelp=fs::exists(outputHelp);
  bool installedSyntax=false;
  bool installedHelp=false;
  try {
    if(hadSyntax) {
      fs::rename(outputSyntax,oldSyntax);
    }
    if(hadHelp) {
      fs::rename(outputHelp,oldHelp);
    }
    fs::rename(newSyntax,outputSyntax);
    installedSyntax=true;
    fs::rename(newHelp,outputHelp);
    installedHelp=true;
  } catch(const std::exception& original) {
    std::error_code ec;
    bool restored=true;
    if(installedHelp) {
      fs::remove_all(outputHelp,ec);
      restored=restored && !ec;
      ec.clear();
    }
    if(installedSyntax) {
      fs::remove_all(outputSyntax,ec);
      restored=restored && !ec;
      ec.clear();
    }
    const bool savedHelp=fs::exists(oldHelp,ec);
    restored=restored && !ec;
    ec.clear();
    if(hadHelp && savedHelp) {
      fs::rename(oldHelp,outputHelp,ec);
      restored=restored && !ec;
      ec.clear();
    }
    const bool savedSyntax=fs::exists(oldSyntax,ec);
    restored=restored && !ec;
    ec.clear();
    if(hadSyntax && savedSyntax) {
      fs::rename(oldSyntax,outputSyntax,ec);
      restored=restored && !ec;
    }
    if(!restored) {
      stage.preserve();
      throw std::runtime_error(std::string("failed to install Vim files and restore backups; backups remain in ")+
                               staging.string()+": "+original.what());
    }
    throw;
  }
}

} // namespace

//+PLUMEDOC TOOLS gen_vim
/*
Generate the PLUMED Vim syntax plugin and per-action help files from the live action registry.

The command writes `syntax/plumed.vim`, `syntax/plumedf.vim`, and `help/<ACTION>.txt` below `--outdir`. The static Vim templates are read from `--template-dir`; by default they are loaded from the installed PLUMED root.

*/
//+ENDPLUMEDOC

class GenVim : public CLTool {
public:
  /// Register command-line options for Vim file generation.
  /// \param keys Keyword collection receiving the generator's options.
  static void registerKeywords(Keywords& keys);
  /// Construct the Vim generator.
  /// \param options Command-line tool construction options.
  explicit GenVim(const CLToolOptions& options);
  /// Generate Vim syntax and help files in the requested output directory.
  /// \return Zero after all generated files have been installed successfully.
  /// \throws std::runtime_error If templates are invalid or files cannot be written.
  /// \throws std::filesystem::filesystem_error If output directories cannot be updated.
  int main(FILE* in, FILE* out, Communicator& pc) override;
  /// Return a brief description of this tool.
  std::string description() const override;
};

PLUMED_REGISTER_CLTOOL(GenVim,"gen_vim")

void GenVim::registerKeywords(Keywords& keys) {
  CLTool::registerKeywords(keys);
  keys.add("optional","--outdir","directory for generated Vim files");
  keys.add("optional","--template-dir","directory containing Vim templates");
}

GenVim::GenVim(const CLToolOptions& options):
  CLTool(options) {
  inputdata=inputType::commandline;
}

std::string GenVim::description() const {
  return "generate Vim syntax and help files from the registered actions";
}

int GenVim::main(FILE*,FILE*,Communicator&) {
  std::string outputDirectory=".";
  std::string templateDirectory=config::getPlumedRoot()+"/vim/templates";
  parse("--outdir",outputDirectory);
  parse("--template-dir",templateDirectory);

  const fs::path output=fs::absolute(outputDirectory);
  fs::create_directories(output);
  StagingDirectory staging(output);
  fs::create_directories(staging.path()/"syntax");
  fs::create_directories(staging.path()/"help");

  const auto actions=actionNames(getDLHandles());
  const fs::path templates(templateDirectory);
  const std::string syntaxTemplate=readFile(templates/"plumed.vim.in");
  const std::string fileSyntax=readFile(templates/"plumedf.vim");
  const std::string syntax=fillDictionaryTemplate(syntaxTemplate,actions,getDLHandles());
  writeFile(staging.path()/"syntax"/"plumed.vim",syntax);
  writeFile(staging.path()/"syntax"/"plumedf.vim",fileSyntax);

  for(const auto& action : actions) {
    Keywords keys;
    actionRegister().getKeywords(getDLHandles(),action,keys);
    const std::string header="****************************************\nShort helpfile for action "+
                             action+"\n****************************************\n";
    writeFile(staging.path()/"help"/(action+".txt"),header+keys.getHelpString()+"\n");
  }

  installStagedDirectories(staging,output);
  return 0;
}

} // namespace cltools
} // namespace PLMD
