/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2019 The plumed team
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

namespace PLMD {

Keywords CLToolOptions::emptyKeys;

CLToolOptions::CLToolOptions(const std::string &name):
  line(1,name),
  keys(emptyKeys)
{
}

CLToolOptions::CLToolOptions(const CLToolOptions& co, const Keywords& k):
  line(co.line),
  keys(k)
{
}

void CLTool::registerKeywords( Keywords& keys ) {
  keys.addFlag("--help/-h",false,"print this help");
}

CLTool::CLTool(const CLToolOptions& co ):
  name(co.line[0]),
  keywords(co.keys),
  inputdata(unset)
{
}

void CLTool::parseFlag( const std::string&key, bool&t ) {
  plumed_massert(keywords.exists(key),"keyword " + key + " has not been registered");
  plumed_massert(keywords.style(key,"flag"),"keyword " + key + " has not been registered as a flag");
  plumed_assert(inputData.count(key)>0);
  if( inputData[key]=="true") t=true;
  else if( inputData[key]=="false") t=false;
  else plumed_error();
}

bool CLTool::readInput( int argc, char**argv, FILE* in, FILE*out ) {
  plumed_massert( inputdata!=unset,"You have not specified where your tool reads its input. "
                  "If it is from the command line (like driver) add inputdata=commandline to the "
                  "tools constructor. If it reads everything from an input file (like simplemd) "
                  "add inputdata=ifile to the tools constructor");
  if(inputdata==commandline) return readCommandLineArgs( argc, argv, out );
  if(inputdata==ifile) return readInputFile( argc, argv, in, out );
  return true;
}

bool CLTool::readCommandLineArgs( int argc, char**argv, FILE*out ) {
  plumed_assert(inputdata==commandline);
  std::string prefix(""), a(""), thiskey;

  // Set all flags to default false
  for(unsigned k=0; k<keywords.size(); ++k) {
    thiskey=keywords.get(k);
    if( keywords.style(thiskey,"flag") ) inputData.insert(std::pair<std::string,std::string>(thiskey,"false"));
  }

  // Read command line arguments
  bool printhelp=false;
  for(int i=1; i<argc; i++) {
    a=prefix+argv[i];
    if(a.length()==0) continue;
    if(a=="-h" || a=="--help") {
      printhelp=true;
    } else {
      bool found=false;
      for(unsigned k=0; k<keywords.size(); ++k) {
        thiskey=keywords.get(k);
        if( keywords.style(thiskey,"flag") ) {
          if( a==thiskey ) { found=true; inputData[thiskey]="true"; }
        } else {
          if( a==thiskey ) {
            prefix=thiskey+"="; found=true;
            inputData.insert(std::pair<std::string,std::string>(thiskey,""));
          } else if(Tools::startWith(a,thiskey+"=")) {
            a.erase(0,a.find("=")+1); prefix=""; found=true;
            if(inputData.count(thiskey)==0) {
              inputData.insert(std::pair<std::string,std::string>(thiskey,a));
            } else {
              inputData[thiskey]=a;
            }
          }
        }
      }
      if(!found) {
        fprintf(stderr,"ERROR in input for command line tool %s : %s option is unknown \n\n", name.c_str(), a.c_str() );
        fprintf(out,"Usage: %s < inputFile \n", name.c_str() );
        fprintf(out,"inputFile should contain one directive per line.  The directives should come from amongst the following\n\n");
        keywords.print( out );
        printhelp=true;
      }
    }
    if(printhelp) break;
  }

  if(!printhelp) setRemainingToDefault(out);

  if(printhelp) {
    fprintf(out,"Usage: %s [options] \n\n", name.c_str() );
    keywords.print( out );
  }

  return !printhelp;
}

void CLTool::setRemainingToDefault(FILE* out) {
  std::string def, thiskey;
  for(unsigned k=0; k<keywords.size(); ++k) {
    thiskey=keywords.get(k);
    if( keywords.style(thiskey,"compulsory") ) {
      if( inputData.count(thiskey)==0 ) {
        if( keywords.getDefaultValue(thiskey,def) ) {
          plumed_assert( def.length()>0 );
          inputData.insert(std::pair<std::string,std::string>(thiskey,def));
        } else {
          fprintf(out,"ERROR : argument %s is compulsory. Use --help option for help\n",thiskey.c_str() );
          plumed_error();
        }
      }
    }
  }
}

bool CLTool::readInputFile( int argc, char**argv, FILE* in, FILE*out ) {
  plumed_assert(inputdata==ifile);

  // Check if use is just asking for help
  std::string a;
  for(int i=1; i<argc; i++) {
    a=argv[i];
    if(a.length()==0) continue;
    if(a=="-h" || a=="--help") {
      fprintf(out,"Usage: %s < inputFile \n", name.c_str() );
      fprintf(out,"inputFile should contain one directive per line.  The directives should come from amongst the following\n\n");
      keywords.print( out );
      return false;
    }
  }

  FILE* mystdin=in;
  if(argc==2) {
    mystdin=fopen(argv[1],"r");
    if(!mystdin) {
      fprintf(stderr,"ERROR: cannot open file %s\n\n",argv[1]);
      fprintf(out,"Usage: %s < inputFile \n", name.c_str() );
      fprintf(out,"inputFile should contain one directive per line.  The directives should come from amongst the following\n\n");
      keywords.print( out );
      return false;
    }
  }

  plumed_assert(mystdin);

  char buffer[256]; std::string line; line.resize(256);
  while(fgets(buffer,256,mystdin)) {
    line=buffer;
    for(unsigned i=0; i<line.length(); ++i) if(line[i]=='#' || line[i]=='\n') line.erase(i);
    Tools::stripLeadingAndTrailingBlanks( line );
    if(line.length()==0) continue;
    sscanf(line.c_str(),"%255s",buffer);
    std::string keyword=buffer; bool found=false;
    for(unsigned i=0; i<keywords.size(); ++i) {
      std::string thiskey=keywords.get(i);
      if(thiskey==keyword) {
        found=true;
        std::size_t keypos=line.find_first_of(keyword)+keyword.length();
        inputData.insert(std::pair<std::string,std::string>(thiskey,line.substr(keypos)));
        Tools::stripLeadingAndTrailingBlanks( inputData[thiskey] );
      }
    }
    if(!found) {
      fprintf(stderr,"ERROR in input for command line tool %s : unknown keyword %s found in input file\n\n",name.c_str(),keyword.c_str());
      fprintf(out,"Usage: %s < inputFile \n", name.c_str() );
      fprintf(out,"inputFile should contain one directive per line.  The directives should come from amongst the following\n\n");
      keywords.print( out );
      fclose(mystdin);
      return false;
    }
  }

  if(argc==2) fclose(mystdin);
  setRemainingToDefault(out);
  return true;
}

void CLTool::error( const std::string& msg ) {
  fprintf(stderr,"ERROR : in input for command line tool %s : %s\n",name.c_str(),msg.c_str());
  plumed_error();
}

}
