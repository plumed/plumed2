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
#include "CLToolRegister.h"
#include "tools/Tools.h"
#include "config/Config.h"
#include "core/ActionRegister.h"
#include "core/ActionWithValue.h"
#include "core/ActionWithVirtualAtom.h"
#include "core/ActionSet.h"
#include "core/PlumedMain.h"
#include "tools/IFile.h"
#include <cstdio>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

namespace PLMD {
namespace cltools {

//+PLUMEDOC TOOLS gen_example
/*
gen_example is a tool that you can use to construct a example for the manual that users can interact with to understand

The example constructed by this action is in html. In all probability you will never need to use this
tool. However, it is used within the scripts that generate the html manual for PLUMED.  If you need to use this
tool outside those scripts the input is specified using the following command line arguments.

\par Examples

The following generates an example based on the contents of the plumed file plumed.dat
\verbatim
plumed gen_example --plumed plumed.dat
\endverbatim


*/
//+ENDPLUMEDOC

class GenExample:
  public CLTool
{
public:
  static void registerKeywords( Keywords& keys );
  explicit GenExample(const CLToolOptions& co );
  int main(FILE* in, FILE*out,Communicator& pc) override;
  string description()const override {
    return "construct an example for the manual that users can interact with";
  }
};

PLUMED_REGISTER_CLTOOL(GenExample,"gen_example")

void GenExample::registerKeywords( Keywords& keys ) {
  CLTool::registerKeywords( keys );
  keys.add("compulsory","--plumed","plumed.dat","convert the input in this file to the html manual");
  keys.add("compulsory","--out","example.html","the file on which to output the example in html");
  keys.add("compulsory","--name","ppp","the name to use for this particular input");
}

GenExample::GenExample(const CLToolOptions& co ):
  CLTool(co)
{
  inputdata=commandline;
}

int GenExample::main(FILE* in, FILE*out,Communicator& pc) {

  PlumedMain myplumed; int rr=sizeof(double), natoms=10000000; 
  myplumed.cmd("setRealPrecision",&rr); myplumed.cmd("setNatoms",&natoms); myplumed.cmd("init");
  std::string fname, egname, outfile; parse("--plumed",fname); parse("--name",egname); parse("--out",outfile);
  IFile ifile; ifile.open(fname); ifile.allowNoEOL(); std::ofstream ofile; ofile.open(outfile);
  ofile<<"<div id=\"value_details_"<<egname<<"\"></div>\n<pre class=\"fragment\">\n"; 
  std::vector<std::string> labellist, words; 
  while( Tools::getParsedLine(ifile, words, false) ) {
     if( words.empty() ) continue; 
     if( words[0]=="#" ) {
         ofile<<"<span style=\"color:blue\">"<<words[0];
         for(unsigned i=1; i<words.size(); ++i) ofile<<" "<<words[i];
         ofile<<"</span>"<<std::endl;;
     } else {
         // Interpret the label if this needs to be done
         std::vector<std::string> interpreted( words ); Tools::interpretLabel(interpreted); std::string lab, myinputline;
         // Now read in the label
         if( Tools::parse(interpreted,"LABEL",lab) ) {
             ofile<<"<b name=\""<<egname<<lab<<"\" onclick=\'showPath(\""<<egname<<"\",\""<<egname<<lab<<"\")\'>"<<lab<<": </b>";
             labellist.push_back(lab); myinputline = lab + ": ";
         }
         // Print the keyword in use in the action
         std::string action = interpreted[0]; myinputline += interpreted[0] + " ";
         Keywords keys; actionRegister().getKeywords( interpreted[0], keys );
         // Handle conversion of action names to links
         std::transform(action.begin(), action.end(), action.begin(), [](unsigned char c){ return std::tolower(c); });
         std::string version="master"; if( config::getVersion()!="2.6") version=config::getVersion(); 
         ofile<<"<a href=\"https://www.plumed.org/doc-"<<version<<"/user-doc/html/"; 
         for(unsigned n=0;;++n) {
             std::size_t und=action.find_first_of("_");
             if( und==std::string::npos ) break;
             std::string first=action.substr(0,und);
             for(auto c : first ) { if( isdigit(c) ) ofile<<c; else ofile<<"_"<<c; }
             ofile<<"_"; action=action.substr(und+1);
         }
         for(auto c : action ) { if( isdigit(c) ) ofile<<c; else ofile<<"_"<<c; }
         ofile<<".html\" style=\"color:green\">"<<interpreted[0]<<"</a> "; 
         // And write out everything else in the input line
         bool trailingcomment=false; 
         for(unsigned i=1; i<interpreted.size(); ++i) {
             if( interpreted[i]=="@newline" && i==1 ) { ofile<<"...\n   "; continue; } 
             else if( interpreted[i]=="@newline" ) { 
                if( trailingcomment ) ofile<<"</span>"; 
                if( interpreted[i+1]=="..." ) ofile<<"\n";
                else ofile<<"\n   "; 
                continue; 
             } 
             if( interpreted[i]=="#" ) { trailingcomment=true; ofile<<"<span style=\"color:blue\">"; }

             if( !trailingcomment ) {
                 std::size_t eq=interpreted[i].find_first_of("=");
                 if( eq!=std::string::npos ) {
                    std::string keyword=interpreted[i].substr(0,eq), rest=interpreted[i].substr(eq+1); 
                    ofile<<"<div class=\"tooltip\">"<<keyword<<"<div class=\"right\">"<<keys.getTooltip(keyword)<<"<i></i></div></div>";
                    if( Tools::getWords(rest).size()>1 ) {
                        ofile<<"={"<<rest<<"}"; myinputline += keyword + "=" + "{" + rest + "} ";
                    } else {
                        std::vector<std::string> args=Tools::getWords(rest,"\t\n ,"); ofile<<"=";
                        for(unsigned i=0;i<args.size();++i) {
                            bool islabel=false; std::string thislab;
                            for(unsigned j=0;j<labellist.size();++j) {
                                std::size_t dot=args[i].find_first_of("."); std::string lll=args[i].substr(0,dot);
                                if( lll==labellist[j] ) { islabel=true; thislab=labellist[j]; break; }
                            }
                            if( islabel ) ofile<<"<b name=\""<<egname<<thislab<<"\">"<<args[i]<<"</b>";
                            else ofile<<args[i];
                            if( i!=args.size()-1 ) ofile<<",";
                        }
                        myinputline += interpreted[i] + " ";
                    }
                    ofile<<" ";
                 } else if( interpreted[i]!="@newline" && interpreted[i]!="..." ) {
                    myinputline += interpreted[i] + " ";
                    ofile<<"<div class=\"tooltip\">"<<interpreted[i]<<"<div class=\"right\">"<<keys.getTooltip(interpreted[i])<<"<i></i></div></div> ";
                 } else if( interpreted[i]=="..." ) ofile<<"...";
             } else ofile<<interpreted[i]<<" ";
         }
         if( trailingcomment ) ofile<<"</span>";
         // This builds the hidden content that tells the user about what is calculated
         if( lab.length()>0 ) {
            ofile<<"<span style=\"display:none;\" id=\""<<egname<<lab<<"\">";
            ofile<<"The "<<interpreted[0]<<" action with label <b>"<<lab<<"</b>";
            myplumed.readInputLine( myinputline );
            ActionWithValue* av=dynamic_cast<ActionWithValue*>( myplumed.getActionSet().selectWithLabel<Action*>(lab) );
            if( av ) {
                if( av->getNumberOfComponents()==1 ){ ofile<<" calculates a single scalar value"; }
                else if( av->getNumberOfComponents()>0 ) {
                   ofile<<" calculates the following quantities:\n";
                   ofile<<"<table  align=\"center\" frame=\"void\" width=\"95%%\" cellpadding=\"5%%\">\n";
                   ofile<<"<tr><td width=\"5%%\"><b> Quantity </b>  </td><td><b> Description </b> </td></tr>\n";
                   unsigned ncomp = av->getNumberOfComponents();
                   for(unsigned k=0; k<ncomp; ++k ) {
                       std::string myname = av->copyOutput(k)->getName(); std::size_t dot=myname.find_first_of(".");
                       ofile<<"<tr><td width=\"5%%\">"<<myname<<"</td><td>"<<keys.getOutputComponentDescription(myname.substr(dot+1))<<"</td></tr>";
                   } 
                   ofile<<"</table>\n";
                }
            } else {
                ActionWithVirtualAtom* avv=dynamic_cast<ActionWithVirtualAtom*>( myplumed.getActionSet().selectWithLabel<Action*>(lab) );
                if( avv ) ofile<<" calculates the position of a virtual atom";
                else if( interpreted[0]=="GROUP" ) ofile<<" defines a group of atoms so that they can be referred to later in the input";
            }
            ofile<<"</span>"<<std::endl; 
         } else ofile<<std::endl;
     }
  }
  ofile<<"</pre>\n"; ofile.close(); return 0;
}

} // End of namespace
}
