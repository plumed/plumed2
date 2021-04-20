/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2019,2020 The plumed team
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
#include "core/ActionShortcut.h"
#include "core/ActionSet.h"
#include "core/PlumedMain.h"
#include "tools/IFile.h"
#include <cstdio>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

namespace PLMD {
namespace cltools {

//+PLUMEDOC TOOLS gen_example
/*
gen_example is a tool that you can use to construct an example for the manual that users can interact with to understand

The example constructed by this action is in html. In all probability you will never need to use this
tool. However, it is used within the scripts that generate the html manual for PLUMED.  If you need to use this
tool outside those scripts the input is specified using the following command line arguments.

\par Examples

The following generates an example based on the contents of the plumed file plumed.dat
\verbatim
plumed gen_example --plumed plumed.dat --status working
\endverbatim


*/
//+ENDPLUMEDOC

class GenExample:
  public CLTool
{
private:
  int multi;
  std::string status, version;
  Communicator intracomm;
  Communicator intercomm;
public:
  static void registerKeywords( Keywords& keys );
  explicit GenExample(const CLToolOptions& co );
  int main(FILE* in, FILE*out,Communicator& pc) override;
  std::string description()const override {
    return "construct an example for the manual that users can interact with";
  }
  void printExampleInput( const std::vector<std::vector<std::string> >& input, const std::string& egname, const std::string& divname, std::ofstream& ofile );
  std::vector<std::vector<std::string> > createLongInput( const std::vector<std::vector<std::string> >& input );
};

PLUMED_REGISTER_CLTOOL(GenExample,"gen_example")

void GenExample::registerKeywords( Keywords& keys ) {
  CLTool::registerKeywords( keys );
  keys.add("compulsory","--plumed","plumed.dat","convert the input in this file to the html manual");
  keys.add("compulsory","--out","example.html","the file on which to output the example in html");
  keys.add("compulsory","--name","ppp","the name to use for this particular input");
  keys.add("compulsory","--status","nobadge","whether or not the input file works");
  keys.add("compulsory","--multi","0","set number of replicas for multi environment (needs MPI)");
}

GenExample::GenExample(const CLToolOptions& co ):
  CLTool(co),
  multi(0),
  status("nobadge"),
  version("master")
{
  inputdata=commandline;
}

int GenExample::main(FILE* in, FILE*out,Communicator& pc) {

// set up for multi replica driver:
  parse("--multi",multi);
  if(multi) {
    int ntot=pc.Get_size(); int nintra=ntot/multi;
    if(multi*nintra!=ntot) error("invalid number of processes for multi environment");
    pc.Split(pc.Get_rank()/nintra,pc.Get_rank(),intracomm);
    pc.Split(pc.Get_rank()%nintra,pc.Get_rank(),intercomm);
  } else {
    intracomm.Set_comm(pc.Get_comm());
  }

  if( config::getVersionLong().find("dev")==std::string::npos ) version="v"+config::getVersion();
  std::string fname, egname, outfile; parse("--plumed",fname);
  parse("--name",egname); parse("--out",outfile); parse("--status",status);

  int r=0;
  if(intracomm.Get_rank()==0) r=intercomm.Get_rank();
  intracomm.Bcast(r,0);
  if(r>0) outfile="/dev/null";

  IFile ifile; ifile.open(fname); ifile.allowNoEOL(); std::ofstream ofile; ofile.open(outfile); std::vector<bool> shortcuts;
  bool hasshortcuts=false, endplumed=false; std::vector<std::vector<std::string> > input; std::vector<std::string> words;
  while( Tools::getParsedLine(ifile, words, false) ) {
    input.push_back( words ); shortcuts.push_back( false );
    if( words.empty() || words[0].find("#")!=std::string::npos || endplumed ) continue;
    std::vector<std::string> interpreted( words ); Tools::interpretLabel(interpreted);
    if( interpreted[0]=="ENDPLUMED" ) { endplumed=true; continue; }
    Keywords keys; actionRegister().getKeywords( interpreted[0], keys );
    if( status=="working" && keys.exists("IS_SHORTCUT") ) hasshortcuts=shortcuts[shortcuts.size()-1]=true;
  }
  ifile.close();
  if( hasshortcuts ) {
    ofile<<"<div style=\"width: 80%; float:left\" id=\"value_details_"<<egname<<"\"> Click on the labels of the actions for more information on what each action computes </div>\n";
    ofile<<"<div style=\"width: 10%; float:left\"><button type=\"button\" id=\""<<egname<<"_button\" onclick=\'swapInput(\""<<egname<<"\")\'>contract shortcuts</button></div>";
  } else {
    ofile<<"<div style=\"width: 90%; float:left\" id=\"value_details_"<<egname<<"\"> Click on the labels of the actions for more information on what each action computes </div>\n";
  }
  ofile<<"<div style=\"width: 10%; float:left\">";
  ofile<<"<img src=\"https://img.shields.io/badge/";
  if(status=="working") ofile<<version<<"-passing-green";
  else if(status=="broken") ofile<<version<<"-failed-red";
  else if(status=="loads") ofile<<"with-LOAD-yellow";
  else if(status=="incomplete") ofile<<version<<"-incomplete-yellow";
  else error("unknown status");
  ofile<<".svg\" alt=\"tested on "<<version<<"\" /></div>";
  ofile.flush();
  if( hasshortcuts ) {
    // Write out the short version of the input
    ofile<<"<div style=\"width: 100%; float:left\" id=\"input_"<<egname<<"\"></div>"<<std::endl;
    // Write an extra pre to make sure the html after the example is put in the right place on the page
    ofile<<"<pre style=\"width: 97%;\" class=\"fragment\"></pre>"<<std::endl;
    ofile<<"<script type=\"text/javascript\">"<<std::endl;
    ofile<<"if (window.addEventListener) { // Mozilla, Netscape, Firefox"<<std::endl;
    ofile<<"    window.addEventListener('load', "<<egname<<"Load, false);"<<std::endl;
    ofile<<"} else if (window.attachEvent) { // IE"<<std::endl;
    ofile<<"    window.attachEvent('onload', "<<egname<<"Load);"<<std::endl;
    ofile<<"}"<<std::endl;
    ofile<<"function "<<egname<<"Load(event) {"<<std::endl;
    ofile<<"       swapInput(\""<<egname<<"\");"<<std::endl;
    ofile<<"}"<<std::endl;
    ofile<<"</script>"<<std::endl;
    ofile<<"<div style=\"display:none;\" id=\""<<egname<<"short\">"<<std::endl;
    printExampleInput( input, egname + "short", egname, ofile );
    ofile<<"</div>"<<std::endl;
    // Write out long version of the input
    ofile<<"<div style=\"display:none;\" id=\""<<egname<<"long\">";
    std::vector<std::vector<std::string> > long_input = createLongInput( input );
    printExampleInput( long_input, egname + "long", egname, ofile );
    ofile<<"</div>"<<std::endl;
  } else printExampleInput( input, egname, egname, ofile );
  ofile.close(); return 0;
}

std::vector<std::vector<std::string> > GenExample::createLongInput( const std::vector<std::vector<std::string> >& input ) {
  std::vector<std::vector<std::string> > long_input; PlumedMain myplumed; int rr=sizeof(double), natoms=10000000; double kt=2.49;
  myplumed.cmd("setRealPrecision",&rr);
  if(Communicator::initialized()) {
    if(multi) {
      if(intracomm.Get_rank()==0) myplumed.cmd("GREX setMPIIntercomm",&intercomm.Get_comm());
      myplumed.cmd("GREX setMPIIntracomm",&intracomm.Get_comm()); myplumed.cmd("GREX init");
    }
    myplumed.cmd("setMPIComm",&intracomm.Get_comm());
  }
  bool endplumed=false; myplumed.cmd("setNatoms",&natoms); myplumed.cmd("setKbT",&kt); myplumed.cmd("init");
  for(unsigned ll=0; ll<input.size(); ++ll) {
    if( input[ll].empty() || endplumed ) { long_input.push_back( input[ll] ); continue; }
    if( input[ll][0].find("#")!=std::string::npos ) { long_input.push_back( input[ll] ); continue; }
    std::vector<std::string> interpreted( input[ll] ); Tools::interpretLabel(interpreted);
    if( interpreted[0]=="ENDPLUMED" ) { endplumed=true; long_input.push_back( input[ll] ); continue; }
    Keywords keys; plumed_assert( actionRegister().check( interpreted[0] ) );
    actionRegister().getKeywords( interpreted[0], keys ); std::string lab, myinputline;
    if( Tools::parse(interpreted, "LABEL", lab ) ) myinputline = lab + ": ";
    myinputline += interpreted[0] + " "; bool trailingcomment=false;
    for(unsigned i=1; i<interpreted.size(); ++i) {
      if( trailingcomment && interpreted[i]=="@newline") { trailingcomment=false; continue; }
      if( interpreted[i].find("#")!=std::string::npos ) { trailingcomment=true; continue; }
      if( interpreted[i]=="@newline" || interpreted[i]=="..." ) continue;
      std::size_t pos = 0;  while ((pos = interpreted[i].find("@newline",pos)) != std::string::npos) { interpreted[i].replace(pos, 8, "\n"); pos++; }
      myinputline += interpreted[i] + " ";
    }
    if( status=="working" && keys.exists("IS_SHORTCUT") ) {
      myplumed.readInputLine( myinputline );
      ActionShortcut* as=dynamic_cast<ActionShortcut*>( myplumed.getActionSet()[myplumed.getActionSet().size()-1].get() );
      plumed_assert( as ); std::vector<std::string> shortcut_commands = as->getSavedInputLines();
      for(unsigned i=0; i<shortcut_commands.size(); ++i) {
        std::vector<std::string> words = Tools::getWords( shortcut_commands[i] ); long_input.push_back( words );
      }
    } else { long_input.push_back( input[ll] ); myplumed.readInputLine( myinputline ); }
  }
  return long_input;
}

void GenExample::printExampleInput( const std::vector<std::vector<std::string> >& input, const std::string& egname, const std::string& divname, std::ofstream& ofile ) {
  PlumedMain myplumed; int rr=sizeof(double), natoms=10000000; double kt=2.49;
  myplumed.cmd("setRealPrecision",&rr);
  if(Communicator::initialized()) {
    if(multi) {
      if(intracomm.Get_rank()==0) myplumed.cmd("GREX setMPIIntercomm",&intercomm.Get_comm());
      myplumed.cmd("GREX setMPIIntracomm",&intracomm.Get_comm()); myplumed.cmd("GREX init");
    }
    myplumed.cmd("setMPIComm",&intracomm.Get_comm());
  }
  myplumed.cmd("setNatoms",&natoms); myplumed.cmd("setKbT",&kt); myplumed.cmd("init");
  std::vector<std::string> labellist; bool endplumed=false;
  ofile<<"<pre style=\"width: 97%;\" class=\"fragment\">"<<std::endl;
  for(unsigned ll=0; ll<input.size(); ++ll) {
    if( input[ll].empty() ) { ofile<<std::endl; continue; }
    if( input[ll][0].find("#")!=std::string::npos || endplumed ) {
      ofile<<"<span style=\"color:blue\">"<<input[ll][0];
      for(unsigned i=1; i<input[ll].size(); ++i) ofile<<" "<<input[ll][i];
      ofile<<"</span>"<<std::endl;;
    } else {
      // Interpret the label if this needs to be done
      std::vector<std::string> interpreted( input[ll] ); Tools::interpretLabel(interpreted); std::string lab, myinputline;
      // Now read in the label
      if( Tools::parse(interpreted,"LABEL",lab) ) {
        ofile<<"<b name=\""<<egname<<lab<<"\" onclick=\'showPath(\""<<divname<<"\",\""<<egname<<lab<<"\")\'>"<<lab<<": </b>";
        labellist.push_back(lab); myinputline = lab + ": ";
      }
      // Print the keyword in use in the action
      std::string action = interpreted[0]; myinputline += interpreted[0] + " ";
      if( action=="ENDPLUMED" ) endplumed=true;
      Keywords keys; actionRegister().getKeywords( interpreted[0], keys );
      // Handle conversion of action names to links
      std::transform(action.begin(), action.end(), action.begin(), [](unsigned char c) { return std::tolower(c); });
      ofile<<"<a href=\"https://www.plumed.org/doc-"<<version<<"/user-doc/html/";
      for(unsigned n=0;; ++n) {
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
        if( interpreted[i]=="@newline" && i==1 ) { ofile<<"..."<<std::endl<<"   "; continue; }
        else if( interpreted[i]=="@newline" ) {
          if( trailingcomment ) { ofile<<"</span>"; trailingcomment=false; }
          if( interpreted[i+1]=="..." ) ofile<<std::endl;
          else ofile<<std::endl<<"   ";
          continue;
        } else if( interpreted[i]=="__FILL__" ) {
          if( status!="incomplete" ) error("found __FILL__ statement but status is " + status);
          ofile<<"<span style=\"background-color:yellow\">__FILL__</span>";
          continue;
        } else if( interpreted[i]==action ) continue;
        if( interpreted[i].find("#")!=std::string::npos ) { trailingcomment=true; ofile<<"<span style=\"color:blue\">"; }

        if( !trailingcomment ) {
          std::size_t eq=interpreted[i].find_first_of("=");
          if( eq!=std::string::npos ) {
            std::string keyword=interpreted[i].substr(0,eq), rest=interpreted[i].substr(eq+1);
            ofile<<"<div class=\"tooltip\">"<<keyword<<"<div class=\"right\">"<<keys.getTooltip(keyword)<<"<i></i></div></div>";
            if( rest=="__FILL__" ) {
              if( status!="incomplete" ) error("found __FILL__ statement but status is " + status);
              ofile<<"=<span style=\"background-color:yellow\">__FILL__</span>";
            } else if( rest.find_first_of("{")!=std::string::npos ) {
              std::size_t pos = 0;  while ((pos = rest.find("@newline",pos)) != std::string::npos) { rest.replace(pos, 8, "\n"); pos++; }
              ofile<<"="<<rest<<" "; myinputline += keyword + "=" + rest + " ";
            } else {
              std::vector<std::string> args=Tools::getWords(rest,"\t\n ,"); ofile<<"=";
              for(unsigned i=0; i<args.size(); ++i) {
                bool islabel=false; std::string thislab;
                for(unsigned j=0; j<labellist.size(); ++j) {
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
      if( status=="working" ) {
        ofile<<"<span style=\"display:none;\" id=\""<<egname<<lab<<"\">";
        ofile<<"The "<<interpreted[0]<<" action with label <b>"<<lab<<"</b>";
        myplumed.readInputLine( myinputline );
        ActionWithValue* av=dynamic_cast<ActionWithValue*>( myplumed.getActionSet().selectWithLabel<Action*>(lab) );
        if( av ) {
          if( av->getNumberOfComponents()==1 ) { ofile<<" calculates a single scalar value"; }
          else if( av->getNumberOfComponents()>0 ) {
            ofile<<" calculates the following quantities:"<<std::endl;
            ofile<<"<table  align=\"center\" frame=\"void\" width=\"95%%\" cellpadding=\"5%%\">"<<std::endl;
            ofile<<"<tr><td width=\"5%%\"><b> Quantity </b>  </td><td><b> Description </b> </td></tr>"<<std::endl;
            unsigned ncomp = av->getNumberOfComponents();
            for(unsigned k=0; k<ncomp; ++k ) {
              std::string myname = av->copyOutput(k)->getName(); std::size_t dot=myname.find_first_of(".");
              std::string tname=myname.substr(dot+1); std::size_t und=tname.find_first_of("_"); std::size_t hyph=tname.find_first_of("-");
              if( und!=std::string::npos && hyph!=std::string::npos ) plumed_merror("cannot use underscore and hyphen in name");
              ofile<<"<tr><td width=\"5%%\">"<<myname<<"</td><td>";
              if( und!=std::string::npos ) {
                ofile<<keys.getOutputComponentDescription(tname.substr(und))<<" This particular component measures this quantity for the input CV named ";
                ofile<<tname.substr(0,und);
              } else if( hyph!=std::string::npos ) {
                ofile<<keys.getOutputComponentDescription(tname.substr(0,hyph))<<"  This is the "<<tname.substr(hyph+1)<<"th of these quantities";
              } else ofile<<keys.getOutputComponentDescription(tname);
              ofile<<"</td></tr>";
            }
            ofile<<"</table>"<<std::endl;
          }
        } else {
          ActionWithVirtualAtom* avv=dynamic_cast<ActionWithVirtualAtom*>( myplumed.getActionSet().selectWithLabel<Action*>(lab) );
          if( avv ) ofile<<" calculates the position of a virtual atom";
          else if( interpreted[0]=="GROUP" ) ofile<<" defines a group of atoms so that they can be referred to later in the input";
        }
        ofile<<"</span>"<<std::endl;
      } else if( status!="working" ) {
        ofile<<"<span style=\"display:none;\" id=\""<<egname<<lab<<"\"> You cannot view the components that are calculated by each action for this input file. Sorry </span>"<<std::endl;
      } else ofile<<std::endl;
    }
    ofile.flush();
  }
  ofile<<"</pre>"<<std::endl;
}

} // End of namespace
}
