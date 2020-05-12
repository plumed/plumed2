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
#include "CLTool.h"
#include "CLToolRegister.h"
#include "tools/Tools.h"
#include "core/Action.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "tools/Communicator.h"
#include "tools/Random.h"
#include <cstdio>
#include <string>
#include <vector>
#include <iostream>
#include "tools/File.h"
#include "core/Value.h"
#include "tools/Matrix.h"

using namespace std;

namespace PLMD {
namespace cltools {

//+PLUMEDOC TOOLS sum_hills
/*
sum_hills is a tool that allows one to to use plumed to post-process an existing hills/colvar file

\par Examples

a typical case is about the integration of a hills file:

\verbatim
plumed sum_hills  --hills PATHTOMYHILLSFILE
\endverbatim

The default name for the output file will be fes.dat
Note that starting from this version plumed will automatically detect the
number of the variables you have and their periodicity.
Additionally, if you use flexible hills (multivariate Gaussian kernels), plumed will understand it from the HILLS file.

The sum_hills tool will also accept multiple files that will be integrated one after the other

\verbatim
plumed sum_hills  --hills PATHTOMYHILLSFILE1,PATHTOMYHILLSFILE2,PATHTOMYHILLSFILE3
\endverbatim

if you want to integrate out some variable you do

\verbatim
plumed sum_hills  --hills PATHTOMYHILLSFILE   --idw t1 --kt 0.6
\endverbatim

where with --idw you define the variables that you want
all the others will be integrated out. --kt defines the temperature of the system in energy units.
(be consistent with the units you have in your hills: plumed will not check this for you)
If you need more variables then you may use a comma separated syntax

\verbatim
plumed sum_hills  --hills PATHTOMYHILLSFILE   --idw t1,t2 --kt 0.6
\endverbatim

You can define the output grid only with the number of bins you want
while min/max will be detected for you

\verbatim
plumed sum_hills --bin 99,99 --hills PATHTOMYHILLSFILE
\endverbatim

or full grid specification

\verbatim
plumed sum_hills --bin 99,99 --min -pi,-pi --max pi,pi --hills PATHTOMYHILLSFILE
\endverbatim

You can of course use numbers instead of -pi/pi.

You can use a --stride keyword to have a dump each bunch of hills you read
\verbatim
plumed sum_hills --stride 300 --hills PATHTOMYHILLSFILE
\endverbatim

You can also have, in case of well tempered metadynamics, only the negative
bias instead of the free energy through the keyword --negbias

\verbatim
plumed sum_hills --negbias --hills PATHTOMYHILLSFILE
\endverbatim

Here the default name will be negativebias.dat

From time to time you might need to use HILLS or a COLVAR file
as it was just a simple set  of points from which you want to build
a free energy by using -(1/beta)log(P)
then you use --histo

\verbatim
plumed sum_hills --histo PATHTOMYCOLVARORHILLSFILE  --sigma 0.2,0.2 --kt 0.6
\endverbatim

in this case you need a --kt to do the reweighting and then you
need also some width (with the --sigma keyword) for the histogram calculation (actually will be done with
Gaussian kernels, so it will be a continuous histogram)
Here the default output will be histo.dat.
Note that also here you can have multiple input files separated by a comma.

Additionally, if you want to do histogram and hills from the same file you can do as this
\verbatim
plumed sum_hills --hills --histo PATHTOMYCOLVARORHILLSFILE  --sigma 0.2,0.2 --kt 0.6
\endverbatim
The two files can be eventually the same

Another interesting thing one can do is monitor the difference in blocks as a metadynamics goes on.
When the bias deposited is constant over the whole domain one can consider to be at convergence.
This can be done with the --nohistory keyword

\verbatim
plumed sum_hills --stride 300 --hills PATHTOMYHILLSFILE  --nohistory
\endverbatim

and similarly one can do the same for an histogram file

\verbatim
plumed sum_hills --histo PATHTOMYCOLVARORHILLSFILE  --sigma 0.2,0.2 --kt 0.6 --nohistory
\endverbatim

just to check the hypothetical free energy calculated in single blocks of time during a simulation
and not in a cumulative way

Output format can be controlled via the --fmt field

\verbatim
plumed sum_hills --hills PATHTOMYHILLSFILE  --fmt %8.3f
\endverbatim

where here we chose a float with length of 8 and 3 digits

The output can be named in a arbitrary way  :

\verbatim
plumed sum_hills --hills PATHTOMYHILLSFILE  --outfile myfes.dat
\endverbatim

will produce a file myfes.dat which contains the free energy.

If you use stride, this keyword is the suffix

\verbatim
plumed sum_hills --hills PATHTOMYHILLSFILE  --outfile myfes_ --stride 100
\endverbatim

will produce myfes_0.dat,  myfes_1.dat, myfes_2.dat etc.

The same is true for the output coming from histogram
\verbatim
plumed sum_hills --histo HILLS --kt 2.5 --sigma 0.01 --outhisto myhisto.dat
\endverbatim

is producing a file myhisto.dat
while, when using stride, this is the suffix

\verbatim
plumed sum_hills --histo HILLS --kt 2.5 --sigma 0.01 --outhisto myhisto_ --stride 100
\endverbatim

that gives  myhisto_0.dat,  myhisto_1.dat,  myhisto_3.dat etc..

*/
//+ENDPLUMEDOC

class CLToolSumHills : public CLTool {
public:
  static void registerKeywords( Keywords& keys );
  explicit CLToolSumHills(const CLToolOptions& co );
  int main(FILE* in,FILE*out,Communicator& pc) override;
  string description()const override;
/// find a list of variables present, if they are periodic and which is the period
/// return false if the file does not exist
  static bool findCvsAndPeriodic(std::string filename, std::vector< std::vector <std::string> > &cvs,std::vector<std::string> &pmin,std::vector<std::string> &pmax, bool &multivariate, string &lowI_, string &uppI_);
};

void CLToolSumHills::registerKeywords( Keywords& keys ) {
  CLTool::registerKeywords( keys );
  keys.addFlag("--help-debug",false,"print special options that can be used to create regtests");
  keys.add("optional","--hills","specify the name of the hills file");
  keys.add("optional","--histo","specify the name of the file for histogram a colvar/hills file is good");
  keys.add("optional","--stride","specify the stride for integrating hills file (default 0=never)");
  keys.add("optional","--min","the lower bounds for the grid");
  keys.add("optional","--max","the upper bounds for the grid");
  keys.add("optional","--bin","the number of bins for the grid");
  keys.add("optional","--spacing","grid spacing, alternative to the number of bins");
  keys.add("optional","--idw","specify the variables to be used for the free-energy/histogram (default is all). With --hills the other variables will be integrated out, with --histo the other variables won't be considered");
  keys.add("optional","--outfile","specify the output file for sumhills");
  keys.add("optional","--outhisto","specify the output file for the histogram");
  keys.add("optional","--kt","specify temperature in energy units for integrating out variables");
  keys.add("optional","--sigma"," a vector that specify the sigma for binning (only needed when doing histogram ");
  keys.addFlag("--negbias",false," print the negative bias instead of the free energy (only needed with well tempered runs and flexible hills) ");
  keys.addFlag("--nohistory",false," to be used with --stride:  it splits the bias/histogram in pieces without previous history ");
  keys.addFlag("--mintozero",false," it translate all the minimum value in bias/histogram to zero (useful to compare results) ");
  keys.add("optional","--fmt","specify the output format");
}

CLToolSumHills::CLToolSumHills(const CLToolOptions& co ):
  CLTool(co)
{
  inputdata=commandline;
}

string CLToolSumHills::description()const { return "sum the hills with  plumed"; }

int CLToolSumHills::main(FILE* in,FILE*out,Communicator& pc) {

// Read the hills input file name
  vector<string> hillsFiles;
  bool dohills;
  dohills=parseVector("--hills",hillsFiles);
// Read the histogram file
  vector<string> histoFiles;
  bool dohisto;
  dohisto=parseVector("--histo",histoFiles);

  plumed_massert(dohisto || dohills,"you should use --histo or/and --hills command");

  vector< vector<string> > vcvs;
  vector<string> vpmin;
  vector<string> vpmax;
  string lowI_, uppI_;
  if(dohills) {
    // parse it as it was a restart
    bool vmultivariate;
    findCvsAndPeriodic(hillsFiles[0], vcvs, vpmin, vpmax, vmultivariate, lowI_, uppI_);
  }

  vector< vector<string> > hcvs;
  vector<string> hpmin;
  vector<string> hpmax;

  vector<std::string> sigma;
  if(dohisto) {
    bool hmultivariate;
    findCvsAndPeriodic(histoFiles[0], hcvs, hpmin, hpmax, hmultivariate, lowI_, uppI_);
    // here need also the vector of sigmas
    parseVector("--sigma",sigma);
    if(sigma.size()==0)plumed_merror("you should define --sigma vector when using histogram");
    lowI_=uppI_="-1.";  // Interval is not use for histograms
  }

  if(dohisto && dohills) {
    plumed_massert(vcvs==hcvs,"variables for histogram and bias should have the same labels");
    plumed_massert(hpmin==vpmin,"variables for histogram and bias should have the same min for periodicity");
    plumed_massert(hpmax==vpmax,"variables for histogram and bias should have the same max for periodicity");
  }

  // now put into a neutral vector

  vector< vector<string> > cvs;
  vector<string> pmin;
  vector<string> pmax;

  if(dohills) {
    cvs=vcvs;
    pmin=vpmin;
    pmax=vpmax;
  }
  if(dohisto) {
    cvs=hcvs;
    pmin=hpmin;
    pmax=hpmax;
  }


  // setup grids
  unsigned grid_check=0;
  vector<std::string> gmin(cvs.size());
  if(parseVector("--min",gmin)) {
    if(gmin.size()!=cvs.size() && gmin.size()!=0) plumed_merror("not enough values for --min");
    grid_check++;
  }
  vector<std::string> gmax(cvs.size() );
  if(parseVector("--max",gmax)) {
    if(gmax.size()!=cvs.size() && gmax.size()!=0) plumed_merror("not enough values for --max");
    grid_check++;
  }
  vector<std::string> gbin(cvs.size());
  bool grid_has_bin; grid_has_bin=false;
  if(parseVector("--bin",gbin)) {
    if(gbin.size()!=cvs.size() && gbin.size()!=0) plumed_merror("not enough values for --bin");
    grid_has_bin=true;
  }
  vector<std::string> gspacing(cvs.size());
  bool grid_has_spacing; grid_has_spacing=false;
  if(parseVector("--spacing",gspacing)) {
    if(gspacing.size()!=cvs.size() && gspacing.size()!=0) plumed_merror("not enough values for --spacing");
    grid_has_spacing=true;
  }
  // allowed: no grids only bin
  // not allowed: partial grid definition
  plumed_massert( gmin.size()==gmax.size() && (gmin.size()==0 ||  gmin.size()==cvs.size() ),"you should specify --min and --max together with same number of components");



  PlumedMain plumed;
  std::string ss;
  unsigned nn=1;
  ss="setNatoms";
  plumed.cmd(ss,&nn);
  if(Communicator::initialized())  plumed.cmd("setMPIComm",&pc.Get_comm());
  plumed.cmd("init",&nn);
  vector <bool> isdone(cvs.size(),false);
  for(unsigned i=0; i<cvs.size(); i++) {
    if(!isdone[i]) {
      isdone[i]=true;
      std::vector<std::string> actioninput;
      std::vector <unsigned> inds;
      actioninput.push_back("FAKE");
      actioninput.push_back("ATOMS=1");
      actioninput.push_back("LABEL="+cvs[i][0]);
      std::vector<std::string> comps, periods;
      if(cvs[i].size()>1) {comps.push_back(cvs[i][1]); inds.push_back(i);}
      periods.push_back(pmin[i]); periods.push_back(pmax[i]);
      for(unsigned j=i+1; j<cvs.size(); j++) {
        if(cvs[i][0]==cvs[j][0] && !isdone[j]) {
          if(cvs[i].size()==1 || cvs[j].size()==1  )plumed_merror("you cannot have twice the same label and no components ");
          if(cvs[j].size()>1) {
            comps.push_back(cvs[j][1]);
            periods.push_back(pmin[j]); periods.push_back(pmax[j]);
            isdone[j]=true; inds.push_back(j);
          }
        }

      }
      // drain all the components
      std::string addme;
      if(comps.size()>0) {
        addme="COMPONENTS=";
        for(unsigned i=0; i<comps.size()-1; i++)addme+=comps[i]+",";
        addme+=comps.back();
        actioninput.push_back(addme);
      }
      // periodicity (always explicit here)
      addme="PERIODIC=";
      for(unsigned j=0; j<periods.size()-1; j++) {
        addme+=periods[j]+",";
      }
      addme+=periods.back();
      actioninput.push_back(addme);
      for(unsigned j=0; j<inds.size(); j++) {
        unsigned jj; jj=inds[j];
        if(grid_check==2) {
          double gm;
          double pm;
          if(pmin[jj]!="none") {
            Tools::convert(gmin[jj],gm);
            Tools::convert(pmin[jj],pm);
            if(  gm<pm  ) {
              plumed_merror("Periodicity issue : GRID_MIN value ( "+gmin[jj]+" ) is less than periodicity in HILLS file in "+cvs[jj][0]+ " ( "+pmin[jj]+" ) ");
            }
          }
          if(pmax[jj]!="none") {
            Tools::convert(gmax[jj],gm);
            Tools::convert(pmax[jj],pm);
            if(  gm>pm ) {
              plumed_merror("Periodicity issue : GRID_MAX value ( "+gmax[jj]+" ) is more than periodicity in HILLS file in "+cvs[jj][0]+ " ( "+pmax[jj]+" ) ");
            }
          }
        }
      }

//  for(unsigned i=0;i< actioninput.size();i++){
//    cerr<<"AA "<<actioninput[i]<<endl;
//  }
      plumed.readInputWords(actioninput);
    }

  }
  unsigned ncv=cvs.size();
  std::vector<std::string> actioninput;
  vector<std::string> idw;
  // check if the variables to be used are correct
  if(parseVector("--idw",idw)) {
    for(unsigned i=0; i<idw.size(); i++) {
      bool found=false;
      for(unsigned j=0; j<cvs.size(); j++) {
        if(cvs[j].size()>1) {
          if(idw[i]==cvs[j][0]+"."+cvs[j][1])found=true;
        } else {
          if(idw[i]==cvs[j][0])found=true;
        }
      }
      if(!found)plumed_merror("variable "+idw[i]+" is not found in the bunch of cvs: revise your --idw option" );
    }
    plumed_massert( idw.size()<=cvs.size(),"the number of variables to be integrated should be at most equal to the total number of cvs  ");
    // in this case you neeed a beta factor!
  }

  std::string kt; kt=std::string("1.");// assign an arbitrary value just in case that idw.size()==cvs.size()
  if ( dohisto || idw.size()!=0  ) {
    plumed_massert(parse("--kt",kt),"if you make a dimensionality reduction (--idw) or a histogram (--histo) then you need to define --kt ");
  }

  std::string addme;

  actioninput.push_back("FUNCSUMHILLS");
  actioninput.push_back("ISCLTOOL");

  // set names
  std::string outfile;
  if(parse("--outfile",outfile)) {
    actioninput.push_back("OUTHILLS="+outfile);
  }
  std::string outhisto;
  if(parse("--outhisto",outhisto)) {
    actioninput.push_back("OUTHISTO="+outhisto);
  }


  addme="ARG=";
  for(unsigned i=0; i<(ncv-1); i++) {
    if(cvs[i].size()==1) {
      addme+=std::string(cvs[i][0])+",";
    } else {
      addme+=std::string(cvs[i][0])+"."+std::string(cvs[i][1])+",";
    }
  }
  if(cvs[ncv-1].size()==1) {
    addme+=std::string(cvs[ncv-1][0]);
  } else {
    addme+=std::string(cvs[ncv-1][0])+"."+std::string(cvs[ncv-1][1]);
  }
  actioninput.push_back(addme);
  //for(unsigned i=0;i< actioninput.size();i++){
  //  cerr<<"AA "<<actioninput[i]<<endl;
  //}
  if(dohills) {
    addme="HILLSFILES="; for(unsigned i=0; i<hillsFiles.size()-1; i++)addme+=hillsFiles[i]+","; addme+=hillsFiles[hillsFiles.size()-1];
    actioninput.push_back(addme);
    // set the grid
  }
  if(grid_check==2) {
    addme="GRID_MAX="; for(unsigned i=0; i<(ncv-1); i++)addme+=gmax[i]+","; addme+=gmax[ncv-1];
    actioninput.push_back(addme);
    addme="GRID_MIN="; for(unsigned i=0; i<(ncv-1); i++)addme+=gmin[i]+","; addme+=gmin[ncv-1];
    actioninput.push_back(addme);
  }
  if(grid_has_bin) {
    addme="GRID_BIN="; for(unsigned i=0; i<(ncv-1); i++)addme+=gbin[i]+","; addme+=gbin[ncv-1];
    actioninput.push_back(addme);
  }
  if(grid_has_spacing) {
    addme="GRID_SPACING="; for(unsigned i=0; i<(ncv-1); i++)addme+=gspacing[i]+","; addme+=gspacing[ncv-1];
    actioninput.push_back(addme);
  }
  std::string  stride; stride="";
  if(parse("--stride",stride)) {
    actioninput.push_back("INITSTRIDE="+stride);
    bool  nohistory;
    parseFlag("--nohistory",nohistory);
    if(nohistory) {
      actioninput.push_back("NOHISTORY");
    }
  }
  bool  mintozero;
  parseFlag("--mintozero",mintozero);
  if(mintozero) {
    actioninput.push_back("MINTOZERO");
  }
  if(idw.size()!=0) {
    addme="PROJ=";
    for(unsigned i=0; i<idw.size()-1; i++) {addme+=idw[i]+",";}
    addme+=idw.back();
    actioninput.push_back(addme);
  }

  if(dohisto) {
    if(idw.size()==0) {
      if(sigma.size()!=hcvs.size()) plumed_merror("you should define as many --sigma vector as the number of collective variable used for the histogram ");
    } else {
      if(idw.size()!=sigma.size()) plumed_merror("you should define as many --sigma vector as the number of collective variable used for the histogram ");
    }
  }

  if(idw.size()!=0 || dohisto) {
    actioninput.push_back("KT="+kt);
  }
  if(dohisto) {
    addme="HISTOFILES="; for(unsigned i=0; i<histoFiles.size()-1; i++) {addme+=histoFiles[i]+",";} addme+=histoFiles[histoFiles.size()-1];
    actioninput.push_back(addme);

    addme="HISTOSIGMA=";
    for(unsigned i=0; i<sigma.size()-1; i++) {addme+=sigma[i]+",";}
    addme+=sigma.back();
    actioninput.push_back(addme);
  }

  bool negbias;
  parseFlag("--negbias",negbias);
  if(negbias) {
    actioninput.push_back("NEGBIAS");
  }

  if(lowI_!=uppI_) {
    addme="INTERVAL="; addme+=lowI_+","; addme+=uppI_;
    actioninput.push_back(addme);
  }

  std::string fmt; fmt="";
  parse("--fmt",fmt);
  if(fmt!="")actioninput.push_back("FMT="+fmt);


//  for(unsigned i=0;i< actioninput.size();i++){
//   cerr<<"AA "<<actioninput[i]<<endl;
//  }
  plumed.readInputWords(actioninput);
  // if not a grid, then set it up automatically
  return 0;
}

bool CLToolSumHills::findCvsAndPeriodic(std::string filename, std::vector< std::vector<std::string>  > &cvs, std::vector<std::string> &pmin,std::vector<std::string> &pmax, bool &multivariate, string &lowI_, string &uppI_) {
  IFile ifile;
  ifile.allowIgnoredFields();
  std::vector<std::string> fields;
  if(ifile.FileExist(filename)) {
    cvs.clear(); pmin.clear(); pmax.clear();
    ifile.open(filename);
    ifile.scanFieldList(fields);
    bool before_sigma=true;
    for(unsigned i=0; i<fields.size(); i++) {
      size_t pos = 0;
      size_t founds,foundm,foundp;
      //found=(fields[i].find("sigma_", pos) || fields[i].find("min_", pos) || fields[i].find("max_", pos) ) ;
      founds=fields[i].find("sigma_", pos)  ;
      foundm=fields[i].find("min_", pos)  ;
      foundp=fields[i].find("max_", pos)  ;
      if (founds!=std::string::npos || foundm!=std::string::npos ||  foundp!=std::string::npos )before_sigma=false;
      // cvs are after time and before sigmas
      size_t  found;
      found=fields[i].find("time", pos);
      if( found==std::string::npos && before_sigma) {
        // separate the components
        size_t dot=fields[i].find_first_of('.');
        std::vector<std::string> ss;
        // this loop does not take into account repetitions
        if(dot!=std::string::npos) {
          std::string a=fields[i].substr(0,dot);
          std::string name=fields[i].substr(dot+1);
          ss.push_back(a);
          ss.push_back(name);
          cvs.push_back(ss);
        } else {
          std::vector<std::string> ss;
          ss.push_back(fields[i]);
          cvs.push_back(ss);
        }
        //std::cerr<<"found variable number  "<<cvs.size()<<" :  "<<cvs.back()[0]<<std::endl;
        //if((cvs.back()).size()!=1){
        //	std::cerr<<"component    "<<(cvs.back()).back()<<std::endl;
        //}
        // get periodicity
        pmin.push_back("none");
        pmax.push_back("none");
        std::string mm; if((cvs.back()).size()>1) {mm=cvs.back()[0]+"."+cvs.back()[1];} else {mm=cvs.back()[0];}
        if(ifile.FieldExist("min_"+mm)) {
          std::string val;
          ifile.scanField("min_"+mm,val);
          pmin[pmin.size()-1]=val;
          // std::cerr<<"found min   :  "<<pmin.back()<<std::endl;
        }
        //std::cerr<<"found min   :  "<<pmin.back()<<std::endl;
        if(ifile.FieldExist("max_"+mm)) {
          std::string val;
          ifile.scanField("max_"+mm,val);
          pmax[pmax.size()-1]=val;
          // std::cerr<<"found max   :  "<<pmax.back()<<std::endl;
        }
        //std::cerr<<"found max   :  "<<pmax.back()<<std::endl;
      }
    }
    // is multivariate ???
    std::string sss;
    multivariate=false;
    if(ifile.FieldExist("multivariate")) {
      ;
      ifile.scanField("multivariate",sss);
      if(sss=="true") { multivariate=true;}
      else if(sss=="false") { multivariate=false;}
    }
    // do interval?
    if(ifile.FieldExist("lower_int")) {
      ifile.scanField("lower_int",lowI_);
      ifile.scanField("upper_int",uppI_);
    } else {
      lowI_="-1.";
      uppI_="-1.";
    }
    ifile.scanField();
    return true;
  } else {
    return false;
  }
}


PLUMED_REGISTER_CLTOOL(CLToolSumHills,"sum_hills")



}
}
