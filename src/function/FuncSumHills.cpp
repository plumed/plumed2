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
#include "ActionRegister.h"
#include "Function.h"
#include "tools/Exception.h"
#include "tools/Communicator.h"
#include "tools/BiasRepresentation.h"
#include "tools/KernelFunctions.h"
#include "tools/File.h"
#include "tools/Tools.h"
#include "tools/Stopwatch.h"
#include "tools/Grid.h"
#include <iostream>
#include <memory>

using namespace std;

namespace PLMD {
namespace function {


//+PLUMEDOC FUNCTION FUNCSUMHILLS
/*
This function is intended to be called by the command line tool sum_hills.  It is meant to integrate a HILLS file or an HILLS file interpreted as a histogram in a variety of ways. It is, therefore, not expected that you use this during your dynamics (it will crash!)

In the future one could implement periodic integration during the metadynamics
or straightforward MD as a tool to check convergence

\par Examples

There are currently no examples for this keyword.

*/
//+ENDPLUMEDOC

class FilesHandler {
  vector <string> filenames;
  vector <std::unique_ptr<IFile>>  ifiles;
  Action *action;
  Log *log;
  bool parallelread;
  unsigned beingread;
  bool isopen;
public:
  FilesHandler(const vector<string> &filenames, const bool &parallelread,  Action &myaction, Log &mylog);
  bool readBunch(BiasRepresentation *br, int stride);
  bool scanOneHill(BiasRepresentation *br, IFile *ifile );
  void getMinMaxBin(vector<Value*> vals, Communicator &cc, vector<double> &vmin, vector<double> &vmax, vector<unsigned> &vbin);
  void getMinMaxBin(vector<Value*> vals, Communicator &cc, vector<double> &vmin, vector<double> &vmax, vector<unsigned> &vbin, vector<double> &histosigma);
};
FilesHandler::FilesHandler(const vector<string> &filenames, const bool &parallelread, Action &action, Log &mylog ):filenames(filenames),log(&mylog),parallelread(parallelread),beingread(0),isopen(false) {
  this->action=&action;
  for(unsigned i=0; i<filenames.size(); i++) {
    std::unique_ptr<IFile> ifile(new IFile());
    ifile->link(action);
    plumed_massert((ifile->FileExist(filenames[i])), "the file "+filenames[i]+" does not exist " );
    ifiles.emplace_back(std::move(ifile));
  }

}

// note that the FileHandler is completely transparent respect to the biasrepresentation
// no check are made at this level
bool FilesHandler::readBunch(BiasRepresentation *br, int stride = -1) {
  bool morefiles; morefiles=true;
  if(parallelread) {
    (*log)<<"  doing parallelread \n";
    plumed_merror("parallelread is not yet implemented !!!");
  } else {
    (*log)<<"  doing serialread \n";
    // read one by one hills
    // is the type defined? if not, assume it is a gaussian
    IFile *ff;
    ff=ifiles[beingread].get();
    if(!isopen) {
      (*log)<<"  opening file "<<filenames[beingread]<<"\n";
      ff->open(filenames[beingread]); isopen=true;
    }
    int n;
    while(true) {
      bool fileisover=true;
      while(scanOneHill(br,ff)) {
        // here do the dump if needed
        n=br->getNumberOfKernels();
        if(stride>0 && n%stride==0 && n!=0  ) {
          (*log)<<"  done with this chunk: now with "<<n<<" kernels  \n";
          fileisover=false;
          break;
        }
      }
      if(fileisover) {
        (*log)<<"  closing file "<<filenames[beingread]<<"\n";
        ff->close();
        isopen=false;
        (*log)<<"  now total "<<br->getNumberOfKernels()<<" kernels \n";
        beingread++;
        if(beingread<ifiles.size()) {
          ff=ifiles[beingread].get(); ff->open(filenames[beingread]);
          (*log)<<"  opening file "<<filenames[beingread]<<"\n";
          isopen=true;
        } else {
          morefiles=false;
          (*log)<<"  final chunk: now with "<<n<<" kernels  \n";
          break;
        }
      }
      // if there are no more files to read and this file is over then quit
      if(fileisover && !morefiles) {break;}
      // if you are in the middle of a file and you are here
      // then means that you read what you need to read
      if(!fileisover ) {break;}
    }
  }
  return morefiles;
}
void FilesHandler::getMinMaxBin(vector<Value*> vals, Communicator &cc, vector<double> &vmin, vector<double> &vmax, vector<unsigned> &vbin) {
  // create the representation (no grid)
  BiasRepresentation br(vals,cc);
  // read all the kernels
  readBunch(&br);
  // loop over the kernels and get the support
  br.getMinMaxBin(vmin,vmax,vbin);
}
void FilesHandler::getMinMaxBin(vector<Value*> vals, Communicator &cc, vector<double> &vmin, vector<double> &vmax, vector<unsigned> &vbin, vector<double> &histosigma) {
  BiasRepresentation br(vals,cc,histosigma);
  // read all the kernels
  readBunch(&br);
  // loop over the kernels and get the support
  br.getMinMaxBin(vmin,vmax,vbin);
  //for(unsigned i=0;i<vals.size();i++){cerr<<"XXX "<<vmin[i]<<" "<<vmax[i]<<" "<<vbin[i]<<"\n";}
}
bool FilesHandler::scanOneHill(BiasRepresentation *br, IFile *ifile ) {
  double dummy;
  if(ifile->scanField("time",dummy)) {
    //(*log)<<"   scanning one hill: "<<dummy<<" \n";
    if(ifile->FieldExist("biasf")) ifile->scanField("biasf",dummy);
    if(ifile->FieldExist("clock")) ifile->scanField("clock",dummy);
    // keep this intermediate function in case you need to parse more data in the future
    br->pushKernel(ifile);
    //(*log)<<"  read hill\n";
    if(br->hasSigmaInInput())ifile->allowIgnoredFields();
    ifile->scanField();
    return true;
  } else {
    return false;
  }
}


double  mylog( double v1 ) {
  return log(v1);
}

double  mylogder( double v1 ) {
  return 1./v1;
}



class FuncSumHills :
  public Function
{
  vector<string> hillsFiles,histoFiles;
  vector<string> proj;
  int initstride;
  bool iscltool,integratehills,integratehisto,parallelread;
  bool negativebias;
  bool nohistory;
  bool minTOzero;
  bool doInt;
  double lowI_;
  double uppI_;
  double beta;
  string outhills,outhisto,fmt;
  std::unique_ptr<BiasRepresentation> biasrep;
  std::unique_ptr<BiasRepresentation> historep;
public:
  explicit FuncSumHills(const ActionOptions&);
  void calculate() override; // this probably is not needed
  bool checkFilesAreExisting(const vector<string> & hills );
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(FuncSumHills,"FUNCSUMHILLS")

void FuncSumHills::registerKeywords(Keywords& keys) {
  Function::registerKeywords(keys);
  keys.use("ARG");
  keys.add("optional","HILLSFILES"," source file for hills creation(may be the same as HILLS)"); // this can be a vector!
  keys.add("optional","HISTOFILES"," source file for histogram creation(may be the same as HILLS)"); // also this can be a vector!
  keys.add("optional","HISTOSIGMA"," sigmas for binning when the histogram correction is needed    ");
  keys.add("optional","PROJ"," only with sumhills: the projection on the CVs");
  keys.add("optional","KT"," only with sumhills: the kt factor when projection on CVs");
  keys.add("optional","GRID_MIN","the lower bounds for the grid");
  keys.add("optional","GRID_MAX","the upper bounds for the grid");
  keys.add("optional","GRID_BIN","the number of bins for the grid");
  keys.add("optional","GRID_SPACING","the approximate grid spacing (to be used as an alternative or together with GRID_BIN)");
  keys.add("optional","INTERVAL","set one dimensional INTERVAL");
  keys.add("optional","OUTHILLS"," output file for hills ");
  keys.add("optional","OUTHISTO"," output file for histogram ");
  keys.add("optional","INITSTRIDE"," stride if you want an initial dump ");
  keys.add("optional","STRIDE"," stride when you do it on the fly ");
  keys.addFlag("ISCLTOOL",true,"use via plumed command line: calculate at read phase and then go");
  keys.addFlag("PARALLELREAD",false,"read parallel HILLS file");
  keys.addFlag("NEGBIAS",false,"dump  negative bias ( -bias )   instead of the free energy: needed in well tempered with flexible hills ");
  keys.addFlag("NOHISTORY",false,"to be used with INITSTRIDE:  it splits the bias/histogram in pieces without previous history  ");
  keys.addFlag("MINTOZERO",false,"translate the resulting bias/histogram to have the minimum to zero  ");
  keys.add("optional","FMT","the format that should be used to output real numbers");
}

FuncSumHills::FuncSumHills(const ActionOptions&ao):
  Action(ao),
  Function(ao),
  initstride(-1),
  iscltool(false),
  integratehills(false),
  integratehisto(false),
  parallelread(false),
  negativebias(false),
  nohistory(false),
  minTOzero(false),
  doInt(false),
  lowI_(-1.),
  uppI_(-1.),
  beta(-1.),
  fmt("%14.9f")
{

  // format
  parse("FMT",fmt);
  log<<"  Output format is "<<fmt<<"\n";
  // here read
  // Grid Stuff
  vector<std::string> gmin;
  parseVector("GRID_MIN",gmin);
  if(gmin.size()!=getNumberOfArguments() && gmin.size()!=0) error("not enough values for GRID_MIN");
  plumed_massert(gmin.size()==getNumberOfArguments() || gmin.size()==0,"need GRID_MIN argument for this") ;
  vector<std::string> gmax;
  parseVector("GRID_MAX",gmax);
  if(gmax.size()!=getNumberOfArguments() && gmax.size()!=0) error("not enough values for GRID_MAX");
  plumed_massert(gmax.size()==getNumberOfArguments() || gmax.size()==0,"need GRID_MAX argument for this") ;
  vector<unsigned> gbin;
  vector<double>   gspacing;
  parseVector("GRID_BIN",gbin);
  plumed_massert(gbin.size()==getNumberOfArguments() || gbin.size()==0,"need GRID_BIN argument for this") ;
  if(gbin.size()!=getNumberOfArguments() && gbin.size()!=0) error("not enough values for GRID_BIN");
  parseVector("GRID_SPACING",gspacing);
  plumed_massert(gspacing.size()==getNumberOfArguments() || gspacing.size()==0,"need GRID_SPACING argument for this") ;
  if(gspacing.size()!=getNumberOfArguments() && gspacing.size()!=0) error("not enough values for GRID_SPACING");
  if(gspacing.size()!=0 && gbin.size()==0) {
    log<<"  The number of bins will be estimated from GRID_SPACING\n";
  } else if(gspacing.size()!=0 && gbin.size()!=0) {
    log<<"  You specified both GRID_BIN and GRID_SPACING\n";
    log<<"  The more conservative (highest) number of bins will be used for each variable\n";
  }
  if(gspacing.size()!=0) for(unsigned i=0; i<getNumberOfArguments(); i++) {
      if(gbin.size()==0) gbin.assign(getNumberOfArguments(),1);
      double a,b;
      Tools::convert(gmin[i],a);
      Tools::convert(gmax[i],b);
      unsigned n=((b-a)/gspacing[i])+1;
      if(gbin[i]<n) gbin[i]=n;
    }

  // Inteval keyword
  vector<double> tmpI(2);
  parseVector("INTERVAL",tmpI);
  if(tmpI.size()!=2&&tmpI.size()!=0) error("both a lower and an upper limits must be provided with INTERVAL");
  else if(tmpI.size()==2) {
    lowI_=tmpI.at(0);
    uppI_=tmpI.at(1);
    if(getNumberOfArguments()!=1) error("INTERVAL limits correction works only for monodimensional metadynamics!");
    if(uppI_<lowI_) error("The Upper limit must be greater than the Lower limit!");
    doInt=true;
  }
  if(doInt) {
    log << "  Upper and Lower limits boundaries for the bias are activated at " << lowI_ << " - " << uppI_<<"\n";
    log << "  Using the same values as boundaries for the grid if not other value was defined (default: 200 bins)\n";
    std::ostringstream strsmin, strsmax;
    strsmin << lowI_;
    strsmax << uppI_;
    if(gmin.size()==0) gmin.push_back(strsmin.str());
    if(gmax.size()==0) gmax.push_back(strsmax.str());
    if(gbin.size()==0) gbin.push_back(200);
  }


  // hills file:
  parseVector("HILLSFILES",hillsFiles);
  if(hillsFiles.size()==0) {
    integratehills=false; // default behaviour
  } else {
    integratehills=true;
    for(unsigned i=0; i<hillsFiles.size(); i++) log<<"  hillsfile  : "<<hillsFiles[i]<<"\n";
  }
  // histo file:
  parseVector("HISTOFILES",histoFiles);
  if(histoFiles.size()==0) {
    integratehisto=false;
  } else {
    integratehisto=true;
    for(unsigned i=0; i<histoFiles.size(); i++) log<<"  histofile  : "<<histoFiles[i]<<"\n";
  }
  vector<double> histoSigma;
  if(integratehisto) {
    parseVector("HISTOSIGMA",histoSigma);
    for(unsigned i=0; i<histoSigma.size(); i++) log<<"  histosigma  : "<<histoSigma[i]<<"\n";
  }

  // needs a projection?
  proj.clear();
  parseVector("PROJ",proj);
  if(integratehills) {
    plumed_massert(proj.size()<getNumberOfArguments()," The number of projection must be less than the full list of arguments ");
  }
  if(integratehisto) {
    plumed_massert(proj.size()<=getNumberOfArguments()," The number of projection must be less or equal to the full list of arguments ");
  }
  if(integratehisto&&proj.size()==0) {
    for(unsigned i=0; i<getNumberOfArguments(); i++) proj.push_back(getPntrToArgument(i)->getName());
  }

  // add some automatic hills width: not in case stride is defined
  // since when you start from zero the automatic size will be zero!
  if(gmin.size()==0 || gmax.size()==0) {
    log<<"   \n";
    log<<"  No boundaries defined: need to do a prescreening of hills \n";
    std::vector<Value*> tmphillsvalues, tmphistovalues;
    if(integratehills) {
      for(unsigned i=0; i<getNumberOfArguments(); i++)tmphillsvalues.push_back( getPntrToArgument(i) );
    }
    if(integratehisto) {
      for(unsigned i=0; i<getNumberOfArguments(); i++) {
        std::string ss = getPntrToArgument(i)->getName();
        for(unsigned j=0; j<proj.size(); j++) {
          if(proj[j]==ss) tmphistovalues.push_back( getPntrToArgument(i) );
        }
      }
    }

    if(integratehills) {
      FilesHandler hillsHandler(hillsFiles,parallelread,*this, log);
      vector<double> vmin,vmax;
      vector<unsigned> vbin;
      hillsHandler.getMinMaxBin(tmphillsvalues,comm,vmin,vmax,vbin);
      log<<"  found boundaries from hillsfile: \n";
      gmin.resize(vmin.size());
      gmax.resize(vmax.size());
      if(gbin.size()==0) {
        gbin=vbin;
      } else {
        log<<"  found nbins in input, this overrides the automatic choice \n";
      }
      for(unsigned i=0; i<getNumberOfArguments(); i++) {
        Tools::convert(vmin[i],gmin[i]);
        Tools::convert(vmax[i],gmax[i]);
        log<<"  variable "<< getPntrToArgument(i)->getName()<<" min: "<<gmin[i]<<" max: "<<gmax[i]<<" nbin: "<<gbin[i]<<"\n";
      }
    }
    // if at this stage bins are not there then do it with histo
    if(gmin.size()==0) {
      FilesHandler histoHandler(histoFiles,parallelread,*this, log);
      vector<double> vmin,vmax;
      vector<unsigned> vbin;
      histoHandler.getMinMaxBin(tmphistovalues,comm,vmin,vmax,vbin,histoSigma);
      log<<"  found boundaries from histofile: \n";
      gmin.resize(vmin.size());
      gmax.resize(vmax.size());
      if(gbin.size()==0) {
        gbin=vbin;
      } else {
        log<<"  found nbins in input, this overrides the automatic choice \n";
      }
      for(unsigned i=0; i<proj.size(); i++) {
        Tools::convert(vmin[i],gmin[i]);
        Tools::convert(vmax[i],gmax[i]);
        log<<"  variable "<< proj[i] <<" min: "<<gmin[i]<<" max: "<<gmax[i]<<" nbin: "<<gbin[i]<<"\n";
      }
    }
    log<<"  done!\n";
    log<<"   \n";
  }


  if( proj.size() != 0 || integratehisto==true  ) {
    parse("KT",beta);
    for(unsigned i=0; i<proj.size(); i++) log<<"  projection "<<i<<" : "<<proj[i]<<"\n";
    // this should be only for projection or free energy from histograms
    plumed_massert(beta>0.,"if you make a projection or a histogram correction then you need KT flag!");
    beta=1./beta;
    log<<"  beta is "<<beta<<"\n";
  }
  // is a cltool: then you start and then die
  parseFlag("ISCLTOOL",iscltool);
  //
  parseFlag("NEGBIAS",negativebias);
  //
  parseFlag("PARALLELREAD",parallelread);
  // stride
  parse("INITSTRIDE",initstride);
  // output suffix or names
  if(initstride<0) {
    log<<"  Doing only one integration: no stride \n";
    outhills="fes.dat"; outhisto="histo.dat";
  }
  else {
    outhills="fes_"; outhisto="histo_";
    log<<"  Doing integration slices every "<<initstride<<" kernels\n";
    parseFlag("NOHISTORY",nohistory);
    if(nohistory)log<<"  nohistory: each stride block has no memory of the previous block\n";
  }
  parseFlag("MINTOZERO",minTOzero);
  if(minTOzero)log<<"  mintozero: bias/histogram will be translated to have the minimum value equal to zero\n";
  //what might it be this?
  // here start
  // want something right now?? do it and return
  // your argument is a set of cvs
  // then you need: a hills / a colvar-like file (to do a histogram)
  // create a bias representation for this
  if(iscltool) {

    std::vector<Value*> tmphillsvalues, tmphistovalues;
    if(integratehills) {
      for(unsigned i=0; i<getNumberOfArguments(); i++) {
        // allocate a new value from the old one: no deriv here
        // if we are summing hills then all the arguments are needed
        tmphillsvalues.push_back( getPntrToArgument(i) );
      }
    }
    if(integratehisto) {
      for(unsigned i=0; i<getNumberOfArguments(); i++) {
        std::string ss = getPntrToArgument(i)->getName();
        for(unsigned j=0; j<proj.size(); j++) {
          if(proj[j]==ss) tmphistovalues.push_back( getPntrToArgument(i) );
        }
      }
    }

    // check if the files exists
    if(integratehills) {
      checkFilesAreExisting(hillsFiles);
      biasrep.reset(new BiasRepresentation(tmphillsvalues,comm, gmin, gmax, gbin, doInt, lowI_, uppI_));
      if(negativebias) {
        biasrep->setRescaledToBias(true);
        log<<"  required the -bias instead of the free energy \n";
        if(initstride<0) {outhills="negativebias.dat";}
        else {outhills="negativebias_";}
      }
    }

    parse("OUTHILLS",outhills);
    parse("OUTHISTO",outhisto);
    if(integratehills)log<<"  output file for fes/bias  is :  "<<outhills<<"\n";
    if(integratehisto)log<<"  output file for histogram is :  "<<outhisto<<"\n";
    checkRead();

    log<<"\n";
    log<<"  Now calculating...\n";
    log<<"\n";

    // here it defines the column to be histogrammed, tmpvalues should be only
    // the list of the collective variable one want to consider
    if(integratehisto) {
      checkFilesAreExisting(histoFiles);
      historep.reset(new BiasRepresentation(tmphistovalues,comm,gmin,gmax,gbin,histoSigma));
    }

    // decide how to source hills ( serial/parallel )
    // here below the input control
    // say how many hills and it will read them from the
    // bunch of files provided, will update the representation
    // of hills (i.e. a list of hills and the associated grid)

    // decide how to source colvars ( serial parallel )
    std::unique_ptr<FilesHandler> hillsHandler;
    std::unique_ptr<FilesHandler> histoHandler;

    if(integratehills)	hillsHandler.reset(new FilesHandler(hillsFiles,parallelread,*this, log));
    if(integratehisto)	histoHandler.reset(new FilesHandler(histoFiles,parallelread,*this, log));

// Stopwatch is logged when it goes out of scope
    Stopwatch sw(log);

// Stopwatch is stopped when swh goes out of scope
    auto swh=sw.startStop("0 Summing hills");

    // read a number of hills and put in the bias representation
    int nfiles=0;
    bool ibias=integratehills; bool ihisto=integratehisto;
    while(true) {
      if(  integratehills  && ibias  ) {
        if(nohistory) {biasrep->clear(); log<<"  clearing history before reading a new block\n";};
        log<<"  reading hills: \n";
        ibias=hillsHandler->readBunch(biasrep.get(),initstride) ; log<<"\n";
      }

      if(  integratehisto  && ihisto ) {
        if(nohistory) {historep->clear(); log<<"  clearing history before reading a new block\n";};
        log<<"  reading histogram: \n";
        ihisto=histoHandler->readBunch(historep.get(),initstride) ;  log<<"\n";
      }

      // dump: need to project?
      if(proj.size()!=0) {

        if(integratehills) {

          log<<"  Bias: Projecting on subgrid... \n";
          BiasWeight Bw(beta);
          Grid biasGrid=*(biasrep->getGridPtr());
          Grid smallGrid=biasGrid.project(proj,&Bw);
          OFile gridfile; gridfile.link(*this);
          std::ostringstream ostr; ostr<<nfiles;
          string myout;
          if(initstride>0) { myout=outhills+ostr.str()+".dat" ;} else {myout=outhills;}
          log<<"  Bias: Writing subgrid on file "<<myout<<" \n";
          gridfile.open(myout);
          if(minTOzero) smallGrid.setMinToZero();
          smallGrid.setOutputFmt(fmt);
          smallGrid.writeToFile(gridfile);
          gridfile.close();
          if(!ibias)integratehills=false;// once you get to the final bunch just give up
        }
        // this should be removed
        if(integratehisto) {

          log<<"  Histo: Projecting on subgrid... \n";
          Grid histoGrid=*(historep->getGridPtr());

          OFile gridfile; gridfile.link(*this);
          std::ostringstream ostr; ostr<<nfiles;
          string myout;
          if(initstride>0) { myout=outhisto+ostr.str()+".dat" ;} else {myout=outhisto;}
          log<<"  Histo: Writing subgrid on file "<<myout<<" \n";
          gridfile.open(myout);

          histoGrid.applyFunctionAllValuesAndDerivatives(&mylog,&mylogder);
          histoGrid.scaleAllValuesAndDerivatives(-1./beta);
          if(minTOzero) histoGrid.setMinToZero();
          histoGrid.setOutputFmt(fmt);
          histoGrid.writeToFile(gridfile);

          if(!ihisto)integratehisto=false;// once you get to the final bunch just give up
        }

      } else {

        if(integratehills) {

          Grid biasGrid=*(biasrep->getGridPtr());
          biasGrid.scaleAllValuesAndDerivatives(-1.);

          OFile gridfile; gridfile.link(*this);
          std::ostringstream ostr; ostr<<nfiles;
          string myout;
          if(initstride>0) { myout=outhills+ostr.str()+".dat" ;} else {myout=outhills;}
          log<<"  Writing full grid on file "<<myout<<" \n";
          gridfile.open(myout);

          if(minTOzero) biasGrid.setMinToZero();
          biasGrid.setOutputFmt(fmt);
          biasGrid.writeToFile(gridfile);
          // rescale back prior to accumulate
          if(!ibias)integratehills=false;// once you get to the final bunch just give up
        }
        if(integratehisto) {

          Grid histoGrid=*(historep->getGridPtr());
          // do this if you want a free energy from a grid, otherwise do not
          histoGrid.applyFunctionAllValuesAndDerivatives(&mylog,&mylogder);
          histoGrid.scaleAllValuesAndDerivatives(-1./beta);

          OFile gridfile; gridfile.link(*this);
          std::ostringstream ostr; ostr<<nfiles;
          string myout;
          if(initstride>0) { myout=outhisto+ostr.str()+".dat" ;} else {myout=outhisto;}
          log<<"  Writing full grid on file "<<myout<<" \n";
          gridfile.open(myout);

          // also this is usefull only for free energy
          if(minTOzero) histoGrid.setMinToZero();
          histoGrid.setOutputFmt(fmt);
          histoGrid.writeToFile(gridfile);

          if(!ihisto)integratehisto=false; // once you get to the final bunch just give up
        }
      }
      if ( !ibias && !ihisto) break; //when both are over then just quit

      nfiles++;
    }

    return;
  }
  // just an initialization but you need to do something on the fly?: need to connect with a metad run and its grid representation
  // your argument is a metad run
  // if the grid does not exist crash and say that you need some data
  // otherwise just link with it

}

void FuncSumHills::calculate() {
  // this should be connected only with a grid representation to metadynamics
  // at regular time just dump it
  plumed_merror("You should have never got here: this stuff is not yet implemented!");
}

bool FuncSumHills::checkFilesAreExisting(const vector<string> & hills ) {
  plumed_massert(hills.size()!=0,"the number of  files provided should be at least one" );
  std::unique_ptr<IFile> ifile(new IFile());
  ifile->link(*this);
  for(unsigned i=0; i< hills.size(); i++) {
    plumed_massert(ifile->FileExist(hills[i]),"missing file "+hills[i]);
  }
  return true;

}

}

}


