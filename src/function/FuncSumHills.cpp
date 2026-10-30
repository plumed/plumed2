/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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
#include "core/ActionRegister.h"
#include "Function.h"
#include "tools/Exception.h"
#include "tools/Communicator.h"
#include "tools/BiasRepresentation.h"
#include "tools/KernelFunctions.h"
#include "tools/File.h"
#include <limits>
#include "tools/Tools.h"
#include "tools/Stopwatch.h"
#include "tools/Grid.h"

namespace PLMD {
namespace function {


//+PLUMEDOC FUNCTION FUNCSUMHILLS
/*
This function is used by the command line tool [sum_hills](sum_hills.md) to sum the hills of a metadynamics simulation
or to compute a histogram. It does all its work when it is created and it cannot be used in a simulation (it stops with an error).

All its options are set by [sum_hills](sum_hills.md): see the documentation of [sum_hills](sum_hills.md) for their meaning.

## Examples

There are currently no examples for this action: use [sum_hills](sum_hills.md) instead.

*/
//+ENDPLUMEDOC

class FilesHandler {
  std::vector <std::string> filenames;
  std::vector <std::unique_ptr<IFile>>  ifiles;
  Action *action;
  Log *log;
  bool parallelread;
  unsigned beingread=0;
  bool isopen=false;
/// key used to merge the files ("time" or "clock"); empty to read the files one after the other
  std::string sortby;
/// a hill read from one of the files and not yet added (used when merging)
  struct PendingHill {
    std::unique_ptr<KernelFunctions> kernel;
    double biasf=1.0;
    double clock=0.0;
    double time=0.0;
    bool valid=false;
  };
  std::vector<PendingHill> pending;
  bool mergestarted=false;
  bool unsortedwarning=false;
/// read the next hill of file i into pending[i] (merge mode)
  void readPending(BiasRepresentation *br, unsigned i);
  bool readBunchMerged(BiasRepresentation *br, int stride);
public:
  FilesHandler(const std::vector<std::string> &filenames,
               bool parallelread,
               Action &myaction,
               Log &mylog,
               const std::string & sortby="");
  bool readBunch(BiasRepresentation *br, int stride=-1);
  bool scanOneHill(BiasRepresentation *br, IFile *ifile );
  void getMinMaxBin(const std::vector<Value*> & vals,
                    Communicator &cc,
                    std::vector<double> &vmin,
                    std::vector<double> &vmax,
                    std::vector<unsigned> &vbin);
  void getMinMaxBin(const std::vector<Value*> & vals,
                    Communicator &cc,
                    std::vector<double> &vmin,
                    std::vector<double> &vmax,
                    std::vector<unsigned> &vbin,
                    const std::vector<double> &histosigma);
};
FilesHandler::FilesHandler(const std::vector<std::string> &fnames,
                           const bool doparallelread,
                           Action &act,
                           Log &mylog,
                           const std::string & sort)
  :filenames(fnames),
   log(&mylog),
   parallelread(doparallelread),
   sortby(sort) {
  plumed_massert(sortby=="" || sortby=="time" || sortby=="clock","files can only be sorted by time or clock");
  action=&act;
  for(unsigned i=0; i<filenames.size(); i++) {
    auto ifile=Tools::make_unique<IFile>();
    ifile->link(act);
    plumed_massert((ifile->FileExist(filenames[i])), "the file "+filenames[i]+" does not exist " );
    ifiles.emplace_back(std::move(ifile));
  }

}

// note that the FileHandler is completely transparent respect to the biasrepresentation
// no check are made at this level
void FilesHandler::readPending(BiasRepresentation *br, unsigned i) {
  IFile *ff=ifiles[i].get();
  PendingHill & p=pending[i];
  double t;
  if(!ff->scanField("time",t)) {
    (*log)<<"  closing file "<<filenames[i]<<"\n";
    ff->close();
    p.valid=false;
    p.kernel.reset();
    return;
  }
  double c=0.0;
  if(ff->FieldExist("clock")) {
    ff->scanField("clock",c);
  } else if(sortby=="clock") {
    plumed_merror("cannot sort by clock: there is no clock field in file "+filenames[i]);
  }
  if(p.valid && !unsortedwarning) {
    const bool back = (sortby=="clock") ? (c<p.clock || (c==p.clock && t<p.time)) : (t<p.time);
    if(back) {
      (*log)<<"  WARNING: file "<<filenames[i]<<" is not ordered by "<<sortby<<": hills are merged in the order in which they appear in each file\n";
      unsortedwarning=true;
    }
  }
  p.time=t;
  p.clock=c;
  p.kernel=br->readKernel(ff,p.biasf);
  if(br->hasSigmaInInput()) {
    ff->allowIgnoredFields();
  }
  ff->scanField();
  p.valid=true;
}

bool FilesHandler::readBunchMerged(BiasRepresentation *br, int stride) {
  if(!mergestarted) {
    (*log)<<"  merging the hills of "<<filenames.size()<<" files in order of "<<sortby<<"\n";
    pending.resize(ifiles.size());
    for(unsigned i=0; i<ifiles.size(); i++) {
      (*log)<<"  opening file "<<filenames[i]<<"\n";
      ifiles[i]->open(filenames[i]);
      readPending(br,i);
    }
    mergestarted=true;
  }
  while(true) {
    // the next hill is the earliest; ties are resolved with time (for clock) and then with the order of the files
    int next=-1;
    for(unsigned i=0; i<pending.size(); i++) {
      if(!pending[i].valid) {
        continue;
      }
      if(next<0) {
        next=i;
        continue;
      }
      const auto & a=pending[i];
      const auto & b=pending[next];
      const bool earlier = (sortby=="clock") ? (a.clock<b.clock || (a.clock==b.clock && a.time<b.time)) : (a.time<b.time);
      if(earlier) {
        next=i;
      }
    }
    if(next<0) {
      (*log)<<"  final chunk: now with "<<br->getNumberOfKernels()<<" kernels  \n";
      return false;
    }
    br->addKernel(std::move(pending[next].kernel),pending[next].biasf,pending[next].time);
    readPending(br,next);
    int n=br->getNumberOfKernels();
    if(stride>0 && n%stride==0) {
      // stop here unless all files are over
      for(const auto & p : pending)
        if(p.valid) {
          (*log)<<"  done with this chunk: now with "<<n<<" kernels  \n";
          return true;
        }
    }
  }
}

bool FilesHandler::readBunch(BiasRepresentation *br, int stride) {
  if(sortby!="" && ifiles.size()>1) {
    return readBunchMerged(br,stride);
  }
  bool morefiles=true;
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
      ff->open(filenames[beingread]);
      isopen=true;
    }
    int n=0;
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
          ff=ifiles[beingread].get();
          ff->open(filenames[beingread]);
          (*log)<<"  opening file "<<filenames[beingread]<<"\n";
          isopen=true;
        } else {
          morefiles=false;
          (*log)<<"  final chunk: now with "<<n<<" kernels  \n";
          break;
        }
      }
      // if there are no more files to read and this file is over then quit
      if(fileisover && !morefiles) {
        break;
      }
      // if you are in the middle of a file and you are here
      // then means that you read what you need to read
      if(!fileisover ) {
        break;
      }
    }
  }
  return morefiles;
}
void FilesHandler::getMinMaxBin(const std::vector<Value*> & vals,
                                Communicator &cc,
                                std::vector<double> &vmin,
                                std::vector<double> &vmax,
                                std::vector<unsigned> &vbin) {
  // create the representation (no grid)
  BiasRepresentation br(vals,cc);
  // read all the kernels
  readBunch(&br);
  // loop over the kernels and get the support
  br.getMinMaxBin(vmin,vmax,vbin);
}
void FilesHandler::getMinMaxBin(const std::vector<Value*> & vals,
                                Communicator &cc,
                                std::vector<double> &vmin,
                                std::vector<double> &vmax,
                                std::vector<unsigned> &vbin,
                                const std::vector<double> &histosigma) {
  BiasRepresentation br(vals,cc,histosigma);
  // read all the kernels
  readBunch(&br);
  // loop over the kernels and get the support
  br.getMinMaxBin(vmin,vmax,vbin);
  //for(unsigned i=0;i<vals.size();i++){cerr<<"XXX "<<vmin[i]<<" "<<vmax[i]<<" "<<vbin[i]<<"\n";}
}
bool FilesHandler::scanOneHill(BiasRepresentation *br, IFile *ifile ) {
  double time;
  if(ifile->scanField("time",time)) {
    // the other fields (e.g. clock) are ignored
    br->pushKernel(ifile,time);
    //(*log)<<"  read hill\n";
    if(br->hasSigmaInInput()) {
      ifile->allowIgnoredFields();
    }
    ifile->scanField();
    return true;
  } else {
    return false;
  }
}


double  mylog( double v1 ) {
  return std::log(v1);
}

double  mylogder( double v1 ) {
  return 1./v1;
}



class FuncSumHills :
  public Function {
  std::vector<std::string> hillsFiles,histoFiles;
  std::vector<std::string> proj;
  int initstride=-1;
  bool iscltool=false;
  bool integratehills=false;
  bool integratehisto=false;
  bool parallelread=false;
  bool negativebias=false;
  bool nohistory=false;
  bool minTOzero=false;
  bool doInt=false;
  double lowI_=-1.0;
  double uppI_=-1.0;
  double beta=-1.0;
  std::string outhills;
  std::string outhisto;
  std::string fmt="%14.9f";
  std::string sortby;
/// time average of the bias
  std::vector<double> average;
/// boxes (min and max for each variable) defining the basins whose free energies are written
  std::vector<std::vector<double>> basinmin,basinmax;
  std::string outbasins;
/// write the free energy of each basin relative to the first one
  void writeBasins(OFile & ofile, const Grid & fes, double time, int nkernels);
  std::unique_ptr<BiasRepresentation> biasrep;
  std::unique_ptr<BiasRepresentation> historep;
public:
  explicit FuncSumHills(const ActionOptions&);
  void calculate() override; // this probably is not needed
  bool checkFilesAreExisting(const std::vector<std::string> & hills );
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(FuncSumHills,"FUNCSUMHILLS")

void FuncSumHills::registerKeywords(Keywords& keys) {
  Function::registerKeywords(keys);
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
  keys.addFlag("ISCLTOOL",false,"use via plumed command line: calculate at read phase and then go");
  keys.addFlag("PARALLELREAD",false,"read parallel HILLS file");
  keys.addFlag("NEGBIAS",false,"dump  negative bias ( -bias )   instead of the free energy: needed in well tempered with flexible hills ");
  keys.addFlag("NOHISTORY",false,"to be used with INITSTRIDE:  it splits the bias/histogram in pieces without previous history  ");
  keys.addFlag("MINTOZERO",false,"translate the resulting bias/histogram to have the minimum to zero  ");
  keys.add("optional","FMT","the format that should be used to output real numbers");
  keys.add("optional","SORT","merge the hills of multiple files in order of time or clock instead of reading the files one after the other");
  keys.add("optional","AVERAGE","time average of the bias: from time t0 (and until t1 if given) the hills are scaled by (t1-t)/(t1-t0), where t1 is by default the time of the last hill");
  keys.add("optional","BASINS","basins (separated by colons) given by the minimum and maximum of each variable of the output, e.g. -3,-1:1,3; their free energies relative to the first basin are written at each output");
  keys.add("optional","OUTBASINS","output file for the free energies of the basins");
  keys.setValueDescription("scalar","a scalar");
}

FuncSumHills::FuncSumHills(const ActionOptions&ao):
  Action(ao),
  Function(ao) {
  // format
  parse("FMT",fmt);
  log<<"  Output format is "<<fmt<<"\n";
  // here read
  // Grid Stuff
  std::vector<std::string> gmin;
  parseVector("GRID_MIN",gmin);
  if(gmin.size()!=getNumberOfArguments() && gmin.size()!=0) {
    error("not enough values for GRID_MIN");
  }
  plumed_massert(gmin.size()==getNumberOfArguments() || gmin.size()==0,"need GRID_MIN argument for this") ;
  std::vector<std::string> gmax;
  parseVector("GRID_MAX",gmax);
  if(gmax.size()!=getNumberOfArguments() && gmax.size()!=0) {
    error("not enough values for GRID_MAX");
  }
  plumed_massert(gmax.size()==getNumberOfArguments() || gmax.size()==0,"need GRID_MAX argument for this") ;
  std::vector<unsigned> gbin;
  std::vector<double>   gspacing;
  parseVector("GRID_BIN",gbin);
  plumed_massert(gbin.size()==getNumberOfArguments() || gbin.size()==0,"need GRID_BIN argument for this") ;
  if(gbin.size()!=getNumberOfArguments() && gbin.size()!=0) {
    error("not enough values for GRID_BIN");
  }
  parseVector("GRID_SPACING",gspacing);
  plumed_massert(gspacing.size()==getNumberOfArguments() || gspacing.size()==0,"need GRID_SPACING argument for this") ;
  if(gspacing.size()!=getNumberOfArguments() && gspacing.size()!=0) {
    error("not enough values for GRID_SPACING");
  }
  if(gspacing.size()!=0 && gbin.size()==0) {
    log<<"  The number of bins will be estimated from GRID_SPACING\n";
  } else if(gspacing.size()!=0 && gbin.size()!=0) {
    log<<"  You specified both GRID_BIN and GRID_SPACING\n";
    log<<"  The more conservative (highest) number of bins will be used for each variable\n";
  }
  // the number of bins from GRID_SPACING is computed once the boundaries are known,
  // possibly after they have been detected from the hills
  const bool userbin=(gbin.size()!=0);
  auto applySpacing=[&]() {
    if(gspacing.size()==0) {
      return;
    }
    if(gbin.size()==0) {
      gbin.assign(gmin.size(),1);
    }
    for(unsigned i=0; i<gmin.size(); i++) {
      double a,b;
      Tools::convert(gmin[i],a);
      Tools::convert(gmax[i],b);
      unsigned n=std::ceil((b-a)/gspacing[i]);
      if(gbin[i]<n) {
        gbin[i]=n;
      }
    }
  };
  if(gmin.size()!=0 && gmax.size()!=0) {
    applySpacing();
  }

  // Inteval keyword
  std::vector<double> tmpI(2);
  parseVector("INTERVAL",tmpI);
  if(tmpI.size()!=2&&tmpI.size()!=0) {
    error("both a lower and an upper limits must be provided with INTERVAL");
  } else if(tmpI.size()==2) {
    lowI_=tmpI.at(0);
    uppI_=tmpI.at(1);
    if(getNumberOfArguments()!=1) {
      error("INTERVAL limits correction works only for monodimensional metadynamics!");
    }
    if(uppI_<lowI_) {
      error("The Upper limit must be greater than the Lower limit!");
    }
    doInt=true;
  }
  if(doInt) {
    log << "  Upper and Lower limits boundaries for the bias are activated at " << lowI_ << " - " << uppI_<<"\n";
    log << "  Using the same values as boundaries for the grid if not other value was defined (default: 200 bins)\n";
    std::ostringstream strsmin, strsmax;
    strsmin << lowI_;
    strsmax << uppI_;
    if(gmin.size()==0) {
      gmin.push_back(strsmin.str());
    }
    if(gmax.size()==0) {
      gmax.push_back(strsmax.str());
    }
    if(gbin.size()==0 && gspacing.size()==0) {
      gbin.push_back(200);
    } else if(gbin.size()==0) {
      applySpacing();
    }
  }


  // hills file:
  parseVector("HILLSFILES",hillsFiles);
  if(hillsFiles.size()==0) {
    integratehills=false; // default behaviour
  } else {
    integratehills=true;
    for(unsigned i=0; i<hillsFiles.size(); i++) {
      log<<"  hillsfile  : "<<hillsFiles[i]<<"\n";
    }
  }
  // histo file:
  parseVector("HISTOFILES",histoFiles);
  if(histoFiles.size()==0) {
    integratehisto=false;
  } else {
    integratehisto=true;
    for(unsigned i=0; i<histoFiles.size(); i++) {
      log<<"  histofile  : "<<histoFiles[i]<<"\n";
    }
  }
  std::vector<double> histoSigma;
  if(integratehisto) {
    parseVector("HISTOSIGMA",histoSigma);
    for(unsigned i=0; i<histoSigma.size(); i++) {
      log<<"  histosigma  : "<<histoSigma[i]<<"\n";
    }
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
    for(unsigned i=0; i<getNumberOfArguments(); i++) {
      proj.push_back(getPntrToArgument(i)->getName());
    }
  }

  // add some automatic hills width: not in case stride is defined
  // since when you start from zero the automatic size will be zero!
  if(gmin.size()==0 || gmax.size()==0) {
    log<<"   \n";
    log<<"  No boundaries defined: need to do a prescreening of hills \n";
    std::vector<Value*> tmphillsvalues, tmphistovalues;
    if(integratehills) {
      for(unsigned i=0; i<getNumberOfArguments(); i++) {
        tmphillsvalues.push_back( getPntrToArgument(i) );
      }
    }
    if(integratehisto) {
      for(unsigned i=0; i<getNumberOfArguments(); i++) {
        std::string ss = getPntrToArgument(i)->getName();
        for(unsigned j=0; j<proj.size(); j++) {
          if(proj[j]==ss) {
            tmphistovalues.push_back( getPntrToArgument(i) );
          }
        }
      }
    }

    if(integratehills) {
      FilesHandler hillsHandler(hillsFiles,parallelread,*this, log);
      std::vector<double> vmin,vmax;
      std::vector<unsigned> vbin;
      hillsHandler.getMinMaxBin(tmphillsvalues,comm,vmin,vmax,vbin);
      log<<"  found boundaries from hillsfile: \n";
      gmin.resize(vmin.size());
      gmax.resize(vmax.size());
      for(unsigned i=0; i<getNumberOfArguments(); i++) {
        Tools::convert(vmin[i],gmin[i]);
        Tools::convert(vmax[i],gmax[i]);
      }
      if(userbin || gspacing.size()!=0) {
        log<<"  found nbins or spacing in input, this overrides the automatic choice \n";
        applySpacing();
      } else {
        gbin=vbin;
      }
      for(unsigned i=0; i<getNumberOfArguments(); i++) {
        log<<"  variable "<< getPntrToArgument(i)->getName()<<" min: "<<gmin[i]<<" max: "<<gmax[i]<<" nbin: "<<gbin[i]<<"\n";
      }
    }
    // if at this stage bins are not there then do it with histo
    if(gmin.size()==0) {
      FilesHandler histoHandler(histoFiles,parallelread,*this, log);
      std::vector<double> vmin,vmax;
      std::vector<unsigned> vbin;
      histoHandler.getMinMaxBin(tmphistovalues,comm,vmin,vmax,vbin,histoSigma);
      log<<"  found boundaries from histofile: \n";
      gmin.resize(vmin.size());
      gmax.resize(vmax.size());
      for(unsigned i=0; i<proj.size(); i++) {
        Tools::convert(vmin[i],gmin[i]);
        Tools::convert(vmax[i],gmax[i]);
      }
      if(userbin || gspacing.size()!=0) {
        log<<"  found nbins or spacing in input, this overrides the automatic choice \n";
        applySpacing();
      } else {
        gbin=vbin;
      }
      for(unsigned i=0; i<proj.size(); i++) {
        log<<"  variable "<< proj[i] <<" min: "<<gmin[i]<<" max: "<<gmax[i]<<" nbin: "<<gbin[i]<<"\n";
      }
    }
    log<<"  done!\n";
    log<<"   \n";
  }


  parseVector("AVERAGE",average);
  if(average.size()!=0) {
    if(average.size()>2) {
      error("AVERAGE needs the initial time and optionally the final time");
    }
    if(!integratehills) {
      error("AVERAGE can only be used with hills files");
    }
    if(average.size()==2 && average[1]<=average[0]) {
      error("the final time of AVERAGE should be larger than the initial time");
    }
    log<<"  time average of the bias from time "<<average[0];
    if(average.size()==2) {
      log<<" to time "<<average[1];
    }
    log<<"\n";
  }
  std::string basins;
  parse("BASINS",basins);
  if(basins.length()>0) {
    // the basins are defined on the variables of the output
    const unsigned nout=(proj.size()>0 ? proj.size() : getNumberOfArguments());
    for(const auto & b : Tools::getWords(basins,":")) {
      std::vector<std::string> words=Tools::getWords(b,",");
      if(words.size()!=2*nout) {
        error("each basin needs a minimum and a maximum for each of the "+std::to_string(nout)+" variables of the output: "+b);
      }
      std::vector<double> bmin(nout),bmax(nout);
      for(unsigned i=0; i<nout; i++) {
        if(!Tools::convertNoexcept(words[2*i],bmin[i]) || !Tools::convertNoexcept(words[2*i+1],bmax[i])) {
          error("cannot read basin "+b);
        }
      }
      basinmin.push_back(bmin);
      basinmax.push_back(bmax);
    }
    if(basinmin.size()<2) {
      error("BASINS needs at least two basins");
    }
    log<<"  free energies of "<<basinmin.size()<<" basins relative to the first one are written at each output\n";
  }
  if( proj.size() != 0 || integratehisto==true || basinmin.size()>0 ) {
    parse("KT",beta);
    for(unsigned i=0; i<proj.size(); i++) {
      log<<"  projection "<<i<<" : "<<proj[i]<<"\n";
    }
    // this should be only for projection or free energy from histograms
    plumed_massert(beta>0.,"if you make a projection, a histogram correction or use basins then you need KT flag!");
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
    outhills="fes.dat";
    outhisto="histo.dat";
  } else {
    outhills="fes_";
    outhisto="histo_";
    log<<"  Doing integration slices every "<<initstride<<" kernels\n";
    parseFlag("NOHISTORY",nohistory);
    if(nohistory) {
      log<<"  nohistory: each stride block has no memory of the previous block\n";
    }
  }
  parse("SORT",sortby);
  if(sortby!="" && sortby!="time" && sortby!="clock") {
    error("SORT can only be time or clock");
  }
  if(sortby!="") {
    log<<"  hills from multiple files are merged in order of "<<sortby<<"\n";
  }
  parseFlag("MINTOZERO",minTOzero);
  if(minTOzero) {
    log<<"  mintozero: bias/histogram will be translated to have the minimum value equal to zero\n";
  }
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
          if(proj[j]==ss) {
            tmphistovalues.push_back( getPntrToArgument(i) );
          }
        }
      }
    }

    // check if the files exists
    if(integratehills) {
      checkFilesAreExisting(hillsFiles);
      biasrep=Tools::make_unique<BiasRepresentation>(tmphillsvalues,comm, gmin, gmax, gbin, doInt, lowI_, uppI_);
      if(average.size()>0) {
        biasrep->setTimeAverage(average[0],average.size()==2,(average.size()==2?average[1]:0.0));
      }
      if(negativebias && basinmin.size()>0) {
        error("BASINS requires a free energy and cannot be used with NEGBIAS");
      }
      if(negativebias) {
        biasrep->setRescaledToBias(true);
        log<<"  required the -bias instead of the free energy \n";
        if(initstride<0) {
          outhills="negativebias.dat";
        } else {
          outhills="negativebias_";
        }
      }
    }

    parse("OUTHILLS",outhills);
    parse("OUTHISTO",outhisto);
    outbasins="basins.dat";
    parse("OUTBASINS",outbasins);
    if(integratehills) {
      log<<"  output file for fes/bias  is :  "<<outhills<<"\n";
    }
    if(integratehisto) {
      log<<"  output file for histogram is :  "<<outhisto<<"\n";
    }
    checkRead();

    log<<"\n";
    log<<"  Now calculating...\n";
    log<<"\n";

    // here it defines the column to be histogrammed, tmpvalues should be only
    // the list of the collective variable one want to consider
    if(integratehisto) {
      checkFilesAreExisting(histoFiles);
      historep=Tools::make_unique<BiasRepresentation>(tmphistovalues,comm,gmin,gmax,gbin,histoSigma);
    }

    // decide how to source hills ( serial/parallel )
    // here below the input control
    // say how many hills and it will read them from the
    // bunch of files provided, will update the representation
    // of hills (i.e. a list of hills and the associated grid)

    // decide how to source colvars ( serial parallel )
    std::unique_ptr<FilesHandler> hillsHandler;
    std::unique_ptr<FilesHandler> histoHandler;

    if(integratehills)	{
      hillsHandler=Tools::make_unique<FilesHandler>(hillsFiles,parallelread,*this, log, sortby);
    }
    if(integratehisto)	{
      histoHandler=Tools::make_unique<FilesHandler>(histoFiles,parallelread,*this, log, sortby);
    }

// Stopwatch is logged when it goes out of scope
    OFile basinsfile, histobasinsfile;
    if(basinmin.size()>0) {
      if(integratehills) {
        basinsfile.link(*this);
        basinsfile.open(outbasins);
        basinsfile.fmtField(fmt);
        log<<"  output file for the free energies of the basins is :  "<<outbasins<<"\n";
      }
      if(integratehisto) {
        histobasinsfile.link(*this);
        histobasinsfile.open("histo_"+outbasins);
        histobasinsfile.fmtField(fmt);
        log<<"  output file for the free energies of the basins from the histogram is :  histo_"<<outbasins<<"\n";
      }
    }

    Stopwatch sw(log);

// Stopwatch is stopped when swh goes out of scope
    auto swh=sw.startStop("0 Summing hills");

    // read a number of hills and put in the bias representation
    int nfiles=0;
    bool ibias=integratehills;
    bool ihisto=integratehisto;
    while(true) {
      if(  integratehills  && ibias  ) {
        if(nohistory) {
          biasrep->clear();
          log<<"  clearing history before reading a new block\n";
        };
        log<<"  reading hills: \n";
        const int before=biasrep->getNumberOfKernels();
        ibias=hillsHandler->readBunch(biasrep.get(),initstride) ;
        log<<"\n";
        // when the number of hills is a multiple of the stride the last block is empty: do not write it again
        if(!ibias && nfiles>0 && biasrep->getNumberOfKernels()==before) {
          integratehills=false;
        }
      }

      if(  integratehisto  && ihisto ) {
        if(nohistory) {
          historep->clear();
          log<<"  clearing history before reading a new block\n";
        };
        log<<"  reading histogram: \n";
        const int before=historep->getNumberOfKernels();
        ihisto=histoHandler->readBunch(historep.get(),initstride) ;
        log<<"\n";
        if(!ihisto && nfiles>0 && historep->getNumberOfKernels()==before) {
          integratehisto=false;
        }
      }

      // dump: need to project?
      if(proj.size()!=0) {

        if(integratehills) {

          log<<"  Bias: Projecting on subgrid... \n";
          BiasWeight Bw(beta);
          Grid biasGrid=*(biasrep->getGridPtr());
          Grid smallGrid=biasGrid.project(proj,&Bw);
          OFile gridfile;
          gridfile.link(*this);
          std::ostringstream ostr;
          ostr<<nfiles;
          std::string myout;
          if(initstride>0) {
            myout=outhills+ostr.str()+".dat" ;
          } else {
            myout=outhills;
          }
          log<<"  Bias: Writing subgrid on file "<<myout<<" \n";
          gridfile.open(myout);
          if(minTOzero) {
            smallGrid.setMinToZero();
          }
          smallGrid.setOutputFmt(fmt);
          smallGrid.writeToFile(gridfile);
          if(basinmin.size()>0) {
            writeBasins(basinsfile,smallGrid,biasrep->getLastTime(),biasrep->getNumberOfKernels());
          }
          gridfile.close();
          if(!ibias) {
            integratehills=false;  // once you get to the final bunch just give up
          }
        }
        // this should be removed
        if(integratehisto) {

          log<<"  Histo: Projecting on subgrid... \n";
          Grid histoGrid=*(historep->getGridPtr());

          OFile gridfile;
          gridfile.link(*this);
          std::ostringstream ostr;
          ostr<<nfiles;
          std::string myout;
          if(initstride>0) {
            myout=outhisto+ostr.str()+".dat" ;
          } else {
            myout=outhisto;
          }
          log<<"  Histo: Writing subgrid on file "<<myout<<" \n";
          gridfile.open(myout);

          histoGrid.applyFunctionAllValuesAndDerivatives(&mylog,&mylogder);
          histoGrid.scaleAllValuesAndDerivatives(-1./beta);
          if(minTOzero) {
            histoGrid.setMinToZero();
          }
          histoGrid.setOutputFmt(fmt);
          histoGrid.writeToFile(gridfile);
          if(basinmin.size()>0) {
            writeBasins(histobasinsfile,histoGrid,historep->getLastTime(),historep->getNumberOfKernels());
          }

          if(!ihisto) {
            integratehisto=false;  // once you get to the final bunch just give up
          }
        }

      } else {

        if(integratehills) {

          Grid biasGrid=*(biasrep->getGridPtr());
          biasGrid.scaleAllValuesAndDerivatives(-1.);

          OFile gridfile;
          gridfile.link(*this);
          std::ostringstream ostr;
          ostr<<nfiles;
          std::string myout;
          if(initstride>0) {
            myout=outhills+ostr.str()+".dat" ;
          } else {
            myout=outhills;
          }
          log<<"  Writing full grid on file "<<myout<<" \n";
          gridfile.open(myout);

          if(minTOzero) {
            biasGrid.setMinToZero();
          }
          biasGrid.setOutputFmt(fmt);
          biasGrid.writeToFile(gridfile);
          if(basinmin.size()>0) {
            writeBasins(basinsfile,biasGrid,biasrep->getLastTime(),biasrep->getNumberOfKernels());
          }
          // rescale back prior to accumulate
          if(!ibias) {
            integratehills=false;  // once you get to the final bunch just give up
          }
        }
        if(integratehisto) {

          Grid histoGrid=*(historep->getGridPtr());
          // do this if you want a free energy from a grid, otherwise do not
          histoGrid.applyFunctionAllValuesAndDerivatives(&mylog,&mylogder);
          histoGrid.scaleAllValuesAndDerivatives(-1./beta);

          OFile gridfile;
          gridfile.link(*this);
          std::ostringstream ostr;
          ostr<<nfiles;
          std::string myout;
          if(initstride>0) {
            myout=outhisto+ostr.str()+".dat" ;
          } else {
            myout=outhisto;
          }
          log<<"  Writing full grid on file "<<myout<<" \n";
          gridfile.open(myout);

          // also this is useful only for free energy
          if(minTOzero) {
            histoGrid.setMinToZero();
          }
          histoGrid.setOutputFmt(fmt);
          histoGrid.writeToFile(gridfile);
          if(basinmin.size()>0) {
            writeBasins(histobasinsfile,histoGrid,historep->getLastTime(),historep->getNumberOfKernels());
          }

          if(!ihisto) {
            integratehisto=false;  // once you get to the final bunch just give up
          }
        }
      }
      if ( !ibias && !ihisto) {
        break;  //when both are over then just quit
      }

      nfiles++;
    }

    return;
  }
  // just an initialization but you need to do something on the fly?: need to connect with a metad run and its grid representation
  // your argument is a metad run
  // if the grid does not exist crash and say that you need some data
  // otherwise just link with it

}

void FuncSumHills::writeBasins(OFile & ofile, const Grid & fes, double time, int nkernels) {
  const unsigned nb=basinmin.size();
  const unsigned nd=fes.getDimension();
  const auto pbc=fes.getIsPeriodic();
  // -kT log sum_s exp(-F(s)/kT) over the points of each basin, computed with a running maximum
  std::vector<double> fmax(nb,-std::numeric_limits<double>::max()), sum(nb,0.0);
  std::vector<double> x(nd);
  for(Grid::index_t i=0; i<fes.getSize(); i++) {
    const double mbf=-beta*fes.getValue(i);
    if(!std::isfinite(mbf)) {
      continue;
    }
    fes.getPoint(i,x);
    for(unsigned b=0; b<nb; b++) {
      bool inside=true;
      for(unsigned d=0; d<nd && inside; d++) {
        const double lo=basinmin[b][d];
        const double hi=basinmax[b][d];
        if(pbc[d] && lo>hi) {
          // along a periodic variable a basin can cross the boundary
          inside=(x[d]>=lo || x[d]<=hi);
        } else {
          inside=(x[d]>=lo && x[d]<=hi);
        }
      }
      if(!inside) {
        continue;
      }
      if(mbf>fmax[b]) {
        sum[b]=sum[b]*std::exp(fmax[b]-mbf)+1.0;
        fmax[b]=mbf;
      } else {
        sum[b]+=std::exp(mbf-fmax[b]);
      }
    }
  }
  std::vector<double> f(nb);
  for(unsigned b=0; b<nb; b++) {
    plumed_massert(sum[b]>0.0,"basin "+std::to_string(b)+" does not contain any point of the grid");
    f[b]=-(std::log(sum[b])+fmax[b])/beta;
  }
  ofile.printField("time",time);
  ofile.printField("nkernels",nkernels);
  for(unsigned b=1; b<nb; b++) {
    ofile.printField("deltaf_"+std::to_string(b),f[b]-f[0]);
  }
  ofile.printField();
}

void FuncSumHills::calculate() {
  // this should be connected only with a grid representation to metadynamics
  // at regular time just dump it
  plumed_merror("You should have never got here: this stuff is not yet implemented!");
}

bool FuncSumHills::checkFilesAreExisting(const std::vector<std::string> & hills ) {
  plumed_massert(hills.size()!=0,"the number of  files provided should be at least one" );
  auto ifile=Tools::make_unique<IFile>();
  ifile->link(*this);
  for(unsigned i=0; i< hills.size(); i++) {
    plumed_massert(ifile->FileExist(hills[i]),"missing file "+hills[i]);
  }
  return true;

}

}

}


