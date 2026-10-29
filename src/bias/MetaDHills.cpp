/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2026 The plumed team
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
#include "MetaDHills.h"
#include "MetaDGrid.h"
#include "core/ActionWithValue.h"
#include "core/Value.h"
#include "tools/Exception.h"
#include "tools/File.h"
#include "tools/Log.h"
#include "tools/Matrix.h"
#include "tools/Tools.h"
#include <cstdio>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <exception>
#include <filesystem>
#include <fstream>
#include <system_error>

namespace PLMD {
namespace bias {

std::string exactString(double x) {
  char buffer[64];
  for(int precision=15; precision<17; precision++) {
    std::snprintf(buffer,sizeof(buffer),"%.*g",precision,x);
    if(std::strtod(buffer,nullptr)==x) {
      return buffer;
    }
  }
  std::snprintf(buffer,sizeof(buffer),"%.17g",x);
  return buffer;
}

void bcastString(Communicator& c, std::string& s) {
  unsigned n=s.size();
  c.Bcast(n,0);
  s.resize(n);
  if(n>0) {
    c.Bcast(s,0);
  }
}

std::string removeIncompleteLastLine(const std::string& path) {
  // compressed files cannot be fixed in place
  if(path.length()>3 && path.substr(path.length()-3)==".gz") {
    return "";
  }
  std::ifstream f(path,std::ios::binary);
  if(!f) {
    return "";
  }
  f.seekg(0,std::ios::end);
  const std::streamoff size=f.tellg();
  if(size<=0) {
    return "";
  }
  char last;
  f.seekg(size-1);
  f.get(last);
  if(last=='\n') {
    return "";
  }
  // look backward for the end of the last complete line
  std::streamoff pos=size-1;
  while(pos>0) {
    char c;
    f.seekg(pos-1);
    f.get(c);
    if(c=='\n') {
      break;
    }
    pos--;
  }
  std::string tail(size-pos,' ');
  f.seekg(pos);
  f.read(&tail[0],size-pos);
  f.close();
  std::string message="  WARNING: the last line of "+path+" is incomplete (the previous run was probably interrupted while writing it) and is removed:\n";
  message+="  "+tail+"\n";
  std::error_code ec;
  std::filesystem::resize_file(path,pos,ec);
  if(ec) {
    message+="  WARNING: the incomplete last line of "+path+" could not be removed ("+ec.message()+")\n";
  }
  return message;
}


MetaDHillsFiles::MetaDHillsFiles(ActionWithValue& action, Communicator& comm, const std::vector<Value*>& args):
  action_(action),
  comm_(comm),
  args_(args) {
}

MetaDHillsFiles::~MetaDHillsFiles() {
}

std::string MetaDHillsFiles::getWalkerFileName(const std::string& name, const std::string& dir, int nwalkers, int walker, bool mpi) {
  const std::string suffix=(nwalkers>1 ? "."+std::to_string(walker) : std::string());
  if(dir!="" && (nwalkers>1 || mpi)) {
    return dir+"/"+name+suffix;
  }
  return name+suffix;
}

unsigned MetaDHillsFiles::addInput(const std::string& name, bool nosuffix) {
  ifiles_.emplace_back(Tools::make_unique<IFile>());
  ifilesnames_.push_back(name);
  nread_.push_back(0);
  IFile* ifile=ifiles_.back().get();
  ifile->link(action_);
  if(nosuffix) {
    ifile->enforceSuffix("");
  }
  // the file operations must not synchronize the processes, as with comm they would do
  ifile->link(selfcomm_);
  const unsigned i=ifiles_.size()-1;
  openIfExists(i);
  return i;
}

bool MetaDHillsFiles::isOpen(unsigned i) {
  int isopen=ifiles_[i]->isOpen();
  if(comm_.Get_size()>1) {
    comm_.Bcast(isopen,0);
  }
  return isopen;
}

std::string MetaDHillsFiles::getPath(unsigned i) const {
  return ifiles_[i]->isOpen() ? ifiles_[i]->getPath() : std::string();
}

bool MetaDHillsFiles::openIfExists(unsigned i) {
  int exists=0;
  if(comm_.Get_rank()==0) {
    exists=ifiles_[i]->FileExist(ifilesnames_[i]);
  }
  if(comm_.Get_size()>1) {
    comm_.Bcast(exists,0);
  }
  if(exists && comm_.Get_rank()==0) {
    ifiles_[i]->open(ifilesnames_[i]);
  }
  return exists;
}

unsigned MetaDHillsFiles::readHills(unsigned i, double* lasttime, bool addtobias, double skipupto, unsigned long long skipcount) {
  Log& log=action_.log;
  // the callers have started a line of the log
  if(bias_.grid) {
    bias_.grid->setMidline(true);
  }
  const std::string path=getPath(i);
  unsigned long long nskipped=0;
  unsigned nhills=0;
  try {
    nhills=read(i,[&](double time,bool multivariate,double height,const std::vector<double>& center,const std::vector<double>& sigma) {
      if(lasttime) {
        *lasttime=time;
      }
      if(lasttime && bias_.diffrestart && bias_.diffrestart->pending) {
        // the own hills file is being read: keep the last hills for the statistics of ADAPTIVE=DIFF
        bias_.diffrestart->add(path,time,multivariate,center,sigma);
      }
      if(!addtobias) {
        return;
      }
      if(nskipped<skipcount || time<=skipupto+bias_.tolerance) {
        nskipped++;
        return;
      }
      bias_.add(MetaDGaussian(multivariate,height*bias_.heightscale,center,sigma));
    });
  } catch(...) {
    if(bias_.grid) {
      bias_.grid->setMidline(false);
    }
    throw;
  }
  if(bias_.grid) {
    bias_.grid->setMidline(false);
  }
  if(!addtobias) {
    log.printf("      %u Gaussians scanned (the bias is read from the grid)\n",nhills);
  } else if(nskipped>0) {
    log.printf("      %u Gaussians read, %llu of them already in the grid\n",nhills,nskipped);
  } else {
    log.printf("      %u Gaussians read\n",nhills);
  }
  return nhills;
}

void MetaDHillsFiles::reset(unsigned i) {
  ifiles_[i]->reset(false);
}

void MetaDHillsFiles::close(unsigned i) {
  ifiles_[i]->close();
}

void MetaDHillsFiles::inputError(const std::string& msg) const {
  // as Action::error, but without writing on the log: the files are read by the first process only, and
  // the log requires the same output on all the processes (read throws the error on all of them)
  plumed_merror("ERROR in input to action "+action_.getName()+" with label "+action_.getLabel()+" : "+msg);
}

void MetaDHillsFiles::noStretchWarning() {
  if(!noStretchWarningDone_) {
    action_.log<<"\nWARNING: you are using a HILLS file with Gaussian kernels, PLUMED 2.8 uses stretched Gaussians by default\n";
  }
  noStretchWarningDone_=true;
}

unsigned MetaDHillsFiles::read(unsigned i, const Callback& f) {
  IFile* ifile=ifiles_[i].get();
  const unsigned ncv=args_.size();
  std::vector<double> center(ncv);
  std::vector<double> sigma(ncv);
  double height;
  unsigned nhills=0;
  bool multivariate=false;

  std::vector<Value> tmpvalues;
  for(unsigned j=0; j<ncv; ++j) {
    tmpvalues.push_back( Value( &action_, args_[j]->getName(), false ) );
  }

  // The hills are read and broadcast in chunks; each one is stored as
  // time, multivariate, height, center (ncv), sigma (ncv or ncv*(ncv+1)/2).
  const unsigned nsigmamax=ncv*(ncv+1)/2;
  const unsigned recordsize=3+ncv+nsigmamax;
  const unsigned chunksize=4096;
  const bool parallel=(comm_.Get_size()>1);
  const bool reader=(comm_.Get_rank()==0);
  std::vector<double> buffer;
  bool more=true;
  while(more) {
    unsigned nread=0;
    if(reader) {
      buffer.clear();
      int failed=0;
      std::string message;
      std::exception_ptr exception;
      try {
        double time;
        while(nread<chunksize && (more=scanOneHill(ifile,tmpvalues,center,sigma,height,multivariate,time))) {
          buffer.push_back(time);
          buffer.push_back(multivariate ? 1.0 : 0.0);
          buffer.push_back(height);
          buffer.insert(buffer.end(),center.begin(),center.end());
          buffer.insert(buffer.end(),sigma.begin(),sigma.end());
          buffer.resize(buffer.size()+nsigmamax-sigma.size(),0.0);
          nread++;
        }
      } catch(const std::exception& e) {
        failed=1;
        message=e.what();
        exception=std::current_exception();
      }
      if(parallel) {
        // the other processes must know if reading failed, so that all of them stop
        comm_.Bcast(failed,0);
        if(failed) {
          bcastString(comm_,message);
          // all the processes throw the same exception: its message can be written on the log,
          // which requires the same output on all the processes
          throw ExceptionError(message);
        }
      }
      if(failed) {
        std::rethrow_exception(exception);
      }
    } else {
      int failed=0;
      comm_.Bcast(failed,0);
      if(failed) {
        std::string message;
        bcastString(comm_,message);
        throw ExceptionError(message);
      }
    }
    // the warning about Gaussian kernels is written on the log after reading, by all the processes
    // (the log requires the same output on all the processes)
    int gaussiankernels=gaussianKernelsFound_;
    if(parallel) {
      int imore=more;
      comm_.Bcast(nread,0);
      comm_.Bcast(imore,0);
      comm_.Bcast(gaussiankernels,0);
      more=imore;
      buffer.resize(nread*recordsize);
      if(nread>0) {
        comm_.Bcast(buffer,0);
      }
    }
    if(gaussiankernels) {
      noStretchWarning();
    }
    for(unsigned h=0; h<nread; h++) {
      const double* record=&buffer[h*recordsize];
      const double time=record[0];
      multivariate=(record[1]!=0.0);
      height=record[2];
      center.assign(record+3,record+3+ncv);
      sigma.assign(record+3+ncv,record+3+ncv+(multivariate ? nsigmamax : ncv));
      nhills++;
      f(time,multivariate,height,center,sigma);
    }
  }
  nread_[i]+=nhills;
  return nhills;
}

namespace {
const std::string gridCommentHeader="# WALKERS_N: hills read from the hills files of the walkers:";
}

std::string MetaDHillsFiles::getGridComment(unsigned own) const {
  std::string line=gridCommentHeader;
  for(unsigned i=0; i<nread_.size(); i++) {
    line+=" "+(i==own ? std::string("-") : std::to_string(nread_[i]));
  }
  return line+"\n";
}

std::vector<long long> MetaDHillsFiles::readGridComment(Communicator& comm, const std::string& path, unsigned nwalkers) {
  std::vector<long long> counts;
  if(comm.Get_rank()==0) {
    counts=readGridComment(path,nwalkers);
  }
  if(comm.Get_size()>1) {
    // (as doubles: they are exact up to 2^53)
    std::vector<double> buffer(counts.begin(),counts.end());
    unsigned n=buffer.size();
    comm.Bcast(n,0);
    buffer.resize(n);
    if(n>0) {
      comm.Bcast(buffer,0);
    }
    counts.assign(buffer.begin(),buffer.end());
  }
  return counts;
}

std::vector<long long> MetaDHillsFiles::readGridComment(const std::string& path, unsigned nwalkers) {
  std::vector<long long> counts;
  // compressed files are not read (the hills of the other walkers are then selected by time)
  if(path.length()>3 && path.substr(path.length()-3)==".gz") {
    return counts;
  }
  std::ifstream f(path);
  std::string line;
  // the comment follows the grid (with STORE_GRIDS, the first grid is the one that is read)
  while(std::getline(f,line)) {
    if(line.compare(0,gridCommentHeader.length(),gridCommentHeader)!=0) {
      continue;
    }
    std::vector<std::string> words=Tools::getWords(line.substr(gridCommentHeader.length()));
    if(words.size()!=nwalkers) {
      return counts;
    }
    for(const auto & w : words) {
      long long n=-1;
      if(w!="-" && !Tools::convertNoexcept(w,n)) {
        counts.clear();
        return counts;
      }
      counts.push_back(n);
    }
    return counts;
  }
  return counts;
}

/// read the next hill: tmpvalues are used to read the CVs (with their periodicity); returns false at the end of the file
bool MetaDHillsFiles::scanOneHill(IFile* ifile, std::vector<Value>& tmpvalues, std::vector<double>& center, std::vector<double>& sigma,
                                  double& height, bool& multivariate, double& time) {
  double dummy;
  multivariate=false;
  if(ifile->scanField("time",time)) {
    unsigned ncv=tmpvalues.size();
    for(unsigned i=0; i<ncv; ++i) {
      ifile->scanField( &tmpvalues[i] );
      if( tmpvalues[i].isPeriodic() && ! args_[i]->isPeriodic() ) {
        inputError("in hills file periodicity for variable " + tmpvalues[i].getName() + " does not match periodicity in input");
      } else if( tmpvalues[i].isPeriodic() ) {
        std::string imin, imax;
        tmpvalues[i].getDomain( imin, imax );
        std::string rmin, rmax;
        args_[i]->getDomain( rmin, rmax );
        if( imin!=rmin || imax!=rmax ) {
          inputError("in hills file periodicity for variable " + tmpvalues[i].getName() + " does not match periodicity in input");
        }
      }
      center[i]=tmpvalues[i].get();
    }
    // scan for kerneltype
    std::string ktype="stretched-gaussian";
    if( ifile->FieldExist("kerneltype") ) {
      ifile->scanField("kerneltype",ktype);
    }
    if( ktype=="gaussian" ) {
      // the warning is written by read
      gaussianKernelsFound_=true;
    } else if( ktype!="stretched-gaussian") {
      inputError("non Gaussian kernels are not supported in MetaD");
    }
    // scan for multivariate label
    std::string sss;
    ifile->scanField("multivariate",sss);
    if(sss=="true") {
      multivariate=true;
    } else if(sss=="false") {
      multivariate=false;
    } else {
      plumed_merror("cannot parse multivariate = "+ sss);
    }
    if(multivariate) {
      sigma.resize(ncv*(ncv+1)/2);
      Matrix<double> upper(ncv,ncv);
      Matrix<double> lower(ncv,ncv);
      for(unsigned i=0; i<ncv; i++) {
        for(unsigned j=0; j<ncv-i; j++) {
          ifile->scanField("sigma_"+args_[j+i]->getName()+"_"+args_[j]->getName(),lower(j+i,j));
          upper(j,j+i)=lower(j+i,j);
        }
      }
      Matrix<double> mymult(ncv,ncv);
      Matrix<double> invmatrix(ncv,ncv);
      mult(lower,upper,mymult);
      // now invert and get the sigmas
      Invert(mymult,invmatrix);
      // put the sigmas in the usual order: upper diagonal (this time in normal form and not in band form)
      unsigned k=0;
      for(unsigned i=0; i<ncv; i++) {
        for(unsigned j=i; j<ncv; j++) {
          sigma[k]=invmatrix(i,j);
          k++;
        }
      }
    } else {
      for(unsigned i=0; i<ncv; ++i) {
        ifile->scanField("sigma_"+args_[i]->getName(),sigma[i]);
      }
    }

    ifile->scanField("height",height);
    ifile->scanField("biasf",dummy);
    if(ifile->FieldExist("clock")) {
      ifile->scanField("clock",dummy);
    }
    if(ifile->FieldExist("lower_int")) {
      ifile->scanField("lower_int",dummy);
    }
    if(ifile->FieldExist("upper_int")) {
      ifile->scanField("upper_int",dummy);
    }
    ifile->scanField();
    return true;
  } else {
    return false;
  }
}

void MetaDHillsFiles::write(OFile& file, double time, const MetaDGaussian& hill, double height, double biasf, bool clock) const {
  const unsigned ncv=args_.size();
  file.printField("time",time);
  for(unsigned i=0; i<ncv; ++i) {
    file.printField(args_[i],hill.center[i]);
  }
  file.printField("kerneltype","stretched-gaussian");
  if(hill.multivariate) {
    file.printField("multivariate","true");
    Matrix<double> mymatrix(ncv,ncv);
    unsigned k=0;
    for(unsigned i=0; i<ncv; i++) {
      for(unsigned j=i; j<ncv; j++) {
        // recompose the full inverse matrix
        mymatrix(i,j)=mymatrix(j,i)=hill.sigma[k];
        k++;
      }
    }
    // invert it
    Matrix<double> invmatrix(ncv,ncv);
    Invert(mymatrix,invmatrix);
    // enforce symmetry
    for(unsigned i=0; i<ncv; i++) {
      for(unsigned j=i; j<ncv; j++) {
        invmatrix(i,j)=invmatrix(j,i);
      }
    }

    // do cholesky so to have a "sigma like" number
    Matrix<double> lower(ncv,ncv);
    cholesky(invmatrix,lower);
    // loop in band form
    for(unsigned i=0; i<ncv; i++) {
      for(unsigned j=0; j<ncv-i; j++) {
        file.printField("sigma_"+args_[j+i]->getName()+"_"+args_[j]->getName(),lower(j+i,j));
      }
    }
  } else {
    file.printField("multivariate","false");
    for(unsigned i=0; i<ncv; ++i) {
      file.printField("sigma_"+args_[i]->getName(),hill.sigma[i]);
    }
  }
  file.printField("height",height).printField("biasf",biasf);
  if(clock) {
    file.printField("clock",int(std::time(0)));
  }
  file.printField();
}

void MetaDDiffusionRestart::add(const std::string& file, double time, bool multivariate, const std::vector<double>& center,
                                const std::vector<double>& sigma) {
  if(!multivariate) {
    // e.g. a previous run without ADAPTIVE=DIFF: these hills cannot be used
    clear();
  } else {
    time_.push_back(time);
    center_.push_back(center);
    sigma_.push_back(sigma);
    if(time_.size()>maxhills) {
      time_.pop_front();
      center_.pop_front();
      sigma_.pop_front();
    }
  }
  file_=file;
}

unsigned MetaDDiffusionRestart::compute(double tol, double now, unsigned blocksize, unsigned offset, double decay, double updatetime,
                                        const std::function<double(unsigned,double,double)>& difference,
                                        const std::function<double(unsigned,double)>& bringBackInPbc,
                                        std::vector<double>& average, std::vector<double>& variance, double& lasttime) const {
  // the hills of this walker
  std::vector<unsigned> own;
  for(unsigned k=0; k<time_.size();) {
    unsigned l=k;
    while(l<time_.size() && std::abs(time_[l]-time_[k])<tol) {
      l++;
    }
    if(blocksize==0) {
      for(unsigned j=k; j<l; j++) {
        own.push_back(j);
      }
    } else if(l-k==blocksize) {
      own.push_back(k+offset);
    }
    k=l;
  }
  // only the hills deposited up to now (when restarting from a grid, or from a checkpoint older than the
  // end of the hills file, the file can go further)
  {
    std::vector<unsigned> before;
    for(const auto j : own) {
      if(time_[j]<=now+tol) {
        before.push_back(j);
      }
    }
    own.swap(before);
  }
  // only the last sequence of increasing times: the file can contain several runs (e.g. a run restarted
  // from an older checkpoint, or a new simulation appended to an old file)
  if(!own.empty()) {
    unsigned first=own.size()-1;
    while(first>0 && time_[own[first-1]]<time_[own[first]]+tol) {
      first--;
    }
    own.erase(own.begin(),own.begin()+first);
  }
  if(own.empty()) {
    return 0;
  }
  // covariance: the one of the last hill, i.e. the running estimate of the previous run at that time
  const unsigned last=own.back();
  const unsigned ncv=center_[last].size();
  Matrix<double> mymatrix(ncv,ncv);
  unsigned k=0;
  for(unsigned i=0; i<ncv; i++) {
    for(unsigned j=i; j<ncv; j++) {
      mymatrix(i,j)=mymatrix(j,i)=sigma_[last][k];
      k++;
    }
  }
  Matrix<double> mycov(ncv,ncv);
  Invert(mymatrix,mycov);
  variance.resize(ncv*(ncv+1)/2);
  k=0;
  for(unsigned i=0; i<ncv; i++) {
    for(unsigned j=i; j<ncv; j++) {
      variance[k]=mycov(i,j);
      k++;
    }
  }
  // running average: the center of a hill is the instantaneous value of the CVs, which differs from
  // the running average by about one standard deviation, so the centers of the last hills are averaged
  // with the same exponential weights used by the running average
  const std::vector<double>& ref=center_[last];
  average.assign(ncv,0.0);
  double wsum=0.0;
  for(const auto j : own) {
    // age in number of updates of the running average
    const double age=std::max(0.0,(time_[last]-time_[j])/updatetime);
    const double w=std::pow(1.0-decay,age);
    for(unsigned i=0; i<ncv; i++) {
      average[i]+=w*difference(i,ref[i],center_[j][i]);
    }
    wsum+=w;
  }
  for(unsigned i=0; i<ncv; i++) {
    average[i]=bringBackInPbc(i,ref[i]+average[i]/wsum);
  }
  lasttime=time_[last];
  return own.size();
}

}
}
