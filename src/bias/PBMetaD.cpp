/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2015 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#include "Bias.h"
#include "ActionRegister.h"
#include "core/ActionSet.h"
#include "tools/Grid.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/Exception.h"
#include "tools/Random.h"
#include <string>
#include <cstring>
#include "tools/File.h"
#include "time.h"
#include <iostream>
#include <limits>

#define DP2CUTOFF 6.25

using namespace std;


namespace PLMD{
namespace bias{

//+PLUMEDOC BIAS PBMETAD 
/*

Put doc here

*/
//+ENDPLUMEDOC

class PBMetaD : public Bias{

private:
  struct Gaussian {
   vector<double> center;
   vector<double> sigma;
   double height;
   Gaussian(const vector<double> & center,const vector<double> & sigma, double height):
     center(center),sigma(sigma),height(height){}
  };
  vector<double> sigma0_;
  vector< vector<Gaussian> > hills_;
  vector<OFile*> hillsOfiles_;
  vector<IFile*> ifiles;
  vector<Grid*> BiasGrids_;
  bool    grid_;
  double  height0_;
  double  biasf_;
  double  kbt_;
  int     stride_;
  bool    welltemp_;
  bool    multiple_w;
  vector<double> uppI_;
  vector<double> lowI_;
  bool doInt_;
  bool isFirstStep;

  void   readGaussians(int iarg, IFile*);
  bool   readChunkOfGaussians(int iarg, IFile *ifile, unsigned n);
  void   writeGaussian(int iarg, const Gaussian&, OFile*);
  void   addGaussian(int iarg, const Gaussian&);
  double getBiasAndDerivatives(int iarg, const vector<double>&,double* der=NULL);
  double evaluateGaussian(int iarg, const vector<double>&, const Gaussian&,double* der=NULL);
  vector<unsigned> getGaussianSupport(int iarg, const Gaussian&);
  bool   scanOneHill(int iarg, IFile *ifile,  vector<Value> &v, vector<double> &center, vector<double>  &sigma, double &height);
  std::string fmt;

public:
  PBMetaD(const ActionOptions&);
  ~PBMetaD();
  void calculate();
  void update();
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(PBMetaD,"PBMETAD")

void PBMetaD::registerKeywords(Keywords& keys){
  Bias::registerKeywords(keys);
  componentsAreNotOptional(keys);
  keys.addOutputComponent("bias","default","the instantaneous value of the bias potential");
  keys.use("ARG");
  keys.add("compulsory","SIGMA","the widths of the Gaussian hills");
  keys.add("compulsory","PACE","the frequency for hill addition");
  keys.add("compulsory","FILE","files in which the lists of added hills are stored");
  keys.add("optional","HEIGHT","the heights of the Gaussian hills. Compulsory unless TAU, TEMP and BIASFACTOR are given");
  keys.add("optional","FMT","specify format for HILLS files (useful for decrease the number of digits in regtests)");
  keys.add("optional","BIASFACTOR","use well tempered metadynamics and use this biasfactor.  Please note you must also specify temp");
  keys.add("optional","TEMP","the system temperature - this is only needed if you are doing well-tempered metadynamics");
  keys.add("optional","TAU","in well tempered metadynamics, sets height to (kb*DeltaT*pace*timestep)/tau");
  keys.add("optional","GRID_MIN","the lower bounds for the grid");
  keys.add("optional","GRID_MAX","the upper bounds for the grid");
  keys.add("optional","GRID_BIN","the number of bins for the grid");
  keys.add("optional","GRID_SPACING","the approximate grid spacing (to be used as an alternative or together with GRID_BIN)");
  keys.addFlag("GRID_SPARSE",false,"use a sparse grid to store hills");
  keys.addFlag("GRID_NOSPLINE",false,"don't use spline interpolation with grids");
  keys.add("optional","INTERVAL_MIN","monodimensional lower limits, outside the limits the system will not feel the biasing force.");
  keys.add("optional","INTERVAL_MAX","monodimensional upper limits, outside the limits the system will not feel the biasing force.");
  keys.addFlag("MULTIPLE_WALKERS",false,"Switch on MPI version of multiple walkers");
}

PBMetaD::~PBMetaD(){
  for(int i=0; i<BiasGrids_.size();   ++i) delete BiasGrids_[i];
  for(int i=0; i<hillsOfiles_.size(); ++i){
   hillsOfiles_[i]->close();
   delete hillsOfiles_[i];
  }
}

PBMetaD::PBMetaD(const ActionOptions& ao):
PLUMED_BIAS_INIT(ao),
height0_(std::numeric_limits<double>::max()), biasf_(1.0), kbt_(0.0),
grid_(false), stride_(0), welltemp_(false),
multiple_w(false), doInt_(false), isFirstStep(true)
{

  parse("FMT",fmt);

  // parse the sigma
  parseVector("SIGMA",sigma0_);
  if( sigma0_.size()!=getNumberOfArguments() ) error("number of arguments does not match number of SIGMA parameters");
 
  // note: HEIGHT is not compulsory, since one could use the TAU keyword, see below
  parse("HEIGHT",height0_);
  parse("PACE",stride_);
  if(stride_<=0) error("frequency for hill addition is nonsensical");
  
  vector<string> hillsfname;
  parseVector("FILE",hillsfname);
  if( hillsfname.size()!=getNumberOfArguments() ) error("number of FILE arguments does not match number of HILLS files");

  parse("BIASFACTOR",biasf_);
  if( biasf_<1.0 ) error("well tempered bias factor is nonsensical");
  double temp=0.0;
  parse("TEMP",temp);
  if(temp>0.0) kbt_=plumed.getAtoms().getKBoltzmann()*temp;
  else kbt_=plumed.getAtoms().getKbT();
  if(biasf_>1.0){
    if(kbt_==0.0) error("Unless the MD engine passes the temperature to plumed, with well-tempered metad you must specify it using TEMP");
    welltemp_=true;
  }
  double tau=0.0;
  parse("TAU",tau);
  if(tau==0.0){
    if(height0_==std::numeric_limits<double>::max()) error("At least one between HEIGHT and TAU should be specified");
    // if tau is not set, we compute it here from the other input parameters
    if(welltemp_) tau=(kbt_*(biasf_-1.0))/height0_*getTimeStep()*stride_;
  } else {
    if(!welltemp_)error("TAU only makes sense in well-tempered metadynamics");
    if(height0_!=std::numeric_limits<double>::max()) error("At most one between HEIGHT and TAU should be specified");
    height0_=(kbt_*(biasf_-1.0))/tau*getTimeStep()*stride_;
  }
  
  // MPI version
  parseFlag("MULTIPLE_WALKERS",multiple_w);
  
  // Grid Stuff
  vector<std::string> gmin(getNumberOfArguments());
  parseVector("GRID_MIN",gmin);
  if(gmin.size()!=getNumberOfArguments() && gmin.size()!=0) error("not enough values for GRID_MIN");
  vector<std::string> gmax(getNumberOfArguments());
  parseVector("GRID_MAX",gmax);
  if(gmax.size()!=getNumberOfArguments() && gmax.size()!=0) error("not enough values for GRID_MAX");
  vector<unsigned> gbin(getNumberOfArguments());
  vector<double>   gspacing;
  parseVector("GRID_BIN",gbin);
  if(gbin.size()!=getNumberOfArguments() && gbin.size()!=0) error("not enough values for GRID_BIN");
  parseVector("GRID_SPACING",gspacing);
  if(gspacing.size()!=getNumberOfArguments() && gspacing.size()!=0) error("not enough values for GRID_SPACING");
  if(gmin.size()!=gmax.size()) error("GRID_MAX and GRID_MIN should be either present or absent");
  if(gspacing.size()!=0 && gmin.size()==0) error("If GRID_SPACING is present also GRID_MIN should be present");
  if(gbin.size()!=0     && gmin.size()==0) error("If GRID_SPACING is present also GRID_MIN should be present");
  if(gmin.size()!=0){
    if(gbin.size()==0 && gspacing.size()==0){
        log<<"  Binsize not specified, 1/5 of sigma will be be used\n";
        plumed_assert(sigma0_.size()==getNumberOfArguments());
        gspacing.resize(getNumberOfArguments());
        for(unsigned i=0;i<gspacing.size();i++) gspacing[i]=0.2*sigma0_[i];
    } else if(gspacing.size()!=0 && gbin.size()==0){
      log<<"  The number of bins will be estimated from GRID_SPACING\n";
    } else if(gspacing.size()!=0 && gbin.size()!=0){
      log<<"  You specified both GRID_BIN and GRID_SPACING\n";
      log<<"  The more conservative (highest) number of bins will be used for each variable\n";
    }
    if(gbin.size()==0) gbin.assign(getNumberOfArguments(),1);
    if(gspacing.size()!=0) for(unsigned i=0;i<getNumberOfArguments();i++){
      double a,b;
      Tools::convert(gmin[i],a);
      Tools::convert(gmax[i],b);
      unsigned n=((b-a)/gspacing[i])+1;
      if(gbin[i]<n) gbin[i]=n;
    }
  }
  bool sparsegrid=false;
  parseFlag("GRID_SPARSE",sparsegrid);
  bool nospline=false;
  parseFlag("GRID_NOSPLINE",nospline);
  bool spline=!nospline;
  if(gbin.size()>0){grid_=true;}
  
  // Interval keyword
  parseVector("INTERVAL_MIN",lowI_);
  parseVector("INTERVAL_MAX",uppI_);
  // various checks
  if(lowI_.size()!=uppI_.size()) error("both a lower and an upper limits must be provided with INTERVAL");
  if(lowI_.size()!=0 && lowI_.size()!=getNumberOfArguments()) error("check number of argument of INTERVAL");
  for(unsigned i=0; i<lowI_.size(); ++i) if(uppI_[i]<lowI_[i]) error("The Upper limit must be greater than the Lower limit!");
  if(lowI_.size()>0) doInt_=true;
  
  checkRead();

  log.printf("  Gaussian width ");
  for(unsigned i=0;i<sigma0_.size();++i) log.printf(" %f",sigma0_[i]);
  log.printf("  Gaussian height %f\n",height0_);
  log.printf("  Gaussian deposition pace %d\n",stride_); 

  log.printf("  Gaussian files ");
  for(unsigned i=0; i<hillsfname.size(); ++i) log.printf("%s ",hillsfname[i].c_str());
  log.printf("\n");
  if(welltemp_){
    log.printf("  Well-Tempered Bias Factor %f\n",biasf_);
    log.printf("  Hills relaxation time (tau) %f\n",tau);
    log.printf("  KbT %f\n",kbt_);
  }
  if(multiple_w) log.printf("  Multiple walkers active using MPI communnication\n");
  if(doInt_) log.printf("  Upper and Lower limits boundaries for the bias are activated\n");
  if(grid_){
   log.printf("  Grid min");
   for(unsigned i=0;i<gmin.size();++i) log.printf(" %s",gmin[i].c_str() );
   log.printf("\n");
   log.printf("  Grid max");
   for(unsigned i=0;i<gmax.size();++i) log.printf(" %s",gmax[i].c_str() );
   log.printf("\n");
   log.printf("  Grid bin");
   for(unsigned i=0;i<gbin.size();++i) log.printf(" %u",gbin[i]);
   log.printf("\n");
   if(spline){log.printf("  Grid uses spline interpolation\n");}
   if(sparsegrid){log.printf("  Grid uses sparse grid\n");}
  }

  addComponentWithDerivatives("bias"); componentIsNotPeriodic("bias");

// initializing vector of hills
  hills_.resize(getNumberOfArguments());

// initializing and checking grid
  if(grid_){
   // check for mesh and sigma size
   for(unsigned i=0;i<getNumberOfArguments();i++) {
     double a,b;
     Tools::convert(gmin[i],a);
     Tools::convert(gmax[i],b);
     double mesh=(b-a)/((double)gbin[i]);
     if(mesh>0.5*sigma0_[i]) log<<"  WARNING: Using a METAD with a Grid Spacing larger than half of the Gaussians width can produce artifacts\n";
   }
   std::string funcl=getLabel() + ".bias";
   for(unsigned i=0; i<getNumberOfArguments(); ++i){
    std::vector<Value*> args(1);
    args[0] = getPntrToArgument(i);
    vector<std::string> gmin_t(1);
    vector<std::string> gmax_t(1);
    vector<unsigned>    gbin_t(1);
    gmin_t[0] = gmin[i];
    gmax_t[0] = gmax[i];
    gbin_t[0] = gbin[i];
    Grid* BiasGrid_;
    if(!sparsegrid){BiasGrid_=new Grid(funcl,args,gmin_t,gmax_t,gbin_t,spline,true);}
    else           {BiasGrid_=new SparseGrid(funcl,args,gmin_t,gmax_t,gbin_t,spline,true);}
    std::vector<std::string> actualmin=BiasGrid_->getMin();
    std::vector<std::string> actualmax=BiasGrid_->getMax();
    if(gmin_t[0]!=actualmin[0]) log<<"  WARNING: GRID_MIN["<<i<<"] has been adjusted to "<<actualmin[0]<<" to fit periodicity\n";
    if(gmax_t[0]!=actualmax[0]) log<<"  WARNING: GRID_MAX["<<i<<"] has been adjusted to "<<actualmax[0]<<" to fit periodicity\n";
    BiasGrids_.push_back(BiasGrid_);
   }
  }

// read Gaussians if restarting
  for(int i=0;i<hillsfname.size();++i){
   IFile *ifile = new IFile();
   ifile->link(*this);
   if(ifile->FileExist(hillsfname[i])){
    ifile->open(hillsfname[i]);
    if(getRestart()){
     log.printf("  Restarting from %s:",hillsfname[i].c_str());                  
     readGaussians(i, ifile);                                                    
    }
    ifile->reset(false);
    ifile->close();
    delete ifile;
   }
  }

// open hills files for writing
 for(int i=0;i<hillsfname.size();++i){
  OFile *ofile = new OFile();
  ofile->link(*this);
  string hillsfname_tmp = hillsfname[i];
  // if multiple walkers, only rank 0 will write to file
  if(multiple_w){
    int r=0;
    if(comm.Get_rank()==0) r=multi_sim_comm.Get_rank();
    comm.Bcast(r,0);
    if(r>0) hillsfname_tmp="/dev/null";
    ofile->enforceSuffix("");
  }
  ofile->open(hillsfname_tmp);
  if(fmt.length()>0) ofile->fmtField(fmt);
  ofile->addConstantField("multivariate");
  if(doInt_) {
    ofile->addConstantField("lower_int").printField("lower_int",lowI_[i]);
    ofile->addConstantField("upper_int").printField("upper_int",uppI_[i]);
  }
  ofile->setHeavyFlush();
  // output periodicities of variables
  ofile->setupPrintValue( getPntrToArgument(i) );
  // push back
  hillsOfiles_.push_back(ofile);
 }

  log<<"  Bibliography "<<plumed.cite("Pfaendtner and Bonomi. J. Chem. Theory Comput. 11, 5062 (2015)");
  if(doInt_) log<<plumed.cite(
     "Baftizadeh, Cossio, Pietrucci, and Laio, Curr. Phys. Chem. 2, 79 (2012)");
  if(multiple_w) log<<plumed.cite(
    "Raiteri, Laio, Gervasio, Micheletti, and Parrinello, J. Phys. Chem. B 110, 3533 (2006)");   
 
  log<<"\n";

  turnOnDerivatives();
}

void PBMetaD::readGaussians(int iarg, IFile *ifile){
 vector<double> center(1);
 vector<double> sigma(1);
 double height;
 int nhills=0; 

 std::vector<Value> tmpvalues;
 tmpvalues.push_back( Value( this, getPntrToArgument(iarg)->getName(), false ) ); 

 while(scanOneHill(iarg,ifile,tmpvalues,center,sigma,height)){;
  nhills++;
  if(welltemp_){height*=(biasf_-1.0)/biasf_;}
  addGaussian(iarg, Gaussian(center,sigma,height));
 }     
 log.printf("      %d Gaussians read\n",nhills);
}

bool PBMetaD::readChunkOfGaussians(int iarg, IFile *ifile, unsigned n){
 vector<double> center(1);
 vector<double> sigma(1);
 double height;
 unsigned nhills=0;
 std::vector<Value> tmpvalues;
 tmpvalues.push_back( Value( this, getPntrToArgument(iarg)->getName(), false ) ); 

 while(scanOneHill(iarg,ifile,tmpvalues,center,sigma,height)){;
  if(welltemp_){height*=(biasf_-1.0)/biasf_;}
  addGaussian(iarg, Gaussian(center,sigma,height));
  if(nhills==n){
      log.printf("      %u Gaussians read\n",nhills);
      return true;
  }
  nhills++;
 }     
 log.printf("      %u Gaussians read\n",nhills);
 return false;
}

void PBMetaD::writeGaussian(int iarg, const Gaussian& hill, OFile *ofile){
  ofile->printField("time",getTimeStep()*getStep());
  ofile->printField(getPntrToArgument(iarg),hill.center[0]);
  ofile->printField("multivariate","false");
  ofile->printField("sigma_"+getPntrToArgument(iarg)->getName(),hill.sigma[0]);
  double height=hill.height;
  if(welltemp_){height *= biasf_/(biasf_-1.0);}
  ofile->printField("height",height);
  ofile->printField("biasf",biasf_);
  ofile->printField();
}

void PBMetaD::addGaussian(int iarg, const Gaussian& hill){
 if(!grid_){hills_[iarg].push_back(hill);} 
 else{
  vector<unsigned> nneighb=getGaussianSupport(iarg, hill);
  vector<Grid::index_t> neighbors=BiasGrids_[iarg]->getNeighbors(hill.center,nneighb);
  vector<double> der(1);
  vector<double> xx(1);
  if(comm.Get_size()==1){
    for(unsigned i=0;i<neighbors.size();++i){
     Grid::index_t ineigh=neighbors[i];
     der[0]=0.0;
     BiasGrids_[iarg]->getPoint(ineigh,xx);
     double bias=evaluateGaussian(iarg,xx,hill,&der[0]);
     BiasGrids_[iarg]->addValueAndDerivatives(ineigh,bias,der);
    } 
  } else {
    unsigned stride=comm.Get_size();
    unsigned rank=comm.Get_rank();
    vector<double> allder(neighbors.size(),0.0);
    vector<double> allbias(neighbors.size(),0.0);
    for(unsigned i=rank;i<neighbors.size();i+=stride){
     Grid::index_t ineigh=neighbors[i];
     BiasGrids_[iarg]->getPoint(ineigh,xx);
     allbias[i]=evaluateGaussian(iarg,xx,hill,&allder[i]);
    }
    comm.Sum(allbias);
    comm.Sum(allder);
    for(unsigned i=0;i<neighbors.size();++i){
     Grid::index_t ineigh=neighbors[i];
     der[0]=allder[i];
     BiasGrids_[iarg]->addValueAndDerivatives(ineigh,allbias[i],der);
    }
  }
 }
}

vector<unsigned> PBMetaD::getGaussianSupport(int iarg, const Gaussian& hill){
 vector<unsigned> nneigh;
 
 if(doInt_){
   double cutoff=sqrt(2.0*DP2CUTOFF)*hill.sigma[0];
   if(hill.center[0]+cutoff > uppI_[iarg] || hill.center[0]-cutoff < lowI_[iarg]) { 
     // in this case, we updated the entire grid to avoid problems
     return BiasGrids_[iarg]->getNbin();
   } else {
     nneigh.push_back( static_cast<unsigned>(ceil(cutoff/BiasGrids_[iarg]->getDx()[0])));
     return nneigh;
   }
 }

 double cutoff=sqrt(2.0*DP2CUTOFF)*hill.sigma[0];
 nneigh.push_back( static_cast<unsigned>(ceil(cutoff/BiasGrids_[iarg]->getDx()[0])) );

 return nneigh;
}

double PBMetaD::getBiasAndDerivatives(int iarg, const vector<double>& cv, double* der)
{
 double bias=0.0;
 if(!grid_){
  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  for(unsigned i=rank;i<hills_[iarg].size();i+=stride){
   bias += evaluateGaussian(iarg,cv,hills_[iarg][i],der);
  }
  comm.Sum(bias);
  if(der) comm.Sum(der,1);
 }else{
  if(der){
   vector<double> vder(1);
   bias = BiasGrids_[iarg]->getValueAndDerivatives(cv,vder);
   der[0] = vder[0];
  }else{
   bias = BiasGrids_[iarg]->getValue(cv);
  }
 }
 return bias;
}

double PBMetaD::evaluateGaussian
 (int iarg, const vector<double>& cv, const Gaussian& hill, double* der)
{
 double bias=0.0;
// I use a pointer here because cv is const (and should be const)
// but when using doInt it is easier to locally replace cv[0] with
// the upper/lower limit in case it is out of range
 const double *pcv=NULL; // pointer to cv
 double tmpcv[1]; // tmp array with cv (to be used with doInt_)
 if(cv.size()>0) pcv=&cv[0];
 if(doInt_){
   plumed_dbg_assert(cv.size()==1);
   pcv=&(tmpcv[0]);
   tmpcv[0]=cv[0];
   if(cv[0]<lowI_[iarg]) tmpcv[0]=lowI_[iarg];
   if(cv[0]>uppI_[iarg]) tmpcv[0]=uppI_[iarg];
 }
 double dp = difference(0,hill.center[0],pcv[0]) / hill.sigma[0];
 double dp2 = 0.5 * dp * dp;
 if(dp2<DP2CUTOFF){
       bias = hill.height*exp(-dp2);
       if(der){der[0]+= -bias * dp / hill.sigma[0];}
 }
 if(doInt_){
   if((cv[0]<lowI_[iarg] || cv[0]>uppI_[iarg]) && der ) der[0] = 0.0;
 }
 return bias;
}

void PBMetaD::calculate()
{
  vector<double> cv(1);
  double* der=new double[1];
  vector<double> bias(getNumberOfArguments());
  vector<double> deriv(getNumberOfArguments());

  double ene = 0.;
  double ncv = (double) getNumberOfArguments();
  for(unsigned i=0; i<getNumberOfArguments(); ++i){
   cv[0] = getArgument(i);
   der[0] = 0.0;
   bias[i] = getBiasAndDerivatives(i, cv, der);
   deriv[i] = der[0];
   ene += exp(-bias[i]/kbt_);
  }
      
  // set Forces 
  for(unsigned i=0; i<getNumberOfArguments(); ++i){
    const double f = - exp(-bias[i]/kbt_) / (ene) * deriv[i];
    setOutputForce(i, f);
    getPntrToComponent("bias")->addDerivative(i,-f);
  }
  delete [] der;

  // set bias
  ene = -kbt_ * (std::log(ene) - std::log(ncv));
  getPntrToComponent("bias")->set(ene);
}

void PBMetaD::update()
{
  // adding hills criteria
  bool nowAddAHill; 
  if(getStep()%stride_==0 && !isFirstStep ){nowAddAHill=true;}else{nowAddAHill=false;isFirstStep=false;}

  if(nowAddAHill){

   // get all CVs value
   vector<double> cv(getNumberOfArguments());
   for(unsigned i=0; i<getNumberOfArguments(); ++i) cv[i] = getArgument(i);

   // get all biases and heights
   vector<double> bias(getNumberOfArguments());
   vector<double> height(getNumberOfArguments());
   vector<double> cv_tmp(1);
   vector<double> sigma_tmp(1);
   double norm = 0.0;
   for(unsigned i=0; i<getNumberOfArguments(); ++i){
     cv_tmp[0] = cv[i];
     bias[i] = getBiasAndDerivatives(i, cv_tmp);
     double h = exp(-bias[i]/kbt_);
     norm += h;
     height[i] = h;
   }
   // normalize and apply welltemp correction
   double ncv = static_cast<double> (getNumberOfArguments());
   for(unsigned i=0; i<getNumberOfArguments(); ++i){
     height[i] *=  ncv * height0_ / norm;
     if(welltemp_) height[i] *= exp(-bias[i]/(kbt_*(biasf_-1.0)));
   }
   
   // Multiple walkers: share hills and add them all
   if(multiple_w){
     int nw = 0;
     int mw = 0;  
     if(comm.Get_rank()==0){
     // Only root of group can communicate with other walkers
       nw = multi_sim_comm.Get_size();
       mw = multi_sim_comm.Get_rank();
     }
     // Communicate to the other members of the same group
     // info about number of walkers and walker index
     comm.Bcast(nw,0);
     comm.Bcast(mw,0);
     // Allocate arrays to store all walkers hills
     std::vector<double> all_cv(nw*cv.size(), 0.0);
     std::vector<double> all_height(nw*height.size(), 0.0);
     if(comm.Get_rank()==0){
     // Communicate (only root)
       multi_sim_comm.Allgather(cv, all_cv);
       multi_sim_comm.Allgather(height, all_height);
     }
     // Share info with group members
     comm.Bcast(all_cv,0);
     comm.Bcast(all_height,0);
     // now add hills one by one
     for(int j=0; j<nw; ++j){
      for(unsigned i=0; i<getNumberOfArguments(); ++i){
       cv_tmp[0]    = all_cv[j*cv.size()+i];
       sigma_tmp[0] = sigma0_[i];
       // new Gaussian
       Gaussian newhill = Gaussian(cv_tmp, sigma_tmp, all_height[j*cv.size()+i]);
       addGaussian(i, newhill);
       // print on HILLS file
       writeGaussian(i, newhill, hillsOfiles_[i]);
      }
     }  
   // just add your own hills  
   }else{
    for(unsigned i=0; i<getNumberOfArguments(); ++i){
      cv_tmp[0]    = cv[i];
      sigma_tmp[0] = sigma0_[i];
      // new Gaussian
      Gaussian newhill = Gaussian(cv_tmp, sigma_tmp, height[i]);
      addGaussian(i, newhill);
      // print on HILLS file
      writeGaussian(i, newhill, hillsOfiles_[i]);
    } 
   }
 }
}

/// takes a pointer to the file and a template string with values v and gives back the next center, sigma and height 
bool PBMetaD::scanOneHill(int iarg, IFile *ifile,  vector<Value> &tmpvalues, vector<double> &center, vector<double>  &sigma, double &height){
  double dummy;
 
   if(ifile->scanField("time",dummy)){
     ifile->scanField( &tmpvalues[0] );
     if( tmpvalues[0].isPeriodic() && ! getPntrToArgument(iarg)->isPeriodic() ){
        error("in hills file periodicity for variable " + tmpvalues[0].getName() + " does not match periodicity in input");
     } else if( tmpvalues[0].isPeriodic() ){
        std::string imin, imax; tmpvalues[0].getDomain( imin, imax );
        std::string rmin, rmax; getPntrToArgument(iarg)->getDomain( rmin, rmax );
        if( imin!=rmin || imax!=rmax ){
          error("in hills file periodicity for variable " + tmpvalues[0].getName() + " does not match periodicity in input");
        }
     }
     center[0]=tmpvalues[0].get();
     // scan for multivariate label: record the actual file position so to eventually rewind 
     std::string sss;
     ifile->scanField("multivariate",sss);
     ifile->scanField("sigma_"+getPntrToArgument(iarg)->getName(),sigma[0]); 
     ifile->scanField("height",height);
     ifile->scanField("biasf",dummy);
     if(ifile->FieldExist("clock")) ifile->scanField("clock",dummy);
     if(ifile->FieldExist("lower_int")) ifile->scanField("lower_int",dummy);
     if(ifile->FieldExist("upper_int")) ifile->scanField("upper_int",dummy);
     ifile->scanField();
     return true;
  }else{ 
    return false; 
  }; 
}

}
}

