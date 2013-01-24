/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#include "tools/BiasRepresentation.h"
#include "tools/File.h"
#include "tools/Tools.h"
#include <iostream>

using namespace std;

namespace PLMD{
namespace function{


//+PLUMEDOC FUNCTION FUNCSUMHILLS 
/*

*/
//+ENDPLUMEDOC

class FilesHandler{
	vector <string> filenames; 
        vector <IFile*>  ifiles;
        Action *action;
        Log *log;
        bool parallelread;
        unsigned beingread;  
        bool isopen;
	public:
		FilesHandler(const vector<string> &filenames, const bool &parallelread ,  Action &myaction , Log &mylog);
		bool readBunch(BiasRepresentation *br, unsigned stride);
		bool scanOneHill(BiasRepresentation *br, IFile *ifile );
		void getMinMaxBin(vector<Value*> vals, Communicator &cc, vector<double> &vmin, vector<double> &vmax, vector<unsigned> &vbin);
		void getMinMaxBin(vector<Value*> vals, Communicator &cc, vector<double> &vmin, vector<double> &vmax, vector<unsigned> &vbin, vector<double> &histosigma);
}; 
FilesHandler::FilesHandler(const vector<string> &filenames, const bool &parallelread , Action &action , Log &mylog ):filenames(filenames),log(&mylog),parallelread(parallelread),beingread(0),isopen(false){
   this->action=&action;
   for(unsigned i=0;i<filenames.size();i++){
      IFile *ifile = new IFile();
      ifile->link(action);
      ifiles.push_back(ifile);
      plumed_massert((ifile->FileExist(filenames[i])), "the file "+filenames[i]+" does not exist " );
   }
   
};

// note that the FileHandler is completely transparent respect to the biasrepresentation 
// no check are made at this level
bool FilesHandler::readBunch(BiasRepresentation *br , unsigned stride = -1){
        bool morefiles; morefiles=true;
	if(parallelread){
		(*log)<<"  doing parallelread \n";
        }else{
		(*log)<<"  doing serialread \n";
        	// read one by one hills      
		// is the type defined? if not, assume it is a gaussian 
                IFile *ff; 
                ff=ifiles[beingread];
                if(!isopen){
                	(*log)<<"  opening file "<<filenames[beingread]<<"\n";
			ff->open(filenames[beingread]);isopen=true;
                }
		int n;
                while(true){
			bool fileisover=true;
			while(scanOneHill(br,ff)){
				// here do the dump if needed 
				n=br->getNumberOfKernels();
				if(stride>0 && n%stride==0 && n!=0  ){
                	        	(*log)<<"  done with this chunk: now with "<<n<<" kernels  \n";
					fileisover=false;
					break;	
                                }
                        }   
                        if(fileisover){
                	        (*log)<<"  closing file "<<filenames[beingread]<<"\n";
                	        ff->close();
				isopen=false;
                	        (*log)<<"  now total "<<br->getNumberOfKernels()<<" kernels \n"; 
				beingread++;
 			        if(beingread<ifiles.size()){
					ff=ifiles[beingread];ff->open(filenames[beingread]);
                			(*log)<<"  opening file "<<filenames[beingread]<<"\n";
					isopen=true;
				}else{morefiles=false; 
                                      (*log)<<"  final chunk: now with "<<n<<" kernels  \n";
				      break;
				}  
                        }
			// if there are no more files to read and this file is over then quit 	
                        if(fileisover && !morefiles){break;} 
                        // if you are in the middle of a file and you are here
			// then means that you read what you need to read
                        if(!fileisover ){break;} 
                } 
        }        
	return morefiles;
};
void FilesHandler::getMinMaxBin(vector<Value*> vals, Communicator &cc, vector<double> &vmin, vector<double> &vmax, vector<unsigned> &vbin){
        // create the representation (no grid)
     	BiasRepresentation br(vals,cc);
        // read all the kernels
        readBunch(&br);
        // loop over the kernels and get the support 
 	br.getMinMaxBin(vmin,vmax,vbin);
};
void FilesHandler::getMinMaxBin(vector<Value*> vals, Communicator &cc, vector<double> &vmin, vector<double> &vmax, vector<unsigned> &vbin, vector<double> &histosigma){
     	BiasRepresentation br(vals,cc,histosigma);
        // read all the kernels
        readBunch(&br);
        // loop over the kernels and get the support 
 	br.getMinMaxBin(vmin,vmax,vbin);
        //for(unsigned i=0;i<vals.size();i++){cerr<<"XXX "<<vmin[i]<<" "<<vmax[i]<<" "<<vbin[i]<<"\n";}	
};
bool FilesHandler::scanOneHill(BiasRepresentation *br, IFile *ifile ){
	double dummy;
	if(ifile->scanField("time",dummy)){
	        //(*log)<<"   scanning one hill: "<<dummy<<" \n";
	        ifile->scanField("biasf",dummy);
	        if(ifile->FieldExist("clock")) ifile->scanField("clock",dummy);
                // keep this intermediate function in case you need to parse more data in the future
                br->pushKernel(ifile);
	 	//(*log)<<"  read hill\n";	
		if(br->hasSigmaInInput())ifile->allowIgnoredFields();
	        ifile->scanField();  
 		return true;
	}else{
		return false;
	}
};


double  mylog( double v1 ){
      return log(v1);
};

double  mylogder( double v1 ){
      return 1./v1;
};



class FuncSumHills :
  public Function
{
  vector<string> hillsFiles,histoFiles; 
  vector<string> proj; 
  unsigned highdim, lowdim;
  int initstride,stride;
  bool iscltool,integratehills,integratehisto,parallelread;
  bool negativebias;
  double beta;
  string outhills,outhisto;
  BiasRepresentation *biasrep;
  BiasRepresentation *historep;
public:
  FuncSumHills(const ActionOptions&);
  ~FuncSumHills();
  void calculate(); // this probably is not needed
  bool checkFilesAreExisting(const vector<string> & hills ); 
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(FuncSumHills,"FUNCSUMHILLS")

void FuncSumHills::registerKeywords(Keywords& keys){
  Function::registerKeywords(keys);
  keys.use("ARG"); 
  keys.add("optional","HILLSFILES"," source file for hills creation(may be the same as HILLS)"); // this can be a vector! 
  keys.add("optional","HISTOFILES"," source file for histogram creation(may be the same as HILLS)"); // also this can be a vector!
  keys.add("optional","HISTOSIGMA"," sigmas for binning when the histogram correction is needed    "); 
  keys.add("optional","PROJ"," only with sumhills: the projection on the cvs");
  keys.add("optional","KT"," only with sumhills: the kt factor when projection on cvs");
  keys.add("optional","GRID_MIN","the lower bounds for the grid");
  keys.add("optional","GRID_MAX","the upper bounds for the grid");
  keys.add("optional","GRID_BIN","the number of bins for the grid"); 
  keys.add("optional","OUTHILLS"," output file for hills ");
  keys.add("optional","OUTHISTO"," output file for histogram ");
  keys.add("optional","INITSTRIDE"," stride if you want an initial dump ");
  keys.add("optional","STRIDE"," stride when you do it on the fly ");
  keys.addFlag("ISCLTOOL",true,"use via plumed commandline: calculate at read phase and then go");
  keys.addFlag("PARALLELREAD",false,"read parallel HILLS file");
  keys.addFlag("NEGBIAS",false,"dump  negative bias ( -bias )   instead of the free energy: needed in welltempered with flexible hills ");
}

FuncSumHills::FuncSumHills(const ActionOptions&ao):
Action(ao),
Function(ao),
initstride(-1),
stride(-1),
iscltool(false),
integratehills(false),
integratehisto(false),
parallelread(false),
negativebias(false),
beta(-1.)
{
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
  parseVector("GRID_BIN",gbin);
  plumed_massert(gbin.size()==getNumberOfArguments() || gbin.size()==0,"need GRID_BIN argument for this"); 
  if(gbin.size()!=getNumberOfArguments() && gbin.size()!=0) error("not enough values for GRID_BIN");
  //plumed_assert(getNumberOfArguments()==gbin.size());

  // hills file: 
  parseVector("HILLSFILES",hillsFiles);
  if(hillsFiles.size()==0){
  	integratehills=false; // default behaviour 
  }else{
  	integratehills=true; 
        for(unsigned i=0;i<hillsFiles.size();i++) log<<"  hillsfile  : "<<hillsFiles[i]<<"\n";
  }
  // histo file: 
  parseVector("HISTOFILES",histoFiles);
  if(histoFiles.size()==0){
  	integratehisto=false;  
  }else{
  	integratehisto=true;  
        for(unsigned i=0;i<histoFiles.size();i++) log<<"  histofile  : "<<histoFiles[i]<<"\n";
  }
  vector<double> histoSigma;
  if(integratehisto){
  	parseVector("HISTOSIGMA",histoSigma);
	plumed_massert(histoSigma.size()==getNumberOfArguments()," The number of sigmas must be the same of the number of arguments ");
        for(unsigned i=0;i<histoSigma.size();i++) log<<"  histosigma  : "<<histoSigma[i]<<"\n";
  }

  // add some automatic hills width: not in case stride is defined  
  // since when you start from zero the automatic size will be zero!
  if(gmin.size()==0 || gmax.size()==0){
	log<<"   \n"; 
	log<<"  No boundaries defined: need to do a prescreening of hills \n"; 
        std::vector<Value*> tmpvalues; 
        for(unsigned i=0;i<getNumberOfArguments();i++)tmpvalues.push_back( getPntrToArgument(i) );
        if(integratehills) {
        	FilesHandler *hillsHandler;
        	hillsHandler=new FilesHandler(hillsFiles,parallelread,*this, log);
		vector<double> vmin,vmax;
        	vector<unsigned> vbin;  
        	hillsHandler->getMinMaxBin(tmpvalues,comm,vmin,vmax,vbin);
		log<<"  found boundaries from hillsfile: \n";
		gmin.resize(vmin.size());
		gmax.resize(vmax.size());
                if(gbin.size()==0){
			gbin=vbin;
                }else{
			log<<"  found nbins in input, this overrides the automatic choice \n"; 
		}
		for(unsigned i=0;i<getNumberOfArguments();i++){
		 	Tools::convert(vmin[i],gmin[i]);
		 	Tools::convert(vmax[i],gmax[i]);
			log<<"  variable "<< getPntrToArgument(i)->getName()<<" min: "<<gmin[i]<<" max: "<<gmax[i]<<" nbin: "<<gbin[i]<<"\n";
		}
        } 
	// if at this stage bins are not there then do it with histo
	if(gmin.size()==0){
    	   	FilesHandler *histoHandler;
	        histoHandler=new FilesHandler(histoFiles,parallelread,*this, log);
		vector<double> vmin,vmax;
        	vector<unsigned> vbin;  
        	histoHandler->getMinMaxBin(tmpvalues,comm,vmin,vmax,vbin,histoSigma);
		log<<"  found boundaries from histofile: \n";
		gmin.resize(vmin.size());
		gmax.resize(vmax.size());
                if(gbin.size()==0){
			gbin=vbin;
                }else{
			log<<"  found nbins in input, this overrides the automatic choice \n"; 
		}
		for(unsigned i=0;i<getNumberOfArguments();i++){
		 	Tools::convert(vmin[i],gmin[i]);
		 	Tools::convert(vmax[i],gmax[i]);
			log<<"  variable "<< getPntrToArgument(i)->getName()<<" min: "<<gmin[i]<<" max: "<<gmax[i]<<" nbin: "<<gbin[i]<<"\n";
		}
        }
	log<<"  done!\n"; 
	log<<"   \n"; 
  }

  // needs a projection? 
  proj.clear();
  parseVector("PROJ",proj);
  plumed_massert(proj.size()<getNumberOfArguments()," The number of projection must be less than the full list of arguments ");

  if( proj.size() != 0 || integratehisto==true  ) {
    parse("KT",beta);
    for(unsigned i=0;i<proj.size();i++) log<<"  projection "<<i<<" : "<<proj[i]<<"\n";
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
  if(initstride<0){ log<<"  Doing only one integration: no stride \n";
  	 outhills="fes.dat";outhisto="correction.dat";}
  else{outhills="fes_";outhisto="correction_";
	log<<"  Doing integration slices every "<<initstride<<" kernels\n";
  }

  //what might it be this? 
  // here start 
  // want something right now?? do it and return
  // your argument is a set of cvs 
  // then you need: a hills / a colvar-like file (to do a histogram) 
  // create a bias representation for this
  if(iscltool){

     std::vector<Value*> tmpvalues; 
    for(unsigned i=0;i<getNumberOfArguments();i++){
        // allocate a new value from the old one: no deriv here
	tmpvalues.push_back( getPntrToArgument(i) );
    }

    // check if the files exists 
    if(integratehills){
         checkFilesAreExisting(hillsFiles); 
         biasrep=new BiasRepresentation(tmpvalues,comm, gmin, gmax, gbin);
	 if(negativebias){
		biasrep->setRescaledToBias(true);
	        log<<"  required the -bias instead of the free energy \n";
		if(initstride<0){outhills="negativebias.dat";}
		else{outhills="negativebias_";}
	 }
    }

    parse("OUTHILLS",outhills);
    parse("OUTHISTO",outhisto);
    if(integratehills)log<<"  output file for fes/bias   is :  "<<outhills<<"\n";   
    if(integratehisto)log<<"  output file for correction is :  "<<outhisto<<"\n";   
    checkRead();

    log<<"\n";
    log<<"  Now calculating...\n";
    log<<"\n";

    if(integratehisto){
         checkFilesAreExisting(histoFiles); 
         historep=new BiasRepresentation(tmpvalues,comm,gmin,gmax,gbin,histoSigma);
    }

    // decide how to source hills ( serial/parallel )
    // here below the input control 
    // say how many hills and it will read them from the 
    // bunch of files provided, will update the representation 
    // of hills (i.e. a list of hills and the associated grid)

    // decide how to source colvars ( serial parallel )
    FilesHandler *hillsHandler;
    FilesHandler *histoHandler;

    if(integratehills)	hillsHandler=new FilesHandler(hillsFiles,parallelread,*this, log);
    if(integratehisto)	histoHandler=new FilesHandler(histoFiles,parallelread,*this, log);

    // read a number of hills and put in the bias representation
    int nfiles=0;
    bool ibias=integratehills; bool ihisto=integratehisto;
    while(true){
        if(  integratehills  && ibias  ){ log<<"  reading hills: \n"; ibias=hillsHandler->readBunch(biasrep,initstride) ; log<<"\n"; }   
        if(  integratehisto  && ihisto ){ log<<"  reading histogram: \n"; ihisto=histoHandler->readBunch(historep,initstride) ;  log<<"\n";  }    
	// dump: need to project?	
        if(proj.size()!=0){

		if(integratehills){

    	      		log<<"  Projecting on subgrid... \n";
              		BiasWeight *Bw=new BiasWeight(beta); 
              		WeightBase *Wb=dynamic_cast<WeightBase*>(Bw); 
             		Grid biasGrid=*(biasrep->getGridPtr());
   	      		Grid smallGrid=biasGrid.project(proj,Wb);
              		OFile gridfile; gridfile.link(*this);
	      		std::ostringstream ostr;ostr<<nfiles;
              		string myout; 
                        if(initstride>0){ myout=outhills+ostr.str()+".dat" ;}else{myout=outhills;}
              		log<<"  Writing subgrid on file "<<myout<<" \n";
              		gridfile.open(myout);	
         
   	      		smallGrid.writeToFile(gridfile);
              		gridfile.close();
                        if(!ibias)integratehills=false;// once you get to the final bunch just give up 
		}
		if(integratehisto){

                        ProbWeight *Pw=new ProbWeight(beta);
                        WeightBase *Wb=dynamic_cast<WeightBase*>(Pw);
             		Grid histoGrid=*(historep->getGridPtr());
   	      		Grid smallGrid=histoGrid.project(proj,Wb);

              		OFile gridfile; gridfile.link(*this);
	      		std::ostringstream ostr;ostr<<nfiles;
              		string myout; 
                        if(initstride>0){ myout=outhisto+ostr.str()+".dat" ;}else{myout=outhisto;}
              		log<<"  Writing subgrid on file "<<myout<<" \n";
              		gridfile.open(myout);	
         
   	      		smallGrid.writeToFile(gridfile);
              		gridfile.close();

                        if(!ihisto)integratehisto=false;// once you get to the final bunch just give up 
                } 

	}else{

		if(integratehills){

	                Grid biasGrid=*(biasrep->getGridPtr());
			biasGrid.scaleAllValuesAndDerivatives(-1.);
	
	                OFile gridfile; gridfile.link(*this);
			std::ostringstream ostr;ostr<<nfiles;
			string myout;
			if(initstride>0){ myout=outhills+ostr.str()+".dat" ;}else{myout=outhills;}
	                log<<"  Writing full grid on file "<<myout<<" \n";
	                gridfile.open(myout);	
	
	                biasGrid.writeToFile(gridfile);
	                gridfile.close();
			// rescale back prior to accumulate
                        if(!ibias)integratehills=false;// once you get to the final bunch just give up 
		}
		if(integratehisto){

	                Grid histoGrid=*(historep->getGridPtr());
                        histoGrid.applyFunctionAllValuesAndDerivatives(&mylog,&mylogder);
                        histoGrid.scaleAllValuesAndDerivatives(-1./beta);	

			OFile gridfile; gridfile.link(*this);
			std::ostringstream ostr;ostr<<nfiles;
			string myout;
			if(initstride>0){ myout=outhisto+ostr.str()+".dat" ;}else{myout=outhisto;}
	                log<<"  Writing full grid on file "<<myout<<" \n";
	                gridfile.open(myout);	
	
	                histoGrid.writeToFile(gridfile);
	                gridfile.close();

                        if(!ihisto)integratehisto=false;// once you get to the final bunch just give up 
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

void FuncSumHills::calculate(){
  // this should be connected only with a grid representation to metadynamics 
  // at regular time just dump it
   plumed_merror("You should have never got here: this stuff is not yet implemented!"); 
}

FuncSumHills::~FuncSumHills(){
}

bool FuncSumHills::checkFilesAreExisting(const vector<string> & hills ){
	plumed_massert(hills.size()!=0,"the number of  files provided should be at least one" );
        IFile *ifile = new IFile();
        ifile->link(*this);
        for(unsigned i=0; i< hills.size();i++){  
          plumed_massert(ifile->FileExist(hills[i]),"missing file "+hills[i]);
        }
        return true; 

}

}

}


