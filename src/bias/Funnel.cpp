/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2014 The plumed team
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
#include "tools/Grid.h"
#include "tools/Exception.h"
#include "tools/File.h"
#include <cstring>
#include "tools/Communicator.h"
#include "core/ActionSet.h"
#include "tools/FileBase.h"
#include <memory>

using namespace std;


namespace PLMD{
namespace bias{

class Funnel : public Bias{

private:
  std::unique_ptr<GridBase> BiasGrid_;

  /////////////////////
  // old 2.3
  //Grid* BiasGrid_;
  /////////////////////
  //Optional parameters
  double NBINS;
  double NBINZ;
  double MINS;
  double KAPPA;
  double RCYL;
  double safety;
  double slope;
  double ALPHA;
  //Compulsory parameters
  double MAXS;
  double ZCC;
  double scale_;


public:
  explicit Funnel(const ActionOptions&);

  // old Funnel-2.3
  // ~Funnel();

  void calculate();
  static void registerKeywords(Keywords& keys);
  void createBIAS(const double& R_cyl, const double& z_cc, const double& alpha, const double& KAPPA,
		  const double& MIN_S, const double& MAX_S, const double& NBIN_S, const double& NBIN_Z,
		  const double& safety, const bool& sphere, const double& slope, const string& funcl,
		  const string& file);
//  void createBIAS3D(const double& R_cyl, const double& z_cc, const double& alpha,
//			const double& KAPPA, const double& MIN_S, const double& MAX_S, const double& NBIN_S,
//			const double& NBIN_Z, const double& safety, const bool& sphere, const double& slope,
//			const string& funcl);
};

PLUMED_REGISTER_ACTION(Funnel,"FUNNEL")

void Funnel::registerKeywords(Keywords& keys){
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.addFlag("NOSPLINE",false,"specifies that no spline interpolation is to be used when calculating the energy and forces due to the external potential");
  keys.addFlag("SPARSE",false,"specifies that the external potential uses a sparse grid");
  keys.addFlag("SPHERE",false, "Shape of your potential");
  keys.add("compulsory","SCALE","1.0","a factor that multiplies the external potential, useful to invert free energies");
// che sia roba vecchia?
  //  componentsAreNotOptional(keys);
  //  keys.addOutputComponent("bias","default","the instantaneous value of the bias potential");

  //Defining optional arguments
  keys.add("optional","NBINS","NBINS");
  keys.add("optional","NBINZ","NBINZ");
  keys.add("optional","MINS","MINS");
  keys.add("optional","KAPPA","KAPPA");
  keys.add("optional","RCYL","RCYL");
  keys.add("optional","SAFETY","SAFETY");
  keys.add("optional","SLOPE","SLOPE");
  keys.add("optional","ALPHA","ALPHA");
  //Defining compulsory arguments
  keys.add("compulsory","MAXS","MAXS");
  keys.add("compulsory","ZCC","ZCC");
  keys.add("compulsory","FILE","Name of the Funnel potential file for unlimited funnel works -cit");
  keys.addFlag("WALKERS_MPI",false,"There is a problem with gromacs walkers");
}

// Old version 2.3
//Funnel::~Funnel(){
//  delete BiasGrid_;
//}

Funnel::Funnel(const ActionOptions& ao):
PLUMED_BIAS_INIT(ao),
// Old version 2.3
// BiasGrid_(NULL),
NBINS(500.0),
NBINZ(500.0),
MINS(0.0),
KAPPA(84.0),
RCYL(0.1),
safety(1.0),
slope(1.0),
ALPHA(1.413)

{
  bool sparsegrid=false;
  parseFlag("SPARSE",sparsegrid);
  bool nospline=false;
  parseFlag("NOSPLINE",nospline);
  bool spline=!nospline;
  bool walkers_mpi=false;
  parseFlag("WALKERS_MPI",walkers_mpi);
//  bool components=false;
//  parseFlag("POINTS",components);
  bool sphere=false;
  parseFlag("SPHERE",sphere);
  parse("SAFETY",safety);
  string file;
  parse("FILE",file);
  if( file.length()==0 ) error("No funnel file name was specified");
  parse("SCALE",scale_);

  //Reading optional arguments
  parse("KAPPA",KAPPA);
  parse("NBINS",NBINS);
  parse("NBINZ",NBINZ);
  parse("MINS",MINS);
  parse("RCYL",RCYL);
  parse("SLOPE",slope);
  parse("ALPHA",ALPHA);
  //Reading compulsory arguments
  parse("MAXS",MAXS);
  parse("ZCC",ZCC);


  checkRead();

  log.printf("  External potential from file %s\n",file.c_str());
  log.printf("  Multiplied by %lf\n",scale_);
  if(spline){log.printf("  External potential uses spline interpolation\n");}
  if(sparsegrid){log.printf("  External potential uses sparse grid\n");}

  // Non più necessario dalla 2.3
//  addComponent("bias"); componentIsNotPeriodic("bias");

  std::string funcl=getLabel() + ".bias";

//  int size = plumed.comm.Get_size();
//  int rank = plumed.comm.Get_rank();
  IFile I_hate_this;
  bool do_exist=I_hate_this.FileExist(file);
  if(comm.Get_rank()==0){
	  if(multi_sim_comm.Get_rank()==0 && walkers_mpi){
		  if(!do_exist){
			  createBIAS(RCYL, ZCC, ALPHA, KAPPA, MINS, MAXS, NBINS, NBINZ, safety, sphere, slope, funcl, file);
		  }
	  } else {
		  if(!do_exist){
			  createBIAS(RCYL, ZCC, ALPHA, KAPPA, MINS, MAXS, NBINS, NBINZ, safety, sphere, slope, funcl, file);
		  }
	  }
	  if(walkers_mpi) multi_sim_comm.Barrier();
  }
  comm.Barrier();

// read grid
  IFile gridfile; gridfile.open(file);
  BiasGrid_=Grid::create(funcl,getArguments(),gridfile,sparsegrid,spline,true);
//not necessary anymore?  gridfile.close();
  if(BiasGrid_->getDimension()!=getNumberOfArguments()) error("mismatch between dimensionality of input grid and number of arguments");
  for(unsigned i=0;i<getNumberOfArguments();++i){
    if( getPntrToArgument(i)->isPeriodic()!=BiasGrid_->getIsPeriodic()[i] ) error("periodicity mismatch between arguments and input bias");
  }
  comm.Barrier();
  if(comm.Get_rank()==0 && walkers_mpi) multi_sim_comm.Barrier();
}


void Funnel::createBIAS(const double& R_cyl, const double& z_cc, const double& alpha,
		const double& KAPPA, const double& MIN_S, const double& MAX_S, const double& NBIN_S,
		const double& NBIN_Z, const double& safety, const bool& sphere, const double& slope,
		const string& funcl, const string& file){
	//R_cyl and z_cc forms the parameters of the cylinder.
	//alpha defines the angle in degrees.

	//PARAMETERS OF THE CONE
	double tg_alpha= tan(alpha);

	//parameters for PROGRESSION
	//parameters for DISTANCE
	double MIN_Z=0;
	double MAX_Z;
	if (sphere==false){
		MAX_Z=R_cyl + tg_alpha * (z_cc - MIN_S);
	}
	else {
		MAX_Z=z_cc+safety;
	}

	//bin size
	double DX_Z = (MAX_Z - MIN_Z) / NBIN_Z;
	double DX_S;
	if (sphere==false){
			DX_S=(MAX_S - MIN_S) / NBIN_S;
		}
		else {
			DX_S=(MAX_S + z_cc + safety)/NBIN_S;
		}

	double SS, Zmax, ZZ, D, d;
	double POT, FZ, FS;

	PLMD::OFile pof;
	pof.open(file);

	//Write the header
	pof.printf("#! FIELDS %s %s %s der_%s der_%s \n", getPntrToArgument(0)->getName().c_str(), getPntrToArgument(1)->getName().c_str(), funcl.c_str(), getPntrToArgument(0)->getName().c_str(), getPntrToArgument(1)->getName().c_str());
	if (sphere==false) pof.printf("#! SET min_%s %f\n", getPntrToArgument(0)->getName().c_str(), MIN_S);
	else pof.printf("#! SET min_%s %f\n", getPntrToArgument(0)->getName().c_str(), -z_cc-safety);
	pof.printf("#! SET max_%s %f\n", getPntrToArgument(0)->getName().c_str(), MAX_S);
	pof.printf("#! SET nbins_%s %f\n", getPntrToArgument(0)->getName().c_str(), NBIN_S);
	pof.printf("#! SET periodic_%s false\n", getPntrToArgument(0)->getName().c_str());
	pof.printf("#! SET min_%s %f\n", getPntrToArgument(1)->getName().c_str(), MIN_Z);
	pof.printf("#! SET max_%s %f\n", getPntrToArgument(1)->getName().c_str(), MAX_Z);
	pof.printf("#! SET nbins_%s %f\n", getPntrToArgument(1)->getName().c_str(), NBIN_Z);
	pof.printf("#! SET periodic_%s false\n", getPntrToArgument(1)->getName().c_str());

	//Calculate and write the GRID
	//Cone or cylinder?

	for(int is=0; is <= NBIN_S; is++){
		if (sphere==false){
			SS = MIN_S + is * DX_S;
		}
		else{
			SS = - z_cc - safety + is * DX_S;
		}
		bool cone = false;
		if (sphere==false){
			if(SS <= z_cc) cone = true;
		}
		else {
			if (SS <= sqrt(pow(z_cc,2)-pow(R_cyl,2))) cone = true;
		}
		//Set wall boundaries properly
		if(cone == true) {
			if(sphere==false){
				Zmax = R_cyl + (z_cc - SS) * tg_alpha;
			}
			else{
				if (SS > -z_cc){
					Zmax = sqrt(pow(z_cc,2) - pow(SS,2));
				}
				else{
					Zmax = 0;
				}
			}
		}
		else Zmax = R_cyl;

		for(int iz=0; iz <= NBIN_Z; iz++){
			ZZ = MIN_Z + iz * DX_Z;

			//Inside or outside?
			bool inside;
			if(ZZ < Zmax) inside = true;
			else inside = false;

			if(inside == true){
				POT = 0;
				FS = 0;
				FZ = 0;
			}
			else {
				if(cone == true){
					if(sphere==false){
						POT = 0.5 * KAPPA * (ZZ - Zmax) * (ZZ - Zmax);
						FZ = - KAPPA * (ZZ - Zmax);
						FS = - KAPPA * (ZZ - Zmax) * tg_alpha;
					}
			        else {
			        	D = sqrt(pow(ZZ,2)+pow(SS,2));
			        	d = D - z_cc;
			        	POT = 0.5 * KAPPA * pow(d,2);
			        	FZ = - KAPPA * d * ZZ / D;
			        	FS = - KAPPA * d * SS / D;
			        }
				}
				else {
					if(sphere==false){
						POT = 0.5 * KAPPA * (ZZ - Zmax) * (ZZ - Zmax);
						FZ = - KAPPA * (ZZ - Zmax);
						FS = 0;
					}
					else{
						D = sqrt(pow(ZZ,2)+pow(SS,2));
						d = D - z_cc;
						if(ZZ>=R_cyl+slope*(SS-z_cc)){
							POT = 0.5 * KAPPA * pow(d,2);
							FZ = - KAPPA * d * ZZ / D;
							FS = - KAPPA * d * SS / D;
						}
						else{
							POT = 0.5 * KAPPA * pow(sqrt(pow((ZZ+slope*z_cc-R_cyl)/slope,2)+pow(ZZ,2))-
									z_cc,2);
							FZ = - KAPPA*(sqrt(pow((ZZ+slope*z_cc-R_cyl)/slope,2)+pow(ZZ,2))-z_cc)*
									ZZ/sqrt(pow((ZZ+slope*z_cc-R_cyl)/slope,2)+pow(ZZ,2));
							FS = 0;
						}
					}
				}
			}
			pof.printf("%13.6lf %13.6lf %13.6lf %13.6lf %13.6lf\n", SS, ZZ, POT, FS, FZ);
		}
		pof.printf("\n");
	}
	pof.close();
}

void Funnel::calculate()
{
  unsigned ncv=getNumberOfArguments();
  vector<double> cv(ncv), der(ncv);

  for(unsigned i=0;i<ncv;++i){cv[i]=getArgument(i);}

//  log.printf("  In Funnel: %13.6lf  %13.6lf\n", cv[0], cv[1]);

  double ene=scale_*BiasGrid_->getValueAndDerivatives(cv,der);

  setBias(ene);

// set Forces
  for(unsigned i=0;i<ncv;++i){
   const double f=-scale_*der[i];
   setOutputForce(i,f);
  }
}

}
}
