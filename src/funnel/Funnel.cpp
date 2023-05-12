/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2019-2020

   This file is part of funnel code module.

   The FM code respects the CC BY-NC license.
   Users are free to download, adapt and use the code as long as it is not for commercial purposes.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "bias/Bias.h"
#include "bias/ActionRegister.h"
#include "tools/Grid.h"
#include "tools/Exception.h"
#include "tools/File.h"
#include <cstring>
#include "tools/Communicator.h"
#include "core/ActionSet.h"
#include "tools/FileBase.h"
#include <memory>
#include "core/PlumedMain.h"

using namespace std;
using namespace PLMD::bias;


namespace PLMD {
namespace funnel {

//+PLUMEDOC FUNNELMOD_BIAS FUNNEL
/*
Calculate a funnel-shape restraint potential that is defined on a grid that is read during the setup.

If the input file is not already present, it will create one with the name specified in the FILE flag.
The potential has a two-dimensional resolution since it has been devised to be used with the two
components of \ref FUNNEL_PS (i.e., fps.lp and fps.ld) and it is divided in two sections, a cone shape
attached to a cylindrical one. The user can customize the shape of both the sections by modifying a
number of flags. In particular the cone section of the funnel is calculated with the following formula:

\f[
MAX_Z=R_{cyl} + tg_{alpha} * (z_{cc} - MIN_S)
\f]

where  \f$ MAX_Z \f$ is the radius of the cone base,  \f$ R_{cyl} \f$ is the radius of the cylinder part,
\f$ tg_{alpha} \f$ is the angle regulating how steep the cone is, \f$ z_{cc} \f$ is the switching point
between cone and cylinder, and \f$ MIN_S \f$ is the lowest possible value assumed by fps.lp of \ref FUNNEL_PS.
As for the cylinder, it starts from the value of \f$ z_{cc} \f$ and stops at the value of \f$ MAX_S \f$
with a section of \f$ pi*r_{cyl}^2 \f$.

There is the option of transforming the cone region into a sphere with the use of the SPHERE flag. In this
case, the new shape will have a radius of \f$ z_{cc} \f$. It might be necessary tuning the SAFETY option
to select how much the potential extends from the sphere.

\par Examples

The following is an input for a calculation with a funnel potential that is defined in the file BIAS
and that acts on the collective variables defined by FUNNEL_PS.
\plumedfile
lig: COM ATOMS=3221,3224,3225,3228,3229,3231,3233,3235,3237
fps: FUNNEL_PS LIGAND=lig REFERENCE=start.pdb ANCHOR=2472 POINTS=4.724,5.369,4.069,4.597,5.721,4.343

FUNNEL ARG=fps.lp,fps.ld ZCC=1.8 ALPHA=0.55 RCYL=0.1 MINS=-0.5 MAXS=3.7 KAPPA=35100 NBINS=500 NBINZ=500 FILE=BIAS LABEL=funnel
\endplumedfile

The BIAS will then look something like this:
\auxfile{BIAS}
#! FIELDS fps.lp fps.ld funnel.bias der_fps.lp der_fps.ld
#! SET min_fps.lp -0.500000
#! SET max_fps.lp 3.700000
#! SET nbins_fps.lp 500.000000
#! SET periodic_fps.lp false
#! SET min_fps.ld 0.000000
#! SET max_fps.ld 1.510142
#! SET nbins_fps.ld 500.000000
#! SET periodic_fps.ld false
    -0.500000      0.000000      0.000000      0.000000      0.000000
    -0.500000      0.003020      0.000000      0.000000      0.000000
    -0.500000      0.006041      0.000000      0.000000      0.000000
    -0.500000      0.009061      0.000000      0.000000      0.000000
    -0.500000      0.012081      0.000000      0.000000      0.000000
    -0.500000      0.015101      0.000000      0.000000      0.000000
\endauxfile

The Funnel potential should always be used in combination with the collective variable  \ref FUNNEL_PS, since it
is constructed to take as inputs fps.lp and fps.ld (the former linepos and linedist of Funnel-Metadynamics
\cite FM).  In the first block of data the value of fps.lp (the value in the first column) is kept fixed
and the value of the function is given at 500 equally spaced values for fps.ld between 0 and 1.51. In
the second block of data fps.lp is fixed at \f$-0.5 + \frac{4.2}{500}\f$ and the value of the function
is given at 500 equally spaced values for fps.ld between 0 and 1.51. In the third block of data the same
is done but fps.lp is fixed at \f$-0.5 + \frac{8.4}{100}\f$ and so on until you get to the five hundredth
block of data where fps.lp is fixed at \f$3.7\f$.

It is possible to switch the shape of the cone region, transforming it in a sphere, with the flag SPHERE.
\plumedfile
lig: COM ATOMS=545,546,547,548,549,550,551,552,553
fps: FUNNEL_PS LIGAND=lig REFERENCE=ref.pdb ANCHOR=52 POINTS=2.793,3.696,3.942,3.607,4.298,3.452

FUNNEL ARG=fps.lp,fps.ld ZCC=4.0 RCYL=0.1 MINS=0.2 MAXS=4.9 KAPPA=100000 NBINS=500 NBINZ=500 SPHERE SAFETY=1.0 FILE=BIAS LABEL=funnel
\endplumedfile

*/
//+ENDPLUMEDOC
class Funnel : public Bias {

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

void Funnel::registerKeywords(Keywords& keys) {
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.addFlag("NOSPLINE",false,"specifies that no spline interpolation is to be used when calculating the energy and forces due to the external potential");
  keys.addFlag("SPARSE",false,"specifies that the external potential uses a sparse grid");
  keys.addFlag("SPHERE",false, "The Funnel potential including the binding site can be spherical instead of a cone");
  keys.add("compulsory","SCALE","1.0","a factor that multiplies the external potential, useful to invert free energies");
// old stuff?
  //  componentsAreNotOptional(keys);
  //  keys.addOutputComponent("bias","default","the instantaneous value of the bias potential");

  //Defining optional arguments
  keys.add("optional","NBINS","number of bins along fps.lp");
  keys.add("optional","NBINZ","number of bins along fps.ld");
  keys.add("optional","MINS","minimum value assumed by fps.lp, if the ligand is able to go beyond this value the simulation will crash");
  keys.add("optional","KAPPA","constant to be used for the funnel-shape restraint potential");
  keys.add("optional","RCYL","radius of the cylindrical section");
  keys.add("optional","SAFETY","To be used in case the SPHERE flag is chosen, it regulates how much the potential extends (in nm)");
  keys.add("optional","SLOPE","Adjust the behavior of the potential outside the funnel, greater values than 1.0 will tend to push the ligand more towards the cylinder and vice versa");
  keys.add("optional","ALPHA","angle to change the width of the cone section");
  //Defining compulsory arguments
  keys.add("compulsory","MAXS","MAXS","maximum value assumed by fps.lp");
  keys.add("compulsory","ZCC","ZCC","switching point between cylinder and cone");
  keys.add("compulsory","FILE","name of the Funnel potential file");
  keys.addFlag("WALKERS_MPI",false,"To be used when gromacs + multiple walkers are used");
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
  if(spline) {
    log.printf("  External potential uses spline interpolation\n");
  }
  if(sparsegrid) {
    log.printf("  External potential uses sparse grid\n");
  }

  // Non più necessario dalla 2.3
//  addComponent("bias"); componentIsNotPeriodic("bias");

  std::string funcl=getLabel() + ".bias";

//  int size = plumed.comm.Get_size();
//  int rank = plumed.comm.Get_rank();
  IFile I_hate_this;
  bool do_exist=I_hate_this.FileExist(file);

  if(walkers_mpi) {
    if(comm.Get_rank()==0 && multi_sim_comm.Get_rank()==0) {
      if(!do_exist) {
        createBIAS(RCYL, ZCC, ALPHA, KAPPA, MINS, MAXS, NBINS, NBINZ, safety, sphere, slope, funcl, file);
      }
    }
    multi_sim_comm.Barrier();
  } else {
    if(comm.Get_rank()==0) {
      if(!do_exist) {
        createBIAS(RCYL, ZCC, ALPHA, KAPPA, MINS, MAXS, NBINS, NBINZ, safety, sphere, slope, funcl, file);
      }
    }
  }

  /*
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
  */
  comm.Barrier();

// read grid
  IFile gridfile;
  gridfile.open(file);
  BiasGrid_=Grid::create(funcl,getArguments(),gridfile,sparsegrid,spline,true);
//not necessary anymore?  gridfile.close();
  if(BiasGrid_->getDimension()!=getNumberOfArguments()) error("mismatch between dimensionality of input grid and number of arguments");
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    if( getPntrToArgument(i)->isPeriodic()!=BiasGrid_->getIsPeriodic()[i] ) error("periodicity mismatch between arguments and input bias");
  }
  comm.Barrier();
  if(comm.Get_rank()==0 && walkers_mpi) multi_sim_comm.Barrier();
  log<<"  Bibliography "<<plumed.cite("Limongelli, Bonomi, and Parrinello, PNAS 110, 6358 (2013)")<<"\n";
}


void Funnel::createBIAS(const double& R_cyl, const double& z_cc, const double& alpha,
                        const double& KAPPA, const double& MIN_S, const double& MAX_S, const double& NBIN_S,
                        const double& NBIN_Z, const double& safety, const bool& sphere, const double& slope,
                        const string& funcl, const string& file) {
  //R_cyl and z_cc forms the parameters of the cylinder.
  //alpha defines the angle in degrees.

  //PARAMETERS OF THE CONE
  double tg_alpha= tan(alpha);

  //parameters for PROGRESSION
  //parameters for DISTANCE
  double MIN_Z=0;
  double MAX_Z;
  if (sphere==false) {
    MAX_Z=R_cyl + tg_alpha * (z_cc - MIN_S);
  }
  else {
    MAX_Z=z_cc+safety;
  }

  //bin size
  double DX_Z = (MAX_Z - MIN_Z) / NBIN_Z;
  double DX_S;
  if (sphere==false) {
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

  for(int is=0; is <= NBIN_S; is++) {
    if (sphere==false) {
      SS = MIN_S + is * DX_S;
    }
    else {
      SS = - z_cc - safety + is * DX_S;
    }
    bool cone = false;
    if (sphere==false) {
      if(SS <= z_cc) cone = true;
    }
    else {
      if (SS <= sqrt(pow(z_cc,2)-pow(R_cyl,2))) cone = true;
    }
    //Set wall boundaries properly
    if(cone == true) {
      if(sphere==false) {
        Zmax = R_cyl + (z_cc - SS) * tg_alpha;
      }
      else {
        if (SS > -z_cc) {
          Zmax = sqrt(pow(z_cc,2) - pow(SS,2));
        }
        else {
          Zmax = 0;
        }
      }
    }
    else Zmax = R_cyl;

    for(int iz=0; iz <= NBIN_Z; iz++) {
      ZZ = MIN_Z + iz * DX_Z;

      //Inside or outside?
      bool inside;
      if(ZZ < Zmax) inside = true;
      else inside = false;

      if(inside == true) {
        POT = 0;
        FS = 0;
        FZ = 0;
      }
      else {
        if(cone == true) {
          if(sphere==false) {
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
          if(sphere==false) {
            POT = 0.5 * KAPPA * (ZZ - Zmax) * (ZZ - Zmax);
            FZ = - KAPPA * (ZZ - Zmax);
            FS = 0;
          }
          else {
            D = sqrt(pow(ZZ,2)+pow(SS,2));
            d = D - z_cc;
            if(ZZ>=R_cyl+slope*(SS-z_cc)) {
              POT = 0.5 * KAPPA * pow(d,2);
              FZ = - KAPPA * d * ZZ / D;
              FS = - KAPPA * d * SS / D;
            }
            else {
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

  for(unsigned i=0; i<ncv; ++i) {
    cv[i]=getArgument(i);
  }

//  log.printf("  In Funnel: %13.6lf  %13.6lf\n", cv[0], cv[1]);

  double ene=scale_*BiasGrid_->getValueAndDerivatives(cv,der);

  setBias(ene);

// set Forces
  for(unsigned i=0; i<ncv; ++i) {
    const double f=-scale_*der[i];
    setOutputForce(i,f);
  }
}

}
}
