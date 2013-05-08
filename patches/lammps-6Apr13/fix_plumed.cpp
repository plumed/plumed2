#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "respa.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "fix_plumed.h"

using namespace LAMMPS_NS;
using namespace PLMD;
using namespace FixConst;

FixPlumed::FixPlumed(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  p(NULL),
  nlocal(0),
  gatindex(NULL),
  masses(NULL),
  charges(NULL)
{
// Not sure this is really necessary:
  if (!atom->tag_enable) error->all(FLERR,"fix plumed requires atom tags");
// Initialize plumed:
  p=new PLMD::Plumed;
  p->cmd("setMPIComm",&world);

// Set up units
// LAMMPS units wrt kj/mol - nm - ps
// Set up units

  if (force->boltz == 1.0){
// LAMMPS units lj
    p->cmd("setNaturalUnits");
  } else {
    double energyUnits=1.0;
    double lengthUnits=1.0;
    double timeUnits=1.0;
    if (force->boltz == 0.0019872067){
// LAMMPS units real :: kcal/mol; angstrom; fs
      energyUnits=4.184;
      lengthUnits=0.1;
      timeUnits=0.001;
    } else if (force->boltz == 8.617343e-5){
// LAMMPS units metal :: eV; angstrom; ps
      energyUnits=96.48530749925792;
      lengthUnits=0.1;
      timeUnits=1.0;
    } else if (force->boltz == 1.3806504e-23){
// LAMMPS units si :: Joule, m; s
      energyUnits=0.001;
      lengthUnits=1.e-9;
      timeUnits=1.e-12;
    } else if (force->boltz == 1.3806504e-16){
// LAMMPS units cgs :: erg; cms;, s
      energyUnits=6.0221418e13;
      lengthUnits=1.e-7;
      timeUnits=1.e-12;
    } else if (force->boltz == 3.16681534e-6){
// LAMMPS units electron :: Hartree, bohr, fs
      energyUnits=2625.5257;
      lengthUnits=0.052917725;
      timeUnits=0.001;
    } else error->all(FLERR,"Odd LAMMPS units, plumed cannot work with that");
    p->cmd("setMDEnergyUnits",&energyUnits);
    p->cmd("setMDLengthUnits",&lengthUnits);
    p->cmd("setMDTimeUnits",&timeUnits);
  }

// Read fix parameters:
  int next=0;
  for(int i=3;i<narg;++i){
    if(!strcmp(arg[i],"outfile")) next=1;
    else if(next==1){
      p->cmd("setLogFile",arg[i]);
      next=0;
    }
    else if(!strcmp(arg[i],"plumedfile"))next=2;
    else if(next==2){
      p->cmd("setPlumedDat",arg[i]);
      next=0;
    }
    else error->all(FLERR,"syntax error in fix plumed - use 'fix name plumed plumedfile plumed.dat outfile plumed.out' ");
  }
  if(next==1) error->all(FLERR,"missing argument for outfile option");
  if(next==2) error->all(FLERR,"missing argument for plumedfile option");

  p->cmd("setMDEngine","LAMMPS");

  int natoms=int(atom->natoms);
  p->cmd("setNatoms",&natoms);

  double dt=update->dt;
  p->cmd("setTimestep",&dt);

  virial_flag=1;

// This is the real initialization:
  p->cmd("init");

}

FixPlumed::~FixPlumed()
{
  delete p;
}

int FixPlumed::setmask()
{
  // set with a bitmask how and when apply the force from plumed 
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

void FixPlumed::init()
{
  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

void FixPlumed::setup(int vflag)
{
  if (strcmp(update->integrate_style,"verlet") == 0)
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

void FixPlumed::min_setup(int vflag)
{
  post_force(vflag);
}

void FixPlumed::post_force(int vflag)
{
  int update_gatindex=0;
// Try to find out if the domain decomposition has been updated:
  if(nlocal!=atom->nlocal){
    if(charges) delete [] charges;
    if(masses) delete [] masses;
    if(gatindex) delete [] gatindex;
    nlocal=atom->nlocal;
    gatindex=new int [nlocal];
    masses=new double [nlocal];
    charges=new double [nlocal];
    update_gatindex=1;
  } else {
    for(int i=0;i<nlocal;i++){
      if(gatindex[i]!=atom->tag[i]-1){
        update_gatindex=1;
        break;
      }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE,&update_gatindex,1,MPI_INT,MPI_SUM,world);

// In case it has been updated, rebuild the local mass/charges array
// and tell plumed about the change:
  if(update_gatindex){
    for(int i=0;i<nlocal;i++){
      gatindex[i]=atom->tag[i]-1;
      masses[i]=atom->mass[atom->type[i]];
      charges[i]=atom->q[atom->type[i]];
    }
    p->cmd("setAtomsNlocal",&nlocal);
    p->cmd("setAtomsGatindex",gatindex);
  }


// set up local virial/box. plumed uses full 3x3 matrices
  double virial[3][3];
  for(int i=0;i<3;i++) for(int j=0;j<3;j++) virial[i][j]=0.0;
  double box[3][3];
  for(int i=0;i<3;i++) for(int j=0;j<3;j++) box[i][j]=0.0;
  box[0][0]=domain->h[0];
  box[1][1]=domain->h[1];
  box[2][2]=domain->h[2];
  box[2][1]=domain->h[3];
  box[2][0]=domain->h[4];
  box[1][0]=domain->h[5];
  
// local variable with timestep:
  int step=update->ntimestep;

// pass all pointers to plumed:
  p->cmd("setStep",&step);
  p->cmd("setPositions",&atom->x[0][0]);
  p->cmd("setBox",&box[0][0]);
  p->cmd("setForces",&atom->f[0][0]);
  p->cmd("setMasses",&masses[0]);
  p->cmd("setCharges",&charges[0]);
  p->cmd("setVirial",&virial[0][0]);

// do the real calculation:
  p->cmd("calc");

// retransform virial to lammps representation:
  Fix::virial[0]=-virial[0][0];
  Fix::virial[1]=-virial[1][1];
  Fix::virial[2]=-virial[2][2];
  Fix::virial[3]=-virial[0][1];
  Fix::virial[4]=-virial[0][2];
  Fix::virial[5]=-virial[1][2];
}

void FixPlumed::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

void FixPlumed::min_post_force(int vflag)
{
  post_force(vflag);
}

