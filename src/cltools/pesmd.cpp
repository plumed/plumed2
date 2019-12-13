/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2019 The plumed team
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
#include "core/PlumedMain.h"
#include "tools/Vector.h"
#include "tools/Random.h"
#include "tools/Communicator.h"
#include <string>
#include <cstdio>
#include <cmath>
#include <vector>
#include <memory>

//+PLUMEDOC TOOLS pesmd
/*
Pesmd allows one to do (biased) Langevin dynamics on a two-dimensional potential energy surface.

The energy landscape that you are moving about on is specified using a plumed input file.
The directives that are available for this command line tool are as follows:

\par Examples

You run a Langevin simulation using pesmd with the following command:
\verbatim
plumed pesmd < input
\endverbatim

The following is an example of an input file for a pesmd simulation. This file
instructs pesmd to do 50 steps of Langevin dynamics on a 2D potential energy surface
at a temperature of 0.722
\verbatim
temperature 0.722
tstep 0.005
friction 1
dimension 2
nstep 50
ipos 0.0 0.0
\endverbatim

If you run the following a description of all the directives that can be used in the
input file will be output.
\verbatim
plumed pesmd --help
\endverbatim

The energy landscape to explore is given within the plumed input file.  For example the following
example input uses \ref MATHEVAL to define a two dimensional potential.

\verbatim
d1: DISTANCE ATOMS=1,2 COMPONENTS
ff: MATHEVAL ARG=d1.x,d1,y PERIODIC=NO FUNC=()
bb: BIASVALUE ARG=ff
\endverbatim

Atom 1 is placed at the origin.  The x and y components on our surface are the
positions of the particle on our two dimensional energy landscape.  By calculating the
vector connecting atom 1 (the origin) to atom 2 (the position of our particle) we are thus
getting the position of the atom on the energy landscape.  This is then inserted into the function
that is calculated on the second line.  The value of this function is then used as a bias.

We can also specify a potential on a grid and look at the dynamics on this function using pesmd.
A plumed input for an example such as this one might look something like this:

\verbatim
d1: DISTANCE ATOMS=1,2 COMPONENTS
bb: EXTERNAL ARG=d1.x,d1,y FILE=fes.dat
\endverbatim

In this way we can use pesmd to do a dynamics on a free energy surface calculated using metadynamics
and sum_hills.  On a final note once we have defined our potential we can use all the biasing functions
within plumed in addition in order to do a biased dynamics on the potential energy landscape of interest.

*/
//+ENDPLUMEDOC

using namespace std;

namespace PLMD {
namespace cltools {

class PesMD  : public PLMD::CLTool {
  string description() const override {
    return "Langevin dynamics on PLUMED energy landscape";
  }
public:
  static void registerKeywords( Keywords& keys ) {
    keys.add("compulsory","nstep","The number of steps of dynamics you want to run");
    keys.add("compulsory","temperature","NVE","the temperature at which you wish to run the simulation in LJ units");
    keys.add("compulsory","friction","off","The friction (in LJ units) for the Langevin thermostat that is used to keep the temperature constant");
    keys.add("compulsory","tstep","0.005","the integration timestep in LJ units");
    keys.add("compulsory","dimension","the dimension of your energy landscape");
    keys.add("compulsory","plumed","plumed.dat","the name of the plumed input file containing the potential");
    keys.add("compulsory","ipos","0.0","the initial position of the system");
    keys.add("compulsory","idum","0","The random number seed");
    keys.addFlag("periodic","false","are your input coordinates periodic");
    keys.add("optional","min","minimum value the coordinates can take for a periodic domain");
    keys.add("optional","max","maximum value the coordinates can take for a periodic domain");
  }

  explicit PesMD( const CLToolOptions& co ) :
    CLTool(co)
  {
    inputdata=ifile;
  }

private:

  void read_input(double& temperature,
                  double& tstep,
                  double& friction,
                  int& dim,
                  std::string& plumedin,
                  std::vector<double>& ipos,
                  int&    nstep,
                  bool&   lperiod,
                  std::vector<double>& periods,
                  int&    idum)
  {
    // Read everything from input file
    std::string tempstr; parse("temperature",tempstr);
    if( tempstr!="NVE" ) Tools::convert(tempstr,temperature);
    parse("tstep",tstep);
    std::string frictionstr; parse("friction",frictionstr);
    if( tempstr!="NVE" ) {
      if(frictionstr=="off") { fprintf(stderr,"Specify friction for thermostat\n"); exit(1); }
      Tools::convert(frictionstr,friction);
    }
    parse("plumed",plumedin); parse("dimension",dim);
    parse("nstep",nstep); parse("idum",idum);
    ipos.resize( dim ); parseVector("ipos",ipos);

    parseFlag("periodic",lperiod);
    if( lperiod ) {
      if( dim>3 ) error("can only do three dimensional periodic functions");
      std::vector<double> min( dim ); parseVector("min",min);
      std::vector<double> max( dim ); parseVector("max",max);
      periods.resize( dim );
      for(int i=0; i<dim; ++i) {
        if( max[i]<min[i] ) error("invalid periods specified max is less than min");
        periods[i]=max[i]-min[i];
      }
    }
  }


public:

  int main( FILE* in, FILE* out, PLMD::Communicator& pc) override {
    std::string plumedin; std::vector<double> ipos;
    double temp, tstep, friction; bool lperiod;
    int dim, nsteps, seed; std::vector<double> periods;
    int plumedWantsToStop;
    Random random;

    read_input( temp, tstep, friction, dim, plumedin, ipos, nsteps, lperiod, periods, seed );
    // Setup random number generator
    random.setSeed(seed);

    // Setup box if we have periodic domain
    std::vector<double> box(9, 0.0);
    if( lperiod && dim==1 ) { box[0]=box[5]=box[9]=periods[0]; }
    else if( lperiod && dim==2 ) { box[0]=periods[0]; box[5]=box[9]=periods[1]; }
    else if( lperiod && dim==3 ) { box[0]=periods[0]; box[5]=periods[1]; box[9]=periods[2]; }
    else if( lperiod ) error("invalid dimension for periodic potential must be 1, 2 or 3");

    // Create plumed object and initialize
    std::unique_ptr<PLMD::PlumedMain> plumed(new PLMD::PlumedMain);
    int s=sizeof(double);
    plumed->cmd("setRealPrecision",&s);
    if(Communicator::initialized()) plumed->cmd("setMPIComm",&pc.Get_comm());
    plumed->cmd("setNoVirial");
    int natoms=( std::floor(dim/3) +  2 );
    plumed->cmd("setNatoms",&natoms);
    plumed->cmd("setMDEngine","pesmd");
    plumed->cmd("setTimestep",&tstep);
    plumed->cmd("setPlumedDat",plumedin.c_str());
    plumed->cmd("init");

    // Now create some fake atoms
    int nat = std::floor( dim/3 ) + 1;
    std::vector<double> masses( 1+nat, 1 );
    std::vector<Vector> velocities( nat ), positions( nat+1 ), forces( nat+1 );
    // Will set these properly eventually
    int k=0; positions[0].zero(); // Atom zero is fixed at origin
    for(int i=0; i<nat; ++i) for(unsigned j=0; j<3; ++j) {
        if( k<dim ) { positions[1+i][j]=ipos[k]; } else { positions[1+i][j]=0;}
        k++;
      }
    // And initialize the velocities
    for(int i=0; i<nat; ++i) for(int j=0; j<3; ++j) velocities[i][j]=random.Gaussian() * sqrt( temp );
    // And calcualte the kinetic energy
    double tke=0;
    for(int i=0; i<nat; ++i) {
      for(int j=0; j<3; ++j) {
        if( 3*i+j>dim-1 ) break;
        tke += 0.5*velocities[i][j]*velocities[i][j];
      }
    }

    // Now call plumed to get initial forces
    int istep=0; double zero=0;
    plumed->cmd("setStep",&istep);
    plumed->cmd("setMasses",&masses[0]);
    for(unsigned i=0; i<forces.size(); ++i) forces[i].zero();
    plumed->cmd("setForces",&forces[0]);
    plumed->cmd("setEnergy",&zero);
    if( lperiod ) plumed->cmd("setBox",&box[0]);
    plumed->cmd("setPositions",&positions[0]);
    plumed->cmd("calc");


    double therm_eng=0;
    FILE* fp=fopen("stats.out","w+");

    for(int istep=0; istep<nsteps; ++istep) {

      if( istep%20==0 && pc.Get_rank()==0 ) printf("Doing step %i\n",istep);

      // Langevin thermostat
      double lscale=exp(-0.5*tstep/friction);
      double lrand=sqrt((1.-lscale*lscale)*temp);
      for(int j=0; j<nat; ++j) {
        for(int k=0; k<3; ++k) {
          if( 3*j+k>dim-1 ) break;
          therm_eng=therm_eng+0.5*velocities[j][k]*velocities[j][k];
          velocities[j][k]=lscale*velocities[j][k]+lrand*random.Gaussian();
          therm_eng=therm_eng-0.5*velocities[j][k]*velocities[0][k];
        }
      }

      // First step of velocity verlet
      for(int j=0; j<nat; ++j) {
        for(int k=0; k<3; ++k) {
          if( 3*j+k>dim-1 ) break;
          velocities[j][k] = velocities[j][k] + 0.5*tstep*forces[1+j][k];
          positions[1+j][k] = positions[1+j][k] + tstep*velocities[j][k];
        }
      }

      int istepplusone=istep+1;
      plumedWantsToStop=0;
      plumed->cmd("setStep",&istepplusone);
      plumed->cmd("setMasses",&masses[0]);
      for(unsigned i=0; i<forces.size(); ++i) forces[i].zero();
      plumed->cmd("setForces",&forces[0]);
      double fenergy=0.0;
      plumed->cmd("setEnergy",&fenergy);
      plumed->cmd("setPositions",&positions[0]);
      plumed->cmd("setStopFlag",&plumedWantsToStop);
      plumed->cmd("calc");
      // if(istep%2000==0) plumed->cmd("writeCheckPointFile");
      if(plumedWantsToStop) nsteps=istep;

      // Second step of velocity verlet
      for(int j=0; j<nat; ++j) {
        for(int k=0; k<3; ++k) {
          if( 3*j+k>dim-1 ) break;
          velocities[j][k] = velocities[j][k] + 0.5*tstep*forces[1+j][k];
        }
      }

      // Langevin thermostat
      lscale=exp(-0.5*tstep/friction);
      lrand=sqrt((1.-lscale*lscale)*temp);
      for(int j=0; j<nat; ++j) {
        for(int k=0; k<3; ++k) {
          if( 3*j+k>dim-1) break;
          therm_eng=therm_eng+0.5*velocities[j][k]*velocities[j][k];
          velocities[j][k]=lscale*velocities[j][k]+lrand*random.Gaussian();
          therm_eng=therm_eng-0.5*velocities[j][k]*velocities[j][k];
        }
      }
      // Calculate total kinetic energy
      tke=0;
      for(int i=0; i<nat; ++i) {
        for(int j=0; j<3; ++j) {
          if( 3*i+j>dim-1 ) break;
          tke += 0.5*velocities[i][j]*velocities[i][j];
        }
      }

      // Print everything
      // conserved = potential+1.5*ttt+therm_eng;
      if( pc.Get_rank()==0 ) fprintf(fp,"%i %f %f %f \n", istep, istep*tstep, tke, therm_eng );
    }

    fclose(fp);

    return 0;
  }
};

PLUMED_REGISTER_CLTOOL(PesMD,"pesmd")

}
}
