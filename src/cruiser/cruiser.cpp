/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
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
#include "cltools/CLTool.h"
#include "cltools/CLToolRegister.h"
#include "wrapper/Plumed.h"
#include "tools/Vector.h"
#include "tools/Random.h"
#include "tools/Communicator.h"
#include <string>
#include <cstdio>
#include <cmath>
#include <vector>

using namespace std;

namespace PLMD{
namespace cruiser{

//+PLUMEDOC TOOLS cruiser
/*

Simple piece of code to do MD on a 3D energy landscape

\par Examples

*/
//+ENDPLUMEDOC

class Cruiser : public PLMD::CLTool {
  string description() const {
    return "dynamics of one atom on energy landscape";
  }
public:
  static void registerKeywords( Keywords& keys ){
    keys.add("compulsory","-nstep","The number of steps of dynamics you want to run");
    keys.add("compulsory","-tstep","0.005","the integration timestep");
    keys.add("compulsory","-temp","temperature units kT/epsilon");
    keys.add("compulsory","-lang","langevin thermostat, give relaxation time");
    keys.add("compulsory","-iseed","2983","value of random number seed");
    keys.addFlag("-plumed",false,"switch on metadynamics");
  }

  Cruiser( const CLToolOptions& co ) :
  CLTool(co)
  {
     inputdata=commandline;
  }

private:

  void read_input(double& temp, double& tstep, double& friction, int& nsteps, int& seed, bool& plumedon ){
     parse("-nstep",nsteps); parse("-tstep",tstep); parse("-temp",temp);
     parse("-lang",friction); parse("-iseed",seed); parseFlag("-plumed",plumedon);
     printf("Doing %d step simulation at temperature %f with friction of %f\n",nsteps,temp,friction);
     if(plumedon) printf("Reading input for plumed from plumed.dat\n" );
  }

public:

  double calc_energy( const std::vector<Vector>& pos, std::vector<Vector>& forces ){
     double x=sin( pos[0][0] ), y=sin( pos[0][1] ), z=sin( pos[0][2] );
     double tmp = exp( 3*( 3 - pow(x,4) - pow(y,4) - pow(z,4) ) );
     forces[0][0] = 12 * ( pow(x,3) ) * cos( pos[0][0] ) * tmp;
     forces[0][1] = 12 * ( pow(y,3) ) * cos( pos[0][1] ) * tmp;
     forces[0][2] = 12 * ( pow(z,3) ) * cos( pos[0][2] ) * tmp;
     return tmp - 1;
  }

  double calc_temp( const std::vector<Vector>& vel ){
     double total_KE=0.0;
     //! Double the total kinetic energy of the system
     for(unsigned j=0;j<3;++j) total_KE+=vel[0][j]*vel[0][j];
     return total_KE/3.;
  }

  int main( FILE* in, FILE* out, PLMD::Communicator& pc){
     bool plumedon;
     double temp, tstep, friction;
     int nsteps, seed; 
     int plumedWantsToStop;
     Random random;

     PLMD::Plumed* plumed=NULL;

     read_input( temp, tstep, friction, nsteps, seed, plumedon );
     // Setup random number generator
     random.setSeed(seed);

     if(plumedon) plumed=new PLMD::Plumed;

     if(plumed){
       int s=sizeof(double);
       plumed->cmd("setRealPrecision",&s);
       plumed->cmd("setMPIComm",&pc.Get_comm());
     }

     if(plumed){
       plumed->cmd("setNoVirial");
       int natoms=1; 
       plumed->cmd("setNatoms",&natoms);
       plumed->cmd("setMDEngine","cruiser");
       plumed->cmd("setTimestep",&tstep);
       plumed->cmd("setPlumedDat","plumed.dat");
       plumed->cmd("init");
     }

     double potential, therm_eng=0; std::vector<double> masses(1,1);
     std::vector<Vector> positions(1), velocities(1), forces(1);
     positions[0][0]=pi/2.; positions[0][1]=pi/2.; positions[0][2]=pi/2.;
     for(unsigned k=0;k<3;++k) velocities[0][k]=random.Gaussian() * sqrt( temp );

     potential=calc_energy(positions,forces); double ttt=calc_temp(velocities);
     double conserved = potential+1.5*ttt+therm_eng; FILE* fp=fopen("stats.out","w+");
     if( pc.Get_rank()==0 ) fprintf(fp,"%d %f %f %f %f %f \n", 0, 0., conserved, ttt, potential, therm_eng );

     for(unsigned istep=0;istep<nsteps;++istep){

        if( istep%20==0 && pc.Get_rank()==0 ) printf("Doing step %d\n",istep);

        // Langevin thermostat
        double lscale=exp(-0.5*tstep/friction);
        double lrand=sqrt((1.-lscale*lscale)*temp);
        for(unsigned k=0;k<3;++k){ 
          therm_eng=therm_eng+0.5*velocities[0][k]*velocities[0][k];
          velocities[0][k]=lscale*velocities[0][k]+lrand*random.Gaussian(); 
          therm_eng=therm_eng-0.5*velocities[0][k]*velocities[0][k];
        }

        // First step of velocity verlet
        for(unsigned k=0;k<3;++k){ 
           velocities[0][k] = velocities[0][k] + 0.5*tstep*forces[0][k];
           positions[0][k] = positions[0][k] + tstep*velocities[0][k];
           // Apply pbc
           if( positions[0][k]>pi ) positions[0][k]-=2*pi;
           if( positions[0][k]<=-pi ) positions[0][k]+=2*pi; 
        }

        potential=calc_energy(positions,forces);

        if(plumed){
          int istepplusone=istep+1;
          plumedWantsToStop=0;
          plumed->cmd("setStep",&istepplusone);
          plumed->cmd("setMasses",&masses[0]);
          plumed->cmd("setForces",&forces[0]);
          plumed->cmd("setEnergy",&potential);
          plumed->cmd("setPositions",&positions[0]);
          plumed->cmd("setStopFlag",&plumedWantsToStop);
          plumed->cmd("calc");
          if(istep%2000==0) plumed->cmd("writeCheckPointFile");
          if(plumedWantsToStop) nsteps=istep;
        }
        
        // Second step of velocity verlet
        for(unsigned k=0;k<3;++k){
           velocities[0][k] = velocities[0][k] + 0.5*tstep*forces[0][k];
        }

        // Langevin thermostat
        lscale=exp(-0.5*tstep/friction);
        lrand=sqrt((1.-lscale*lscale)*temp);
        for(unsigned k=0;k<3;++k){
          therm_eng=therm_eng+0.5*velocities[0][k]*velocities[0][k];
          velocities[0][k]=lscale*velocities[0][k]+lrand*random.Gaussian(); 
          therm_eng=therm_eng-0.5*velocities[0][k]*velocities[0][k];
        }

        // Print everything
        ttt = calc_temp( velocities ); conserved = potential+1.5*ttt+therm_eng;
        if( pc.Get_rank()==0 ) fprintf(fp,"%d %f %f %f %f %f \n", istep, istep*tstep, conserved, ttt, potential, therm_eng ); 
     }

     if(plumed) delete plumed;
     fclose(fp);

     return 0; 
  }
};

PLUMED_REGISTER_CLTOOL(Cruiser,"cruiser")

}
}
