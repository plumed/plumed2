/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015,2016 The plumed team
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
#include "Bias.h"
#include "ActionRegister.h"
#include "tools/Random.h"
#include "tools/File.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"

#include <iostream>


using namespace std;


namespace PLMD{
namespace bias{

//+PLUMEDOC BIAS EDS
/*
Add a linear bias on a set of observables.

This force is the same as the linear part of the bias in \ref RESTRAINT, 
but this bias has the ability to compute prefactors
adaptively using the scheme of White and Voth (JCTC 2014) in order to match 
target observable values for a set of CVs.

The addition to the potential is of the form 
\f[
  \sum_i \frac{\alpha}_i}{s_i}*x_i
\f]

where for CV \f$x_i\f$, a coupling \f${\alpha}_i\f$ is determined
adaptively or set by the user to match a target value for
\f$x_i\f$. \f$s_i$\f is a scale parameter, which by default is set to
the target value. It may also be set separately. 

\warning
If the target observable value is set to zero while using the adaptive scheme this will cause a divide-by-zero error

\par Examples

The following input for a harmonic oscillator of two beads will adaptively find a linear bias to change the mean and variance to the target values. The PRINT line shows how to access the value of the coupling constants. 

\verbatim 
dist: DISTANCE ATOMS=1,2
# this is the squared of the distance
dist2: COMBINE ARG=dist POWERS=2 PERIODIC=NO

#bias mean and variance
eds: EDS ARG=dist,dist2 CENTER=2.0,1.0 PERIOD=50000 TEMP=1.0 
PRINT ARG=dist,dist2,eds.dist_coupling,eds.dist2_coupling,eds.bias,eds.force2 FILE=colvars.dat STRIDE=100
\endverbatim

Rather than trying to find the coupling constants adaptively, can ramp up to a constant value.
\verbatim
#ramp couplings from 0,0 to -1,1 over 50000 steps
eds: EDS ARG=dist,dist2 CENTER=2.0,1.0 FIXED=-1,1 RAMP PERIOD=50000 TEMP=1.0

#same as above, except starting at -0.5,0.5 rather than default of 0,0
eds: EDS ARG=dist,dist2 CENTER=2.0,1.0 FIXED=-1,1 INIT=-0.5,0.5 RAMP PERIOD=50000 TEMP=1.0
\endverbatim

A restart file can be added to dump information needed to restart/continue simulation using these parameters every STRIDE.
\verbatim 
#add the option to write to a restart file
eds: EDS ARG=dist,dist2 CENTER=2.0,1.0 PERIOD=50000 TEMP=1.0 ORESTARTFILE=restart.dat

#add the option to read in a previous restart file
eds: EDS ARG=dist,dist2 CENTER=2.0,1.0 PERIOD=50000 TEMP=1.0 IRESTARTFILE=restart.dat EDSRESTART

#add the option to read in a previous restart file and freeze the bias at the final level from the previous simulation
eds: EDS ARG=dist,dist2 CENTER=2.0,1.0 PERIOD=50000 TEMP=1.0 IRESTARTFILE=restart.dat EDSRESTART FREEZE
\endverbatim


*/
//+ENDPLUMEDOC

class EDS : public Bias{
  std::vector<double> center;
  std::vector<double> scale;
  std::vector<double> current_coupling;
  std::vector<double> set_coupling;
  std::vector<double> target_coupling;
  std::vector<double> max_coupling_range;
  std::vector<double> max_coupling_grad;
  std::vector<double> coupling_rate;
  std::vector<double> coupling_accum;
  std::vector<double> means;
  std::vector<double> ssds;
  std::vector<Value*> outCoupling;
  std::string _irestartfilename;
  std::string _orestartfilename;
  OFile orestartfile_;
  IFile irestartfile_;
  bool adaptive;
  bool freeze;
  bool equilibration;
  bool ramp;
  bool restart;
  bool b_write_restart;
  int seed;
  int update_calls;
  int avg_coupling_count;
  int update_period;
  double kbt;
  double coupling_range_increase_factor;
  int b_hard_coupling_range;
  Random rand;
  Value* valueForce2;

public:
  explicit EDS(const ActionOptions&);
  void calculate();
  void update();
  void turnOnDerivatives();
  void read_irestart();
  void setup_orestart();
  void write_orestart();
  static void registerKeywords(Keywords& keys);
  ~EDS();
};

PLUMED_REGISTER_ACTION(EDS,"EDS")

void EDS::registerKeywords(Keywords& keys){
   Bias::registerKeywords(keys);
   keys.use("ARG");
   keys.add("compulsory","CENTER","The desired centers (equilibrium values) which will be sought during the adaptive linear biasing.");
   keys.add("compulsory","PERIOD","Steps over which to adjust bias");

   keys.add("compulsory","RANGE","3.0","The largest magnitude of the force constant which one expects (in kBT) for each CV based");
   keys.add("compulsory","SEED","0","Seed for random order of changing bias");
   keys.add("compulsory","INIT","0","Starting value for coupling coefficients");
   keys.add("compulsory","FIXED","0","Fixed target values for bias factors (not adaptive)");
   keys.add("optional","BIAS_SCALE","A divisor to set the units of the bias. If not set, this will be the experimental value by default (as is done in White and Voth 2014).");
   keys.add("optional","TEMP","The system temperature. If not provided will be taken from MD code (if available)");

   keys.add("optional","ORESTARTFILE","Output file for all information needed to continue EDS simulation");
   keys.add("optional","IRESTARTFILE","Read this file to continue an EDS simulation (if same as above, will be overwritten)");

   keys.addFlag("RAMP",false,"Slowly increase bias constant to a fixed value");
   keys.addFlag("FREEZE",false,"Fix bias at current level (only used for restarting). Can also set PERIOD=0 if not using EDSRESTART.");
   keys.addFlag("EDSRESTART",false,"Get settings from IRESTARTFILE");

   //keys.addOutputComponent("bias","default","the instantaneous value of the bias potential");
   keys.addOutputComponent("force2","default","squared value of force from the bias");
   keys.addOutputComponent("_coupling","default","For each named CV biased, there will be a corresponding output CV_coupling storing the current linear bias prefactor.");
}

EDS::EDS(const ActionOptions&ao):
PLUMED_BIAS_INIT(ao),
center(getNumberOfArguments(),1.0),
current_coupling(getNumberOfArguments(),0.0),
set_coupling(getNumberOfArguments(),0.0),
target_coupling(getNumberOfArguments(),0.0),
max_coupling_range(getNumberOfArguments(),3.0),
max_coupling_grad(getNumberOfArguments(),0.0),
coupling_rate(getNumberOfArguments(),1.0),
coupling_accum(getNumberOfArguments(),1.0),
means(getNumberOfArguments(),0.0),
ssds(getNumberOfArguments(),0.0),
outCoupling(getNumberOfArguments(),NULL),
_irestartfilename(""),
_orestartfilename(""),
update_period(0),
kbt(0.0),
coupling_range_increase_factor(1.25),
b_hard_coupling_range(0),
adaptive(true),
equilibration(true),
ramp(false),
freeze(false),
restart(false),
b_write_restart(false),
seed(0),
update_calls(0),
avg_coupling_count(1),
valueForce2(NULL)
{
  double temp=-1.0;

  addComponent("force2");
  componentIsNotPeriodic("force2");
  valueForce2=getPntrToComponent("force2");

  for(unsigned i=0;i<getNumberOfArguments();i++){
    std::string comp=getPntrToArgument(i)->getName()+"_coupling";
    addComponent(comp);
    componentIsNotPeriodic(comp);
    outCoupling[i]=getPntrToComponent(comp);
  }
    
  parseVector("CENTER",center);
  parseVector("BIAS_SCALE", scale);
  parseVector("RANGE",max_coupling_range);
  parseVector("FIXED",target_coupling);
  parseVector("INIT",set_coupling);
  parse("PERIOD",update_period);
  parse("TEMP",temp);
  parse("SEED",seed);
  parse("ORESTARTFILE",_orestartfilename);
  parseFlag("RAMP",ramp);
  parseFlag("FREEZE",freeze);
  parseFlag("EDSRESTART",restart);
  parse("IRESTARTFILE",_irestartfilename);
  parse("ORESTARTFILE",_orestartfilename);
  checkRead();  

  if(restart){
      log.printf("  reading all simulation information from file: %s\n",_irestartfilename.c_str());
      read_irestart();
  }
  else{
    
      if(temp>=0.0) kbt=plumed.getAtoms().getKBoltzmann()*temp;
      else kbt=plumed.getAtoms().getKbT();
    
      log.printf("  with kBT = %f\n",kbt);
      log.printf("  Updating every %i steps\n",update_period);
    
      log.printf("  with centers:");
      for(unsigned i=0;i<center.size();i++){
          log.printf(" %f",center[i]);
      }

      log.printf("\n  and scaling:");
      if(scale.size() > 0  && scale.size() < getNumberOfArguments()) {
	error("the number of BIAS_SCALE values be the same as number of CVs");
      } else if(scale.size() == 0) {

	log.printf("(default) ");

	scale.resize(center.size());
	for(unsigned int i = 0; i < scale.size(); i++) {
          if(center[i]==0) 
	    error("BIAS_SCALE parameter has been set to CENTER value of 0 (as is default). This will divide by 0, so giving up. See doc for EDS bias");
	  scale[i] = center[i];
	}
      }

      for(unsigned int i = 0; i < scale.size(); i++) 
	log.printf(" %f",center[i]);
      
      log.printf("\n  with initial ranges / rates:\n");
      for(unsigned i=0;i<max_coupling_range.size();i++) {
          //this is just an empirical guess. Bigger range, bigger grads. Less frequent updates, bigger changes
          max_coupling_range[i]*=kbt;
          max_coupling_grad[i] = max_coupling_range[i]*update_period/100.;
          log.printf("    %f / %f\n",max_coupling_range[i],max_coupling_grad[i]);
      }
    
      if(seed>0){
         log.printf("  setting random seed = %i",seed);
         rand.setSeed(seed);
      }
    
      for(unsigned i=0;i<getNumberOfArguments();++i) if(target_coupling[i]!=0.0) adaptive=false;
    
      if(!adaptive){
        if(ramp>0) {
            log.printf("  ramping up coupling constants over %i steps\n",update_period);
        }
    
        log.printf("  with starting coupling constants");
        for(unsigned i=0;i<set_coupling.size();i++) log.printf(" %f",set_coupling[i]);
        log.printf("\n");
        log.printf("  and final coupling constants");
        for(unsigned i=0;i<target_coupling.size();i++) log.printf(" %f",target_coupling[i]);
        log.printf("\n");
      }
    
      //now do setup
      if(ramp){
          update_period*=-1;
      }
    
      for(unsigned i=0;i<set_coupling.size();i++) current_coupling[i] = set_coupling[i];
    
      // if adaptive, then first half will be used for equilibrating and second half for statistics
      if(update_period>0){
          update_period/=2;
      }


    }

    if(freeze){
        adaptive=false;
        log.printf("  freezing bias at current level\n");
    }

    if(_orestartfilename.length()>0) {
        log.printf("  writing restart information every %i steps to file: %s\n",abs(update_period),_orestartfilename.c_str());
        b_write_restart = true;
        setup_orestart();
    }

    log<<"  Bibliography "<<plumed.cite("White and Voth, J. Chem. Theory Comput. 10 (8), 3023-3030 (2014)")<<"\n";
}

void EDS::read_irestart(){
    int adaptive_i;


    irestartfile_.open(_irestartfilename);

    //some sample code to get the field names:
    /*
    std::vector<std::string> fields;
    irestartfile_.scanFieldList(fields);
    log.printf("field");
    for(unsigned i=0;i<fields.size();i++) log.printf(" %s",fields[i].c_str());
    log.printf("\n");
    */

    if(irestartfile_.FieldExist("kbt")){
        irestartfile_.scanField("kbt",kbt);
    }else{ error("No field 'kbt' in restart file"); }
    log.printf("  with kBT = %f\n",kbt);

    if(irestartfile_.FieldExist("update_period")){
        irestartfile_.scanField("update_period",update_period);
    }else{ error("No field 'update_period' in restart file"); }
    log.printf("  Updating every %i steps\n",update_period);

    if(irestartfile_.FieldExist("adaptive")){
        //note, no version of scanField for boolean
        irestartfile_.scanField("adaptive",adaptive_i);
    }else{ error("No field 'adaptive' in restart file"); }
    adaptive = bool(adaptive_i);

    if(irestartfile_.FieldExist("seed")){
        irestartfile_.scanField("seed",seed);
    }else{ error("No field 'seed' in restart file"); }
    if(seed>0){
       log.printf("  setting random seed = %i",seed);
       rand.setSeed(seed);
    }

    double time;
    std::string cv_name;

    while(irestartfile_.scanField("time",time)){

        for(unsigned i=0;i<getNumberOfArguments();++i) {
   	    cv_name = getPntrToArgument(i)->getName();
            irestartfile_.scanField(cv_name + +"_center",center[i]);
            irestartfile_.scanField(cv_name + +"_scale",scale[i]);
            irestartfile_.scanField(cv_name + "_init", set_coupling[i]);
            irestartfile_.scanField(cv_name + "_target",target_coupling[i]);
            irestartfile_.scanField(cv_name + "_coupling",current_coupling[i]);
            irestartfile_.scanField(cv_name + "_maxrange",max_coupling_range[i]);
            irestartfile_.scanField(cv_name + "_maxgrad",max_coupling_grad[i]);
       }

       irestartfile_.scanField();
   }

   log.printf("  with centers:");
   for(unsigned i=0;i<center.size();i++) {
       log.printf(" %f",center[i]);
   }
   log.printf("\n  and scaling:");
   for(unsigned i=0;i<scale.size();i++) {
       log.printf(" %f",scale[i]);
   }

   log.printf("\n  with initial ranges / rates:\n");
   for(unsigned i=0;i<max_coupling_range.size();i++) {
      log.printf("    %f / %f\n",max_coupling_range[i],max_coupling_grad[i]);
   }
    
   if(!adaptive && update_period<0){
       log.printf("  ramping up coupling constants over %i steps\n",-update_period);
   }

   log.printf("  with current coupling constants:\n    ");
   for(unsigned i=0;i<current_coupling.size();i++) log.printf(" %f",current_coupling[i]);
   log.printf("\n");
   log.printf("  with initial coupling constants:\n    ");
   for(unsigned i=0;i<set_coupling.size();i++) log.printf(" %f",set_coupling[i]);
   log.printf("\n");
   log.printf("  and final coupling constants:\n    ");
   for(unsigned i=0;i<target_coupling.size();i++) log.printf(" %f",target_coupling[i]);
   log.printf("\n");
    
   irestartfile_.close();
}

void EDS::setup_orestart(){
    orestartfile_.link(*this);
    orestartfile_.open(_orestartfilename);
    orestartfile_.setHeavyFlush();

    orestartfile_.addConstantField("adaptive").printField("adaptive",adaptive);
    orestartfile_.addConstantField("update_period").printField("update_period",update_period);
    orestartfile_.addConstantField("seed").printField("seed",seed);
    orestartfile_.addConstantField("kbt").printField("kbt",kbt);

    write_orestart();
}

void EDS::write_orestart(){
    std::string cv_name;
    orestartfile_.printField("time",getTimeStep()*getStep());

    for(unsigned i=0;i<getNumberOfArguments();++i) {
        cv_name = getPntrToArgument(i)->getName();

        orestartfile_.printField(cv_name + "_center",center[i]);
        orestartfile_.printField(cv_name + "_scale",center[i]);
        orestartfile_.printField(cv_name + "_init",set_coupling[i]);
        orestartfile_.printField(cv_name + "_target",target_coupling[i]);
        orestartfile_.printField(cv_name + "_coupling",current_coupling[i]);
        orestartfile_.printField(cv_name + "_maxrange",max_coupling_range[i]);
        orestartfile_.printField(cv_name + "_maxgrad",max_coupling_grad[i]);
    }
    orestartfile_.printField();
}

void EDS::calculate(){
  int ncvs = getNumberOfArguments();

  //Compute linear force as in "restraint"
  double ene=0.0;
  double totf2=0.0;
  for(unsigned i=0;i<ncvs;++i){
    const double cv=difference(i,center[i],getArgument(i));
    const double m=current_coupling[i];
    const double f=-m;
    ene+=m*cv;
    setOutputForce(i,f);
    totf2+=f*f;
  };
  setBias(ene);
  valueForce2->set(totf2);
  
  //adjust parameters according to EDS recipe
  update_calls++;

  if(b_write_restart && update_calls%abs(update_period)==0){
     write_orestart();
  }

  int b_finished_equil_flag = 1;
  double delta;

  //assume forces already applied and saved
  
  for(unsigned i=0;i<ncvs;++i){
      //are we ramping to a constant value and not done equilibrating
      if(update_period<0){
          if(update_calls<fabs(update_period) && !freeze){
              current_coupling[i] += (target_coupling[i]-set_coupling[i])/fabs(update_period);
          }
          //make sure we don't reset update calls
          b_finished_equil_flag = 0;
          continue;
      } 
      //not updating
      else if(update_period==0){
          continue;
      }

      //if we aren't wating for the bias to equilibrate, collect data
      if(!equilibration){
          //Welford, West, and Hanso online variance method
          delta = difference(i,means[i],getArgument(i));
          means[i]+=delta/update_calls;
          ssds[i]+=delta*difference(i,means[i],getArgument(i));
      }
      // equilibrating
      else {
          //check if we've reached the setpoint
          if(coupling_rate[i]==0 || pow(current_coupling[i]-set_coupling[i],2)<pow(coupling_rate[i],2)) {
             b_finished_equil_flag &= 1;
          }
          else{
              current_coupling[i]+=coupling_rate[i];
              b_finished_equil_flag = 0;
          }
      }
      //Update max coupling range if allowed
      if(!b_hard_coupling_range && fabs(current_coupling[i])>max_coupling_range[i]) {
         max_coupling_range[i]*=coupling_range_increase_factor; 
         max_coupling_grad[i]*=coupling_range_increase_factor; 
      }
  }

  //reduce all the flags 
  if(equilibration && b_finished_equil_flag) {
    equilibration = false;
    update_calls = 0;
  }

  //Now we update coupling constant, if necessary
  if(!equilibration && update_period > 0 && update_calls == update_period && !freeze) {
    double step_size = 0;
    double tmp;
    for(unsigned i=0;i<ncvs;++i){
       //calulcate step size
       //uses scale here, which by default is center
       tmp = 2. * (means[i]/scale[i] - 1) * ssds[i] / (update_calls - 1);
       step_size = tmp / kbt;

       //check if the step_size exceeds maximum possible gradient
       step_size = copysign(fmin(fabs(step_size), max_coupling_grad[i]), step_size);

       //reset means/vars
       means[i] = 0;
       ssds[i] = 0;

      //multidimesional stochastic step
      if(ncvs == 1 || (rand.RandU01() < (1. / ncvs) ) ) {
	    coupling_accum[i] += step_size * step_size;
            //equation 5 in White and Voth, JCTC 2014
	    //no negative sign because it's in step_size
            set_coupling[i] += max_coupling_range[i]/sqrt(coupling_accum[i])*step_size;
            coupling_rate[i] = (set_coupling[i]-current_coupling[i])/update_period;
//            coupling_rate[i] = copysign( fmin(fabs(coupling_rate[i]),
//                                            max_coupling_rate[i])
//                                          ,coupling_rate[i]);
        } else {
            //do not change the bias
            coupling_rate[i] = 0;
        }


   }

    update_calls = 0;
    avg_coupling_count++;
    equilibration = true; //back to equilibration now
  } //close update if

  //pass couplings out so they are accessible
  for(unsigned i=0;i<ncvs;++i){
      outCoupling[i]->set(current_coupling[i]);
  }
}

void EDS::update(){
}

EDS::~EDS(){
    orestartfile_.close();
}

void EDS::turnOnDerivatives(){
  // do nothing
  // this is to avoid errors triggered when a bias is used as a CV
  // (This is done in ExtendedLagrangian.cpp)
}


}


}
