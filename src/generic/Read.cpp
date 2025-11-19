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
#include "core/ActionPilot.h"
#include "core/ActionWithValue.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "tools/IFile.h"
#include <memory>

namespace PLMD {
namespace generic {

//+PLUMEDOC GENERIC READ
/*
Read quantities from a colvar file.

This Action can be used with driver to read in a colvar file that was generated during
an MD simulation. The following example shows how this works.

```plumed
#SETTINGS INPUTFILES=regtest/basic/rt-fametad-1/Input-COLVAR
sum_abs: READ VALUES=sum_abs FILE=regtest/basic/rt-fametad-1/Input-COLVAR IGNORE_FORCES
PRINT ARG=sum_abs STRIDE=1 FILE=colvar
```

The input file `Input-COLVAR` is a colvar file that was generated using a [PRINT](PRINT.md)
command.  The instruction `VALUES=sum_abs` tells PLUMED that we want to read the column headed
`sum_abs` from that file.  The value outputted from a read command is thus always a scalar.

The `IGNORE_FORCES` flag tells PLUMED that any forces that are applied on
this action can be safely ignored.  If you try to run a simulation in which a bias acts upon a
quantity that is read in from a file using a READ command PLUMED will crash with an error.

## Dealing with components

The following example provides another example where the READ command is used:

```plumed
#SETTINGS INPUTFILES=regtest/basic/rt41/input_colvar
r1: READ VALUES=p2.X  FILE=regtest/basic/rt41/input_colvar
r2: READ VALUES=p3.* FILE=regtest/basic/rt41/input_colvar
PRINT ARG=r1.X,r2.* FILE=colvar
```

Notice that the READ command preseves the part of the name that appears after the ".".  Consequently,
when we read the value `p2.X` from the input file the output value the action with label `r1` calls the output
value that contains this information `r1.X`.  The READ command is implemented this way so that you can use the
wildcard syntax that is illustrated in the command labelled `r2` in the above input.

## Reading COLVAR files and trajectories

The example input below indicates a case where a trajectory is being post processed and where some COLVAR files
that were generated when the MD simulation was run are read in.

```plumed
#SETTINGS INPUTFILES=regtest/basic/rt19/input_colvar2
# The distance here is being calculated from the trajectory
d: DISTANCE ATOMS=1,2
# This CV is being read in from a file that was output with the same frequency as frames
# were output from the trajectory. Notice that you can read from zip files
r1: READ VALUES=rmsd0  FILE=regtest/basic/rt19/input_colvar.gz
# This CV is being read in from a file that was output with twice as frequency as frames
r2: READ VALUES=rmsd0  FILE=regtest/basic/rt19/input_colvar2 EVERY=2 IGNORE_TIME
# We start reading this CV from this file after we have read the first 100 ps of the trajectory
# and stop reading from it after we have read the first 800 ps of the trajectory
r3: READ VALUES=rmsd0 FILE=regtest/basic/rt19/input_colvar2 UPDATE_FROM=100 UPDATE_UNTIL=800 IGNORE_TIME
# This outputs our the three quantities that are determined for every step
PRINT ARG=d,r1,r2 FILE=colvar
# This outputs the data that we only have for the 700 ps of the trajectory after the first 100 ps
PRINT ARG=d,r1,r2,r3 FILE=colvar2 UPDATE_FROM=100 UPDATE_UNTIL=800
```

The driver command to run this script as follows:

````
plumed driver --plumed plumed.dat --trajectory-stride 10 --timestep 0.005 --ixyz trajectory.xyz
````

When you are doing analyses like these, which involve mixing using READ command and analysing a trajectory, PLUMED forces you
to take care to ensure that everything in the input file was generated from the same step in the input trajectory.  You must
use the `--trajectory-stride` and `--timestep` commands when using driver above so PLUMED can correctly calculate the simulation
time and compare it with the time stamps in any colvar files that are being read in using READ commands.

If you want to turn off these checks either because you are confident that you have set things up correctly or if you are not
mixing variables that have been calculated from a trajectory with variables that are being read from a file you can use the
`IGNORE_TIME` flag.  Notice also that you can use the `EVERY` flag to tell PLUMED to ignore parts of the COLVAR file if data
has been output to that file more frequently that data has been output to the trajectory.

*/
//+ENDPLUMEDOC

class Read :
  public ActionPilot,
  public ActionWithValue {
private:
  bool ignore_time;
  bool ignore_forces;
  bool cloned_file;
  unsigned nlinesPerStep;
  std::string filename;
/// Unique pointer with the same scope as ifile.
  std::unique_ptr<IFile> ifile_ptr;
/// Pointer to input file.
/// It is either pointing to the content of ifile_ptr
/// or to the file it is cloned from.
  IFile* ifile;
  std::vector<std::unique_ptr<Value>> readvals;
public:
  static void registerKeywords( Keywords& keys );
  explicit Read(const ActionOptions&);
  std::string getOutputComponentDescription( const std::string& cname, const Keywords& keys ) const override ;
  void prepare() override;
  void apply() override {}
  void calculate() override;
  void update() override;
  std::string getFilename() const;
  IFile* getFile();
  unsigned getNumberOfDerivatives() override;
  void turnOnDerivatives() override;
};

PLUMED_REGISTER_ACTION(Read,"READ")

void Read::registerKeywords(Keywords& keys) {
  Action::registerKeywords(keys);
  ActionPilot::registerKeywords(keys);
  ActionWithValue::registerKeywords(keys);
  keys.add("compulsory","STRIDE","1","the frequency with which the file should be read.");
  keys.add("compulsory","EVERY","1","only read every nth line of the colvar file. This should be used if the colvar was written more frequently than the trajectory.");
  keys.add("compulsory","VALUES","the values to read from the file");
  keys.add("compulsory","FILE","the name of the file from which to read these quantities");
  keys.addFlag("IGNORE_TIME",false,"ignore the time in the colvar file. When this flag is not present read will be quite strict "
               "about the start time of the simulation and the stride between frames");
  keys.addFlag("IGNORE_FORCES",false,"use this flag if the forces added by any bias can be safely ignored.  As an example forces can be "
               "safely ignored if you are doing post processing that does not involve outputting forces");
  keys.remove("NUMERICAL_DERIVATIVES");
  keys.use("UPDATE_FROM");
  keys.use("UPDATE_UNTIL");
  ActionWithValue::useCustomisableComponents(keys);
}

Read::Read(const ActionOptions&ao):
  Action(ao),
  ActionPilot(ao),
  ActionWithValue(ao),
  ignore_time(false),
  ignore_forces(false),
  nlinesPerStep(1) {
  // Read the file name from the input line
  parse("FILE",filename);
  // Check if time is to be ignored
  parseFlag("IGNORE_TIME",ignore_time);
  // Check if forces are to be ignored
  parseFlag("IGNORE_FORCES",ignore_forces);
  // Open the file if it is not already opened
  cloned_file=false;
  std::vector<Read*> other_reads=plumed.getActionSet().select<Read*>();
  for(unsigned i=0; i<other_reads.size(); ++i) {
    if( other_reads[i]->getFilename()==filename ) {
      ifile=other_reads[i]->getFile();
      cloned_file=true;
    }
  }
  if( !cloned_file ) {
    ifile_ptr=Tools::make_unique<IFile>();
    ifile=ifile_ptr.get();
    if( !ifile->FileExist(filename) ) {
      error("could not find file named " + filename);
    }
    ifile->link(*this);
    ifile->open(filename);
    ifile->allowIgnoredFields();
  }
  parse("EVERY",nlinesPerStep);
  if(nlinesPerStep>1) {
    log.printf("  only reading every %uth line of file %s\n",nlinesPerStep,filename.c_str() );
  } else {
    log.printf("  reading data from file %s\n",filename.c_str() );
  }
  // Find out what we are reading
  std::vector<std::string> valread;
  parseVector("VALUES",valread);

  if(nlinesPerStep>1 && cloned_file) {
    error("Opening a file multiple times and using EVERY is not allowed");
  }

  std::size_t dot=valread[0].find_first_of('.');
  if( valread[0].find(".")!=std::string::npos ) {
    std::string labelVal=valread[0].substr(0,dot);
    std::string nameVal=valread[0].substr(dot+1);
    if( nameVal=="*" ) {
      if( valread.size()>1 ) {
        error("all values must be from the same Action when using READ");
      }
      std::vector<std::string> fieldnames;
      ifile->scanFieldList( fieldnames );
      for(unsigned i=0; i<fieldnames.size(); ++i) {
        if( fieldnames[i].substr(0,dot)==labelVal ) {
          readvals.emplace_back(Tools::make_unique<Value>(this, fieldnames[i], false) );
          addComponentWithDerivatives( fieldnames[i].substr(dot+1) );
          if( ifile->FieldExist("min_" + fieldnames[i]) ) {
            componentIsPeriodic( fieldnames[i].substr(dot+1), "-pi","pi" );
          } else {
            componentIsNotPeriodic( fieldnames[i].substr(dot+1) );
          }
        }
      }
    } else {
      readvals.emplace_back(Tools::make_unique<Value>(this, valread[0], false) );
      addComponentWithDerivatives( nameVal );
      if( ifile->FieldExist("min_" + valread[0]) ) {
        componentIsPeriodic( valread[0].substr(dot+1), "-pi", "pi" );
      } else {
        componentIsNotPeriodic( valread[0].substr(dot+1) );
      }
      for(unsigned i=1; i<valread.size(); ++i) {
        if( valread[i].substr(0,dot)!=labelVal ) {
          error("all values must be from the same Action when using READ");
        };
        readvals.emplace_back(Tools::make_unique<Value>(this, valread[i], false) );
        addComponentWithDerivatives( valread[i].substr(dot+1) );
        if( ifile->FieldExist("min_" + valread[i]) ) {
          componentIsPeriodic( valread[i].substr(dot+1), "-pi", "pi" );
        } else {
          componentIsNotPeriodic( valread[i].substr(dot+1) );
        }
      }
    }
  } else {
    if( valread.size()!=1 ) {
      error("all values must be from the same Action when using READ");
    }
    readvals.emplace_back(Tools::make_unique<Value>(this, valread[0], false) );
    addValueWithDerivatives();
    if( ifile->FieldExist("min_" + valread[0]) ) {
      setPeriodic( "-pi", "pi" );
    } else {
      setNotPeriodic();
    }
    log.printf("  reading value %s and storing as %s\n",valread[0].c_str(),getLabel().c_str() );
  }
  checkRead();
}

std::string Read::getOutputComponentDescription( const std::string& cname, const Keywords& keys ) const {
  plumed_assert( !exists( getLabel() ) );
  for(unsigned i=0; i<readvals.size(); ++i) {
    if( readvals[i]->getName().find( cname )!=std::string::npos ) {
      return "values from the column labelled " + readvals[i]->getName() + " in the file named " + filename;
    }
  }
  plumed_error();
  return "";
}

std::string Read::getFilename() const {
  return filename;
}

IFile* Read::getFile() {
  return ifile;
}

unsigned Read::getNumberOfDerivatives() {
  return 0;
}

void Read::turnOnDerivatives() {
  if( !ignore_forces )
    error("cannot calculate derivatives for colvars that are read in from a file.  If you are postprocessing and "
          "these forces do not matter add the flag IGNORE_FORCES to all READ actions");
}

void Read::prepare() {
  if( !cloned_file ) {
    double du_time;
    if( !ifile->scanField("time",du_time) ) {
      error("Reached end of file " + filename + " before end of trajectory");
    } else if( std::abs( du_time-getTime() )>getTimeStep() && !ignore_time ) {
      std::string str_dutime,str_ptime;
      Tools::convert(du_time,str_dutime);
      Tools::convert(getTime(),str_ptime);
      error("mismatched times in colvar files : colvar time=" + str_dutime + " plumed time=" + str_ptime + ". Add IGNORE_TIME to ignore error.");
    }
  }
}

void Read::calculate() {
  std::string smin, smax;
  for(unsigned i=0; i<readvals.size(); ++i) {
// .get  returns the raw pointer
// ->get calls the Value::get() method
    ifile->scanField( readvals[i].get() );
    getPntrToComponent(i)->set( readvals[i]->get() );
    if( readvals[i]->isPeriodic() ) {
      readvals[i]->getDomain( smin, smax );
      getPntrToComponent(i)->setDomain( smin, smax );
    }
  }
}

void Read::update() {
  if( !cloned_file ) {
    for(unsigned i=0; i<nlinesPerStep; ++i) {
      ifile->scanField();
      double du_time;
      if( !ifile->scanField("time",du_time) && !plumed.inputsAreActive() ) {
        plumed.stop();
      }
    }
  }
}

}
}
