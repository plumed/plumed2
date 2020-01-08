/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2020 The plumed team
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
#include "core/Atoms.h"
#include "tools/IFile.h"
#include <memory>

namespace PLMD {
namespace generic {

//+PLUMEDOC GENERIC READ
/*
Read quantities from a colvar file.

This Action can be used with driver to read in a colvar file that was generated during
an MD simulation

\par Description of components

The READ command will read those fields that are labelled with the text string given to the
VALUE keyword.  It will also read in any fields that are labeled with the text string
given to the VALUE keyword followed by a dot and a further string. If a single Value is read in
this value can be referenced using the label of the Action.  Alternatively, if multiple quantities
are read in, they can be referenced elsewhere in the input by using the label for the Action
followed by a dot and the character string that appeared after the dot in the title of the field.

\par Examples

This input reads in data from a file called input_colvar.data that was generated
in a calculation that involved PLUMED.  The first command reads in the data from the
column headed phi1 while the second reads in the data from the column headed phi2.

\plumedfile
rphi1:       READ FILE=input_colvar.data  VALUES=phi1
rphi2:       READ FILE=input_colvar.data  VALUES=phi2
PRINT ARG=rphi1,rphi2 STRIDE=500  FILE=output_colvar.data
\endplumedfile

The file input_colvar.data is just a normal colvar file as shown below

\auxfile{input_colvar.data}
#! FIELDS time phi psi metad.bias metad.rbias metad.rct
#! SET min_phi -pi
#! SET max_phi pi
#! SET min_psi -pi
#! SET max_psi pi
 0.000000  -1.2379   0.8942   0.0000   0.0000   0.0000
 1.000000  -1.4839   1.0482   0.0000   0.0000   0.0089
 2.000000  -1.3243   0.6055   0.0753   0.0664   0.0184
\endauxfile

*/
//+ENDPLUMEDOC

class Read :
  public ActionPilot,
  public ActionWithValue
{
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
  keys.add("compulsory","EVERY","1","only read every \\f$n\\f$th line of the colvar file. This should be used if the colvar was written more frequently than the trajectory.");
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
  nlinesPerStep(1)
{
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
    ifile_ptr.reset(new IFile());
    ifile=ifile_ptr.get();
    if( !ifile->FileExist(filename) ) error("could not find file named " + filename);
    ifile->link(*this);
    ifile->open(filename);
    ifile->allowIgnoredFields();
  }
  parse("EVERY",nlinesPerStep);
  if(nlinesPerStep>1) log.printf("  only reading every %uth line of file %s\n",nlinesPerStep,filename.c_str() );
  else log.printf("  reading data from file %s\n",filename.c_str() );
  // Find out what we are reading
  std::vector<std::string> valread; parseVector("VALUES",valread);

  std::size_t dot=valread[0].find_first_of('.');
  if( valread[0].find(".")!=std::string::npos ) {
    std::string label=valread[0].substr(0,dot);
    std::string name=valread[0].substr(dot+1);
    if( name=="*" ) {
      if( valread.size()>1 ) error("all values must be from the same Action when using READ");
      std::vector<std::string> fieldnames;
      ifile->scanFieldList( fieldnames );
      for(unsigned i=0; i<fieldnames.size(); ++i) {
        if( fieldnames[i].substr(0,dot)==label ) {
          readvals.emplace_back(new Value(this, fieldnames[i], false) ); addComponentWithDerivatives( fieldnames[i].substr(dot+1) );
          if( ifile->FieldExist("min_" + fieldnames[i]) ) componentIsPeriodic( fieldnames[i].substr(dot+1), "-pi","pi" );
          else componentIsNotPeriodic( fieldnames[i].substr(dot+1) );
        }
      }
    } else {
      readvals.emplace_back(new Value(this, valread[0], false) ); addComponentWithDerivatives( name );
      if( ifile->FieldExist("min_" + valread[0]) ) componentIsPeriodic( valread[0].substr(dot+1), "-pi", "pi" );
      else componentIsNotPeriodic( valread[0].substr(dot+1) );
      for(unsigned i=1; i<valread.size(); ++i) {
        if( valread[i].substr(0,dot)!=label ) error("all values must be from the same Action when using READ");;
        readvals.emplace_back(new Value(this, valread[i], false) ); addComponentWithDerivatives( valread[i].substr(dot+1) );
        if( ifile->FieldExist("min_" + valread[i]) ) componentIsPeriodic( valread[i].substr(dot+1), "-pi", "pi" );
        else componentIsNotPeriodic( valread[i].substr(dot+1) );
      }
    }
  } else {
    if( valread.size()!=1 ) error("all values must be from the same Action when using READ");
    readvals.emplace_back(new Value(this, valread[0], false) ); addValueWithDerivatives();
    if( ifile->FieldExist("min_" + valread[0]) ) setPeriodic( "-pi", "pi" );
    else setNotPeriodic();
    log.printf("  reading value %s and storing as %s\n",valread[0].c_str(),getLabel().c_str() );
  }
  checkRead();
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
  if( !ignore_forces ) error("cannot calculate derivatives for colvars that are read in from a file.  If you are postprocessing and "
                               "these forces do not matter add the flag IGNORE_FORCES to all READ actions");
}

void Read::prepare() {
  if( !cloned_file ) {
    double du_time;
    if( !ifile->scanField("time",du_time) ) {
      error("Reached end of file " + filename + " before end of trajectory");
    } else if( fabs( du_time-getTime() )>plumed.getAtoms().getTimeStep() && !ignore_time ) {
      std::string str_dutime,str_ptime; Tools::convert(du_time,str_dutime); Tools::convert(getTime(),str_ptime);
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
      ifile->scanField(); double du_time;
      if( plumed.getAtoms().getNatoms()==0 && !ifile->scanField("time",du_time) ) plumed.stop();
    }
  }
}

}
}
