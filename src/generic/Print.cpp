/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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
#include "core/ActionWithArguments.h"
#include "core/ActionWithValue.h"
#include "core/ActionAtomistic.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/AverageBase.h"
//#include "core/SetupMolInfo.h"
#include "core/ActionSet.h"
#include "core/ActionSetup.h"
#include "tools/h36.h"

namespace PLMD {
namespace generic {

//+PLUMEDOC PRINTANALYSIS PRINT
/*
Print quantities to a file.

This directive can be used multiple times
in the input so you can print files with different strides or print different quantities
to different files.  You can control the buffering of output using the \subpage FLUSH keyword.
Output file is either appended or backed up depending on the presence of the \ref RESTART action.
A per-action `RESTART` keyword can be used as well.

Notice that printing happens in the so-called "update" phase. This implies that printing
is affected by the presence of \ref UPDATE_IF actions. In addition, one might decide to start
and stop printing at preassigned values of time using the `UPDATE_FROM` and `UPDATE_UNTIL` keywords.
Keep into account that even on steps when the action is not updated (and thus the file is not printed)
the argument will be activated. In other words, if you use `UPDATE_FROM` to start printing at a given time,
the collective variables this PRINT statement depends on will be computed also before that time.

\par Examples

The following input instructs plumed to print the distance between atoms 3 and 5 on a file
called COLVAR every 10 steps, and the distance and total energy on a file called COLVAR_ALL
every 1000 steps.
\plumedfile
# compute distance:
distance: DISTANCE ATOMS=2,5
# compute total energy (potential)
energy: ENERGY
# print distance on a file
PRINT ARG=distance          STRIDE=10   FILE=COLVAR
# print both variables on another file
PRINT ARG=distance,energy   STRIDE=1000 FILE=COLVAR_ALL
\endplumedfile

Notice that \ref DISTANCE and \ref ENERGY are computed respectively every 10 and 1000 steps, that is
only when required.

*/
//+ENDPLUMEDOC

class Print :
  public ActionPilot,
  public ActionWithArguments,
  public ActionAtomistic
{
  std::string tstyle;
  std::string file;
  OFile ofile;
  std::string description;
  std::vector<std::string> pdb_arg_names;
  std::vector<unsigned> pdb_atom_indices;
  std::string fmt;
  bool hasorigin;
  double lenunit;
  std::vector<std::string> names;
  bool gridinput;
// small internal utility
/////////////////////////////////////////
// these are crazy things just for debug:
// they allow to change regularly the
// printed argument
  int rotate;
  int rotateCountdown;
  int rotateLast;
  std::vector<Value*> rotateArguments;
  std::vector<double> lower, upper;
/////////////////////////////////////////
  bool timeseries;
  double dot_connection_cutoff;
  bool isInTargetRange( const std::vector<double>& argvals ) const ;
private:
  void printAtom( OFile& opdbf, const unsigned& anum, const Vector& pos, const double& m, const double& q ) const ;
public:
  void calculate() override {}
  void prepare() override;
  explicit Print(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
  void apply() override {}
  void update() override;
  void runFinalJobs() override;
  void unlockRequests() override { ActionWithArguments::unlockRequests(); ActionAtomistic::unlockRequests(); }
  void lockRequests() override { ActionWithArguments::lockRequests(); ActionAtomistic::lockRequests(); }
  void calculateNumericalDerivatives( ActionWithValue* a=NULL ) override { plumed_error(); }
  ~Print();
};

PLUMED_REGISTER_ACTION(Print,"PRINT")

void Print::registerKeywords(Keywords& keys) {
  Action::registerKeywords(keys);
  ActionPilot::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys);
  ActionAtomistic::registerKeywords(keys);
  keys.use("ARG");
  keys.add("atoms","ATOMS","the atoms that you would like to you output - only required if using xyz");
  keys.add("compulsory","UNITS","PLUMED","the length units you would like to use when outputting atoms in you xyz file");
  keys.add("compulsory","STRIDE","0","the frequency with which the quantities of interest should be output");
  keys.add("compulsory","CONNECTION_TOL","epsilson","if value of matrix element between i and j is greater than this value they are not connected");
  keys.add("optional","DESCRIPTION","the title to use for your PDB output");
  keys.add("optional","ATOM_INDICES","the indices of the atoms in your PDB output");
  keys.add("optional","ARG_NAMES","the names to use for the arguments in the output PDB");
  keys.add("optional","FILE","the name of the file on which to output these quantities");
  keys.add("optional","FMT","the format that should be used to output real numbers");
  keys.add("atoms","ORIGIN","You can use this keyword to specify the position of an atom as an origin. The positions output will then be displayed relative to that origin");
  keys.add("optional","LESS_THAN_OR_EQUAL","when printing with arguments that are vectors only print components of vectors have a value less than or equal to this value");
  keys.add("optional","GREATER_THAN_OR_EQUAL","when printing with arguments that are vectors only print components of vectors have a value greater than or equal to this value");
  keys.add("hidden","_ROTATE","some funky thing implemented by GBussi");
  keys.use("RESTART");
  keys.use("UPDATE_FROM");
  keys.use("UPDATE_UNTIL");
}

Print::Print(const ActionOptions&ao):
  Action(ao),
  ActionPilot(ao),
  ActionWithArguments(ao),
  ActionAtomistic(ao),
  tstyle("colvar"),
  description(""),
  fmt("%f"),
  hasorigin(false),
  lenunit(1.0),
  gridinput(false),
  rotate(0),
  timeseries(false),
  dot_connection_cutoff(0.)
{
  parse("FILE",file);
  // This checks if we are printing a stored time series
  if( getNumberOfArguments()>0 ) {
      for(unsigned i=0; i<getNumberOfArguments(); ++i) {
          if( getPntrToArgument(i)->getRank()==1 && !getPntrToArgument(i)->hasDerivatives() ) { 
              if( getPntrToArgument(i)->isTimeSeries() ) { timeseries=true; break; }
              ActionSetup* as = dynamic_cast<ActionSetup*>( getPntrToArgument(i)->getPntrToAction() );
              if(as) { timeseries=true; break; }
          }
      }
      if( timeseries ) {
          std::vector<Value*> myvals;
          for(unsigned i=0; i<getNumberOfArguments(); ++i) {
              if( getPntrToArgument(i)->getRank()==1 && !getPntrToArgument(i)->hasDerivatives() ) myvals.push_back( getPntrToArgument(i) );
          }
          unsigned nv = myvals[0]->getNumberOfValues();
          for(unsigned i=0; i<myvals.size();++i) {
              if( myvals[i]->getNumberOfValues()!=nv ) error("for printing of time series all arguments must have same number of values");
          }
          log.printf("  input is printed as time series\n"); requestArguments( myvals, false ); 
      }
  }
  if(file.length()>0) {
    tstyle = Tools::extension( file ); 
    if( tstyle!="xyz" && tstyle!="ndx" && tstyle!="grid" && tstyle!="cube" && tstyle!="dot" && tstyle!="matx" && tstyle!="pdb" ) tstyle="colvar";
    log.printf("  on file %s\n",file.c_str());
    if( !timeseries && tstyle!="grid" && tstyle!="cube" && tstyle!="pdb" && tstyle!="matx" ) { ofile.link(*this); ofile.open(file); }
  } else {
    log.printf("  on plumed log file\n");
    ofile.link(log);
  }
  parse("FMT",fmt);
  if( tstyle=="cube" ) fmt=fmt+" ";
  else fmt=" "+fmt;
  log.printf("  with format %s\n",fmt.c_str());
  if( tstyle=="colvar" ) {
    for(unsigned i=0; i<getNumberOfArguments(); ++i) { 
      if( !timeseries ) getPntrToArgument(i)->buildDataStore( getLabel() );
      if( getPntrToArgument(i)->isPeriodic() ) {
        ofile.addConstantField("min_" + getPntrToArgument(i)->getName() );
        ofile.addConstantField("max_" + getPntrToArgument(i)->getName() );
      }
    }
    if( getStride()==0 ) { 
        if( timeseries ) { setStride(0); log.printf("  printing time series at end of calculation \n"); }
        else { setStride(1); log.printf("  with stride %d\n",getStride()); }
    }
/////////////////////////////////////////
// these are crazy things just for debug:
// they allow to change regularly the
// printed argument
    parse("_ROTATE",rotate);
    if(rotate>0) {
      rotateCountdown=rotate;
      for(unsigned i=0; i<getNumberOfArguments(); ++i) rotateArguments.push_back( getPntrToArgument(i) );
      std::vector<Value*> a(1,rotateArguments[0]);
      requestArguments(std::vector<Value*>(1,rotateArguments[0]),false);
      rotateLast=0;
    }
  } else if( tstyle=="xyz" || tstyle=="ndx" ) {
    unsigned nper=getNumberOfArguments();
    for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      if( getPntrToArgument(i)->getRank()>0 && getPntrToArgument(i)->hasDerivatives() ) { gridinput=true; break; }
      if( getPntrToArgument(i)->getRank()!=1 ) error("problem outputting " + getPntrToArgument(i)->getName() + " can only output vectors in xyz/ndx output" );
      unsigned nt = getPntrToArgument(i)->getNumberOfValues();;
      if( i==0 ) { nper=nt; }
      else if( nt!=nper ) error("mismatched number of values in matrices input in input");
    }
    if( gridinput ) {
      if( getStride()==0 ) {
        setStride(0); log.printf("  printing final grid only \n");
      }
      if( tstyle=="ndx" ) error("grids should be printed to xyz, grid or cube files only");
      if( getNumberOfArguments()!=1 ) error("can only print one grid at a time");
      log.printf("  converting input grid to a set of coordinates and printing \n");
      std::string unitname; parse("UNITS",unitname);
      if(unitname!="PLUMED") {
        Units myunit; myunit.setLength(unitname);
        lenunit=plumed.getAtoms().getUnits().getLength()/myunit.getLength();
      }
    } else {
      if( getStride()==0 ) { setStride(1); log.printf("  with stride %d\n",getStride()); }
      std::vector<std::string> str_upper, str_lower;
      parseVector("LESS_THAN_OR_EQUAL",str_upper); parseVector("GREATER_THAN_OR_EQUAL",str_lower);
      if( str_upper.size()!=getNumberOfArguments() && str_upper.size()>0 ) error("wrong number of arguments for LESS_THAN_OR_EQUAL keyword");
      if( str_lower.size()!=getNumberOfArguments() && str_lower.size()>0 ) error("wrong number of arguments for GREATER_THAN_OR_EQUAL keyword");
      if( str_upper.size()>0 && str_lower.size()>0 ) {
        lower.resize( str_lower.size() ); upper.resize( str_upper.size() );
        for(unsigned i=0; i<upper.size(); ++i) {
          if( str_lower[i]=="none" ) lower[i] = -std::numeric_limits<double>::max();
          else Tools::convert( str_lower[i], lower[i] );
          if( str_upper[i]=="none" ) upper[i] = std::numeric_limits<double>::max();
          else Tools::convert( str_upper[i], upper[i] );
        }
        log.printf("  only printing positions/indices of atoms that have %f <= %s <= %f ", lower[0], getPntrToArgument(0)->getName().c_str(), upper[0] );
        for(unsigned i=1; i<upper.size(); ++i) log.printf("and %f <= %s <= %f ", lower[i], getPntrToArgument(i)->getName().c_str(), upper[i] );
        log.printf("\n");
      } else if( str_upper.size()>0 ) {
        upper.resize( str_upper.size() );
        for(unsigned i=0; i<upper.size(); ++i) {
          if( str_upper[i]=="none" ) upper[i] = std::numeric_limits<double>::max();
          else Tools::convert( str_upper[i], upper[i] );
        }
        log.printf("  only printing positions/indices of atoms that have %s <= %f ", getPntrToArgument(0)->getName().c_str(), upper[0] );
        for(unsigned i=1; i<upper.size(); ++i) log.printf("and %s <= %f ", getPntrToArgument(i)->getName().c_str(), upper[i] );
        log.printf("\n");
      } else if( str_lower.size()>0 ) {
        lower.resize( str_lower.size() );
        for(unsigned i=0; i<lower.size(); ++i) {
          if( str_lower[i]=="none" ) lower[i] = -std::numeric_limits<double>::max();
          else Tools::convert( str_lower[i], lower[i] );
        }
        log.printf("  only printing positions/indices of atoms that have %f <= %s ", lower[0], getPntrToArgument(0)->getName().c_str()  );
        for(unsigned i=1; i<upper.size(); ++i) log.printf("and %f <= %s ", lower[i], getPntrToArgument(i)->getName().c_str() );
        log.printf("\n");
      }

      std::vector<AtomNumber> atoms; parseAtomList("ATOMS",atoms);
      if( atoms.size()!=0 && atoms.size()!=nper ) error("number of atoms should match number of colvars");
      std::vector<AtomNumber> origin; parseAtomList("ORIGIN",origin);
      if( origin.size()==1 ) {
        hasorigin=true; atoms.push_back( origin[0] ); 
        log.printf("  printing atom positions relative to atom %d \n", origin[0].serial() );
      } else if( origin.size()>0 ) error("should only specify one atom for origin");

      if( tstyle=="xyz" ) {
        std::string unitname; parse("UNITS",unitname);
        if(unitname!="PLUMED") {
          Units myunit; myunit.setLength(unitname);
          lenunit=plumed.getAtoms().getUnits().getLength()/myunit.getLength();
        }
        log.printf("  printing xyz file containing poisitions of atoms in columns 1, 2 and 3\n");
        for(unsigned i=0; i<getNumberOfArguments(); ++i) {
          log.printf("  column %d contains components of vector %s \n", 4+i, getPntrToArgument(i)->getName().c_str() );
        }
        log.printf("  atom positions printed are : ");
      } else if( tstyle=="ndx" ) {
        log.printf("  printing ndx file containing indices of atoms that have symmetry functions in ranges prescribed above \n");
        log.printf("  full set of atom indices investigated are : ");
      }
      for(unsigned int i=0; i<atoms.size(); ++i) {
        if ( (i+1) % 25 == 0 ) log.printf("  \n");
        log.printf("  %d", atoms[i].serial());
      }
      log.printf("\n"); 
      std::vector<Value*> args( getArguments() ); requestAtoms( atoms ); requestArguments( args, false );
    }
  } else if( tstyle=="grid" ) {
    if( getStride()==0 ) {
      setStride(0); log.printf("  printing final grid only \n");
    }
    if( getNumberOfArguments()!=1 ) error("when printing a grid you should only have one argument in input");
    if( getPntrToArgument(0)->getRank()==0 || !getPntrToArgument(0)->hasDerivatives() ) error("input argument is not a grid");
    log.printf("  printing function labelled %s at points on a grid in a PLUMED grid file \n", getPntrToArgument(0)->getName().c_str() );
  } else if( tstyle=="cube" ) {
    if( getStride()==0 ) {
      setStride(0); log.printf("  printing final grid only \n");
    }
    if( getNumberOfArguments()!=1 ) error("when printing a grid you should only have one argument in input");
    if( getPntrToArgument(0)->getRank()!=3 || !getPntrToArgument(0)->hasDerivatives() ) error("input argument is not a 3D grid");
    log.printf("  printing function labelled %s at points on a grid in a cube file \n", getPntrToArgument(0)->getName().c_str() );
  } else if( tstyle=="matx" ) {
    if( getNumberOfArguments()!=1 ) error("when printing a matrix to do a file you should only have one argument in input");
    if( getPntrToArgument(0)->getRank()!=2 || getPntrToArgument(0)->hasDerivatives() ) error("input argument is not a matrix");
    if( getStride()==0 ) {
      setStride(0); log.printf("  printing final matrix only \n");
    }
    log.printf("  printing matrix labelled %s to file \n", getPntrToArgument(0)->getName().c_str() );
  } else if( tstyle=="dot" ) {
    if( getNumberOfArguments()!=1 ) error("when printing a matrix to do a dot file you should only have one argument in input");
    if( getPntrToArgument(0)->getRank()!=2 || getPntrToArgument(0)->hasDerivatives() ) error("input argument is not a matrix");
    if( getPntrToArgument(0)->getShape()[0]!=getPntrToArgument(0)->getShape()[1] ) error("should not print non square matrices to dot file");
    if( getStride()==0 ) {
      setStride(0); log.printf("  printing final matrix only \n");
    }
    log.printf("  printing matrix labelled %s to a dot file \n", getPntrToArgument(0)->getName().c_str() );
    std::string ctol; parse("CONNECTION_TOL",ctol);
    if( ctol=="epsilon" ) dot_connection_cutoff = epsilon;
    else Tools::convert( ctol, dot_connection_cutoff );
    log.printf("  elements in graph are shown connected if matrix element is greater than %f \n", dot_connection_cutoff );
  } else if( tstyle=="pdb" ) {
    if( getStride()==0 ) {
      setStride(0); log.printf("  printing pdb at end of calculation \n");
    }
    parseVector("ATOM_INDICES",pdb_atom_indices); parseVector("ARG_NAMES",pdb_arg_names);
    parse("DESCRIPTION",description); log.printf("  printing configurations to a pdb file \n");
    if( getPntrToArgument(0)->getRank()==0 || getPntrToArgument(0)->hasDerivatives() ) error("argument for printing of PDB should be vector or matrix");
    unsigned nv=getPntrToArgument(0)->getShape()[0]; 
    for(unsigned i=0; i<getNumberOfArguments(); ++i) {
        if( getPntrToArgument(i)->getRank()==0 || getPntrToArgument(i)->hasDerivatives() ) error("argument for printing of PDB should be vector or matrix"); 
        if( getPntrToArgument(i)->getShape()[0]!=nv ) error("for printing to pdb files all arguments must have same number of values");
    }
  } else {
    error("expected output does not exist");
  }
/////////////////////////////////////////
  checkRead();
}

void Print::prepare() {
/////////////////////////////////////////
// these are crazy things just for debug:
// they allow to change regularly the
// printed argument
  if(rotate>0) {
    rotateCountdown--;
    if(rotateCountdown==0) {
      rotateCountdown=rotate;
      rotateLast++;
      rotateLast%=rotateArguments.size();
      requestArguments(std::vector<Value*>(1,rotateArguments[rotateLast]), false);
    }
  }
/////////////////////////////////////////
}

bool Print::isInTargetRange( const std::vector<double>& argvals ) const {
  bool printthis=true;
  for(unsigned j=0; j<argvals.size(); ++j) {
    if( upper.size()>0 ) {
      if( argvals[j]>upper[j] ) { printthis=false; break; }
    }
    if( lower.size()>0 ) {
      if( argvals[j]<lower[j] ) { printthis=false; break; }
    }
  }
  return printthis;
}

void Print::update() {
  if( getStep()==0 ) {
    bool dontprint=getNumberOfArguments()>0;
    for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      AverageBase* av = dynamic_cast<AverageBase*>( getPntrToArgument(i)->getPntrToAction() );
      if( !av ) { dontprint=false; break; }
    }
    if( dontprint ) return;  // If everything is an average don't print on first step
  }
  plumed_assert( getStride()>0 );

  if( !timeseries && tstyle=="colvar" ) {
    ofile.fmtField(" %f");
    ofile.printField("time",getTime());
    if( getNumberOfArguments()>0 ) {
        if( getNumberOfArguments()>1 || getPntrToArgument(0)->getRank()==0 ) {
          for(unsigned i=0; i<getNumberOfArguments(); i++) {
            ofile.fmtField(fmt); getPntrToArgument(i)->print( getLabel(), ofile );
          }
        } else {
          for(unsigned i=0; i<getNumberOfArguments(); i++) { ofile.fmtField(fmt); getPntrToArgument(i)->print( getLabel(), ofile ); }
        }
    }
    ofile.printField();
  } else if( tstyle=="colvar" ) {
    OFile ogfile; ogfile.link(*this); 
    ogfile.setBackupString("analysis"); ogfile.open( file );
    unsigned nv = getPntrToArgument(0)->getNumberOfValues();
    std::vector<std::string> arg_names( getNumberOfArguments() );
    for(unsigned j=0;j<getNumberOfArguments();++j) {
        arg_names[j] = getPntrToArgument(j)->getName();
        AverageBase* myav = dynamic_cast<AverageBase*>( getPntrToArgument(j)->getPntrToAction() );
        if( myav ) { std::size_t dot=arg_names[j].find_first_of("."); arg_names[j] = arg_names[j].substr(dot+1); }
    }
    std::vector<double> time(1), time2(1); std::vector<unsigned> ind(1); ind[0]=nv;
    for(unsigned i=0; i<nv; ++i) {
        if( getPntrToArgument(0)->isConstant() ) time[0]=i;
        else (getPntrToArgument(0)->getPntrToAction())->getGridPointIndicesAndCoordinates( i, ind, time ); 
        ogfile.fmtField(" %f"); ogfile.printField("time",time[0]);  
        for(unsigned j=0;j<getNumberOfArguments();++j) {
            if( !getPntrToArgument(0)->isConstant() ) { 
                (getPntrToArgument(0)->getPntrToAction())->getGridPointIndicesAndCoordinates( i, ind, time2 );
                if( time[0]!=time2[0] ) error("mismatched times in printed time series"); 
            }
            ogfile.fmtField(fmt); 
            if( getPntrToArgument(j)->isPeriodic() ) { 
                std::string str_min, str_max; getPntrToArgument(j)->getDomain( str_min, str_max );
                ogfile.printField( "min_" + arg_names[j], str_min ); ogfile.printField("max_" + arg_names[j], str_max ); 
            }
            ogfile.printField( arg_names[j], getPntrToArgument(j)->get(i) );
        }
        ogfile.printField();
    }
    ogfile.close();
  } else if( tstyle=="xyz") {
    if( getNumberOfAtoms()>0 ) {
      unsigned natoms=0, ntatoms=getNumberOfAtoms(); if( hasorigin ) ntatoms = ntatoms - 1;
      std::vector<double> argvals( getNumberOfArguments() );
      for(unsigned i=0; i<ntatoms; ++i) {
        for(unsigned j=0;j<argvals.size();++j) argvals[j] = getPntrToArgument(j)->get(i); 
        if( isInTargetRange( argvals ) ) natoms++;
      }
      ofile.printf("%d\n",natoms);
      const Tensor & t(getPbc().getBox());
      if(getPbc().isOrthorombic()) {
        ofile.printf((" "+fmt+" "+fmt+" "+fmt+"\n").c_str(),lenunit*t(0,0),lenunit*t(1,1),lenunit*t(2,2));
      } else {
        ofile.printf((" "+fmt+" "+fmt+" "+fmt+" "+fmt+" "+fmt+" "+fmt+" "+fmt+" "+fmt+" "+fmt+"\n").c_str(),
                     lenunit*t(0,0),lenunit*t(0,1),lenunit*t(0,2),
                     lenunit*t(1,0),lenunit*t(1,1),lenunit*t(1,2),
                     lenunit*t(2,0),lenunit*t(2,1),lenunit*t(2,2)
                    );
      }
      for(unsigned i=0; i<ntatoms; ++i) {
        const char* defname="X"; const char* name=defname;
        if(names.size()>0) if(names[i].length()>0) name=names[i].c_str();
        for(unsigned j=0;j<argvals.size();++j) argvals[j] = getPntrToArgument(j)->get(i);
        if( isInTargetRange( argvals ) ) {
          if( hasorigin ) {
            Vector fpos=pbcDistance( getPosition(ntatoms), getPosition(i) );
            ofile.printf(("%s "+fmt+" "+fmt+" "+fmt).c_str(),name,lenunit*fpos[0],lenunit*fpos[1],lenunit*fpos[2]);
          } else {
            ofile.printf(("%s "+fmt+" "+fmt+" "+fmt).c_str(),name,lenunit*getPosition(i)[0],lenunit*getPosition(i)[1],lenunit*getPosition(i)[2]);
          }
          for(unsigned j=0; j<argvals.size(); ++j) ofile.printf((" " + fmt).c_str(), argvals[j] );
          ofile.printf("\n");
        }
      }
    } else if( gridinput ) {
      unsigned ngrid = getPntrToArgument(0)->getNumberOfValues();
      ActionWithValue* myaction = getPntrToArgument(0)->getPntrToAction();
      ofile.printf("%d\n",ngrid); ofile.printf("\n"); std::vector<double> pos;
      for(unsigned i=0; i<ngrid; ++i) {
        const char* defname="X"; const char* name=defname; ofile.printf("%s", name);
        myaction->getGridPointAsCoordinate( i, true, pos );
        for(unsigned j=0; j<pos.size(); ++j) ofile.printf((" " + fmt).c_str(), lenunit*pos[j] );
        ofile.printf("\n");
      }
    } else {
      std::set<AtomNumber> tasks ( getPntrToArgument(0)->getPntrToAction()->getTaskSet() );
      std::vector<double> argvals( getNumberOfArguments() );
      ofile.printf("%d\n",tasks.size()); ofile.printf("\n");
      for(const auto & t : tasks) {
        const char* defname="X"; const char* name=defname; ofile.printf("%s", name);
        for(unsigned j=0;j<argvals.size();++j) argvals[j] = getPntrToArgument(j)->get(t.index()); 
        if( isInTargetRange( argvals ) ) {
          for(unsigned j=0; j<argvals.size(); ++j) ofile.printf((" " + fmt).c_str(), argvals[j] );
          ofile.printf("\n");
        }
      }
    }
  } else if( tstyle=="ndx" ) {
    unsigned n=0; std::vector<double> argvals( getNumberOfArguments() );
    ofile.printf("[ %s step %d ] \n", getLabel().c_str(), getStep() );
    for(unsigned i=0; i<getNumberOfAtoms(); ++i) {
      for(unsigned j=0;j<argvals.size();++j) argvals[j] = getPntrToArgument(j)->get(i); 
      if( isInTargetRange( argvals ) ) {
        ofile.printf("%6d", getAbsoluteIndexes()[i].serial() ); n++;
        if( n%15==0 ) ofile.printf("\n");
      }
    }
    if( n%15!=0 ) ofile.printf("\n");
  } else if( tstyle=="grid" ) {
    OFile ogfile; ogfile.link(*this);
    ogfile.setBackupString("analysis");
    ogfile.open( file ); ogfile.addConstantField("normalisation");
    Value* gval=getPntrToArgument(0); ActionWithValue* act=gval->getPntrToAction();
    std::vector<unsigned> ind( gval->getRank() ), nbin( gval->getRank() ); std::string gtype;
    std::vector<double> spacing( gval->getRank() ), xx( gval->getRank() ); std::vector<bool> pbc( gval->getRank() );
    std::vector<std::string> argn( gval->getRank() ), min( gval->getRank() ), max( gval->getRank() );
    act->getInfoForGridHeader( gtype, argn, min, max, nbin, spacing, pbc, false );
    if( gtype=="fibonacci" ) {
      ogfile.addConstantField("nbins"); 
    } else {
      plumed_assert( gtype=="flat" );
      for(unsigned i=0; i<gval->getRank(); ++i) {
        ogfile.addConstantField("min_" + argn[i] );
        ogfile.addConstantField("max_" + argn[i] );
        ogfile.addConstantField("nbins_" + argn[i] );
        ogfile.addConstantField("periodic_" + argn[i] );
      }
    }

    for(unsigned i=0; i<gval->getNumberOfValues(); ++i) {
      // Retrieve and print the grid coordinates
      act->getGridPointIndicesAndCoordinates( i, ind, xx );
      if(i>0 && gval->getRank()==2 && ind[gval->getRank()-2]==0) ogfile.printf("\n");
      ogfile.fmtField(fmt); ogfile.printField("normalisation", gval->getNorm() );
      if( gtype=="fibonacci" ) {
        ogfile.printField("nbins", static_cast<int>(nbin[0]) );
      } else {
        for(unsigned j=0; j<gval->getRank(); ++j) {
          ogfile.printField("min_" + argn[j], min[j] );
          ogfile.printField("max_" + argn[j], max[j] );
          ogfile.printField("nbins_" + argn[j], static_cast<int>(nbin[j]) );
          if( pbc[j] ) ogfile.printField("periodic_" + argn[j], "true" );
          else         ogfile.printField("periodic_" + argn[j], "false" );
        }
      }
      // Print the grid coordinates
      for(unsigned j=0; j<gval->getRank(); ++j) { ogfile.fmtField(fmt); ogfile.printField(argn[j],xx[j]); }
      // Print value
      ogfile.fmtField(fmt); ogfile.printField( gval->getName(), gval->get(i) );
      // Print the derivatives
      for(unsigned j=0; j<gval->getRank(); ++j) { ogfile.fmtField(fmt); ogfile.printField( "d" + gval->getName() + "_" + argn[j], gval->getGridDerivative(i,j) ); }
      ogfile.printField();
    }
    ogfile.close();
  } else if( tstyle=="cube" ) {
    OFile ogfile; ogfile.link(*this);
    ogfile.setBackupString("analysis");
    ogfile.open( file ); Value* gval=getPntrToArgument(0); ActionWithValue* act=gval->getPntrToAction();
    std::vector<unsigned> nbin( 3 ), pp( 3 );
    std::vector<double> xx( 3 ), spacing( 3 ), extent( 3 ); std::vector<bool> pbc( 3 );
    std::vector<std::string> argn( 3 ), min( 3 ), max( 3 ); std::string gtype;
    act->getInfoForGridHeader( gtype, argn, min, max, nbin, spacing, pbc, true );
    if( gtype=="fibonacci" ) error("cannot print fibonacci grids out to cube files");
    for(unsigned j=0; j<3; ++j) {
      double mind, maxd;
      Tools::convert( min[j], mind );
      Tools::convert( max[j], maxd );
      if( pbc[j] ) extent[j]=maxd-mind;
      else { extent[j]=maxd-mind+spacing[j]; nbin[j]++; }
    }
    ogfile.printf("PLUMED CUBE FILE\n");
    ogfile.printf("OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n");
    // Number of atoms followed by position of origin (origin set so that center of grid is in center of cell)
    std::string ostr = "%d " + fmt + fmt + fmt + "\n";
    ogfile.printf(ostr.c_str(),1,-0.5*extent[0],-0.5*extent[1],-0.5*extent[2] );
    ogfile.printf(ostr.c_str(),nbin[0],spacing[0],0.0,0.0);  // Number of bins in each direction followed by
    ogfile.printf(ostr.c_str(),nbin[1],0.0,spacing[1],0.0);  // shape of voxel
    ogfile.printf(ostr.c_str(),nbin[2],0.0,0.0,spacing[2]);
    ogfile.printf(ostr.c_str(),1,0.0,0.0,0.0); // Fake atom otherwise VMD doesn't work
    for(pp[0]=0; pp[0]<nbin[0]; ++pp[0]) {
      for(pp[1]=0; pp[1]<nbin[1]; ++pp[1]) {
        for(pp[2]=0; pp[2]<nbin[2]; ++pp[2]) {
          unsigned ival=pp[pp.size()-1];
          for(unsigned i=pp.size()-1; i>0; --i) ival=ival*nbin[i-1]+pp[i-1];
          ogfile.printf(fmt.c_str(), gval->get(ival) );
          if(pp[2]%6==5) ogfile.printf("\n");
        }
        ogfile.printf("\n");
      }
    }
    ogfile.close();
  } else if( tstyle=="matx" ) {
    OFile ogfile; ogfile.link(*this);
    ogfile.setBackupString("analysis");
    ogfile.open( file ); Value* gval=getPntrToArgument(0);
    // Now print connections
    unsigned nrows = gval->getShape()[0], ncols = gval->getShape()[1];
    for(unsigned i=0; i<nrows; ++i) {
       for(unsigned j=0; j<ncols; ++j) ogfile.printf( fmt.c_str(), gval->get( i*ncols + j ) );
       ogfile.printf("\n");
    }
    ogfile.close();
  } else if( tstyle=="dot" ) {
    OFile ogfile; ogfile.link(*this);
    ogfile.setBackupString("analysis");
    ogfile.open( file ); Value* gval=getPntrToArgument(0);
    ogfile.printf("graph %s { \n", gval->getName().c_str() );
    // Print all nodes
    for(unsigned i=0; i<gval->getShape()[0]; ++i) ogfile.printf("%u [label=\"%u\"];\n",i,i);
    // Now print connections
    unsigned nrows = gval->getShape()[0];
    for(unsigned i=1; i<nrows; ++i) {
      for(unsigned j=0; j<i; ++j) {
        if( fabs(gval->get(i*nrows+j)-gval->get(j*nrows+i))>epsilon ) error("to print undirected graph matrix should be symmetric");
        if( gval->get(i*nrows+j)>dot_connection_cutoff )  ogfile.printf("%u -- %u \n", i, j );
      }
    }
    ogfile.printf("} \n"); ogfile.close();
  } else if( tstyle=="pdb" ) {
    OFile opdbf; opdbf.link(*this);
    opdbf.setBackupString("analysis");
    opdbf.open( file ); unsigned nn=0; 
    std::size_t psign=fmt.find("%"); plumed_assert( psign!=std::string::npos );  
    std::string descr2="%s=%-" + fmt.substr(psign+1) + " "; 
    // Add suitable code in here to print frames for paths here.  Gareth !!!!!
    std::vector<unsigned> argnums, posnums, matnums; bool use_real_arg_names = pdb_arg_names.size()==0;
    for(unsigned k=0;k<getNumberOfArguments();++k) {
        if( getPntrToArgument(k)->getRank()==2 ) {
             if( matnums.size()>0 ) error("can only output one matrix at a time");
             matnums.push_back(k); 
             if( pdb_atom_indices.size()==0 ) { pdb_atom_indices.resize( getPntrToArgument(k)->getShape()[1] / 3 ); for(unsigned i=0;i<pdb_atom_indices.size();++i) pdb_atom_indices[i]=i; }
             plumed_assert( pdb_atom_indices.size()==getPntrToArgument(k)->getShape()[1] / 3 );
        } else if( getPntrToArgument(k)->getName().find(".pos")!=std::string::npos ) {
             if( posnums.size()%3==0 && getPntrToArgument(k)->getName().find(".posx-")==std::string::npos ) error("x coordinate of input positions in wrong place");
             if( posnums.size()%3==1 && getPntrToArgument(k)->getName().find(".posy-")==std::string::npos ) error("y coordinate of input positions in wrong place");
             if( posnums.size()%3==2 && getPntrToArgument(k)->getName().find(".posz-")==std::string::npos ) error("z coordinate of input positions in wrong place");
             posnums.push_back(k);
        } else {
            if( use_real_arg_names ) pdb_arg_names.push_back( getPntrToArgument(k)->getName() );
            argnums.push_back( k );
        }
    }
    if( posnums.size()%3!=0 ) error("found misleading number of stored positions for output");
    if( pdb_atom_indices.size()==0 ) { pdb_atom_indices.resize( posnums.size() / 3 ); for(unsigned i=0;i<pdb_atom_indices.size();++i) pdb_atom_indices[i]=i; }

    std::string argstr;
    if( argnums.size()>0 ) {
        argstr = "ARG=" + pdb_arg_names[0];
        for(unsigned k=1;k<argnums.size();++k) argstr += "," + pdb_arg_names[k]; 
    }
    if( description.length()>0 ) opdbf.printf("# %s AT STEP %d TIME %f \n", description.c_str(), getStep(), getTime() );
    unsigned nvals = getPntrToArgument(0)->getShape()[0]; Vector pos;
    for(unsigned j=0;j<nvals;++j) {
        if( argnums.size()>0 ) { 
            opdbf.printf("REMARK %s \n", argstr.c_str() ); opdbf.printf("REMARK ");
            for(unsigned k=0;k<argnums.size();++k) {
                Value* thisarg=getPntrToArgument(argnums[k]); opdbf.printf( descr2.c_str(), pdb_arg_names[k].c_str(), thisarg->get(j) ); 
            }
            opdbf.printf("\n"); 
        }
        if( posnums.size()==0 && matnums.size()==1 ) {
            unsigned npos = getPntrToArgument(matnums[0])->getShape()[1] / 3;
            for(unsigned k=0;k<npos;++k) {
                pos[0]=getPntrToArgument(matnums[0])->get(npos*(3*j+0) + k);
                pos[1]=getPntrToArgument(matnums[0])->get(npos*(3*j+1) + k);
                pos[2]=getPntrToArgument(matnums[0])->get(npos*(3*j+2) + k);
                printAtom( opdbf, pdb_atom_indices[k], pos, 1.0, 1.0 );
            }
        } else {
            unsigned npos = posnums.size() / 3;
            for(unsigned k=0;k<npos;++k) {
                pos[0]=getPntrToArgument(posnums[3*k+0])->get(j);
                pos[1]=getPntrToArgument(posnums[3*k+1])->get(j);
                pos[2]=getPntrToArgument(posnums[3*k+2])->get(j);
                printAtom( opdbf, pdb_atom_indices[k], pos, 1.0, 1.0 );
            }
        }
        opdbf.printf("END\n");
    }
    opdbf.close();
  }
}

void Print::runFinalJobs() {
  if( getStride()>0 ) return ;
  // Stride is set to one here otherwise update crashes
  setStride(1); retrieveAtoms(); update();
}

void Print::printAtom( OFile& opdbf, const unsigned& anum, const Vector& pos, const double& m, const double& q ) const {
  std::array<char,6> at; const char* msg = h36::hy36encode(5,anum,&at[0]);
  plumed_assert(msg==nullptr) << msg; at[5]=0; double lunits = atoms.getUnits().getLength()/0.1;
  opdbf.printf("ATOM  %s  X   RES  %4u    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
               &at[0], anum, lunits*pos[0], lunits*pos[1], lunits*pos[2], m, q );
}

Print::~Print() {
}

}


}
