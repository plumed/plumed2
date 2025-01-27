// In order to correctly catch the thrown C++11 exceptions,
// we notify the Plumed wrapper that those exceptions are recognized by the compiler.
#define __PLUMED_WRAPPER_LIBCXX17 1
#include "plumed/tools/Stopwatch.h"
#include "plumed/wrapper/Plumed.h"
#include <fstream>
#include <iostream>
#include <vector>

using namespace PLMD;

void test_line(std::ostream & ofs,Plumed & p,const std::string & arg){
  std::string cmd="readInputLine";
  ofs<<cmd<<" "<<arg<<std::endl;
  try{
    p.cmd(cmd.c_str(),arg.c_str());
    ofs<<"+++ !!!! uncatched !!!!"<<std::endl;
  } catch(Plumed::Exception&e) {
    ofs<<"+++ catched"<<std::endl;
  }
}

void test_this(std::ostream & ofs,Plumed & p,const std::string & cmd,const void*arg){
  ofs<<cmd;
  if(!arg) ofs<<" NULL";
  ofs<<std::endl;
  try{
    p.cmd(cmd.c_str(),arg);
    ofs<<"+++ !!!! uncatched !!!!"<<std::endl;
  } catch(Plumed::Exception&e) {
    ofs<<"+++ catched"<<std::endl;
  }
}

int main(){

  std::ofstream ofs("output");

  {
// test a mistake in timer
    Stopwatch sw;
    ofs<<"pause"<<std::endl;
    try{
      sw.pause();
      ofs<<"+++ !!!! uncatched !!!!"<<std::endl;
// this is not of type PLMD::Plumed::Exception since it is thrown within the library
    } catch(std::exception & e) {
      ofs<<"+++ catched"<<std::endl;
    }
  }

  Plumed plumed;

// first try to wrongly set real precision
  unsigned i=127;
  test_this(ofs,plumed,"setRealPrecision",&i);

  int natoms=10;

  std::vector<double> positions(3*natoms,0.0);
  for(int i=0;i<3*natoms;i++) positions[i]=i;
  std::vector<double> masses(natoms,1.0);
  std::vector<double> forces(3*natoms,0.0);
  std::vector<double> box(9,0.0);
  std::vector<double> virial(9,0.0);

  plumed.cmd("setNatoms",&natoms);
  plumed.cmd("setLogFile","test.log");
  plumed.cmd("init");

// I try many mistaken lines.
// Each of them will raise an exception
// Notice that name "d" will not be reserved and it will be possible
// to use it later
  test_line(ofs,plumed,"d: DISTANCE ATOMS=1,2,3");
  test_line(ofs,plumed,"d:DISTANCE ATOMS=1,2");
  test_line(ofs,plumed,"d: DIST ANCE ATOMS=1,2");
  test_line(ofs,plumed,"d: DISTANCE ATOMS=1,2 COMPONENTS SCALED_COMPONENTS");
  test_line(ofs,plumed,"d: GYRATION ATOMS=");
  test_line(ofs,plumed,"d: GYRATION ATOMS=1-4 TYPE=WHAT");
  test_line(ofs,plumed,"d: POSITION ATOM=1,2");
  test_line(ofs,plumed,"d: PUCKERING ATOMS=1-4");
  test_line(ofs,plumed,"d: ANGLE ATOMS=1,2,3,4,5");
  test_line(ofs,plumed,"d: TORSION ATOMS=1,2,3,4,5");
  test_line(ofs,plumed,"d: TORSION ATOMS=1,2,3,4 VECTOR1=1,2 VECTOR2=2,3 AXIS=3,4");
  test_line(ofs,plumed,"d: COORDINATION GROUPA=1 GROUPB=2 R_0=0.5 NN=1.5");
  test_line(ofs,plumed,"d: COORDINATION GROUPA=1 GROUPB=2 R_0=0.5 NN=1.5");
  test_line(ofs,plumed,"d: CENTER ATOMS=2,3,4,5 WEIGHTS=1,2,3");
  test_line(ofs,plumed,"d: CENTER ATOMS=2,3,4,5 MASS WEIGHTS=1,2,3,4");
  test_line(ofs,plumed,"LOAD FILE=nonexist.cpp");
  test_line(ofs,plumed,"d: COORDINATION GROUPA=1 GROUPB=2 SWITCH={WRONGNAME R_0=1.0}");
  test_line(ofs,plumed,"d: RMSD REFERENCE=missing.pdb");
  test_line(ofs,plumed,"d: RMSD REFERENCE=test-too-many-atoms.pdb");
  test_line(ofs,plumed,"d: DRMSD REFERENCE=missing.pdb LOWER_CUTOFF=0.0 UPPER_CUTOFF=15.0");

// these should not fail
  plumed.cmd("readInputLine","d: DISTANCE ATOMS=1,2");
  plumed.cmd("readInputLine","d1: DISTANCE ATOMS={1 2}"); // check if braces are parsed correctly
  plumed.cmd("readInputLine","t: TORSION ATOMS=1,2,3,4");
  plumed.cmd("readInputLine","RESTRAINT ARG=d AT=0 KAPPA=1");

// check regexp not surrounded by parentheses (bug in 2.3.3).
  test_line(ofs,plumed,"RESTRAINT ARG=x(d) KAPPA=5 AT=0");
  test_line(ofs,plumed,"RESTRAINT ARG=(d)x KAPPA=5 AT=0");

// check error in regular expression
  test_line(ofs,plumed,"RESTRAINT ARG=([a) KAPPA=5 AT=0");

// cannot use simultaneously x2 and x
  test_line(ofs,plumed,"COORDINATION GROUPA=1 GROUPB=2 SWITCH={CUSTOM FUNC=x2+x R_0=1.0}");

// cannot use variables other than x2 or x
  test_line(ofs,plumed,"COORDINATION GROUPA=1 GROUPB=2 SWITCH={CUSTOM FUNC=c R_0=1.0}");

  test_line(ofs,plumed,"EXTERNAL ARG=d FILE=potential LABEL=ext");
  test_line(ofs,plumed,"METAD ARG=d PACE=1 SIGMA=1 HEIGHT=0 FILE=H1 RESTART=WHAT");
  test_line(ofs,plumed,"METAD ARG=d PACE=1 SIGMA=1 TAU=5");
  test_line(ofs,plumed,"METAD ARG=d ADAPTIVE=UNKNOWN PACE=1 SIGMA=1 HEIGHT=5");
  test_line(ofs,plumed,"METAD ARG=d,d1 ADAPTIVE=GEOM PACE=1 SIGMA=1,2 HEIGHT=5");
  test_line(ofs,plumed,"METAD ARG=d,d1 ADAPTIVE=DIFF PACE=1.5 SIGMA=1 HEIGHT=5");
  test_line(ofs,plumed,"METAD ARG=d,d1 ADAPTIVE=GEOM PACE=1 SIGMA=1 HEIGHT=5 SIGMA_MIN=3");
  test_line(ofs,plumed,"METAD ARG=d,d1 ADAPTIVE=GEOM PACE=1 SIGMA=1 HEIGHT=5 SIGMA_MAX=4");
  test_line(ofs,plumed,"PIECEWISE ARG=t POINT0=1.2,10 POINT1=1.3,0 POINT2=1.4,5");
  test_line(ofs,plumed,"SORT ARG=t,d");
  test_line(ofs,plumed,"COMBINE ARG=d,d1 COEFFICIENTS=3");
  test_line(ofs,plumed,"COMBINE ARG=d,d1 COEFFICIENTS=3,3 PARAMETERS=1");
  test_line(ofs,plumed,"COMBINE ARG=d,d1 COEFFICIENTS=3,3 PARAMETERS=1,2 POWERS=4");
  test_line(ofs,plumed,"METAD ARG=d PACE=1 SIGMA=5 HEIGHT=1 GRID_MIN=bla GRID_MAX=100");

// these should not fail
  plumed.cmd("readInputLine","m1: METAD ARG=d PACE=1 SIGMA=5 HEIGHT=1 FILE=H1 FMT=%9.5f");
  plumed.cmd("readInputLine","m2: METAD ARG=d PACE=2 SIGMA=5 HEIGHT=1 FILE=H2 FMT=%9.5f");
  plumed.cmd("readInputLine","PRINT ARG=d,d1,m1.bias FILE=COLVAR FMT=%9.5f");

  test_this(ofs,plumed,"something random here",NULL);
  for(int step=0;step<3;step++){

// this should fail
    test_this(ofs,plumed,"setStep",NULL);

    plumed.cmd("setStep",&step);
    plumed.cmd("setPositions",&positions[0]);
    plumed.cmd("setBox",&box[0]);
    plumed.cmd("setForces",&forces[0]);
    plumed.cmd("setVirial",&virial[0]);
    plumed.cmd("setMasses",&masses[0]);
// set positions after having passed the pointer. They should be accessed here (at "calc").
    for(int i=0;i<3*natoms;i++) positions[i]=i*step;
    plumed.cmd("calc");

// this should fail
    test_this(ofs,plumed,"setMasses",&masses[0]);
  }

  {
    PLMD::Plumed p;

    p.cmd("setNestedExceptions",1);

#define TEST_STD(type) try { p.cmd("throw", #type " msg"); } catch (type & e ) { plumed_assert(std::string(e.what())== #type " msg"); }
    TEST_STD(std::logic_error);
    TEST_STD(std::invalid_argument);
    TEST_STD(std::domain_error);
    TEST_STD(std::length_error);
    TEST_STD(std::out_of_range);
    TEST_STD(std::runtime_error);
    TEST_STD(std::range_error);
    TEST_STD(std::overflow_error);
    TEST_STD(std::underflow_error);

#define TEST_STD_NOMSG(type) try { p.cmd("throw", #type);} catch (type & e ) { }
    TEST_STD_NOMSG(std::bad_cast);
    TEST_STD_NOMSG(std::bad_weak_ptr);
    TEST_STD_NOMSG(std::bad_function_call);
    TEST_STD_NOMSG(std::bad_typeid);
    TEST_STD_NOMSG(std::bad_alloc);
    TEST_STD_NOMSG(std::bad_array_new_length);
    TEST_STD_NOMSG(std::bad_exception);

    TEST_STD_NOMSG(std::bad_variant_access);
    TEST_STD_NOMSG(std::bad_optional_access);
    TEST_STD_NOMSG(std::bad_any_cast);


#define TEST_REGEX(type) try { p.cmd("throw", "std::regex_error std::regex_constants::error_" #type);} catch (std::regex_error & e) { plumed_assert(e.code()==std::regex_constants::error_ ##type); }
    TEST_REGEX(collate);
    TEST_REGEX(ctype);
    TEST_REGEX(escape);
    TEST_REGEX(backref);
    TEST_REGEX(brack);
    TEST_REGEX(paren);
    TEST_REGEX(brace);
    TEST_REGEX(badbrace);
    TEST_REGEX(range);
    TEST_REGEX(space);
    TEST_REGEX(badrepeat);
    TEST_REGEX(complexity);
    TEST_REGEX(stack);

#define TEST_FUTURE(type) try { p.cmd("throw", "std::future_error std::future_errc::" #type);} catch (std::future_error & e) { plumed_assert(e.code()==std::make_error_code(std::future_errc::type)); }

    TEST_FUTURE(broken_promise);
    TEST_FUTURE(future_already_retrieved);
    TEST_FUTURE(promise_already_satisfied);
    TEST_FUTURE(no_state);

    try { p.cmd("throw","PLMD::Exception msg"); } catch (PLMD::Plumed::Exception &e) {
    }
    try { p.cmd("throw","PLMD::ExceptionError msg"); } catch (PLMD::Plumed::ExceptionError &e) {
    }
    try { p.cmd("throw","PLMD::ExceptionDebug msg"); } catch (PLMD::Plumed::ExceptionDebug &e) {
    }
    try { p.cmd("throw","PLMD::lepton::Exception msg"); } catch (PLMD::Plumed::LeptonException &e) {
      plumed_assert(std::string(e.what())=="PLMD::lepton::Exception msg");
    }
    try { p.cmd("throw","std::system_error std::generic_category 100"); } catch (std::system_error & e) {
      plumed_assert(e.code().value()==100)<<" value="<<e.code().value();
      plumed_assert(e.code().category()==std::generic_category());
    }
    try { p.cmd("throw","std::system_error std::system_category 200"); } catch (std::system_error & e) {
      plumed_assert(e.code().value()==200);
      plumed_assert(e.code().category()==std::system_category());
    }
    try { p.cmd("throw","std::system_error std::iostream_category 300"); } catch (std::system_error & e) {
      plumed_assert(e.code().value()==300);
      plumed_assert(e.code().category()==std::iostream_category());
    }
    try { p.cmd("throw","std::system_error std::future_category 400"); } catch (std::system_error & e) {
      plumed_assert(e.code().value()==400);
      plumed_assert(e.code().category()==std::future_category());
    }
    try { p.cmd("throw","std::ios_base::failure"); } catch (std::ios_base::failure & e) {
    }

    try { p.cmd("throw","std::filesystem::filesystem_error std::generic_category 100 a/b/c x/y/z"); } catch (std::filesystem::filesystem_error & e) {
      plumed_assert(e.code().value()==100);
      plumed_assert(e.code().category()==std::generic_category());
      plumed_assert(!e.path1().empty());
      plumed_assert(!e.path2().empty());
      plumed_assert(e.path1().native()=="a/b/c");
      plumed_assert(e.path2().native()=="x/y/z");
    }
    try { p.cmd("throw","std::filesystem::filesystem_error std::generic_category 200 a/b/c"); } catch (std::filesystem::filesystem_error & e) {
      plumed_assert(e.code().value()==200);
      plumed_assert(e.code().category()==std::generic_category());
      plumed_assert(!e.path1().empty());
      plumed_assert(e.path2().empty());
      plumed_assert(e.path1().native()=="a/b/c");
    }
    try { p.cmd("throw","std::filesystem::filesystem_error std::generic_category 300"); } catch (std::filesystem::filesystem_error & e) {
      plumed_assert(e.code().value()==300);
      plumed_assert(e.code().category()==std::generic_category());
      plumed_assert(e.path1().empty());
      plumed_assert(e.path2().empty());
    }

    try { p.cmd("throw","unknown_name"); } catch (PLMD::Plumed::Exception &e) {
    }

    
    {
      bool ok=true;
      try {
        p.cmd("throw","int");
        ok=false;
      } catch(const std::bad_exception &e) {
        std::rethrow_if_nested(e); // this checks that it's not nested anymore
      }
      plumed_assert(ok)<<"should not arrive here";
    }

    {
      bool ok=true;
      try {
        p.cmd("throw","test_nested1");
        ok=false;
      } catch(const PLMD::Plumed::Exception &e) {
        plumed_assert(!std::strcmp(e.what(),"\nouter test_nested1"))<<e.what();
        bool ok=true;
        try {
          std::rethrow_if_nested(e);
          ok=false;
        } catch(const PLMD::Plumed::Exception & e) {
          plumed_assert(!std::strcmp(e.what(),"\nmiddle test_nested1"))<<e.what();
          bool ok=true;
          try {
            std::rethrow_if_nested(e);
            ok=false;
          } catch(const PLMD::Plumed::Exception & e) {
            plumed_assert(!std::strcmp(e.what(),"\ninner test_nested1"))<<e.what();
            std::rethrow_if_nested(e); // this checks that it's not nested anymore
          }
          plumed_assert(ok)<<"should not arrive here";
        }
        plumed_assert(ok)<<"should not arrive here";
      }
      plumed_assert(ok)<<"should not arrive here";
    }

    {
      bool ok=true;
      try {
        p.cmd("throw","test_nested2");
        ok=false;
      } catch(const PLMD::Plumed::Exception &e) {
        plumed_assert(!std::strcmp(e.what(),"\nouter test_nested2"))<<e.what();
        bool ok=true;
        try {
          std::rethrow_if_nested(e);
          ok=false;
        } catch(const PLMD::Plumed::Exception & e) {
          plumed_assert(!std::strcmp(e.what(),"\nmiddle test_nested2"))<<e.what();
          bool ok=true;
          try {
            std::rethrow_if_nested(e);
            ok=false;
          } catch(const std::bad_alloc & e) {
            std::rethrow_if_nested(e); // this checks that it's not nested anymore
          }
          plumed_assert(ok)<<"should not arrive here";
        }
        plumed_assert(ok)<<"should not arrive here";
      }
      plumed_assert(ok)<<"should not arrive here";
    }

    {
      bool ok=true;
      try {
        p.cmd("throw","test_nested3");
        ok=false;
      } catch(const PLMD::Plumed::Exception &e) {
        plumed_assert(!std::strcmp(e.what(),"\nouter test_nested3"))<<e.what();
        bool ok=true;
        try {
          std::rethrow_if_nested(e);
          ok=false;
        } catch(const PLMD::Plumed::Exception & e) {
          plumed_assert(!std::strcmp(e.what(),"\nmiddle test_nested3"))<<e.what();
          bool ok=true;
          try {
            std::rethrow_if_nested(e);
            ok=false;
          } catch(const PLMD::Plumed::Exception & e) {
            plumed_error()<<"should not arrive here";
          } catch(const std::exception & e) {
            plumed_assert(!std::strcmp(e.what(),"inner"))<<e.what();
            std::rethrow_if_nested(e); // this checks that it's not nested anymore
          }
          plumed_assert(ok)<<"should not arrive here";
        }
        plumed_assert(ok)<<"should not arrive here";
      }
      plumed_assert(ok)<<"should not arrive here";
    }

    p.cmd("setNestedExceptions",0);

    {
      bool ok=true;
      try {
        p.cmd("throw","test_nested1");
        ok=false;
      } catch(const PLMD::Plumed::Exception &e) {
         plumed_assert(!std::strcmp(e.what(),
           "\ninner test_nested1\n\nThe above exception was the direct cause of the following exception:\n"
           "\nmiddle test_nested1\n\nThe above exception was the direct cause of the following exception:\n"
           "\nouter test_nested1"
         ))<<e.what();
        std::rethrow_if_nested(e); // this checks that it's not nested anymore
      }
      plumed_assert(ok)<<"should not arrive here";
    }

    {
      bool ok=true;
      try {
        p.cmd("throw","test_nested3");
        ok=false;
      } catch(const PLMD::Plumed::Exception &e) {
         plumed_assert(!std::strcmp(e.what(),
           "inner\n\nThe above exception was the direct cause of the following exception:\n"
           "\nmiddle test_nested3\n\nThe above exception was the direct cause of the following exception:\n"
           "\nouter test_nested3"
         ))<<e.what();
        std::rethrow_if_nested(e); // this checks that it's not nested anymore
      }
      plumed_assert(ok)<<"should not arrive here";
    }
  }

  return 0;

}
