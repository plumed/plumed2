#include <functional>
#include "core/ActionRegister.h"
#include "core/ActionWithValue.h"
#include "tools/ConjugateGradient.h" 

namespace PLMD {

class myfunc {
public:
  double getEng( const double& x );
};

class my2dfunc {
public:
  double getEng( const std::vector<double>& pp, std::vector<double>& der );
};

class TestMin : 
  public ActionWithValue {
private:
public:
  static void registerKeywords( Keywords& keys );
  TestMin(const ActionOptions&);
  unsigned getNumberOfDerivatives(){ return 0; }
  void calculate();
  void apply(){} 
};

PLUMED_REGISTER_ACTION(TestMin,"TESTMIN")

void TestMin::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
}

TestMin::TestMin(const ActionOptions& ao):
Action(ao),
ActionWithValue(ao)
{
  addComponent("1d"); componentIsNotPeriodic("1d");
  addComponent("2d-1"); componentIsNotPeriodic("2d-1");
  addComponent("2d-2"); componentIsNotPeriodic("2d-2"); 
}

void TestMin::calculate(){

  myfunc ff;
  Minimise1DBrent<myfunc> bb(ff);

  bb.bracket( 5.0, 4.0, &myfunc::getEng );
  double xmin=bb.minimise( &myfunc::getEng );
  getPntrToComponent(0)->set( xmin );

  my2dfunc* ff2 = new my2dfunc(); 
  ConjugateGradient<my2dfunc> bb2( ff2 );
  double tol; std::vector<double> pp(2);
  pp[0]=5; pp[1]=4;
  bb2.minimise( tol, pp, &my2dfunc::getEng );
  getPntrToComponent(1)->set( pp[0] );
  getPntrToComponent(2)->set( pp[1] );
  delete ff2;
}


double myfunc::getEng( const double& x ){
  return (x-10) *(x-10);
}

double my2dfunc::getEng( const std::vector<double>& pp, std::vector<double>& der ){
  double eng=0;
  std::vector<double> min(2); min[0]=3; min[1]=2.5;
  for(unsigned i=0;i<pp.size();++i){
      der[i]=pp[i]-min[i];
      eng+=der[i]*der[i];
  }
  return eng;
}

}
