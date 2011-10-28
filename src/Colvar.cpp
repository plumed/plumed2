#include "Colvar.h"
#include "PlumedMain.h"

using namespace std;
using namespace PLMD;

Colvar::Colvar(const ActionOptions&ao) :
ActionAtomistic(ao),
doall(true),
domin(false),
domax(false),
dototal(false),
domean(false),
dolt(false),
domt(false),
dohist(false),
isCSphereF(false)
{
  forbidKeyword("STRIDE");
  registerKeyword(0, "MIN", "calculate the minimum for the defined colvars");
  registerKeyword(0, "MAX", "calculate the maximum for the defined colvars");
  registerKeyword(0, "BETA", "(default=50) value used to create smooth derivatives in calculations of MIN/MAX"); 
  registerKeyword(0, "SUM", "calculate the sum of all the defined colvars");
  registerKeyword(0, "AVERAGE", "compute the average of the defined colvars");
  registerKeyword(0, "LESS_THAN", "compute the number of colvars that that are less than than a particular value using a smooth switching function");
  registerKeyword(0, "MORE_THAN", "compute the number of colvars that that are more than than a particular value using a smooth switching function");
  registerKeyword(0, "LOGIC_NN", "(default=6) value of NN in switching functions for MORE_THAN/LESS_THAN");
  registerKeyword(0, "LOGIC_MM", "(default=12) value of MM in switching functions for MORE_THAN/LESS_THAN");
  registerKeyword(0, "HISTOGRAM", "calculate the histogram of a distribution of collective coordinates within a certain range");
  registerKeyword(0, "NBINS", "number of bins to use in the calculation of the histogram");
  registerKeyword(0, "SMEAR", "(default=0.5) the ammount to smear the values, relative to bin size, in the calculation of the histogram");
  registerKeyword(3, "WITHIN_RANGE", "calculate the number of cvs from the distribution that are within a particular range");
}

void Colvar::readActionColvar( int natoms, const std::vector<double>& domain ){
  if(isCSphereF){
    int natoms=0; unsigned ngrp=static_cast<unsigned>(2);
    readActionAtomistic( natoms, ngrp ); 
  } else { 
    unsigned ngrp=static_cast<unsigned>(natoms); 
    readActionAtomistic( natoms, ngrp ); 
  }
  readActionWithExternalArguments( 3*getNumberOfAtoms()+9, domain );

  // Setup everything for calculation of individual colvars
  skipto.resize( getNumberOfColvars() ); 
  for(unsigned i=0;i<getNumberOfColvars();++i) skipto[i]=i+1;

  // Resize stuff for applying forces
  f.resize( getNumberOfAtoms() ); forces.resize( 3*getNumberOfAtoms()+9 );

  // Read in everything that tells us what sort of calculation we are doing
  parseFlag("MIN",domin);
  if (domin){
     doall=false; 
     beta=50.0/plumed.getAtoms().getUnits().length; parse("BETA",beta);
     addValue("min", true, true);
     log.printf("  using minimum value \n");
  }
  parseFlag("MAX",domax);
  if (domax){
     doall=false;
     error("Dont know how to implement this yet");
     addValue("max", true, true);
     log.printf("  using maximum value \n");
  }
  parseFlag("SUM",dototal);
  if( dototal ){
     doall=false;
     addValue("sum", true, true);
     log.printf("  using sum of all values \n");
  }
  parseFlag("AVERAGE",domean);
  if( domean ){
     doall=false;
     addValue("average", true, true);
     log.printf("  using average value \n");
  }

  std::vector<double> r_0;
  parseVector("LESS_THAN",r_0);
  if( r_0.size()!=0 ){
     doall=false; dolt=true;
     if ( r_0.size()!=1 ) error("input for LESS_THAN makes no sense");
     int nn=6; parse("LOGIC_NN",nn);
     int mm=12; parse("LOGIC_MM",mm);
     ltswitch.set(nn, mm, r_0[0], 0.0);
     log.printf("  number of values less than %f.  Switching function paramers are %d %d \n", r_0[0], nn, mm );
     std::string lm; Tools::convert( r_0[0], lm );
     addValue("less_than" + lm, true, true);
  }

  std::vector<double> r_1;
  parseVector("MORE_THAN",r_1);
  if ( r_1.size()!=0 ) {
     doall=false; domt=true;
     if ( r_1.size()!=1 ) error("Input for MORE_THAN makes no sense");
     int nn=6; parse("LOGIC_NN",nn);
     int mm=12; parse("LOGIC_MM",mm);
     mtswitch.set(nn, mm, r_1[0], 0.0);
     log.printf("  number of values greater than %f.  Switching function paramers are %d %d \n", r_1[0], nn, mm );
     std::string lm; Tools::convert( r_1[0], lm );
     addValue("more_than" + lm, true, true);
  }

  std::vector<double> hrange;
  parseVector("HISTOGRAM",hrange);
  if( hrange.size()!=0 ){
      if( !doall ) error("you cannot mix HISTOGRAM calculation with other COLVAR modifiers");
      if( hrange.size()!=2 || hrange[0]>=hrange[1] ) error("range specified for histogram makes no sense");
      doall=false; dohist=true;
      int nbins=-1; parse("NBINS",nbins);
      if (nbins<0) error("number of bins is unspecified for histogram");
      double smear=0.5; parse("SMEAR",smear);
      double delr = ( hrange[1] - hrange[0] ) / static_cast<double>(nbins);
      histogram.resize( nbins );
      for (int i=0; i<nbins; ++i){ 
        histogram[i].set( hrange[0]+i*delr, hrange[0]+(i+1)*delr, smear*delr );
        std::string lb, ub; Tools::convert( hrange[0]+i*delr, lb); Tools::convert( hrange[0]+(i+1)*delr, ub); 
        addValue("between" + lb + "&" + ub, true, true );
      }
      hvalues.resize(nbins);
      log.printf("  calculating a histogram using %d bins \n",histogram.size() );
      for (unsigned i=0; i<histogram.size(); ++i) {
         log.printf("  bin %d counts values between %f and %f \n", i+1, histogram[i].getlowb(), histogram[i].getbigb() );
      }
  }

  if( testForKey("WITHIN_RANGE") ){
     std::vector<double> range; HistogramBead tmpbead;
     if( !doall ) error("you cannot mix WITHIN_RANGE calculation with other COLVAR modifiers");
     doall=false; dohist=true;
     double smear=0.5; parse("SMEAR",smear);
     if( testForNumberedKeys("WITHIN_RANGE") ){
         for(int i=1;; ++i ){
            range.clear();
            if( !parseNumberedVector( "WITHIN_RANGE", i, range ) ) break;
            if( range.size()!=2 || range[0]>=range[1] ) {
                std::string num; Tools::convert(i,num);
                error("range specified for WITHIN_RANGE" + num + " makes no sense");
            }
            tmpbead.set( range[0], range[1], smear*( range[1]-range[0] ) );
            histogram.push_back( tmpbead );
            std::string lb, ub; Tools::convert( range[0], lb); Tools::convert( range[1], ub);
            addValue("between" + lb + "&" + ub, true, true ); 
         }
         hvalues.resize( histogram.size() );
     } else {
         parseVector("WITHIN_RANGE", range); 
         if( range.size()!=2 || range[0]>=range[1] ) {
             error("range specified for WITHIN_RANGE makes no sense");
         }
         tmpbead.set( range[0], range[1], smear*( range[1]-range[0] ) );
         histogram.push_back( tmpbead );
         std::string lb, ub; Tools::convert( range[0], lb); Tools::convert( range[1], ub);
         addValue("between" + lb + "&" + ub, true, true );
         hvalues.resize(1);
     }
     log.printf("  calculating a histogram using %d bins \n",histogram.size() );
     for (unsigned i=0; i<histogram.size(); ++i) {
         log.printf("  bin %d counts values between %f and %f \n", i+1, histogram[i].getlowb(), histogram[i].getbigb() );
     }
  }

  if( doall ){
     std::string n;
     if ( !isCSphereF && updateIsOn() ){
        error("UPDATE keyword cannot be used when calculating a set of distinguishable colvars");
     } else if( usingDynamicGroups() ){
        error("using dynamic groups is incompatible with calculating all colvars consider MIN/LESS_THAN/etc");
     }

     for(unsigned i=0;i<getNumberOfColvars();++i){
        Tools::convert(i,n); addValue("value" + n, false, true );
     }
  }
}

void Colvar::calculate(){
  double df, tmp, value, mintotal, ttotal, atotal, maxtotal, lttotal, mttotal;
  mintotal=maxtotal=ttotal=atotal=lttotal=mttotal=0.0; 
  hvalues.assign( hvalues.size(), 0.0 );

  if( updateTime() ) updateDynamicAtoms();

  std::string mtstring, ltstring;
  if( dolt ) ltstring="less_than" + ltswitch.get_r0_string();
  if( domt ) mtstring="more_than" + mtswitch.get_r0_string(); 

  for (unsigned i=0; i<skipto.size(); i=skipto[i] ) {
     value=calcFunction( i );
     if (doall) {
        mergeFunctions( i, i, 1.0 );
        setValue( i, value, 1.0 );
     } else if (dohist) {
        for (unsigned j=0; j<histogram.size(); ++j) {
          hvalues[j]+=histogram[j].calculate( value, df );
          mergeFunctions( j, i, df );
        }
     } else {
        if ( domin ) {
           // Minimum
           tmp=exp( beta/value ); mintotal+=tmp; df=tmp/(value*value); 
           mergeFunctions("min", i, df );
        }
        if ( domax ) {
           // Maximum
           assert(false);
        } 
        if ( dototal ) {
           ttotal+=value;
           mergeFunctions("sum", i, 1.0 );
        } 
        if ( domean ) {
           // The mean
           atotal+=value;
           mergeFunctions("average", i, 1.0 );
        } 
        if ( dolt ) {
           // Less than
           lttotal+=ltswitch.calculate(value, df);
           mergeFunctions(ltstring, i, df*value );
        } 
        if ( domt ) {
           // More than
           mttotal+=1.0 - mtswitch.calculate(value, df);
           mergeFunctions(mtstring, i, -df*value );
        }
     }
  }
  if ( dohist ){ 
    for (unsigned j=0; j<histogram.size(); ++j) setValue( j, hvalues[j], 1.0 );
  } else {
    if ( domin ){ double dist=beta/std::log(mintotal); setValue("min", dist, dist*dist/mintotal ); }
    if ( domax ) { }
    if ( dototal ) { setValue("sum", ttotal, 1.0 ); }
    if ( domean ) { setValue("average", atotal/skipto.size(), 1.0/skipto.size() ); }
    if ( dolt ) {  setValue(ltstring, lttotal, 1.0 ); }
    if ( domt ) { setValue(mtstring, mttotal, 1.0 ); }
  }

  // Do derivatives for dynamic group
  if( usingDynamicGroups() ) addGroupDerivatives(); 
}

void Colvar::apply(){

  const unsigned nat=f.size();
  for(unsigned i=0;i<nat;i++){
    f[i][0]=0.0; f[i][1]=0.0; f[i][2]=0.0;
  }

  Tensor v; v.clear();
  for(int i=0;i<getNumberOfValues();++i){
    if( getForces( i, forces ) ){
       for(unsigned j=0;j<nat;++j){
          f[j][0]+=forces[3*j+0];
          f[j][0]+=forces[3*j+1];
          f[j][0]+=forces[3*j+2];
       }
       v(0,0)+=forces[3*nat+0];
       v(0,1)+=forces[3*nat+1];
       v(0,2)+=forces[3*nat+2];
       v(1,0)+=forces[3*nat+3];
       v(1,1)+=forces[3*nat+4];
       v(1,2)+=forces[3*nat+5];
       v(2,0)+=forces[3*nat+6];
       v(2,1)+=forces[3*nat+7];
       v(2,2)+=forces[3*nat+8]; 
    }
  }
  applyForces( f, v );
}





