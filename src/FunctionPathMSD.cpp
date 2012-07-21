#include <cmath>
// this is an hack for cmath disambiguation

double (*cmathLog)(double) = log; 

#include "Function.h"
#include "ActionRegister.h"

#include <string>
#include <cstring>
#include <cassert>
//#include "ColvarRMSD.h"
#include <iostream>

using namespace std;

namespace PLMD{

//+PLUMEDOC FUNCTION PATHMSD
/*

This is the Path Collective Variables implementation 
( see \cite brand07 ).
This variable computes the progress along a given set of frames that is provided  
in input ("s" component) and the distance from them ("z" component). 
It is a function of MSD that are obtained by the joint use of MSD variable and SQUARED flag 
(see below).

\par Examples

Here below is a case where you have defined three frames and you want to  
calculate the progress alng the path and the distance from it in p1

\verbatim
t1: RMSD REFERENCE=frame_1.dat TYPE=OPTIMAL SQUARED
t2: RMSD REFERENCE=frame_21.dat TYPE=OPTIMAL SQUARED
t3: RMSD REFERENCE=frame_42.dat TYPE=OPTIMAL SQUARED
p1: PATHMSD ARG=t1,t2,t3 LAMBDA=500.0 
PRINT ARG=t1,t2,t3,p1.s,p1.z STRIDE=1 FILE=colvar FMT=%8.4f
\endverbatim

*/
//+ENDPLUMEDOC
   
class FunctionPathMSD : public Function {
  double lambda;
  bool pbc;
  int neigh_size;
  double neigh_stride;
  vector< pair<unsigned,double> > neighpair;
public:
  FunctionPathMSD(const ActionOptions&);
// active methods:
  virtual void calculate();
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(FunctionPathMSD,"PATHMSD")

void FunctionPathMSD::registerKeywords(Keywords& keys){
  Function::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","LAMBDA","all compulsory keywords should be added like this with a description here");
  keys.add("optional","NEIGH_SIZE","all optional keywords that have input should be added like a description here");
  keys.add("optional","NEIGH_STRIDE","all optional keywords that have input should be added like a description here");
}
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// this below is useful when one wants to sort a vector of double and have back the order 
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// create a custom sorter
typedef vector<double>::const_iterator myiter;
struct ordering {
       bool operator ()(pair<unsigned , myiter> const& a, pair<unsigned , myiter> const& b) {
           return *(a.second) < *(b.second);
       };
};
// sorting utility
vector<int> increasingOrder( vector<double> &v){
   // make a pair
   vector< pair<unsigned , myiter> > order(v.size());
   unsigned n = 0;
   for (myiter it = v.begin(); it != v.end(); ++it, ++n){
       order[n] = make_pair(n, it); // note: heere i do not put the values but the addresses that point to the value 
   };
   // now sort according the second value
   sort(order.begin(), order.end(), ordering());
   typedef vector< pair<unsigned , myiter> >::const_iterator pairiter;
   vector<int> vv(v.size());n=0;
   for ( pairiter it = order.begin(); it != order.end(); ++it ){
       vv[n]=(*it).first;n++; 
   }
   return vv;
};
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

struct pairordering {
       bool operator ()(pair<unsigned , double> const& a, pair<unsigned , double> const& b) {
           return (a).second > (b).second;
       };
};

FunctionPathMSD::FunctionPathMSD(const ActionOptions&ao):
Action(ao),
Function(ao),
pbc(true),
neigh_size(-1),
neigh_stride(-1.)
{
//  bool nopbc=!pbc;
//  parseFlag("NOPBC",nopbc);
//  pbc=!nopbc;

  parse("LAMBDA",lambda);
  parse("NEIGH_SIZE",neigh_size);
  parse("NEIGH_STRIDE",neigh_stride);
  checkRead();
  log.printf("  lambda is %f\n",lambda);
  // list the action involved and check the type 
  for(unsigned i=0;i<getNumberOfArguments();i++){
       // for each value get the name and the label of the corresponding action
       std::string mylabel=getPntrToArgument(i)->getPntrToAction()->getLabel(); 
       std::string myname=getPntrToArgument(i)->getPntrToAction()->getName(); 
       if(myname!="RMSD")plumed_merror("This argument is not of RMSD type!!!");
  }   
  log.printf("  Consistency check completed! Your path cvs look good!\n"); 
  // do some neighbor printout
  if(neigh_stride>0. || neigh_size>0){
           if(neigh_size>getNumberOfArguments()){
           	log.printf(" List size required ( %d ) is too large: resizing to the maximum number of arg required: %d  \n",neigh_size,getNumberOfArguments());
 		neigh_size=getNumberOfArguments();
           }
           log.printf("  Neighbor list enabled: \n");
           log.printf("                size   :  %d elements\n",neigh_size);
           log.printf("                stride :  %f time \n",neigh_stride);
  }else{
           log.printf("  Neighbor list NOT enabled \n");
  }

  addComponentWithDerivatives("s"); componentIsNotPeriodic("s");
  addComponentWithDerivatives("z"); componentIsNotPeriodic("z");
}
// calculator
void FunctionPathMSD::calculate(){
  double s_path=0.;
  double partition=0.;
  double tmp;
  if(neighpair.empty()){
       neighpair.resize(getNumberOfArguments());  
       for(unsigned i=0;i<getNumberOfArguments();i++)neighpair[i].first=i; 
  }

  Value* val_s_path=getPntrToComponent("s");
  Value* val_z_path=getPntrToComponent("z");
  typedef  vector< pair<unsigned,double> >::iterator pairiter;
  for(pairiter it=neighpair.begin();it!=neighpair.end();++it){ 
    unsigned n=(*it).first;
    (*it).second=exp(-lambda*getArgument(n));
    s_path+=(n+1)*(*it).second;
    partition+=(*it).second;
  }
  s_path/=partition;
  val_s_path->set(s_path);
  val_z_path->set(-(1./lambda)*cmathLog(partition));
  for(pairiter it=neighpair.begin();it!=neighpair.end();++it){ 
    unsigned n=(*it).first;
    double expval=(*it).second;
    tmp=lambda*expval*(s_path-(n+1))/partition;
    setDerivative(val_s_path,n,tmp);
    setDerivative(val_z_path,n,expval/partition);
  }

  // neighbor list: rank and activate the chain for the next step 

  // neighbor list: if neigh_size<0 never sort and keep the full vector
  // neighbor list: if neigh_size>0  
  //                if the size is full -> sort the vector and decide the dependencies for next step 
  //                if the size is not full -> check if next step will need the full dependency otherwise keep this dependencies 

  vector<int> actions_to_be_suspended; 
  if (neigh_size>0){
     if(neighpair.size()==getNumberOfArguments()){ // I just did the complete round: need to sort, shorten and give it a go
                // sort the values  
		sort(neighpair.begin(),neighpair.end(),pairordering());
                // create a list of the actions to be suspended
                for(unsigned i=neigh_size;i<getNumberOfArguments();i++)actions_to_be_suspended.push_back(neighpair[i].first);
                // resize the effective list
                neighpair.resize(neigh_size);
     }else{
        if( int(getStep())%int(neigh_stride/getTimeStep())==0 ){
                 log.printf(" Time %f : recalculating full neighlist \n",getStep()*getTimeStep());
     		 neighpair.resize(getNumberOfArguments());  
     		 for(unsigned i=0;i<getNumberOfArguments();i++)neighpair[i].first=i; 
        }
     } 
  }
  // TODO prepare dependencies for next step n 
  // make a backup of the pointers to the actions (remove instantaneously the dependency) 
  //for(Dependencies::iterator p=after.begin();p!=after.end();++p) (*p)->activate(); 

}

}



