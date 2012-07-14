#include <cmath>
// this is an hack for cmath disambiguation

double (*cmathLog)(double) = log; 

#include "Function.h"
#include "ActionRegister.h"

#include <string>
#include <cstring>
#include <cassert>
#include "ColvarRMSD.h"
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
  int neigh_stride;
  vector<int> neighlist;
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
   //    cerr<<(*it).first<<" "<<(*(*it).second)<<endl;
       vv[n]=(*it).first;n++; 
   }
   return vv;
};

FunctionPathMSD::FunctionPathMSD(const ActionOptions&ao):
Action(ao),
Function(ao),
pbc(true),
neigh_stride(-1),
neigh_size(-1)
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
       // check of the SQUARED flag is set in all the dependent variables
       ColvarRMSD* ptr=dynamic_cast<ColvarRMSD*>(getPntrToArgument(i)->getPntrToAction());
       if(ptr){
            log.printf("  The cv type for %s is ColvarRMSD: good! \n",mylabel.c_str());
            if((*ptr).squared){
               log.printf("  You enabled the SQUARED option in RMSD! it is the original flavour ! \n");
            }else{
               log.printf("  BEWARE: In ARG %s  You did not enable the SQUARED option in RMSD! it is not the original flavour!\n",mylabel.c_str());
               plumed_merror("There are problems in the pathcv setup. Check the log!!!");
            }
       }else{
            log.printf("  Hey, the CV %s used for the path has wrong type. Must be RMSD with added SQUARED flag!  \n",mylabel.c_str());
            plumed_merror("  There are problems in the pathcv setup. Check the log!!!");
       }
       // check structural consistency: indexing must be the same 
       if(i!=0){
       		ColvarRMSD* ptr0=dynamic_cast<ColvarRMSD*>(getPntrToArgument(i)->getPntrToAction());
                if( ptr->getNumberOfAtoms() != ptr0->getNumberOfAtoms()  ){
	            log.printf("  Hey, the CV %s used for the path has wrong number of atoms compared to the other frames involved in the path cv!! Check it out!  \n",mylabel.c_str());
       		 	    plumed_merror("There are problems in the pathcv setup. Check the log!!!");
		}
                for( unsigned ii=0; ii< ptr->getNumberOfAtoms();ii++ ){
                     if(  ptr->getAbsoluteIndex(ii)!=ptr0->getAbsoluteIndex(ii) ){
		            log.printf("  Hey, the CV %s used for the path has wrong number of atoms compared to the other frames involved in the path cv!! Check it out!  \n",mylabel.c_str());
       		 	    plumed_merror("There are problems in the pathcv setup. Check the log!!!");
                     }
                } 
       }
  }   
  log.printf("  Consistency check completed! Your path cvs look good!\n"); 
  // do some neighbor printout
  if(neigh_stride!=-1 || neigh_size!=-1){
           log.printf("  Neighbor list enabled: \n");
           log.printf("                size   :  %d elements\n",neigh_size);
           log.printf("                stride :  %d steps\n",neigh_stride);
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
  // no neighbor list
  if(neigh_size<0 ){
    neighlist.resize(getNumberOfArguments());  
    for(unsigned i=0;i<getNumberOfArguments();i++)neighlist[i]=i; 
  }else{
  // if on step, recalculate
     
  // use the neighbor list 
  }
  Value* val_s_path=getPntrToComponent("s");
  Value* val_z_path=getPntrToComponent("z");
  vector<double> zvec(getNumberOfArguments());
  for(unsigned i=0;i<getNumberOfArguments();i++){ 
    zvec[i]=exp(-lambda*getArgument(i));
    s_path+=(i+1)*zvec[i];
    partition+=zvec[i];
  }
  s_path/=partition;
  val_s_path->set(s_path);
  val_z_path->set(-(1./lambda)*cmathLog(partition));
  for(unsigned i=0;i<getNumberOfArguments();i++){ 
    tmp=lambda*zvec[i]*(s_path-(i+1))/partition;
    setDerivative(val_s_path,i,tmp);
    setDerivative(val_z_path,i,zvec[i]/partition);
  }
  // neighbor list: rank and activate the chain for the next step 
  vector<int> order=increasingOrder(zvec); // this is the order of the exponential, therefore we should take from the back 
     

}

}



