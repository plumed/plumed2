/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2020 The plumed team
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
#include "ActionForInterface.h"
#include "ActionToPutData.h"
#include "ActionAtomistic.h"
#include "ActionRegister.h"
#include "PlumedMain.h"
#include "ActionSet.h"
#include "tools/Communicator.h"

namespace PLMD {

class DomainDecomposition : public ActionForInterface {
private:
  class DomainComms : public Communicator {
  public:
    bool on;
    bool async;

    std::vector<Communicator::Request> mpi_request_positions;
    std::vector<Communicator::Request> mpi_request_index;

    std::vector<double> positionsToBeSent;
    std::vector<double> positionsToBeReceived;
    std::vector<int>    indexToBeSent;
    std::vector<int>    indexToBeReceived;
    operator bool() const {return on;}
    DomainComms(): on(false), async(false) {}
    void enable(Communicator& c);
  };
  DomainComms dd;
  long int ddStep;

  std::vector<int>    gatindex;
/// Map global indexes to local indexes
/// E.g. g2l[i] is the position of atom i in the array passed from the MD engine.
/// Called "global to local" since originally it was used to map global indexes to local
/// ones used in domain decomposition. However, it is now also used for the NAMD-like
/// interface, where only a small number of atoms is passed to plumed.
  std::vector<int> g2l;

  unsigned shuffledAtoms;

  bool asyncSent;
/// This holds the list of unique atoms 
  std::vector<unsigned> uniq_index;
/// This holds the list of actions that are set from this action
  std::vector<ActionToPutData*> inputs;
/// This holds all the actions that read atoms
  std::vector<ActionAtomistic*> actions;
/// The list that holds all the atom indexes we need
  std::vector<int> fullList;
/// This actually does the sharing of the data across the domains
  void share(const bool& getallatoms, const std::set<AtomNumber>& unique);
public:
  static void registerKeywords(Keywords& keys);
  explicit DomainDecomposition(const ActionOptions&ao);
  bool retrieveRequiredAtoms( const std::vector<int>& g2l );
  void setStride( const std::string& name, const unsigned& sss) override ;
  void resetForStepStart() override ; 
  bool setValuePointer( const std::string& name, void* val ) override ;
  bool setForcePointer( const std::string& name, void* val ) override ;
  unsigned getNumberOfForcesToRescale() const override ;
  void share() override ;
  void shareAll() override ;
  void wait() override ;
  void writeBinary(std::ostream&o);
  void readBinary(std::istream&i);
  void apply() override ;
  void setAtomsNlocal(int n) override ;
  void setAtomsGatindex(int*g,bool fortran) override ;
  void setAtomsContiguous(int start) override ;
  void Set_comm(Communicator& comm) override ;
  void broadcastToDomains( Value* val ) override ;
  void sumOverDomains( Value* val ) override ;
  const long int& getDdStep() const override ;
  const std::vector<int>& getGatindex() const override ;
  bool hasFullList() const override { return true; }
  void createFullList(int*) override ;
  void getFullList(int**) override ;
  void clearFullList() override ;
  bool onStep() const override { return getPntrToOutput(0)->getShape()[0]>0; }
};

PLUMED_REGISTER_ACTION(DomainDecomposition,"DOMAIN_DECOMPOSITION")

void DomainDecomposition::DomainComms::enable(Communicator& c) {
  on=true;
  Set_comm(c.Get_comm());
  async=Get_size()<10;
  if(std::getenv("PLUMED_ASYNC_SHARE")) {
    std::string s(std::getenv("PLUMED_ASYNC_SHARE"));
    if(s=="yes") async=true;
    else if(s=="no") async=false;
    else plumed_merror("PLUMED_ASYNC_SHARE variable is set to " + s + "; should be yes or no");
  }
}

void DomainDecomposition::registerKeywords(Keywords& keys) {
  ActionForInterface::registerKeywords( keys ); keys.remove("ROLE");
  keys.add("compulsory","NATOMS","the number of atoms across all the domains");
  keys.add("numbered","VALUE","value to create for the domain decomposition");
  keys.add("numbered","UNIT","unit of the ith value that is being passed through the domain decomposition");
  keys.add("numbered","CONSTANT","is the ith value that is being passed through the domain decomposition constant");
  keys.add("numbered","PERIODIC","if the value being passed to plumed is periodic then you should specify the periodicity of the function.  If the value "
                                   "is not periodic you must state this using PERIODIC=NO.  Positions are passed with PERIODIC=NO even though special methods are used "
                                    "to deal with pbc");
  keys.add("numbered","ROLE","Get the role this value plays in the code can be x/y/z/m/q to signify that this is x, y, z positions of atoms or masses or charges of atoms");
}

DomainDecomposition::DomainDecomposition(const ActionOptions&ao):
Action(ao),
ActionForInterface(ao),
ddStep(0),
shuffledAtoms(0),
asyncSent(false)
{
  // Read in the number of atoms
  int natoms; parse("NATOMS",natoms);
  std::string str_natoms; Tools::convert( natoms, str_natoms );
  // Create a value to hold something
  std::vector<unsigned> shape(1); shape[0]=natoms; addValue( shape ); setNotPeriodic();
  // Setup the gat index 
  gatindex.resize(natoms); for(unsigned i=0; i<gatindex.size(); i++) gatindex[i]=i;
  // Turn on the domain decomposition
  if( Communicator::initialized() ) Set_comm(comm);
  // Now read in the values we are creating here
  for(unsigned i=1;;++i) {
      std::string valname;
      if( !parseNumbered("VALUE",i,valname) ) break;
      std::string unit; parseNumbered("UNIT",i,unit);
      std::string constant; parseNumbered("CONSTANT",i,constant);
      std::string period; parseNumbered("PERIODIC",i,period);
      std::string role; parseNumbered("ROLE",i,role);
      if( constant=="True") plumed.readInputLine( valname + ": PUT CONSTANT SHAPE=" + str_natoms + " UNIT=" + unit + " PERIODIC=" + period + " ROLE=" + role ); 
      else if( constant=="False") plumed.readInputLine( valname + ": PUT SHAPE=" + str_natoms + " UNIT=" + unit + " PERIODIC=" + period + " ROLE=" + role );
      else plumed_merror("missing information on whether value is constant");
      // And save the list of values that are set from here
      ActionToPutData* ap=plumed.getActionSet().selectWithLabel<ActionToPutData*>(valname); ap->addDependency( this ); inputs.push_back( ap );
  }
  plumed.readInputLine("Box: PBC");
}

void DomainDecomposition::Set_comm(Communicator& comm) {
  dd.enable(comm); setAtomsNlocal(getPntrToOutput(0)->getShape()[0]); setAtomsContiguous(0);
  if( dd.Get_rank()!=0 ) {
      ActionToPutData* ap=plumed.getActionSet().selectWithLabel<ActionToPutData*>("Box"); ap->noforce=true; 
  } 
}

void DomainDecomposition::resetForStepStart() { 
  for(const auto & pp : inputs) pp->resetForStepStart(); 
}

void DomainDecomposition::setStride( const std::string& name, const unsigned& sss) {
  for(const auto & pp : inputs) {
      if( pp->getLabel()==name ) { pp->setStride(name, sss); return; }
  } 
  plumed_error();
}

bool DomainDecomposition::setValuePointer( const std::string& name, void* val ) {
  wasset=true;  // Once the domain decomposition stuff is transferred moved the setting of this to where the g2l vector is setup
  for(const auto & pp : inputs) {
      if( pp->setValuePointer( name, val ) ) return true; 
  }
  return false;
}

bool DomainDecomposition::setForcePointer( const std::string& name, void* val ) {
  for(const auto & pp : inputs) {
      if( pp->setForcePointer( name, val ) ) return true;
  }
  return false;
}

void DomainDecomposition::setAtomsNlocal(int n) {
  gatindex.resize(n);
  g2l.resize(getPntrToOutput(0)->getShape()[0],-1);
  if(dd) {
// Since these vectors are sent with MPI by using e.g.
// &dd.positionsToBeSent[0]
// we make sure they are non-zero-sized so as to
// avoid errors when doing boundary check
    if(n==0) n++;
    unsigned nvals = inputs.size(), natoms = getPntrToValue()->getShape()[0];
    dd.positionsToBeSent.resize(n*nvals,0.0);
    dd.positionsToBeReceived.resize(natoms*nvals,0.0);
    dd.indexToBeSent.resize(n,0);
    dd.indexToBeReceived.resize(natoms,0);
  }
}

void DomainDecomposition::setAtomsGatindex(int*g,bool fortran) {
  plumed_massert( g || gatindex.size()==0, "NULL gatindex pointer with non-zero local atoms");
  ddStep=getStep();
  if(fortran) {
    for(unsigned i=0; i<gatindex.size(); i++) gatindex[i]=g[i]-1;
  } else {
    for(unsigned i=0; i<gatindex.size(); i++) gatindex[i]=g[i];
  }
  for(unsigned i=0; i<g2l.size(); i++) g2l[i]=-1;
  if( gatindex.size()==getPntrToOutput(0)->getShape()[0] ) {
    shuffledAtoms=0;
    for(unsigned i=0; i<gatindex.size(); i++) {
      if( gatindex[i]!=i ) { shuffledAtoms=1; break; }
    }
  } else {
    shuffledAtoms=1;
  }
  if(dd) dd.Sum(shuffledAtoms);
  for(unsigned i=0; i<gatindex.size(); i++) g2l[gatindex[i]]=i;
}

void DomainDecomposition::setAtomsContiguous(int start) {
  ddStep=plumed.getStep();
  for(unsigned i=0; i<gatindex.size(); i++) gatindex[i]=start+i;
  for(unsigned i=0; i<g2l.size(); i++) g2l[i]=-1;
  for(unsigned i=0; i<gatindex.size(); i++) g2l[gatindex[i]]=i;
  if(gatindex.size()<getPntrToOutput(0)->getShape()[0]) shuffledAtoms=1;
}

bool DomainDecomposition::retrieveRequiredAtoms( const std::vector<int>& g2l ) {
  clearTaskLists(); Value* val=getPntrToOutput(0);
  if(!(int(gatindex.size())==val->getShape()[0] && shuffledAtoms==0)) {
     for(const auto & p : actions ) {
         if( !p->isActive() || p->getUnique().empty() ) continue ;
         if( g2l.size()>0 ) {
             for(auto pp=p->getUnique().begin(); pp!=p->getUnique().end(); ++pp) {
                 if( g2l[pp->index()]>=0 ) val->addTaskToCurrentList( (*pp) );
             }
         } else val->addTasksToCurrentList( p->getUnique() );
     }
     return true;
  } else {
     for(const auto & p : actions ) {
         if( !p->isActive() || p->getUnique().empty() ) continue ;
         setupCurrentTaskList(); return true ; 
     }
  }
}

void DomainDecomposition::shareAll() {
  if( dd && shuffledAtoms>0 ) {
      clearTaskLists(); Value* val=copyOutput(0); int natoms = val->getShape()[0];
      for(int i=0; i<natoms; ++i) if( g2l[i]>=0 ) val->addTaskToCurrentList( AtomNumber::index(i) );
      share( false, val->getTaskList() );
  } else { 
      setupCurrentTaskList(); 
      std::set<AtomNumber> allatoms; unsigned natoms=getPntrToOutput(0)->getShape()[0];
      for(unsigned i=0; i<natoms; ++i) allatoms.insert(AtomNumber::index(i));
      share( false, allatoms ); 
  }
}

void DomainDecomposition::share() {
/// We can no longer set the pointers after the share
  bool atomsNeeded=false; for(const auto & pp : inputs) pp->share();
// At first step I scatter all the atoms so as to store their mass and charge
// Notice that this works with the assumption that charges and masses are
// not changing during the simulation!
  if( firststep ) { 
      actions = plumed.getActionSet().select<ActionAtomistic*>(); shareAll();
  } else { atomsNeeded = retrieveRequiredAtoms( g2l ); }

  // Now we retrieve the atom numbers we need
  if( atomsNeeded ) share( getPntrToOutput(0)->performAllTasks(), getPntrToOutput(0)->getTaskList() );
}

void DomainDecomposition::share(const bool& getallatoms, const std::set<AtomNumber>& unique) {
  if( !getallatoms && unique.empty() ) return ;

  // This retrieves what values we need to get
  int ndata=0; std::vector<Value*> values_to_get;
  if(int(gatindex.size())==getPntrToOutput(0)->getShape()[0] && shuffledAtoms==0) {
// faster version, which retrieves all atoms
    for(const auto & ip : inputs) {
      if( (!ip->fixed || firststep) && ip->wasset ) { (ip->mydata)->share_data( 0, getPntrToOutput(0)->getShape()[0], ip->copyOutput(0) ); values_to_get.push_back(ip->copyOutput(0)); ndata++; }
    }
  } else {
    uniq_index.clear();
    uniq_index.reserve(unique.size());
    if(shuffledAtoms>0) {
      for(const auto & p : unique) uniq_index.push_back(g2l[p.index()]);
    }
    for(const auto & ip : inputs) {
      if( (!ip->fixed || firststep) && ip->wasset ) { (ip->mydata)->share_data( unique, uniq_index, ip->copyOutput(0) ); values_to_get.push_back(ip->copyOutput(0)); ndata++; }
    }
  }

  if(dd && shuffledAtoms>0) {
    if(dd.async) {
      for(unsigned i=0; i<dd.mpi_request_positions.size(); i++) dd.mpi_request_positions[i].wait();
      for(unsigned i=0; i<dd.mpi_request_index.size(); i++)     dd.mpi_request_index[i].wait();
    }

    int count=0;
    for(const auto & p : unique) {
      dd.indexToBeSent[count]=p.index(); int dpoint=0;
      for(unsigned i=0;i<values_to_get.size();++i) {
          dd.positionsToBeSent[ndata*count+dpoint]=values_to_get[i]->get(p.index());
          dpoint++;
      }
      count++;
    }

    if(dd.async) {
      asyncSent=true;
      dd.mpi_request_positions.resize(dd.Get_size());
      dd.mpi_request_index.resize(dd.Get_size());
      for(int i=0; i<dd.Get_size(); i++) {
        dd.mpi_request_index[i]=dd.Isend(&dd.indexToBeSent[0],count,i,666);
        dd.mpi_request_positions[i]=dd.Isend(&dd.positionsToBeSent[0],ndata*count,i,667);
      }
    } else {
      const int n=(dd.Get_size());
      std::vector<int> counts(n);
      std::vector<int> displ(n);
      std::vector<int> counts5(n);
      std::vector<int> displ5(n);
      dd.Allgather(count,counts);
      displ[0]=0;
      for(int i=1; i<n; ++i) displ[i]=displ[i-1]+counts[i-1];
      for(int i=0; i<n; ++i) counts5[i]=counts[i]*ndata;
      for(int i=0; i<n; ++i) displ5[i]=displ[i]*ndata;
      dd.Allgatherv(&dd.indexToBeSent[0],count,&dd.indexToBeReceived[0],&counts[0],&displ[0]);
      dd.Allgatherv(&dd.positionsToBeSent[0],ndata*count,&dd.positionsToBeReceived[0],&counts5[0],&displ5[0]);
      int tot=displ[n-1]+counts[n-1];
      for(int i=0; i<tot; i++) {
        int dpoint=0;
        for(unsigned j=0;j<values_to_get.size();++j) {
            values_to_get[j]->set( dd.indexToBeReceived[i], dd.positionsToBeReceived[ndata*i+dpoint] );
            dpoint++;
        }
      }
    }
  }
}

void DomainDecomposition::wait() {
  for(const auto & ip : inputs) ip->dataCanBeSet=false; 

  if(dd && shuffledAtoms>0) {
    int ndata=0; std::vector<Value*> values_to_set;
    for(const auto & ip : inputs) {
        if( (!ip->fixed || firststep) && ip->wasset ) { values_to_set.push_back(ip->copyOutput(0)); ndata++; }
    }

// receive toBeReceived
    if(asyncSent) {
      Communicator::Status status;
      std::size_t count=0;
      for(int i=0; i<dd.Get_size(); i++) {
        dd.Recv(&dd.indexToBeReceived[count],dd.indexToBeReceived.size()-count,i,666,status);
        int c=status.Get_count<int>();
        dd.Recv(&dd.positionsToBeReceived[ndata*count],dd.positionsToBeReceived.size()-ndata*count,i,667);
        count+=c;
      }
      for(int i=0; i<count; i++) {
        int dpoint=0;
        for(unsigned j=0;j<values_to_set.size();++j) {
            values_to_set[j]->set(dd.indexToBeReceived[i], dd.positionsToBeReceived[ndata*i+dpoint] );
            dpoint++;
        }
      }
      asyncSent=false;
    }
  }
}

unsigned DomainDecomposition::getNumberOfForcesToRescale() const {
  return gatindex.size(); 
}

void DomainDecomposition::apply() {
  for(const auto & ip : inputs) {
      if( !(ip->getPntrToValue())->forcesWereAdded() || ip->noforce ) {
          continue; 
      } else if( ip->wasscaled || (int(gatindex.size())==getPntrToOutput(0)->getShape()[0] && shuffledAtoms==0) ) {
          (ip->mydata)->add_force( gatindex, ip->getPntrToValue() );
      } else { (ip->mydata)->add_force( getPntrToOutput(0)->taskList, uniq_index, ip->getPntrToValue() ); }
  }
}

void DomainDecomposition::writeBinary(std::ostream&o) {
  for(const auto & ip : inputs) ip->writeBinary(o); 
}

void DomainDecomposition::readBinary(std::istream&i) {
  for(const auto & ip : inputs) ip->readBinary(i);
}

void DomainDecomposition::broadcastToDomains( Value* val ) {
  if( dd ) dd.Bcast( val->data, 0 );
}

void DomainDecomposition::sumOverDomains( Value* val ) {
  if( dd ) dd.Sum( val->data );
}

const long int& DomainDecomposition::getDdStep() const {
  return ddStep;
}
  
const std::vector<int>& DomainDecomposition::getGatindex() const {
  return gatindex;
}

void DomainDecomposition::createFullList(int* n) {
  if( firststep ) {
    int natoms = getPntrToOutput(0)->getShape()[0];
    *n=natoms; fullList.resize(natoms);
    for(unsigned i=0; i<natoms; i++) fullList[i]=i;
  } else {
     std::vector<int> fake; fake.resize(0); 
     shuffledAtoms=1; retrieveRequiredAtoms( fake ); shuffledAtoms=0;
     fullList.resize(0); fullList.reserve( getPntrToOutput(0)->taskList.size() );
     for(const auto & p : getPntrToOutput(0)->taskList) fullList.push_back(p.index());
     *n=fullList.size();
  }
} 

void DomainDecomposition::getFullList(int**x) {
  if(!fullList.empty()) *x=&fullList[0];
  else *x=NULL;
}

void DomainDecomposition::clearFullList() {  
  fullList.resize(0);
}

}
