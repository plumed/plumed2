/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2023 The plumed team
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
#include "DomainDecomposition.h"
#include "ActionToPutData.h"
#include "ActionAtomistic.h"
#include "ActionRegister.h"
#include "PlumedMain.h"
#include "ActionSet.h"

#include "small_vector/small_vector.h"
#include "tools/MergeVectorTools.h"

//+PLUMEDOC ANALYSIS DOMAIN_DECOMPOSITION
/*
Pass domain decomposed properties of atoms to PLUMED

Data is passed from the MD code to PLUMED by creating [PUT](PUT.md) actions.
These PUT actions take the data from a void pointer that is passed to PLUMED from the MD code and transfer it into a
PLMD::Value object. Passing a void pointer and using a PUT action to convert the data within it
to a PLMD::Value is also used when the atomic positions, masses and charges are sent to PLUMED. However,
there are some other subtleties for these quantities because MD codes use a domain decomposition and scatter the properties of atoms across
multiple domains. We thus need to use the action DOMAIN_DECOMPOSITION when passing these quantities to make sense of the
information in the void pointers that the MD code passes.

## Creating a DOMAIN_DECOMPOSITION action

A DOMAIN_DECOMPOSITION can be created by using a call to plumed.cmd as follows:

```c++
plumed.cmd("readInputLine dd: DOMAIN_DECOMPOSITION NATOMS=20 VALUE1=vv UNIT1=length PERIODIC1=NO CONSTANT1=False");
```

The DOMAIN_DECOMPOSTION command above creates a PUT action with the label vv. The pointer to the data that needs to be transferred to the PLMD::Value
object that is created by the PUT action is then set by using the command below:

```c++
plumed.cmd("setInputValue vv, &val);
```

Meanwhile, the pointer to the forces that should be modified is passed as follows:

```c++
plumed.cmd("setValueForces vv", force);
```

In other words, pointers to values and forces in the MD code are passed to PUT actions that are created by the DOMAIN_DECOMPOSION in
[the way you pass data to other PUT actions](PUT.md).

The PLMD::Value objects created by a DOMAIN_DECOMPOSITION action are always vectors with the same number of components as atoms in the system. Furthermore, you can create multiple PUT
actions from a single DOMAIN_DECOMPOSITION action. To see why this is useful, consider the following DOMAIN_DECOMPOSITION command:

```plumed
gromacs: DOMAIN_DECOMPOSITION ...
   NATOMS=2000
   VALUE1=myposx UNIT1=length PERIODIC1=NO CONSTANT1=False ROLE1=x
   VALUE2=myposy UNIT2=length PERIODIC2=NO CONSTANT2=False ROLE2=y
   VALUE3=myposz UNIT3=length PERIODIC3=NO CONSTANT3=False ROLE3=z
   VALUE4=myMasses UNIT4=mass PERIODIC4=NO CONSTANT4=True ROLE4=m
   VALUE5=myCharges UNIT5=charge PERIODIC5=NO CONSTANT5=True ROLE5=q
   PBCLABEL=mybox
...
```

This action is created when you call `plumed.cmd("setNatoms",&natoms)` from gromacs. It makes 5 PLMD::Value called posx, posy, posz, Masses and Charges.
These PLMD::Value objects then hold the x, y and z positions of the atoms and the masses and charges of the atoms. It is important to note that this command will
also, create a PBC_ACTION to hold the cell.

The ROLE keywords above are only needed because the five quantities passed by the command above play a unique role within PLUMED. If you pass
some other quantities, this instruction is not required. PLMD::ActionAtomistic searches for atomic positions, masses and charges by looking for PUT Actions
that have these five particular roles and for ActionWithVirtualAtom objects.

## Differences from regular PUT actions

PUT actions created from a DOMAIN_DECOMPOSITION action behave differently from other PUT actions. In particular:

* If a DOMAIN_DECOMPOSITION action creates a PUT action, then the PUT action depends on the DOMAIN_DECOMPOSITION action. ActionToPutData::apply thus does nothing for these PUT actions.
* Similarly, when DOMAIN_DECOMPOSITION actions create PUT actions, data is transferred from the input pointer to the PLMD::Value object by the DOMAIN_DECOMPOSITION action. When ActionToPutData::wait is called for these PUT Actions `wasset=true`, ActionToPutData::wait does nothing.
* Lastly, if a constant PUT action is created by DOMAIN_DECOMPOSITION, the values in the vector are set during the first step of MD rather than during initialisation.

These differences are necessary because PUT actions cannot understand the information in the pointers that are passed from the MD code. These pointers are only understandable if you know
which atoms are on each processor. This information is only passed to the DOMAIN_DECOMPOSITION action. DOMAIN_DECOMPOSITION must translate the information passed from the MD code before it is
passed back on through PLMD::Value objects created by the PUT actions. DOMAIN_DECOMPOSITION thus keeps pointers to all the PUT actions that it creates. It sets the data in these action's PLMD::Value objects
within `DomainDecomposition::share(const std::vector<AtomNumber>& unique)` and `DomainDecomposition::wait().`  The forces on the PUT actions created by the DOMAIN_DECOMPOSITION action are added in `DomainDecomposition::apply()`.

*/
//+ENDPLUMEDOC

namespace PLMD {

namespace {

enum class Option {automatic, no, yes };

Option interpretEnvString(const char* env,const char* str) {
  if(!str) {
    return Option::automatic;
  }
  if(!std::strcmp(str,"yes")) {
    return Option::yes;
  }
  if(!std::strcmp(str,"no")) {
    return Option::no;
  }
  if(!std::strcmp(str,"auto")) {
    return Option::automatic;
  }
  plumed_error()<<"Cannot understand env var "<<env<<"\nPossible values: yes/no/auto\nActual value: "<<str;
}

/// Use unique list of atoms to manipulate forces and positions.
/// A unique list of atoms is used to manipulate forces and positions in MPI parallel runs.
/// In serial runs, this is done if convenient. The code currently contain
/// some heuristic to decide if the unique list should be used or not.
/// An env var can be used to override this decision.
/// export PLUMED_FORCE_UNIQUE=yes  # enforce using the unique list in serial runs
/// export PLUMED_FORCE_UNIQUE=no   # enforce not using the unique list in serial runs
/// export PLUMED_FORCE_UNIQUE=auto # choose heuristically
/// default: auto
Option getenvForceUnique() {
  static const char* name="PLUMED_FORCE_UNIQUE";
  static const auto opt = interpretEnvString(name,std::getenv(name));
  return opt;
}

}

PLUMED_REGISTER_ACTION(DomainDecomposition,"DOMAIN_DECOMPOSITION")

void DomainDecomposition::DomainComms::enable(Communicator& c) {
  on=true;
  Set_comm(c.Get_comm());
  async=Get_size()<10;
  if(std::getenv("PLUMED_ASYNC_SHARE")) {
    std::string s(std::getenv("PLUMED_ASYNC_SHARE"));
    if(s=="yes") {
      async=true;
    } else if(s=="no") {
      async=false;
    } else {
      plumed_merror("PLUMED_ASYNC_SHARE variable is set to " + s + "; should be yes or no");
    }
  }
}

void DomainDecomposition::registerKeywords(Keywords& keys) {
  ActionForInterface::registerKeywords( keys );
  keys.remove("ROLE");
  keys.add("compulsory","NATOMS","the number of atoms across all the domains");
  keys.add("numbered","VALUE","value to create for the domain decomposition");
  keys.add("numbered","UNIT","unit of the ith value that is being passed through the domain decomposition");
  keys.add("numbered","CONSTANT","is the ith value that is being passed through the domain decomposition constant");
  keys.add("numbered","PERIODIC","if the value being passed to plumed is periodic then you should specify the periodicity of the function.  If the value "
           "is not periodic you must state this using PERIODIC=NO.  Positions are passed with PERIODIC=NO even though special methods are used "
           "to deal with pbc");
  keys.add("numbered","ROLE","Get the role this value plays in the code can be x/y/z/m/q to signify that this is x, y, z positions of atoms or masses or charges of atoms");
  keys.add("compulsory","PBCLABEL","Box","the label to use for the PBC action that will be created");
  keys.remove("NUMERICAL_DERIVATIVES");
  keys.setValueDescription("vector","the domain that each atom is within");
}

DomainDecomposition::DomainDecomposition(const ActionOptions&ao):
  Action(ao),
  ActionForInterface(ao),
  ddStep(0),
  shuffledAtoms(0),
  asyncSent(false),
  unique_serial(false) {
  // Read in the number of atoms
  int natoms;
  parse("NATOMS",natoms);
  std::string str_natoms;
  Tools::convert( natoms, str_natoms );
  // Setup the gat index
  gatindex.resize(natoms);
  for(unsigned i=0; i<gatindex.size(); i++) {
    gatindex[i]=i;
  }
  // Now read in the values we are creating here
  for(unsigned i=1;; ++i) {
    std::string valname;
    if( !parseNumbered("VALUE",i,valname) ) {
      break;
    }
    std::string unit;
    parseNumbered("UNIT",i,unit);
    std::string constant;
    parseNumbered("CONSTANT",i,constant);
    std::string period;
    parseNumbered("PERIODIC",i,period);
    std::string ddrole;
    parseNumbered("ROLE",i,ddrole);
    if( constant=="True") {
      plumed.readInputLine( valname + ": PUT FROM_DOMAINS CONSTANT SHAPE=" + str_natoms + " UNIT=" + unit + " PERIODIC=" + period + " ROLE=" + ddrole, !plumed.hasBeenInitialized() );
    } else if( constant=="False") {
      plumed.readInputLine( valname + ": PUT FROM_DOMAINS SHAPE=" + str_natoms + " UNIT=" + unit + " PERIODIC=" + period + " ROLE=" + ddrole, !plumed.hasBeenInitialized() );
    } else {
      plumed_merror("missing information on whether value is constant");
    }
    // And save the list of values that are set from here
    ActionToPutData* ap=plumed.getActionSet().selectWithLabel<ActionToPutData*>(valname);
    ap->addDependency( this );
    inputs.push_back( ap );
  }
  std::string pbclabel;
  parse("PBCLABEL",pbclabel);
  plumed.readInputLine(pbclabel + ": PBC",true);
  // Turn on the domain decomposition
  if( Communicator::initialized() ) {
    Set_comm(comm);
  }
}

void DomainDecomposition::Set_comm(Communicator& newcomm) {
  dd.enable(newcomm);
  setAtomsNlocal(getNumberOfAtoms());
  setAtomsContiguous(0);
  if( dd.Get_rank()!=0 ) {
    ActionToPutData* ap=plumed.getActionSet().selectWithLabel<ActionToPutData*>("Box");
    ap->noforce=true;
  }
}

unsigned DomainDecomposition::getNumberOfAtoms() const {
  if( inputs.size()==0 ) {
    return 0;
  }
  return (inputs[0]->getPntrToValue())->getShape()[0];
}

void DomainDecomposition::resetForStepStart() {
  for(const auto & pp : inputs) {
    pp->resetForStepStart();
  }
}

void DomainDecomposition::setStart( const std::string& argName, const unsigned& sss) {
  for(const auto & pp : inputs) {
    if( pp->getLabel()==argName ) {
      pp->setStart(argName, sss);
      return;
    }
  }
  plumed_error();
}

void DomainDecomposition::setStride( const std::string& argName, const unsigned& sss) {
  for(const auto & pp : inputs) {
    if( pp->getLabel()==argName ) {
      pp->setStride(argName, sss);
      return;
    }
  }
  plumed_error();
}

bool DomainDecomposition::setValuePointer( const std::string& argName, const TypesafePtr & val ) {
  wasset=true;  // Once the domain decomposition stuff is transferred moved the setting of this to where the g2l vector is setup
  for(const auto & pp : inputs) {
    if( pp->setValuePointer( argName, val ) ) {
      return true;
    }
  }
  return false;
}

bool DomainDecomposition::setForcePointer( const std::string& argName, const TypesafePtr & val ) {
  for(const auto & pp : inputs) {
    if( pp->setForcePointer( argName, val ) ) {
      return true;
    }
  }
  return false;
}

void DomainDecomposition::setAtomsNlocal(int n) {
  gatindex.resize(n);
  g2l.resize(getNumberOfAtoms(),-1);
  if(dd) {
// Since these vectors are sent with MPI by using e.g.
// &dd.positionsToBeSent[0]
// we make sure they are non-zero-sized so as to
// avoid errors when doing boundary check
    if(n==0) {
      n++;
    }
    std::size_t nvals = inputs.size(), natoms = getNumberOfAtoms();
    dd.positionsToBeSent.resize(n*nvals,0.0);
    dd.positionsToBeReceived.resize(natoms*nvals,0.0);
    dd.indexToBeSent.resize(n,0);
    dd.indexToBeReceived.resize(natoms,0);
  }
}

void DomainDecomposition::setAtomsGatindex(const TypesafePtr & g,bool fortran) {
  plumed_massert( g || gatindex.size()==0, "NULL gatindex pointer with non-zero local atoms");
  auto gg=g.get<const int*>({gatindex.size()});
  ddStep=getStep();
  if(fortran) {
    for(unsigned i=0; i<gatindex.size(); i++) {
      gatindex[i]=gg[i]-1;
    }
  } else {
    for(unsigned i=0; i<gatindex.size(); i++) {
      gatindex[i]=gg[i];
    }
  }
  for(unsigned i=0; i<g2l.size(); i++) {
    g2l[i]=-1;
  }
  if( gatindex.size()==getNumberOfAtoms() ) {
    shuffledAtoms=0;
    for(unsigned i=0; i<gatindex.size(); i++) {
      if( gatindex[i]!=static_cast<int>(i) ) {
        shuffledAtoms=1;
        break;
      }
    }
  } else {
    shuffledAtoms=1;
  }
  if(dd) {
    dd.Sum(shuffledAtoms);
  }
  for(unsigned i=0; i<gatindex.size(); i++) {
    g2l[gatindex[i]]=i;
  }
  // keep in unique only those atoms that are local
  for(unsigned i=0; i<actions.size(); i++) {
    actions[i]->unique_local_needs_update=true;
  }
  unique.clear();
  forced_unique.clear();
}

void DomainDecomposition::setAtomsContiguous(int start) {
  ddStep=plumed.getStep();
  for(unsigned i=0; i<gatindex.size(); i++) {
    gatindex[i]=start+i;
  }
  for(unsigned i=0; i<g2l.size(); i++) {
    g2l[i]=-1;
  }
  for(unsigned i=0; i<gatindex.size(); i++) {
    g2l[gatindex[i]]=i;
  }
  if(gatindex.size()<getNumberOfAtoms()) {
    shuffledAtoms=1;
  }
  // keep in unique only those atoms that are local
  for(unsigned i=0; i<actions.size(); i++) {
    actions[i]->unique_local_needs_update=true;
  }
  unique.clear();
  forced_unique.clear();
}

void DomainDecomposition::shareAll() {
  unique.clear();
  forced_unique.clear();
  int natoms = getNumberOfAtoms();
  if( dd && shuffledAtoms>0 ) {
    for(int i=0; i<natoms; ++i)
      if( g2l[i]>=0 ) {
        unique.push_back( AtomNumber::index(i) );
      }
  } else {
    unique.resize(natoms);
    for(int i=0; i<natoms; i++) {
      unique[i]=AtomNumber::index(i);
    }
  }
  forced_unique.resize( unique.size() );
  for(unsigned i=0; i<unique.size(); ++i) {
    forced_unique[i] = unique[i];
  }
  share(unique);
}

void DomainDecomposition::share() {
  // We can no longer set the pointers after the share
  bool atomsNeeded=false;
  for(const auto & pp : inputs) {
    pp->share();
  }
  // At first step I scatter all the atoms so as to store their mass and charge
  // Notice that this works with the assumption that charges and masses are
  // not changing during the simulation!
  if( firststep ) {
    actions = plumed.getActionSet().select<ActionAtomistic*>();
    shareAll();
    return;
  }

  if(getenvForceUnique()==Option::automatic) {
    unsigned largest=0;
    for(unsigned i=0; i<actions.size(); i++) {
      if(actions[i]->isActive()) {
        auto l=actions[i]->getUnique().size();
        if(l>largest) {
          largest=l;
        }
      }
    }
    if(largest*2<getNumberOfAtoms()) {
      unique_serial=true;
    } else {
      unique_serial=false;
    }
  } else if(getenvForceUnique()==Option::yes) {
    unique_serial=true;
  } else if(getenvForceUnique()==Option::no) {
    unique_serial=false;
  } else {
    plumed_error();
  }

  if(unique_serial || !(gatindex.size()==getNumberOfAtoms() && shuffledAtoms==0)) {
    for(unsigned i=0; i<actions.size(); i++) {
      if( actions[i]->unique_local_needs_update ) {
        actions[i]->updateUniqueLocal( !(dd && shuffledAtoms>0), g2l );
      }
    }
    // Now reset unique for the new step
    gch::small_vector<const std::vector<AtomNumber>*,32> forced_vectors;
    gch::small_vector<const std::vector<AtomNumber>*,32> nonforced_vectors;
    forced_vectors.reserve(actions.size());
    nonforced_vectors.reserve(actions.size());
    for(unsigned i=0; i<actions.size(); i++) {
      if(actions[i]->isActive()) {
        if(!actions[i]->getUnique().empty()) {
          // unique are the local atoms
          if( actions[i]->actionHasForces() ) {
            forced_vectors.push_back(&actions[i]->getUniqueLocal());
          } else {
            nonforced_vectors.push_back(&actions[i]->getUniqueLocal());
          }
        }
      }
    }
    if( !(forced_vectors.empty() && nonforced_vectors.empty()) ) {
      atomsNeeded=true;
    }
    // Merge the atoms from the atoms that have a force on
    unique.clear();
    forced_unique.clear();
    mergeVectorTools::mergeSortedVectors(forced_vectors,forced_unique);
    // Merge all the atoms
    nonforced_vectors.push_back( &forced_unique );
    mergeVectorTools::mergeSortedVectors(nonforced_vectors,unique);
  } else {
    for(unsigned i=0; i<actions.size(); i++) {
      if(actions[i]->isActive()) {
        if(!actions[i]->getUnique().empty()) {
          atomsNeeded=true;
        }
      }
    }
  }

  // Now we retrieve the atom numbers we need
  if( atomsNeeded ) {
    share( unique );
  }
}

void DomainDecomposition::share(const std::vector<AtomNumber>& uniqueIdxIn) {
  // This retrieves what values we need to get
  int ndata=0;
  std::vector<Value*> values_to_get;
  if(!(gatindex.size()==getNumberOfAtoms() && shuffledAtoms==0)) {
    uniq_index.resize(uniqueIdxIn.size());
    for(unsigned i=0; i<uniqueIdxIn.size(); i++) {
      uniq_index[i]=g2l[uniqueIdxIn[i].index()];
    }
    for(const auto & ip : inputs) {
      if( (!ip->fixed || firststep) && ip->wasset ) {
        (ip->mydata)->share_data( uniqueIdxIn, uniq_index, ip->copyOutput(0) );
        values_to_get.push_back(ip->copyOutput(0));
        ndata++;
      }
    }
  } else if( unique_serial) {
    uniq_index.resize(uniqueIdxIn.size());
    for(unsigned i=0; i<uniqueIdxIn.size(); i++) {
      uniq_index[i]=uniqueIdxIn[i].index();
    }
    for(const auto & ip : inputs) {
      if( (!ip->fixed || firststep) && ip->wasset ) {
        (ip->mydata)->share_data( uniqueIdxIn, uniq_index, ip->copyOutput(0) );
        values_to_get.push_back(ip->copyOutput(0));
        ndata++;
      }
    }
  } else {
// faster version, which retrieves all atoms
    for(const auto & ip : inputs) {
      if( (!ip->fixed || firststep) && ip->wasset ) {
        (ip->mydata)->share_data( 0, getNumberOfAtoms(), ip->copyOutput(0) );
        values_to_get.push_back(ip->copyOutput(0));
        ndata++;
      }
    }
  }

  if(dd && shuffledAtoms>0) {
    if(dd.async) {
      for(unsigned i=0; i<dd.mpi_request_positions.size(); i++) {
        dd.mpi_request_positions[i].wait();
      }
      for(unsigned i=0; i<dd.mpi_request_index.size(); i++) {
        dd.mpi_request_index[i].wait();
      }
    }

    int count=0;
    for(const auto & p : uniqueIdxIn) {
      dd.indexToBeSent[count]=p.index();
      int dpoint=0;
      for(unsigned i=0; i<values_to_get.size(); ++i) {
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
      for(int i=1; i<n; ++i) {
        displ[i]=displ[i-1]+counts[i-1];
      }
      for(int i=0; i<n; ++i) {
        counts5[i]=counts[i]*ndata;
      }
      for(int i=0; i<n; ++i) {
        displ5[i]=displ[i]*ndata;
      }
      dd.Allgatherv(&dd.indexToBeSent[0],count,&dd.indexToBeReceived[0],&counts[0],&displ[0]);
      dd.Allgatherv(&dd.positionsToBeSent[0],ndata*count,&dd.positionsToBeReceived[0],&counts5[0],&displ5[0]);
      int tot=displ[n-1]+counts[n-1];
      for(int i=0; i<tot; i++) {
        int dpoint=0;
        for(unsigned j=0; j<values_to_get.size(); ++j) {
          values_to_get[j]->data[dd.indexToBeReceived[i]] = dd.positionsToBeReceived[ndata*i+dpoint];
          dpoint++;
        }
      }
    }
  }
}

void DomainDecomposition::wait() {
  for(const auto & ip : inputs) {
    ip->dataCanBeSet=false;
  }

  if(dd && shuffledAtoms>0) {
    int ndata=0;
    std::vector<Value*> values_to_set;
    for(const auto & ip : inputs) {
      if( (!ip->fixed || firststep) && ip->wasset ) {
        values_to_set.push_back(ip->copyOutput(0));
        ndata++;
      }
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
      for(unsigned i=0; i<count; i++) {
        int dpoint=0;
        for(unsigned j=0; j<values_to_set.size(); ++j) {
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
  std::vector<unsigned> forced_uniq_index(forced_unique.size());
  if(!(gatindex.size()==getNumberOfAtoms() && shuffledAtoms==0)) {
    for(unsigned i=0; i<forced_unique.size(); i++) {
      forced_uniq_index[i]=g2l[forced_unique[i].index()];
    }
  } else {
    for(unsigned i=0; i<forced_unique.size(); i++) {
      forced_uniq_index[i]=forced_unique[i].index();
    }
  }
  for(const auto & ip : inputs) {
    if( !(ip->getPntrToValue())->forcesWereAdded() || ip->noforce ) {
      continue;
    } else if( ip->wasscaled || (!unique_serial && gatindex.size()==getNumberOfAtoms() && shuffledAtoms==0) ) {
      (ip->mydata)->add_force( gatindex, ip->getPntrToValue() );
    } else {
      (ip->mydata)->add_force( forced_unique, forced_uniq_index, ip->getPntrToValue() );
    }
  }
}

void DomainDecomposition::reset() {
  if( !unique_serial && gatindex.size()==getNumberOfAtoms() && shuffledAtoms==0 ) {
    return;
  }
  // This is an optimisation to ensure that we don't call std::fill over the whole forces
  // array if there are a small number of atoms passed between the MD code and PLUMED
  if( dd && shuffledAtoms>0 ) {
    getAllActiveAtoms( unique );
  }
  for(const auto & ip : inputs) {
    (ip->copyOutput(0))->clearInputForce( unique );
  }
}

void DomainDecomposition::writeBinary(std::ostream&o) {
  for(const auto & ip : inputs) {
    ip->writeBinary(o);
  }
}

void DomainDecomposition::readBinary(std::istream&i) {
  for(const auto & ip : inputs) {
    ip->readBinary(i);
  }
}

void DomainDecomposition::broadcastToDomains( Value* val ) {
  if( dd ) {
    dd.Bcast( val->data, 0 );
  }
}

void DomainDecomposition::sumOverDomains( Value* val ) {
  if( dd && shuffledAtoms>0 ) {
    dd.Sum( val->data );
  }
}

const long int& DomainDecomposition::getDdStep() const {
  return ddStep;
}

const std::vector<int>& DomainDecomposition::getGatindex() const {
  return gatindex;
}

void DomainDecomposition::getAllActiveAtoms( std::vector<AtomNumber>& u ) {
  gch::small_vector<const std::vector<AtomNumber>*,32> vectors;
  vectors.reserve(actions.size());
  for(unsigned i=0; i<actions.size(); i++) {
    if(actions[i]->isActive()) {
      if(!actions[i]->getUnique().empty()) {
        // unique are the local atoms
        vectors.push_back(&actions[i]->getUnique());
      }
    }
  }
  u.clear();
  mergeVectorTools::mergeSortedVectors(vectors,u);
}

void DomainDecomposition::createFullList(const TypesafePtr & n) {
  if( firststep ) {
    int natoms = getNumberOfAtoms();
    n.set(natoms);
    fullList.resize(natoms);
    for(int i=0; i<natoms; i++) {
      fullList[i]=i;
    }
  } else {
// We update here the unique list defined at Atoms::unique.
// This is not very clear, and probably should be coded differently.
// Hopefully this fix the longstanding issue with NAMD.
    getAllActiveAtoms( unique );
    fullList.clear();
    fullList.reserve(unique.size());
    for(const auto & p : unique) {
      fullList.push_back(p.index());
    }
    n.set(int(fullList.size()));
  }
}

void DomainDecomposition::getFullList(const TypesafePtr & x) {
  auto xx=x.template get<const int**>();
  if(!fullList.empty()) {
    *xx=&fullList[0];
  } else {
    *xx=NULL;
  }
}

void DomainDecomposition::clearFullList() {
  fullList.resize(0);
}

}
