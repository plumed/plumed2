/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#include "PlumedMain.h"
#include "ActionAtomistic.h"
#include "ActionPilot.h"
#include "ActionForInterface.h"
#include "ActionRegister.h"
#include "ActionSet.h"
#include "ActionWithValue.h"
#include "ActionWithVirtualAtom.h"
#include "ActionToGetData.h"
#include "ActionToPutData.h"
#include "CLToolMain.h"
#include "ExchangePatterns.h"
#include "GREX.h"
#include "DomainDecomposition.h"
#include "config/Config.h"
#include "tools/Citations.h"
#include "tools/Communicator.h"
#include "tools/DLLoader.h"
#include "tools/Exception.h"
#include "tools/IFile.h"
#include "tools/Log.h"
#include "tools/OpenMP.h"
#include "tools/Tools.h"
#include "tools/Stopwatch.h"
#include "tools/TypesafePtr.h"
#include "lepton/Exception.h"
#include "DataPassingTools.h"
#include "small_vector/small_vector.h"
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <set>
#include <exception>
#include <stdexcept>
#include <ios>
#include <new>
#include <typeinfo>
#include <iostream>
#include <algorithm>
#include <system_error>
#include <future>
#include <memory>
#include <functional>
#include <regex>
#include <any>
#include <optional>
#include <variant>
#include <filesystem>

namespace PLMD {

/// Small utility just used in this file to throw arbitrary exceptions
[[noreturn]] static void testThrow(const char* what) {
  auto words=Tools::getWords(what);
  plumed_assert(words.size()>0);
#define __PLUMED_THROW_NOMSG(type) if(words[0]==#type) throw type()
#define __PLUMED_THROW_MSG(type) if(words[0]==#type) throw type(what)
  __PLUMED_THROW_MSG(PLMD::ExceptionError);
  __PLUMED_THROW_MSG(PLMD::ExceptionDebug);
  __PLUMED_THROW_MSG(PLMD::Exception);
  __PLUMED_THROW_MSG(PLMD::lepton::Exception);
  __PLUMED_THROW_NOMSG(std::bad_exception);
  __PLUMED_THROW_NOMSG(std::bad_array_new_length);
  __PLUMED_THROW_NOMSG(std::bad_alloc);
  __PLUMED_THROW_NOMSG(std::bad_function_call);
  __PLUMED_THROW_NOMSG(std::bad_weak_ptr);
  __PLUMED_THROW_NOMSG(std::bad_cast);
  __PLUMED_THROW_NOMSG(std::bad_typeid);
  __PLUMED_THROW_NOMSG(std::bad_variant_access);
  __PLUMED_THROW_NOMSG(std::bad_optional_access);
  __PLUMED_THROW_NOMSG(std::bad_any_cast);
  __PLUMED_THROW_MSG(std::underflow_error);
  __PLUMED_THROW_MSG(std::overflow_error);
  __PLUMED_THROW_MSG(std::range_error);
  __PLUMED_THROW_MSG(std::runtime_error);
  __PLUMED_THROW_MSG(std::out_of_range);
  __PLUMED_THROW_MSG(std::length_error);
  __PLUMED_THROW_MSG(std::domain_error);
  __PLUMED_THROW_MSG(std::invalid_argument);
  __PLUMED_THROW_MSG(std::logic_error);



  if(words[0]=="std::system_error") {
    plumed_assert(words.size()>2);
    int error_code;
    Tools::convert(words[2],error_code);
    if(words[1]=="std::generic_category") {
      throw std::system_error(error_code,std::generic_category(),what);
    }
    if(words[1]=="std::system_category") {
      throw std::system_error(error_code,std::system_category(),what);
    }
    if(words[1]=="std::iostream_category") {
      throw std::system_error(error_code,std::iostream_category(),what);
    }
    if(words[1]=="std::future_category") {
      throw std::system_error(error_code,std::future_category(),what);
    }
  }

  if(words[0]=="std::filesystem::filesystem_error") {
    int error_code;
    plumed_assert(words.size()>2);
    Tools::convert(words[2],error_code);
    std::error_code x_error_code;
    if(words[1]=="std::generic_category") {
      x_error_code=::std::error_code(error_code,::std::generic_category());
    }
    if(words[1]=="std::system_category") {
      x_error_code=::std::error_code(error_code,::std::system_category());
    }
    if(words[1]=="std::iostream_category") {
      x_error_code=::std::error_code(error_code,::std::iostream_category());
    }
    if(words[1]=="std::future_category") {
      x_error_code=::std::error_code(error_code,::std::future_category());
    }

    if(words.size()<4) {
      throw std::filesystem::filesystem_error(what,x_error_code);
    }
    if(words.size()<5) {
      throw std::filesystem::filesystem_error(what,std::filesystem::path(words[3]),x_error_code);
    }
    throw std::filesystem::filesystem_error(what,std::filesystem::path(words[3]),std::filesystem::path(words[4]),x_error_code);
  }

#define __PLUMED_THROW_REGEX(name) if(words[1]=="std::regex_constants::error_" #name) throw std::regex_error(std::regex_constants::error_ ##name)
  if(words[0]=="std::regex_error") {
    plumed_assert(words.size()>1);
    __PLUMED_THROW_REGEX(collate);
    __PLUMED_THROW_REGEX(ctype);
    __PLUMED_THROW_REGEX(escape);
    __PLUMED_THROW_REGEX(backref);
    __PLUMED_THROW_REGEX(brack);
    __PLUMED_THROW_REGEX(paren);
    __PLUMED_THROW_REGEX(brace);
    __PLUMED_THROW_REGEX(badbrace);
    __PLUMED_THROW_REGEX(range);
    __PLUMED_THROW_REGEX(space);
    __PLUMED_THROW_REGEX(badrepeat);
    __PLUMED_THROW_REGEX(complexity);
    __PLUMED_THROW_REGEX(stack);
  }

#define __PLUMED_THROW_FUTURE(name) if(words[1]=="std::future_errc::" #name) throw std::future_error(::std::future_errc::name)
  if(words[0]=="std::future_error") {
    plumed_assert(words.size()>1);
    __PLUMED_THROW_FUTURE(broken_promise);
    __PLUMED_THROW_FUTURE(future_already_retrieved);
    __PLUMED_THROW_FUTURE(promise_already_satisfied);
    __PLUMED_THROW_FUTURE(no_state);
  }

  if(words[0]=="std::ios_base::failure") {
    int error_code=0;
    if(words.size()>2) {
      Tools::convert(words[2],error_code);
    }
    if(words.size()>1 && words[1]=="std::generic_category") {
      throw std::ios_base::failure(what,std::error_code(error_code,std::generic_category()));
    }
    if(words.size()>1 && words[1]=="std::system_category") {
      throw std::ios_base::failure(what,std::error_code(error_code,std::system_category()));
    }
    if(words.size()>1 && words[1]=="std::iostream_category") {
      throw std::ios_base::failure(what,std::error_code(error_code,std::iostream_category()));
    }
    if(words.size()>1 && words[1]=="std::future_category") {
      throw std::ios_base::failure(what,std::error_code(error_code,std::future_category()));
    }
    throw std::ios_base::failure(what);
  }

  if(words[0]=="int") {
    int value=0;
    if(words.size()>1) {
      Tools::convert(words[1],value);
    }
    throw value;
  }

  if(words[0]=="test_nested1") {
    try {
      throw Exception(std::string("inner ")+what);
    } catch(...) {
      try {
        std::throw_with_nested(Exception(std::string("middle ")+what));
      } catch(...) {
        std::throw_with_nested(Exception(std::string("outer ")+what));
      }
    }
  }

  if(words[0]=="test_nested2") {
    try {
      throw std::bad_alloc();
    } catch(...) {
      try {
        std::throw_with_nested(Exception(std::string("middle ")+what));
      } catch(...) {
        std::throw_with_nested(Exception(std::string("outer ")+what));
      }
    }
  }

  if(words[0]=="test_nested3") {
    try {
      throw "inner";
    } catch(...) {
      try {
        std::throw_with_nested(Exception(std::string("middle ")+what));
      } catch(...) {
        std::throw_with_nested(Exception(std::string("outer ")+what));
      }
    }
  }

  plumed_error() << "unknown exception " << what;
}

namespace {
/// This is an internal tool used to count how many PlumedMain objects have been created
/// and if they were correctly destroyed.
/// When using debug options, it leads to a crash
/// Otherwise, it just prints a message
class CountInstances {
  std::atomic<int> counter{};
  // private constructor to avoid direct usage
  CountInstances() noexcept {}
  ~CountInstances() {
    if(counter!=0) {
      std::cerr<<"WARNING: internal inconsistency in allocated PlumedMain instances (" <<counter<< ")\n";
      std::cerr<<"Might be a consequence of incorrectly paired plumed_create/plumed_finalize in the C interface\n";
      std::cerr<<"Or it could be due to incorrect calls to std::exit, without properly destroying all PlumedMain objects\n";
#ifndef NDEBUG
      std::cerr<<"This is a debug build, so the warning will make PLUMED abort\n";
      std::abort();
#endif
    }
  }
  static CountInstances & instance() {
    static CountInstances counter;
    return counter;
  }
public:
  /// Only access through these static functions
  /// The first call to increase() ensures the instance is constructed
  /// This should provide the correct construction and destruction order
  /// also in cases where the PlumedMain object is constructed in the
  /// constructor of a static object
  static void increase() noexcept {
    ++instance().counter;
  }
  /// See increase()
  static void decrease() noexcept {
    --instance().counter;
  }
};

}


PlumedMain::PlumedMain():
  datoms_fwd(*this),
// automatically write on log in destructor
  stopwatch_fwd(log),
  actionSet_fwd(*this),
  passtools(DataPassingTools::create(sizeof(double))) {
  passtools->usingNaturalUnits=false;
  increaseReferenceCounter();
  log.link(comm);
  log.setLinePrefix("PLUMED: ");
  // this is at last so as to avoid inconsistencies if an exception is thrown
  CountInstances::increase(); // noexcept
}

// destructor needed to delete forward declarated objects
PlumedMain::~PlumedMain() {
  CountInstances::decrease();
}

/////////////////////////////////////////////////////////////
//  MAIN INTERPRETER

#define CHECK_INIT(ini,word) plumed_assert(ini)<<"cmd(\"" << word << "\") should be only used after plumed initialization"
#define CHECK_NOTINIT(ini,word) plumed_assert(!(ini))<<"cmd(\"" << word << "\") should be only used before plumed initialization"
#define CHECK_NOTNULL(val,word) plumed_assert(val)<<"NULL pointer received in cmd(\"" << word << "\")"


void PlumedMain::cmd(std::string_view word,const TypesafePtr & val) {

// Enumerate all possible commands:
  enum {
#include "PlumedMainEnum.inc"
  };

// Static object (initialized once) containing the map of commands:
  const static Tools::FastStringUnorderedMap<int> word_map = {
#include "PlumedMainMap.inc"
  };

  try {

    auto ss=stopwatch.startPause();

    gch::small_vector<std::string_view> words;
    Tools::getWordsSimple(words,word);

    unsigned nw=words.size();
    if(nw==0) {
      // do nothing
    } else {
      int iword=-1;
      const auto it=word_map.find(words[0]);
      if(it!=word_map.end()) {
        iword=it->second;
      }

      switch(iword) {
      case cmd_setBox:
        CHECK_INIT(initialized,word);
        CHECK_NOTNULL(val,word);
        setInputValue( "Box", 0, 1, val );
        break;
      case cmd_setPositions:
        CHECK_INIT(initialized,word);
        setInputValue("posx", 0, 3, val );
        setInputValue("posy", 1, 3, val );
        setInputValue("posz", 2, 3, val );
        break;
      case cmd_setMasses:
        CHECK_INIT(initialized,word);
        setInputValue("Masses", 0, 1, val );
        break;
      case cmd_setCharges:
        CHECK_INIT(initialized,word);
        setInputValue("Charges", 0, 1, val);
        break;
      case cmd_setPositionsX:
        CHECK_INIT(initialized,word);
        setInputValue("posx", 0, 1, val);
        break;
      case cmd_setPositionsY:
        CHECK_INIT(initialized,word);
        setInputValue("posy", 0, 1, val);
        break;
      case cmd_setPositionsZ:
        CHECK_INIT(initialized,word);
        setInputValue("posz", 0, 1, val);
        break;
      case cmd_setVirial:
        CHECK_INIT(initialized,word);
        CHECK_NOTNULL(val,word);
        setInputForce("Box",val);
        break;
      case cmd_setEnergy:
        CHECK_INIT(initialized,word);
        CHECK_NOTNULL(val,word);
        if( name_of_energy!="" ) {
          setInputValue( name_of_energy, 0, 1,  val );
        }
        break;
      case cmd_setForces:
        CHECK_INIT(initialized,word);
        setInputForce("posx",val);
        setInputForce("posy",val);
        setInputForce("posz",val);
        break;
      case cmd_setForcesX:
        CHECK_INIT(initialized,word);
        setInputForce("posx",val);
        break;
      case cmd_setForcesY:
        CHECK_INIT(initialized,word);
        setInputForce("posy",val);
        break;
      case cmd_setForcesZ:
        CHECK_INIT(initialized,word);
        setInputForce("posz",val);
        break;
      case cmd_calc:
        CHECK_INIT(initialized,word);
        calc();
        break;
      case cmd_prepareDependencies:
        CHECK_INIT(initialized,word);
        prepareDependencies();
        break;
      case cmd_shareData:
        CHECK_INIT(initialized,word);
        shareData();
        break;
      case cmd_prepareCalc:
        CHECK_INIT(initialized,word);
        prepareCalc();
        break;
      case cmd_performCalc:
        CHECK_INIT(initialized,word);
        performCalc();
        break;
      case cmd_performCalcNoUpdate:
        CHECK_INIT(initialized,word);
        performCalcNoUpdate();
        break;
      case cmd_performCalcNoForces:
        CHECK_INIT(initialized,word);
        performCalcNoForces();
        break;
      case cmd_update:
        CHECK_INIT(initialized,word);
        update();
        break;
      case cmd_setStep:
        CHECK_INIT(initialized,word);
        CHECK_NOTNULL(val,word);
        step=val.get<int>();
        startStep();
        break;
      case cmd_setStepLong:
        CHECK_INIT(initialized,word);
        CHECK_NOTNULL(val,word);
        step=val.get<long int>();
        startStep();
        break;
      case cmd_setStepLongLong:
        CHECK_INIT(initialized,word);
        CHECK_NOTNULL(val,word);
        step=val.get<long long int>();
        startStep();
        break;
      case cmd_setValue: {
        CHECK_INIT(initialized,words[0]);
        plumed_assert(nw==2);
        setInputValue( std::string(words[1]), 0, 1, val );
      }
      break;
      /* ADDED WITH API=7 */
      case cmd_setValueForces: {
        CHECK_INIT(initialized,words[0]);
        plumed_assert(nw==2);
        setInputForce( std::string(words[1]), val );
      }
      break;
      // words used less frequently:
      case cmd_setAtomsNlocal:
        CHECK_INIT(initialized,word);
        CHECK_NOTNULL(val,word);
        for(const auto & pp : inputs ) {
          plumed_assert(pp);
          DomainDecomposition* dd=pp->castToDomainDecomposition();
          if( dd ) {
            dd->setAtomsNlocal(val.get<int>());
          }
        }
        break;
      case cmd_setAtomsGatindex:
        CHECK_INIT(initialized,word);
        for(const auto & pp : inputs ) {
          plumed_assert(pp);
          DomainDecomposition* dd=pp->castToDomainDecomposition();
          if( dd ) {
            dd->setAtomsGatindex(val,false);
          }
        }
        break;
      case cmd_setAtomsFGatindex:
        CHECK_INIT(initialized,word);
        for(const auto & pp : inputs ) {
          plumed_assert(pp);
          DomainDecomposition* dd=pp->castToDomainDecomposition();
          if( dd ) {
            dd->setAtomsGatindex(val,true);
          }
        }
        break;
      case cmd_setAtomsContiguous:
        CHECK_INIT(initialized,word);
        CHECK_NOTNULL(val,word);
        for(const auto & pp : inputs ) {
          plumed_assert(pp);
          DomainDecomposition* dd=pp->castToDomainDecomposition();
          if( dd ) {
            dd->setAtomsContiguous(val.get<int>());
          }
        }
        break;
      case cmd_createFullList:
        CHECK_INIT(initialized,word);
        CHECK_NOTNULL(val,word);
        for(const auto & pp : inputs ) {
          plumed_assert(pp);
          DomainDecomposition* dd=pp->castToDomainDecomposition();
          if( dd ) {
            dd->createFullList(val);
          }
        }
        break;
      case cmd_getFullList: {
        CHECK_INIT(initialized,word);
        CHECK_NOTNULL(val,word);
        unsigned nlists=0;
        for(const auto & pp : inputs ) {
          plumed_assert(pp);
          DomainDecomposition* dd=pp->castToDomainDecomposition();
          if( dd ) {
            dd->getFullList(val);
            nlists++;
          }
        }
        plumed_assert( nlists==1 );
      }
      break;
      case cmd_clearFullList:
        CHECK_INIT(initialized,word);
        for(const auto & pp : inputs ) {
          plumed_assert(pp);
          DomainDecomposition* dd=pp->castToDomainDecomposition();
          if( dd ) {
            dd->clearFullList();
          }
        }
        break;
      /* ADDED WITH API==6 */
      case cmd_getDataRank: {
        CHECK_INIT(initialized,words[0]);
        plumed_assert(nw==2 || nw==3);
        std::string vtype="";
        if( nw==3 ) {
          vtype=" TYPE="+std::string(words[2]);
        }
        readInputLine( "grab_" + std::string(words[1]) + ": GET ARG=" + std::string(words[1]) + vtype );
        ActionToGetData* as=actionSet.selectWithLabel<ActionToGetData*>("grab_"+std::string(words[1]));
        plumed_assert( as );
        as->get_rank( val );
      }
      break;
      /* ADDED WITH API==6 */
      case cmd_getDataShape: {
        CHECK_INIT(initialized,std::string(words[0]));
        ActionToGetData* as1=actionSet.selectWithLabel<ActionToGetData*>("grab_"+std::string(words[1]));
        plumed_assert( as1 );
        as1->get_shape( val );
      }
      break;
      /* ADDED WITH API==6 */
      case cmd_setMemoryForData: {
        CHECK_INIT(initialized,words[0]);
        plumed_assert(nw==2 || nw==3);
        ActionToGetData* as2=actionSet.selectWithLabel<ActionToGetData*>("grab_"+std::string(words[1]));
        plumed_assert( as2 );
        as2->set_memory( val );
      }
      break;
      /* ADDED WITH API==6 */
      case cmd_setErrorHandler: {
        if(val) {
          error_handler=*static_cast<const plumed_error_handler*>(val.get<const void*>());
        } else {
          error_handler.handler=NULL;
        }
      }
      break;
      case cmd_read:
        CHECK_INIT(initialized,word);
        if(val) {
          readInputFile(val.getCString());
        } else {
          readInputFile("plumed.dat");
        }
        break;
      case cmd_readInputLine:
        CHECK_INIT(initialized,word);
        CHECK_NOTNULL(val,word);
        readInputLine(val.getCString());
        break;
      case cmd_readInputLines:
        CHECK_INIT(initialized,word);
        CHECK_NOTNULL(val,word);
        readInputLines(val.getCString());
        break;
      case cmd_clear: {
        CHECK_INIT(initialized,word);
        std::vector<int> natoms;
        for(const auto & pp : inputs ) {
          plumed_assert(pp);
          DomainDecomposition* dd=pp->castToDomainDecomposition();
          if ( dd ) {
            natoms.push_back( dd->getNumberOfAtoms() );
          }
        }
        actionSet.clearDelete();
        inputs.clear();
        for(unsigned i=0; i<natoms.size(); ++i) {
          std::string str_natoms;
          Tools::convert( natoms[i], str_natoms );
          readInputLine( MDEngine + ": DOMAIN_DECOMPOSITION NATOMS=" + str_natoms +
                         " VALUE1=posx UNIT1=length PERIODIC1=NO CONSTANT1=False ROLE1=x" +
                         " VALUE2=posy UNIT2=length PERIODIC2=NO CONSTANT2=False ROLE2=y" +
                         " VALUE3=posz UNIT3=length PERIODIC3=NO CONSTANT3=False ROLE3=z" +
                         " VALUE4=Masses UNIT4=mass PERIODIC4=NO CONSTANT4=True ROLE4=m" +
                         " VALUE5=Charges UNIT5=charge PERIODIC5=NO CONSTANT5=True ROLE5=q");

        }
        setUnits( passtools->usingNaturalUnits, passtools->units );
      }
      break;
      case cmd_getApiVersion:
        CHECK_NOTNULL(val,word);
        val.set(int(10));
        break;
      // commands which can be used only before initialization:
      case cmd_init:
        CHECK_NOTINIT(initialized,word);
        init();
        break;
      case cmd_setRealPrecision:
        CHECK_NOTINIT(initialized,word);
        CHECK_NOTNULL(val,word);
        passtools=DataPassingTools::create(val.get<int>());
        passtools->usingNaturalUnits=false;
        break;
      case cmd_setMDLengthUnits:
        CHECK_NOTINIT(initialized,word);
        CHECK_NOTNULL(val,word);
        passtools->MDUnits.setLength(passtools->MD2double(val));
        break;
      case cmd_setMDChargeUnits:
        CHECK_NOTINIT(initialized,word);
        CHECK_NOTNULL(val,word);
        passtools->MDUnits.setCharge(passtools->MD2double(val));
        break;
      case cmd_setMDMassUnits:
        CHECK_NOTINIT(initialized,word);
        CHECK_NOTNULL(val,word);
        passtools->MDUnits.setMass(passtools->MD2double(val));
        break;
      case cmd_setMDEnergyUnits:
        CHECK_NOTINIT(initialized,word);
        CHECK_NOTNULL(val,word);
        passtools->MDUnits.setEnergy(passtools->MD2double(val));
        break;
      case cmd_setMDTimeUnits:
        CHECK_NOTINIT(initialized,word);
        CHECK_NOTNULL(val,word);
        passtools->MDUnits.setTime(passtools->MD2double(val));
        break;
      case cmd_setNaturalUnits:
        // set the boltzman constant for MD in natural units (kb=1)
        // only needed in LJ codes if the MD is passing temperatures to plumed (so, not yet...)
        // use as cmd("setNaturalUnits")
        CHECK_NOTINIT(initialized,word);
        passtools->usingNaturalUnits=true;
        break;
      case cmd_setNoVirial: {
        CHECK_NOTINIT(initialized,word);
        ActionToPutData* ap=actionSet.selectWithLabel<ActionToPutData*>("Box");
        if( ap ) {
          ap->noforce=true;
        } else {
          ActionForInterface* af = actionSet.selectWithLabel<ActionForInterface*>(MDEngine);
          if( af ) {
            plumed_merror("setNoVirial should be called after number of atoms have been set");
          }
        }
      }
      break;
      case cmd_setPlumedDat:
        CHECK_NOTINIT(initialized,word);
        CHECK_NOTNULL(val,word);
        plumedDat=val.getCString();
        break;
      case cmd_setMPIComm:
        CHECK_NOTINIT(initialized,word);
        comm.Set_comm(val);
        for(const auto & pp : inputs ) {
          pp->Set_comm(comm);
        }
        break;
      case cmd_setMPIFComm:
        CHECK_NOTINIT(initialized,word);
        comm.Set_fcomm(val);
        for(const auto & pp : inputs ) {
          pp->Set_comm(comm);
        }
        break;
      case cmd_setMPImultiSimComm:
        CHECK_NOTINIT(initialized,word);
        multi_sim_comm.Set_comm(val);
        break;
      case cmd_setNatoms: {
        CHECK_NOTINIT(initialized,word);
        CHECK_NOTNULL(val,word);
        int natoms = val.get<int>();
        std::string str_natoms;
        Tools::convert( natoms, str_natoms );
        ActionForInterface* dd=actionSet.selectWithLabel<ActionForInterface*>(MDEngine);
        if( !dd && natoms>0 )
          readInputLine( MDEngine + ": DOMAIN_DECOMPOSITION NATOMS=" + str_natoms +  +
                         " VALUE1=posx UNIT1=length PERIODIC1=NO CONSTANT1=False ROLE1=x" +
                         " VALUE2=posy UNIT2=length PERIODIC2=NO CONSTANT2=False ROLE2=y" +
                         " VALUE3=posz UNIT3=length PERIODIC3=NO CONSTANT3=False ROLE3=z" +
                         " VALUE4=Masses UNIT4=mass PERIODIC4=NO CONSTANT4=True ROLE4=m" +
                         " VALUE5=Charges UNIT5=charge PERIODIC5=NO CONSTANT5=True ROLE5=q", true );
      }
      break;
      case cmd_setTimestep: {
        CHECK_NOTINIT(initialized,word);
        CHECK_NOTNULL(val,word);
        ActionToPutData* ts = actionSet.selectWithLabel<ActionToPutData*>("timestep");
        if( !ts ) {
          readInputLine("timestep: PUT UNIT=time PERIODIC=NO CONSTANT", true);
          ts = actionSet.selectWithLabel<ActionToPutData*>("timestep");
        }
        if( !ts->setValuePointer("timestep", val ) ) {
          plumed_error();
        }
        // The following is to avoid extra digits in case the MD code uses floats
        // e.g.: float f=0.002 when converted to double becomes 0.002000000094995
        // To avoid this, we keep only up to 6 significant digits after first one
        if( getRealPrecision()<=4 ) {
          Value* tstepv = ts->copyOutput(0);
          double magnitude=std::pow(10,std::floor(std::log10(tstepv->get())));
          tstepv->set( std::round(tstepv->get()/magnitude*1e6)/1e6*magnitude );
        }
        ts->updateUnits( passtools.get() );
      }
      break;
      /* ADDED WITH API==2 */
      case cmd_setKbT: {
        CHECK_NOTINIT(initialized,word);
        CHECK_NOTNULL(val,word);
        readInputLine("kBT: PUT CONSTANT PERIODIC=NO UNIT=energy", true);
        ActionToPutData* kb = actionSet.selectWithLabel<ActionToPutData*>("kBT");
        if( !kb->setValuePointer("kBT", val ) ) {
          plumed_error();
        }
        kb->updateUnits( passtools.get() );
      }
      break;
      /* ADDED WITH API==3 */
      case cmd_setRestart:
        CHECK_NOTINIT(initialized,word);
        CHECK_NOTNULL(val,word);
        if(val.get<int>()!=0) {
          restart=true;
        }
        break;
      /* ADDED WITH API==4 */
      case cmd_doCheckPoint:
        CHECK_INIT(initialized,word);
        CHECK_NOTNULL(val,word);
        doCheckPoint = false;
        if(val.get<int>()!=0) {
          doCheckPoint = true;
        }
        break;
      /* ADDED WITH API==6 */
      case cmd_setNumOMPthreads:
        CHECK_NOTNULL(val,word);
        {
          auto nt=val.get<unsigned>();
          if(nt==0) {
            nt=1;
          }
          OpenMP::setNumThreads(nt);
        }
        break;
      /* ADDED WITH API==10 */
      case cmd_setGpuDeviceId:
        CHECK_NOTNULL(val,word);
        {
          auto id=val.get<int>();
          if(id>=0) {
            gpuDeviceId=id;
          }
        }
        break;
      /* ADDED WITH API==6 */
      /* only used for testing */
      case cmd_throw:
        CHECK_NOTNULL(val,word);
        testThrow(val.getCString());
      /* ADDED WITH API==10 */
      case cmd_setNestedExceptions:
        CHECK_NOTNULL(val,word);
        if(val.get<int>()!=0) {
          nestedExceptions=true;
        } else {
          nestedExceptions=false;
        }
        break;
      /* STOP API */
      case cmd_setMDEngine:
        CHECK_NOTINIT(initialized,word);
        CHECK_NOTNULL(val,word);
        MDEngine=val.getCString();
        break;
      case cmd_setLog:
        CHECK_NOTINIT(initialized,word);
        log.link(val.get<FILE*>());
        break;
      case cmd_setLogFile:
        CHECK_NOTINIT(initialized,word);
        CHECK_NOTNULL(val,word);
        log.open(val.getCString());
        break;
      // other commands that should be used after initialization:
      case cmd_setStopFlag:
        CHECK_INIT(initialized,word);
        CHECK_NOTNULL(val,word);
        val.get<int*>(); // just check type and discard pointer
        stopFlag=val.copy();
        break;
      case cmd_getExchangesFlag:
        CHECK_INIT(initialized,word);
        CHECK_NOTNULL(val,word);
        exchangePatterns.getFlag(*val.get<int*>()); // note: getFlag changes the value of the reference!
        break;
      case cmd_setExchangesSeed:
        CHECK_INIT(initialized,word);
        CHECK_NOTNULL(val,word);
        exchangePatterns.setSeed(val.get<int>());
        break;
      case cmd_setNumberOfReplicas:
        CHECK_INIT(initialized,word);
        CHECK_NOTNULL(val,word);
        exchangePatterns.setNofR(val.get<int>());
        break;
      case cmd_getExchangesList:
        CHECK_INIT(initialized,word);
        CHECK_NOTNULL(val,word);
        exchangePatterns.getList(val);
        break;
      case cmd_runFinalJobs:
        CHECK_INIT(initialized,word);
        runJobsAtEndOfCalculation();
        break;
      case cmd_isEnergyNeeded: {
        CHECK_INIT(initialized,word);
        CHECK_NOTNULL(val,word);
        if( name_of_energy =="" ) {
          val.set(int(0));
        } else {
          ActionToPutData* ap=actionSet.selectWithLabel<ActionToPutData*>(name_of_energy);
          if(ap->isActive()) {
            val.set(int(1));
          } else {
            val.set(int(0));
          }
        }
      }
      break;
      case cmd_getBias:
        CHECK_INIT(initialized,word);
        CHECK_NOTNULL(val,word);
        plumedQuantityToMD( "energy", getBias(), val );
        break;
      case cmd_checkAction:
        CHECK_NOTNULL(val,word);
        plumed_assert(nw==2);
        val.set(int(actionRegister().check(dlloader.getHandles(), std::string(words[1])) ? 1:0));
        break;
      case cmd_setExtraCV: {
        CHECK_NOTNULL(val,word);
        plumed_assert(nw==2);
        if( valueExists(std::string(words[1])) ) {
          setInputValue( std::string(words[1]), 0, 1, val );
        }
      }
      break;
      case cmd_setExtraCVForce: {
        CHECK_NOTNULL(val,word);
        plumed_assert(nw==2);
        if( valueExists(std::string(words[1])) ) {
          setInputForce( std::string(words[1]), val );
        }
      }
      break;
      /* ADDED WITH API==10 */
      case cmd_isExtraCVNeeded:
        CHECK_NOTNULL(val,word);
        plumed_assert(nw==2);
        val.set(int(0));
        for(const auto & p : inputs) {
          if( p->getLabel()==words[1] && p->isActive() ) {
            val.set(int(1));
            break;
          }
        }
        break;
      case cmd_GREX:
        if(!grex) {
          grex=Tools::make_unique<GREX>(*this);
        }
        plumed_massert(grex,"error allocating grex");
        {
          std::string kk=std::string(words[1]);
          for(unsigned i=2; i<words.size(); i++) {
            kk+=" "+std::string(words[i]);
          }
          grex->cmd(kk.c_str(),val);
        }
        break;
      case cmd_CLTool:
        CHECK_NOTINIT(initialized,word);
        if(!cltool) {
          cltool=Tools::make_unique<CLToolMain>();
        }
        {
          std::string kk(words[1]);
          for(unsigned i=2; i<words.size(); i++) {
            kk+=" "+std::string(words[i]);
          }
          cltool->cmd(kk.c_str(),val);
        }
        break;
        break;
      /* ADDED WITH API==7 */
      case cmd_convert: {
        double v;
        plumed_assert(words.size()==2);
        if(Tools::convertNoexcept(std::string(words[1]),v)) {
          passtools->double2MD(v,val);
        }
      }
      break;
      default:
        plumed_error() << "cannot interpret cmd(\"" << word << "\"). check plumed developers manual to see the available commands.";
        break;
      }
    }

  } catch (...) {
    if(log.isOpen()) {
      try {
        log<<"\n################################################################################\n";
        log<<Tools::concatenateExceptionMessages();
        log<<"\n################################################################################\n";
        log.flush();
      } catch(...) {
        // ignore errors here.
        // in any case, we are rethrowing this below
      }
    }
    throw;
  }
}

////////////////////////////////////////////////////////////////////////

void PlumedMain::init() {
// check that initialization just happens once
  initialized=true;
  if(!log.isOpen()) {
    log.link(stdout);
  }
  log<<"PLUMED is starting\n";
  log<<"Version: "<<config::getVersionLong()<<" (git: "<<config::getVersionGit()<<") "
     <<"compiled on " <<config::getCompilationDate() << " at " << config::getCompilationTime() << "\n";
  log<<"Please cite these papers when using PLUMED ";
  log<<cite("The PLUMED consortium, Nat. Methods 16, 670 (2019)");
  log<<cite("Tribello, Bonomi, Branduardi, Camilloni, and Bussi, Comput. Phys. Commun. 185, 604 (2014)");
  log<<"\n";
  log<<"For further information see the PLUMED web page at http://www.plumed.org\n";
  log<<"Root: "<<config::getPlumedRoot()<<"\n";
  log<<"LibraryPath: "<<config::getLibraryPath()<<"\n";
  log<<"For installed feature, see "<<config::getPlumedRoot() + "/src/config/config.txt\n";
  log.printf("Molecular dynamics engine: %s\n",MDEngine.c_str());
  log.printf("Precision of reals: %d\n",passtools->getRealPrecision());
  log.printf("Running over %d %s\n",comm.Get_size(),(comm.Get_size()>1?"nodes":"node"));
  log<<"Number of threads: "<<OpenMP::getNumThreads()<<"\n";
  log<<"Cache line size: "<<OpenMP::getCachelineSize()<<"\n";
  for(const auto & pp : inputs ) {
    plumed_assert(pp);
    DomainDecomposition* dd=pp->castToDomainDecomposition();
    if ( dd ) {
      log.printf("Number of atoms: %d\n",dd->getNumberOfAtoms());
    }
  }
  if(grex) {
    log.printf("GROMACS-like replica exchange is on\n");
  }
  log.printf("File suffix: %s\n",getSuffix().c_str());
  if(plumedDat.length()>0) {
    readInputFile(plumedDat);
    plumedDat="";
  }
  setUnits( passtools->usingNaturalUnits, passtools->units );
  ActionToPutData* ts = actionSet.selectWithLabel<ActionToPutData*>("timestep");
  if(ts) {
    log.printf("Timestep: %f\n",(ts->copyOutput(0))->get());
  }
  ActionToPutData* kb = actionSet.selectWithLabel<ActionToPutData*>("kBT");
  if(kb) {
    log.printf("KbT: %f\n",(kb->copyOutput(0))->get());
  } else {
    log.printf("KbT has not been set by the MD engine\n");
    log.printf("It should be set by hand where needed\n");
  }
  log<<"Relevant bibliography:\n";
  log<<citations;
  log<<"Please read and cite where appropriate!\n";
  log<<"Finished setup\n";
}

void PlumedMain::setupInterfaceActions() {
  inputs.clear();
  std::vector<ActionForInterface*> ap=actionSet.select<ActionForInterface*>();
  for(unsigned i=0; i<ap.size(); ++i) {
    if( ap[i]->getName()=="ENERGY" || ap[i]->getDependencies().size()==0 ) {
      inputs.push_back( ap[i] );
    }
  }
}

void PlumedMain::readInputFile(const std::string & str) {
  plumed_assert(initialized);
  log<<"FILE: "<<str<<"\n";
  IFile ifile;
  ifile.link(*this);
  ifile.open(str);
  ifile.allowNoEOL();
  readInputFile(ifile);
  log<<"END FILE: "<<str<<"\n";
  log.flush();

}

void PlumedMain::readInputFile(IFile & ifile) {
  std::vector<std::string> words;
  while(Tools::getParsedLine(ifile,words) && !endPlumed) {
    readInputWords(words,false);
  }
  endPlumed=false;
  pilots=actionSet.select<ActionPilot*>();
  setupInterfaceActions();
}

void PlumedMain::readInputLine(const std::string & str, const bool& before_init) {
  if( !before_init ) {
    plumed_assert(initialized);
  }
  if(str.empty()) {
    return;
  }
  std::vector<std::string> words=Tools::getWords(str);
  if( before_init ) {
    plumed_assert( citations.empty() );
  }
  citations.clear();
  readInputWords(words,before_init);
  if(!citations.empty()) {
    log<<"Relevant bibliography:\n";
    log<<citations;
    log<<"Please read and cite where appropriate!\n";
  }
}

void PlumedMain::readInputLines(const std::string & str) {
  plumed_assert(initialized);
  if(str.empty()) {
    return;
  }

  log<<"FILE: (temporary)\n";

  // Open a temporary file
  auto fp=std::tmpfile();
  plumed_assert(fp);

  // make sure file is closed (and thus deleted) also if an exception occurs
  auto deleter=[](auto FP) {
    std::fclose(FP);
  };
  std::unique_ptr<FILE,decltype(deleter)> fp_deleter(fp,deleter);

  auto ret=std::fputs(str.c_str(),fp);
  plumed_assert(ret!=EOF);

  std::rewind(fp);

  IFile ifile;
  ifile.link(*this);
  ifile.link(fp);
  ifile.allowNoEOL();

  readInputFile(ifile);
  log<<"END FILE: (temporary)\n";
}

void PlumedMain::readInputWords(const std::vector<std::string> & words, const bool& before_init) {
  if( !before_init ) {
    plumed_assert(initialized);
  }
  if(words.empty()) {
    return;
  } else if(words[0]=="_SET_SUFFIX") {
    plumed_assert(words.size()==2);
    setSuffix(words[1]);
  } else {
    std::vector<std::string> interpreted(words);
    Tools::interpretLabel(interpreted);
    auto action=actionRegister().create(dlloader.getHandles(),ActionOptions(*this,interpreted));
    if(!action) {
      std::string msg;
      msg ="ERROR\nI cannot understand line:";
      for(unsigned i=0; i<interpreted.size(); ++i) {
        msg+=" "+interpreted[i];
      }
      msg+="\nMaybe a missing space or a typo?";
      log << msg;
      log.flush();
      plumed_merror(msg);
    }
    action->checkRead();
    actionSet.emplace_back(std::move(action));
  };

  pilots=actionSet.select<ActionPilot*>();
  setupInterfaceActions();
}

////////////////////////////////////////////////////////////////////////

void PlumedMain::exit(int c) {
  comm.Abort(c);
}

Log& PlumedMain::getLog() {
  return log;
}

void PlumedMain::calc() {
  prepareCalc();
  performCalc();
}

void PlumedMain::prepareCalc() {
  prepareDependencies();
  shareData();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// here we have the main steps in "calc()"
// they can be called individually, but the standard thing is to
// traverse them in this order:
void PlumedMain::prepareDependencies() {

// Stopwatch is stopped when sw goes out of scope
  auto sw=stopwatch.startStop("1 Prepare dependencies");

// activate all the actions which are on step
// activation is recursive and enables also the dependencies
// before doing that, the prepare() method is called to see if there is some
// new/changed dependency (up to now, only useful for dependences on virtual atoms,
// which can be dynamically changed).

// First switch off all actions
  for(const auto & p : actionSet) {
    p->deactivate();
  }

// for optimization, an "active" flag remains false if no action at all is active
  active=false;
  for(unsigned i=0; i<pilots.size(); ++i) {
    if(pilots[i]->onStep()) {
      pilots[i]->activate();
      active=true;
    }
  };

// This stops the driver calculation if there is not a read action
  if( !active && !inputsAreActive() ) {
    stopFlag.set(int(1));
  }

// also, if one of them is the total energy, tell to atoms that energy should be collected
  for(const auto & p : actionSet) {
    if(p->isActive()) {
      if(p->checkNeedsGradients()) {
        p->setOption("GRADIENTS");
      }
    }
  }

}

bool PlumedMain::inputsAreActive() const {
  for(const auto & ip : inputs) {
    if( ip->onStep() ) {
      return true;
    }
  }
  return false;
}

void PlumedMain::shareAll() {
  for(const auto & ip : inputs) {
    ip->shareAll();
  }
}

void PlumedMain::shareData() {
// atom positions are shared (but only if there is something to do)
  if(!active) {
    return;
  }
// Stopwatch is stopped when sw goes out of scope
  auto sw=stopwatch.startStop("2 Sharing data");
  for(const auto & ip : inputs) {
    ip->share();
  }
}

void PlumedMain::performCalcNoUpdate() {
  waitData();
  justCalculate();
  backwardPropagate();
  resetInputs();
}

void PlumedMain::performCalcNoForces() {
  waitData();
  justCalculate();
}

void PlumedMain::performCalc() {
  waitData();
  justCalculate();
  backwardPropagate();
  update();
  resetInputs();
}

void PlumedMain::waitData() {
  if(!active) {
    return;
  }
// Stopwatch is stopped when sw goes out of scope
  auto sw=stopwatch.startStop("3 Waiting for data");
  for(const auto & ip : inputs) {
    if( ip->isActive() && ip->hasBeenSet() ) {
      ip->wait();
    } else if( ip->isActive() ) {
      ip->warning("input requested but this quantity has not been set");
    }
  }
}

void PlumedMain::justCalculate() {
  if(!active) {
    return;
  }
// Stopwatch is stopped when sw goes out of scope
  auto sw=stopwatch.startStop("4 Calculating (forward loop)");
  bias=0.0;
  work=0.0;

  // Check the input actions to determine if we need to calculate constants that
  // depend on masses and charges
  bool firststep=false;
  for(const auto & ip : inputs) {
    if( ip->firststep ) {
      firststep=true;
    }
  }
  if( firststep ) {
    for(const auto & ip : inputs) {
      ip->firststep=false;
    }
  }

  int iaction=0;
// calculate the active actions in order (assuming *backward* dependence)
  for(const auto & pp : actionSet) {
    Action* p(pp.get());
    plumed_assert(p);
    try {
      if(p->isActive()) {
// Stopwatch is stopped when sw goes out of scope.
// We explicitly declare a Stopwatch::Handler here to allow for conditional initialization.
        Stopwatch::Handler swh;
        if(detailedTimers) {
          auto actionNumberLabel=std::to_string(iaction);
          const unsigned m=actionSet.size();
          unsigned k=0;
          unsigned n=1;
          while(n<m) {
            n*=10;
            k++;
          }
          auto spaces=std::string(k-actionNumberLabel.length(),' ');
          swh=stopwatch.startStop("4A " + spaces + actionNumberLabel+" "+p->getLabel());
        }
        ActionWithValue*av=p->castToActionWithValue();
        ActionAtomistic*aa=p->castToActionAtomistic();
        {
          if(av) {
            av->clearInputForces();
          }
          if(av) {
            av->clearDerivatives();
          }
          if( av && av->calculateOnUpdate() ) {
            continue ;
          }
        }
        {
          if(aa)
            if(aa->isActive()) {
              aa->retrieveAtoms();
            }
        }
        if(p->checkNumericalDerivatives()) {
          p->calculateNumericalDerivatives();
        } else {
          p->calculate();
        }
        // This retrieves components called bias
        if(av) {
          bias+=av->getOutputQuantity("bias");
          work+=av->getOutputQuantity("work");
          av->setGradientsIfNeeded();
        }
        // This makes all values that depend on the (fixed) masses and charges constant
        if( firststep ) {
          p->setupConstantValues( true );
        }
        ActionWithVirtualAtom*avv=p->castToActionWithVirtualAtom();
        if(avv) {
          avv->setGradientsIfNeeded();
        }
      }
    } catch(...) {
      plumed_error_nested() << "An error happened while calculating " << p->getLabel();
    }
    iaction++;
  }
}

void PlumedMain::justApply() {
  backwardPropagate();
  update();
  resetInputs();
}

void PlumedMain::backwardPropagate() {
  if(!active) {
    return;
  }
  int iaction=0;
// Stopwatch is stopped when sw goes out of scope
  auto sw=stopwatch.startStop("5 Applying (backward loop)");
// apply them in reverse order
  for(auto pp=actionSet.rbegin(); pp!=actionSet.rend(); ++pp) {
    const auto & p(pp->get());
    if(p->isActive()) {

// Stopwatch is stopped when sw goes out of scope.
// We explicitly declare a Stopwatch::Handler here to allow for conditional initialization.
      Stopwatch::Handler swh;
      if(detailedTimers) {
        auto actionNumberLabel=std::to_string(iaction);
        const unsigned m=actionSet.size();
        unsigned k=0;
        unsigned n=1;
        while(n<m) {
          n*=10;
          k++;
        }
        auto spaces=std::string(k-actionNumberLabel.length(),' ');
        swh=stopwatch.startStop("5A " + spaces + actionNumberLabel+" "+p->getLabel());
      }

      p->apply();
    }
    iaction++;
  }

// Stopwatch is stopped when sw goes out of scope.
// We explicitly declare a Stopwatch::Handler here to allow for conditional initialization.
  Stopwatch::Handler sw1;
  if(detailedTimers) {
    sw1=stopwatch.startStop("5B Update forces");
  }
}

void PlumedMain::update() {
  if(!active) {
    return;
  }

// Stopwatch is stopped when sw goes out of scope
  auto sw=stopwatch.startStop("6 Update");

// update step (for statistics, etc)
  updateFlags.push(true);
  for(const auto & p : actionSet) {
    p->beforeUpdate();
    if(p->isActive() && p->checkUpdate() && updateFlagsTop()) {
      ActionWithValue* av=dynamic_cast<ActionWithValue*>(p.get());
      if( av && av->calculateOnUpdate() ) {
        p->prepare();
        p->calculate();
      } else {
        p->update();
      }
    }
  }
  while(!updateFlags.empty()) {
    updateFlags.pop();
  }
  if(!updateFlags.empty()) {
    plumed_merror("non matching changes in the update flags");
  }
// Check that no action has told the calculation to stop
  if(stopNow) {
    if(stopFlag) {
      stopFlag.set(int(1));
    } else {
      plumed_merror("your md code cannot handle plumed stop events - add a call to plumed.comm(stopFlag,stopCondition)");
    }
  }

// flush by default every 10000 steps
// hopefully will not affect performance
// also if receive checkpointing signal
  if(step%10000==0||doCheckPoint) {
    fflush();
    log.flush();
    for(const auto & p : actionSet) {
      p->fflush();
    }
  }
}

void PlumedMain::load(const std::string& fileName) {
  if(DLLoader::installed()) {
    std::string libName=fileName;
    size_t n=libName.find_last_of(".");
    std::string extension="";
    std::string base=libName;
    if(n!=std::string::npos && n<libName.length()-1) {
      extension=libName.substr(n+1);
    }
    if(n!=std::string::npos && n<libName.length()) {
      base=libName.substr(0,n);
    }

    if(extension=="cpp") {
      libName="./"+base+"."+config::getVersionLong()+"."+config::getSoExt();
// full path command, including environment setup
// this will work even if plumed is not in the execution path or if it has been
// installed with a name different from "plumed"
      std::string cmd=config::getEnvCommand()+" \""+config::getPlumedRoot()+"\"/scripts/mklib.sh -n -o "+libName+" "+fileName;

      if(std::getenv("PLUMED_LOAD_ACTION_DEBUG")) {
        log<<"Executing: "<<cmd;
      } else {
        log<<"Compiling: "<<fileName<<" to "<<libName;
      }

      if(comm.Get_size()>0) {
        log<<" (only on master node)";
      }
      log<<"\n";

      // On MPI process (intracomm), we use Get_rank to make sure a single process does the compilation
      // Processes from multiple replicas might simultaneously do the compilation.
      if(comm.Get_rank()==0) {
        static Tools::CriticalSectionWithKey<std::string> section;
        // This is only locking commands that are running with identical arguments.
        // It is not necessary for correctness (a second mklib would just result in a no op since
        // the library is already there, even if running simultaneously).
        // It however decreases the system load if many threads are used.
        auto s=section.startStop(cmd);
        int ret=std::system(cmd.c_str());
        if(ret!=0) {
          plumed_error() <<"An error happened while executing command "<<cmd<<"\n";
        }
      }
      comm.Barrier();
    } else {
      libName=base+"."+config::getSoExt();
    }

    // If we have multiple threads (each holding a Plumed object), each of them
    // will load the library, but each of them will only see actions registered
    // from the owned library
    auto *p=dlloader.load(libName);
    log<<"Loading shared library "<<libName.c_str()<<" at "<<p<<"\n";
    log<<"Here is the list of new actions\n";
    log<<"\n";
    for(const auto & a : actionRegister().getKeysWithDLHandle(p)) {
      log<<a<<"\n";
    }
    log<<"\n";
  } else {
    plumed_error()<<"While loading library "<< fileName << " loading was not enabled, please check if dlopen was found at configure time";
  }
}

void PlumedMain::resetInputs() {
  for(const auto & ip : inputs) {
    if( ip->isActive() && ip->hasBeenSet() ) {
      ip->reset();
    }
  }
}

double PlumedMain::getBias() const {
  return bias;
}

double PlumedMain::getWork() const {
  return work;
}

FILE* PlumedMain::fopen(const char *path, const char *mode) {
  std::string mmode(mode);
  std::string ppath(path);
  std::string ppathsuf=ppath+getSuffix();
  FILE*fp=std::fopen(const_cast<char*>(ppathsuf.c_str()),const_cast<char*>(mmode.c_str()));
  if(!fp) {
    fp=std::fopen(const_cast<char*>(ppath.c_str()),const_cast<char*>(mmode.c_str()));
  }
  plumed_massert(fp,"file " + ppath + " cannot be found");
  return fp;
}

int PlumedMain::fclose(FILE*fp) {
  return std::fclose(fp);
}

std::string PlumedMain::cite(const std::string&item) {
  return citations.cite(item);
}

void PlumedMain::fflush() {
  for(const auto  & p : files) {
    p->flush();
  }
}

void PlumedMain::insertFile(FileBase&f) {
  files.insert(&f);
}

void PlumedMain::eraseFile(FileBase&f) {
  files.erase(&f);
}

void PlumedMain::stop() {
  stopNow=true;
}

void PlumedMain::runJobsAtEndOfCalculation() {
  for(const auto & p : actionSet) {
    ActionWithValue* av=dynamic_cast<ActionWithValue*>(p.get());
    if( av && av->calculateOnUpdate() ) {
      p->activate();
    }
  }
  for(const auto & p : actionSet) {
    ActionPilot* ap=dynamic_cast<ActionPilot*>(p.get());
    ActionWithValue* av=dynamic_cast<ActionWithValue*>(p.get());
    if( av && av->calculateOnUpdate() ) {
      p->calculate();
    } else if( ap && !av && ap->getStride()==0 ) {
      p->update();
    } else {
      p->runFinalJobs();
    }
  }
}

unsigned PlumedMain::increaseReferenceCounter() noexcept {
  return ++referenceCounter;
}

unsigned PlumedMain::decreaseReferenceCounter() noexcept {
  return --referenceCounter;
}

unsigned PlumedMain::useCountReferenceCounter() const noexcept {
  return referenceCounter;
}

bool PlumedMain::valueExists( const std::string& name ) const {
  for(const auto & p : inputs) {
    if( p->getLabel()==name ) {
      return true;
    }
  }
  return false;
}

void PlumedMain::setInputValue( const std::string& name, const unsigned& start, const unsigned& stride, const TypesafePtr & val ) {
  bool found=false;
  for(const auto & pp : inputs) {
    if( pp->setValuePointer( name, val ) ) {
      pp->setStart(name, start);
      pp->setStride(name, stride);
      found=true;
      break;
    }
  }
  plumed_massert( found, "found no action to set named " + name );
}

void PlumedMain::setInputForce( const std::string& name, const TypesafePtr & val ) {
  bool found=false;
  for(const auto & pp : inputs) {
    if( pp->setForcePointer( name, val ) ) {
      found=true;
      break;
    }
  }
  plumed_massert( found, "found no action to set named " + name );
}

void PlumedMain::setUnits( const bool& natural, const Units& u ) {
  passtools->usingNaturalUnits = natural;
  passtools->units=u;
  std::vector<ActionToPutData*> idata = actionSet.select<ActionToPutData*>();
  for(const auto & ip : idata) {
    ip->updateUnits( passtools.get() );
  }
  for(const auto & p : actionSet ) {
    p->resetStoredTimestep();
  }
}

void PlumedMain::startStep() {
  for(const auto & ip : inputs) {
    ip->resetForStepStart();
  }
}

void PlumedMain::writeBinary(std::ostream&o)const {
  for(const auto & ip : inputs) {
    ip->writeBinary(o);
  }
}

void PlumedMain::readBinary(std::istream&i) {
  for(const auto & ip : inputs) {
    ip->readBinary(i);
  }
}

void PlumedMain::setEnergyValue( const std::string& name ) {
  name_of_energy = name;
}

int PlumedMain::getRealPrecision() const {
  return passtools->getRealPrecision();
}

bool PlumedMain::usingNaturalUnits() const {
  return passtools->usingNaturalUnits;
}

const Units& PlumedMain::getUnits() {
  return passtools->units;
}

PlumedMain::DeprecatedAtoms& PlumedMain::getAtoms() {
  return datoms;
}

void PlumedMain::plumedQuantityToMD( const std::string& unit, const double& eng, const TypesafePtr & m) const {
  passtools->double2MD( eng/passtools->getUnitConversion(unit),m );
}

double PlumedMain::MDQuantityToPLUMED( const std::string& unit, const TypesafePtr & m) const {
  double x=passtools->MD2double(m);
  return x*passtools->getUnitConversion(unit);
}

double PlumedMain::DeprecatedAtoms::getKBoltzmann() const {
  if( plumed.usingNaturalUnits() ) {
    return 1.0;
  }
  return kBoltzmann/plumed.getUnits().getEnergy();
}

double PlumedMain::DeprecatedAtoms::getKbT() const {
  ActionForInterface* kb=plumed.getActionSet().selectWithLabel<ActionForInterface*>("kBT");
  if( kb ) {
    return (kb->copyOutput(0))->get();
  }
  return 0.0;
}

int PlumedMain::DeprecatedAtoms::getNatoms() const {
  std::vector<ActionToPutData*> atpd=plumed.getActionSet().select<ActionToPutData*>();
  for(const auto & pp : atpd ) {
    if( pp->getRole()=="x" ) {
      return (pp->copyOutput(0))->getShape()[0];
    }
  }
  return 0;
}

bool PlumedMain::DeprecatedAtoms::usingNaturalUnits() const {
  return plumed.usingNaturalUnits();
}

void PlumedMain::DeprecatedAtoms::setCollectEnergy(bool b) const {
  plumed.readInputWords( Tools::getWords(plumed.MDEngine + "_energy: ENERGY"), false );
  plumed.setEnergyValue( plumed.MDEngine + "_energy" );
}

double PlumedMain::DeprecatedAtoms::getEnergy() const {
  ActionToPutData* av = plumed.getActionSet().selectWithLabel<ActionToPutData*>( plumed.MDEngine + "_energy" );
  return (av->copyOutput(0))->get();
}

void PlumedMain::activateParseOnlyMode() {
  doParseOnly=true;
}

bool PlumedMain::parseOnlyMode() const {
  return doParseOnly;
}

void PlumedMain::getKeywordsForAction( const std::string& action, Keywords& keys ) const {
  actionRegister().getKeywords( dlloader.getHandles(), action, keys );
}

#ifdef __PLUMED_HAS_PYTHON
// This is here to stop cppcheck throwing an error
#endif

#ifdef __PLUMED_HAS_DLADDR
// This is here to stop cppcheck throwing an error
#endif

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
