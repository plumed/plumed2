#include "plumed/core/ModuleMap.h"
#include "plumed/core/ActionRegister.h"

#include <fstream>


int main(){
  const auto actionList=PLMD::actionRegister().getActionNames();
  //just checking that all the default-registered actions are in the module map;
  //output should be empty, the diff will be clear
  //this test may break in case of adding a new LOAD mechanism on startup
  std::ofstream file("output");
  for (const auto & name: actionList){
    if(PLMD::getModuleMap().count(name)==0){
      file << name << " is not in the ModuleMap\n";
    }
  }
  return 0;
}
