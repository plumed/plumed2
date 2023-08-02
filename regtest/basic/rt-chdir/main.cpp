#include "plumed/tools/Tools.h"
#include <algorithm>
#include <fstream>
#include <iostream>

void list(std::ostream & os,const std::string & path){
  auto ls=PLMD::Tools::ls(path);
  std::sort(ls.begin(),ls.end());
  for(const auto & l : ls) {
    os<<l<<"\n";
  }
}

int main(){
  std::ofstream os("test");

  list(os,"aa");

  // chdir to existing directory
  {
    auto dc=PLMD::Tools::DirectoryChanger("aa");
    list(os,".");
  }

  // chdir to non existing directory (should fail)
  try {
    auto dc=PLMD::Tools::DirectoryChanger("non_existing");
    os<<"!not catched\n";
  } catch(const std::exception & e) { // generic std::exception as this might change across versions
    os<<"catched\n";
  }

  // chdir to empty directory (should do nothing)
  {
    auto dc=PLMD::Tools::DirectoryChanger("");
    list(os,"aa");
  }

  list(os,"aa");
  return 0;
}
