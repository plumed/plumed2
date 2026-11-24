#include "plumed/wrapper/Plumed.h"
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <exception>
#include <vector>

using PLMD::Plumed;

void init(Plumed& p,const std::string& log) {
  p.cmd("setLogFile",log.c_str());
  p.cmd("init");
}
void line(Plumed& p,const std::string& input) {
  p.cmd("readInputLine",input.c_str());
}
const std::string source="LOAD FILE={source dir/Extension.cpp}";
const std::string extra="special ' $(touch injected)";
void check(Plumed& p,const std::string& dir,int value) {
  line(p,source+" INCLUDE_PATH={"+dir+":"+extra+"}");
  line(p,"CHECK_INCLUDE EXPECTED="+std::to_string(value));
}
void reject(Plumed& p,const std::string& input) {
  try {
    line(p,input);
  } catch(const std::exception&) {
    return;
  }
  throw std::runtime_error("Input should have failed: "+input);
}
std::string environment(const char* name) {
  const auto* value=std::getenv(name);
  return value ? std::string("set:")+value : "unset";
}
int main() {
  try {
    const std::vector<const char*> variables{
      "CPLUS_INCLUDE_PATH","PLUMED_MKLIB_CFLAGS","PLUMED_ROOT","PLUMED_VERSION",
      "PLUMED_HTMLDIR","PLUMED_INCLUDEDIR","PLUMED_PROGRAM_NAME","PLUMED_IS_INSTALLED"
    };
    std::vector<std::string> before;
    for(const auto* name : variables) {
      before.push_back(environment(name));
    }
    Plumed first,second;
    init(first,"first.log");
    init(second,"second.log");
    line(first,"LOAD FILE=./prebuilt." __PLUMED_SOEXT);
    line(first,"CHECK_INCLUDE EXPECTED=13");
    check(first,"first headers",13);
    check(second,"second headers",23);
    // Both images remain alive while another build uses an absolute include path.
    const char* cwd=std::getenv("PWD");
    if(!cwd) {
      throw std::runtime_error("PWD is required by the test");
    }
    check(first,std::string(cwd)+"/second headers",23);
    // A header edit must be visible on the next opted-in load.
    std::ofstream("first headers/value.h") << "#define INCLUDE_VALUE 30\n";
    check(second,"first headers",33);
    check(first,"first headers:second headers",33);
    check(second,"second headers:first headers",23);
    std::exception_ptr failures[2];
    std::thread threads[2];
    for(unsigned i=0; i<2; ++i) {
      threads[i]=std::thread([i,&failures]() {
        try {
          Plumed p;
          init(p,"thread"+std::to_string(i)+".log");
          check(p,i==0 ? "first headers" : "second headers",i==0 ? 33 : 23);
        } catch(...) {
          failures[i]=std::current_exception();
        }
      });
    }
    for(auto& t : threads) {
      t.join();
    }
    for(const auto& failure : failures)
      if(failure) {
        std::rethrow_exception(failure);
      }
    Plumed invalid;
    init(invalid,"invalid.log");
    reject(invalid,source+" INCLUDE_PATH={}");
    reject(invalid,source+" INCLUDE_PATH={first headers:}");
    reject(invalid,"LOAD FILE=missing.so INCLUDE_PATH=first");
    reject(invalid,source+" INCLUDE_PATH=missing");
    // A fresh object cannot inherit include settings from any earlier LOAD.
    reject(invalid,source);
    for(unsigned i=0; i<variables.size(); ++i) {
      if(environment(variables[i])!=before[i]) {
        throw std::runtime_error("LOAD changed the process environment");
      }
    }
    if(std::ifstream("injected") || std::ifstream("injected_root")) {
      throw std::runtime_error("Shell interpreted a path");
    }
    std::ofstream("result") << "OK\n";
  } catch(const std::exception& e) {
    std::ofstream("result") << e.what() << '\n';
    return 1;
  }
}
