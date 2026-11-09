#include "plumed/wrapper/Plumed.h"
#include <cassert>
#include <memory>
#include <string>

using PLMD::Plumed;
std::unique_ptr<Plumed> object() {
  auto p=std::make_unique<Plumed>();
  int atoms=0;
  p->cmd("setNatoms",&atoms);
  p->cmd("setLogFile","/dev/null");
  p->cmd("init");
  return p;
}
bool module(Plumed& p,const char* name) {
  int available=0;
  p.cmd(std::string("checkModule ")+name,&available);
  return available;
}
bool action(Plumed& p,const char* name) {
  int available=0;
  p.cmd(std::string("checkAction ")+name,&available);
  return available;
}
void check(bool dependencyFirst) {
  auto a=object(), b=object(), untouched=object();
  auto loadA=[&] { a->cmd("readInputLine","LOAD MODULE=pamm"); };
  auto loadB=[&] { b->cmd("readInputLine","LOAD MODULE=adjmat"); };
  if(dependencyFirst) {
    loadB();
    loadA();
  } else {
    loadA();
    loadB();
  }
  assert(module(*a,"pamm") && module(*a,"adjmat"));
  assert(module(*b,"adjmat") && !module(*b,"pamm"));
  assert(!module(*untouched,"pamm") && !module(*untouched,"adjmat"));
  assert(action(*a,"PAMM"));
  assert(action(*a,"CONTACT_MATRIX"));
  assert(action(*b,"CONTACT_MATRIX"));
  assert(!action(*b,"PAMM"));
  assert(!action(*untouched,"PAMM"));
  assert(!action(*untouched,"CONTACT_MATRIX"));
  loadA();
  loadB();
  a.reset();
  assert(action(*b,"CONTACT_MATRIX"));
  assert(!action(*b,"PAMM"));
  b.reset();
  assert(!action(*untouched,"PAMM"));
}
int main() {
  check(false);
  check(true);
}
