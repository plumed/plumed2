#include "PlumedException.h"
#include <cstdio>
#include <cstdlib>

using namespace std;
using namespace PLMD;

std::string PlumedException::format(const std::string&msg,const std::string&file,unsigned line,const std::string&function){
  std::string message;
  message="\n+++ Internal PLUMED error";
  if(file.length()>0){
    char cline[1000];
    sprintf(cline,"%u",line);
    message += "\n+++ file "+file+", line "+cline;
    if(function.length()>0) message +=", function "+function;
  }
  if(msg.length()>0) message +="\n+++ message: "+msg;
  return message;
}


PlumedException::PlumedException()
{
  this->msg=format("","",0,"");
  abortIfExceptionsAreDisabled();
}

PlumedException::PlumedException(const std::string&msg)
{
  this->msg=format(msg,"",0,"");
  abortIfExceptionsAreDisabled();
}

PlumedException::PlumedException(const std::string&msg,const std::string&file,unsigned line,const std::string&function)
{
  this->msg=format(msg,file,line,function);
  abortIfExceptionsAreDisabled();
}

void PlumedException::abortIfExceptionsAreDisabled(){
#if ! defined(__PLUMED_EXCEPTIONS)
  fprintf(stderr,"%s",what());
  fprintf(stderr,"\n");
  std::abort();
#endif
}


