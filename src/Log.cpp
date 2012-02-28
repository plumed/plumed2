#include <cstdarg>
#include <cstring>
#include "Log.h"
#include "PlumedCommunicator.h"
#include "PlumedException.h"

using namespace std;
using namespace PLMD;

int Log::printf(const char*fmt,...){
  if(comm.Get_rank()>0)return 0;
  int pointer=strlen(buffer);
  va_list arg;
  va_start(arg, fmt);
  int r=vsnprintf(&buffer[pointer],buflen-pointer,fmt,arg);
  va_end(arg);
  plumed_massert(r>-1 && r<buflen-pointer,"error using fmt string " + std::string(fmt));

// Line is buffered until newline, then written with a PLUMED: prefix
  char*p1=buffer;
  char*p2;
  while((p2=strchr(p1,'\n'))){
    *p2='\0';
    fprintf(fp,"PLUMED: %s\n",p1);
    p1=p2+1;
  };
  memmove(buffer,p1,strlen(p1)+1);
  return r;
}

Log::Log(PlumedCommunicator &comm):
  fp(stdout),
  toBeClosed(false),
  comm(comm){
  buffer=new char[buflen];
  for(int i=0;i<buflen;i++) buffer[i]='\0';
}

Log::~Log(){
  if(comm.Get_rank()>0)return;
  if(!fp && toBeClosed)fclose(fp);
  delete [] buffer;
}

void Log::setFile(string str){
  if(comm.Get_rank()>0)return;
  fp=fopen(str.c_str(),"w");
  toBeClosed=true;
}

void Log::set(FILE*f){
  if(comm.Get_rank()>0)return;
  fp=f;
}

void Log::flush(){
  fflush(fp);
}

void Log::printKeyword( const std::string& key, const std::string& documentation ){
  unsigned nlines; nlines=floor( double(documentation.length() / 60) );
  if ( nlines==0 ){
     printf("%23s - %-60s \n", key.c_str(), documentation.c_str() );
  } else {
     std::vector<std::string> w=Tools::getWords( documentation );
     std::vector<unsigned> lens( nlines + 1 );
     unsigned ll=1, nl=0; 
     for(unsigned i=0;i<w.size();++i){
        nl+=w[i].length() + 1;
        if( nl>=ll*60 ){ lens[ll]=nl-1-w[i].length(); ll++; }
     }

     printf("%23s - %-60s \n", key.c_str(), documentation.substr(0,lens[1]).c_str() );
     std::string blank=" ";
     for(unsigned i=1;i<nlines;++i){
        printf("%23s   %-60s \n", blank.c_str(), documentation.substr(lens[i],lens[i+1]-lens[i]).c_str() );
     }
     printf("%23s   %-60s  \n", blank.c_str(), documentation.substr(lens[nlines]).c_str() );
  }
}

