#include "plumed/tools/File.h"
#include "plumed/tools/Exception.h"
#include <cstdio>
#include <fstream>
#include <string>

using namespace PLMD;

// Checks IFile::getline: unix and dos newlines, very long lines, a last line
// without newline (allowNoEOL or not), and an incomplete line that is read again
// once it has been completed (as for the hills files of multiple walkers).
void readAll(std::ofstream & out,const std::string & name,bool noeol) {
  IFile in;
  if(noeol) {
    in.allowNoEOL();
  }
  in.open(name);
  std::string line;
  unsigned n=0;
  while(in.getline(line)) {
    n++;
    if(line.length()>40) {
      out<<"  line "<<n<<": "<<line.length()<<" characters, starting with "<<line.substr(0,10)
         <<" and ending with "<<line.substr(line.length()-10)<<"\n";
    } else {
      out<<"  line "<<n<<": ["<<line<<"]\n";
    }
  }
  in.close();
}

int main() {
  std::ofstream out("output");
  {
    FILE* f=std::fopen("unix.txt","w");
    std::fprintf(f,"first line\n\nthird line\n");
    std::string longline;
    for(unsigned i=0; i<20000; i++) {
      longline+=char('a'+i%26);
    }
    std::fprintf(f,"%s\n",longline.c_str());
    std::fprintf(f,"last line without newline");
    std::fclose(f);
  }
  out<<"unix file\n";
  readAll(out,"unix.txt",false);
  out<<"unix file, allowNoEOL\n";
  readAll(out,"unix.txt",true);
  {
    FILE* f=std::fopen("dos.txt","wb");
    std::fprintf(f,"first line\r\nsecond line\r\n");
    std::fclose(f);
  }
  out<<"dos file\n";
  readAll(out,"dos.txt",false);
  {
    // an incomplete line is not returned, and it is read again when completed
    FILE* f=std::fopen("growing.txt","w");
    std::fprintf(f,"complete line 1\nincomplete l");
    std::fflush(f);
    IFile in;
    in.open("growing.txt");
    std::string line;
    out<<"growing file\n";
    while(in.getline(line)) {
      out<<"  ["<<line<<"]\n";
    }
    out<<"  end of file reached\n";
    std::fprintf(f,"ine 2\ncomplete line 3\n");
    std::fflush(f);
    in.reset(false);
    while(in.getline(line)) {
      out<<"  ["<<line<<"]\n";
    }
    out<<"  end of file reached\n";
    in.close();
    std::fclose(f);
  }
  {
    // the same, with gzipped files
    OFile o;
    o.open("compressed.gz");
    o.printf("line 1\nline 2\n");
    o.close();
    out<<"gzipped file\n";
    readAll(out,"compressed.gz",false);
  }
  return 0;
}
