#include "plumed/wrapper/Plumed.h"
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using File = std::unique_ptr<FILE,decltype(&std::fclose)>;

void require(bool ok) {
  if(!ok) {
    throw std::runtime_error("borrowed CL-tool stream check failed");
  }
}

void run(const std::vector<const char*>& args,FILE* input,FILE* output) {
  PLMD::Plumed p;
  int argc=args.size();
  int result=-1;
  p.cmd("CLTool setArgc",&argc);
  p.cmd("CLTool setArgv",args.data());
  p.cmd("CLTool setIn",input);
  p.cmd("CLTool setOut",output);
  p.cmd("CLTool run",&result);
  require(result==0);
}

std::string contents(FILE* stream) {
  require(std::fflush(stream)==0);
  std::rewind(stream);
  std::string result;
  char buffer[256];
  while(std::fgets(buffer,sizeof(buffer),stream)) {
    result+=buffer;
  }
  require(!std::ferror(stream));
  return result;
}

int main() {
  File input(std::tmpfile(),&std::fclose);
  File output(std::tmpfile(),&std::fclose);
  require(bool(input) && bool(output));
  require(std::fputs("#! FIELDS x\n1.25",input.get())>=0);
  std::rewind(input.get());
  run({"plumed","h5convert","-i","-","-o","borrowed.h5"},input.get(),output.get());
  require(contents(output.get()).empty());
  // The input stream must still be open after the Plumed object is destroyed.
  std::rewind(input.get());
  require(contents(input.get())=="#! FIELDS x\n1.25");
  run({"plumed","h5convert","--input","borrowed.h5","--output","-"},input.get(),output.get());
  require(contents(output.get())=="#! FIELDS x\n1.25\n");
  // Reuse the same caller-owned output with a separate invocation.
  std::fseek(output.get(),0,SEEK_END);
  run({"plumed","h5cat","borrowed.h5","borrowed.h5"},input.get(),output.get());
  require(contents(output.get())=="#! FIELDS x\n1.25\n#! FIELDS x\n1.25\n#! FIELDS x\n1.25\n");
  require(std::fseek(output.get(),0,SEEK_END)==0);
  require(std::fputs("still open\n",output.get())>=0);
  require(contents(output.get()).find("still open\n")!=std::string::npos);
  File report(std::fopen("streams","w"),&std::fclose);
  require(bool(report));
  require(std::fputs("caller-owned input and output streams: OK\n",report.get())>=0);
}
