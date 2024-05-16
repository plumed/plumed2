#include "plumed/tools/Subprocess.h"
#include <fstream>
#include <unistd.h>

using namespace PLMD;

int main(){
  // if processes are not properly terminated and cleaned
  // only a finite number of them will be allowed

  std::ofstream ofs("should_be_empty");

  // the number of iteration required to trigger the bugs fixed in
  // 2d45954bbb26471eaf48c72d9a9f0a2dcd8536cc
  // are, on MacOS:
  // ~ 250 for the missing close call
  // ~ 2500 for the missing waitpid call

  // notice that we are not spawning simultaneous processes.
  // they are called in a sequential way and killed (by destructor)
  for(unsigned i=0;i<10000;i++) {
    try {
      Subprocess sp("yes");
    } catch(...) {
      ofs<<"failed after "<<i<<"\n";
      break;
    }
  }
  return 0;
}

