#include "plumed/tools/H5Table.h"
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace PLMD::h5table;

Value real(double x) {
  Value value;
  value.type=Value::Double;
  value.d=x;
  return value;
}

Value integer(long long x) {
  Value value;
  value.type=Value::Int;
  value.i=x;
  return value;
}

Value text(const std::string& x) {
  Value value;
  value.s=x;
  return value;
}

void first(Writer& writer) {
  writer.addRow({{"x",real(1.25)},{"count",integer(9007199254740993LL)},
    {"label",text("first")}}, {{"description",text("constant with spaces")}});
}

void second(Writer& writer,unsigned n=1) {
  for(unsigned i=0; i<n; ++i) {
    writer.addRow({{"y",real(-2.5)},{"label",text("different-schema")}},
    {{"description",text("different constant")}});
  }
}

void gate(const char* stage) {
  std::puts(stage);
  std::fflush(stdout);
  if(std::getchar()==EOF) {
    throw std::runtime_error("test controller disappeared");
  }
}

int main(int argc,char** argv) {
  if(argc==1) {
    return 0;
  }
  try {
    if(argc!=3) {
      return 2;
    }
    const std::string mode=argv[1];
    const std::string path=argv[2];
    if(mode=="read" || mode=="hold-reader") {
      Reader reader(path);
      if(mode=="hold-reader") {
        gate("opened");
      }
      std::string line;
      while(reader.getline(line)) {
        std::puts(line.c_str());
      }
    } else if(mode=="trim") {
      truncate(path, {1},path+".removed.h5");
    } else if(mode=="worker" || mode=="uncommitted" || mode=="interrupt") {
      auto writer=std::make_unique<Writer>(path,false);
      first(*writer);
      if(mode=="uncommitted") {
        gate("uncommitted");
        return 0;
      }
      writer->flush();
      gate("committed1");
      second(*writer,mode=="interrupt"?4000:1);
      gate("buffered2");
      writer->flush();
      gate("committed2");
      writer.reset();
      gate("closed");
    } else if(mode=="append") {
      Writer writer(path,true);
      second(writer);
      writer.flush();
    } else if(mode=="once") {
      Writer writer(path,false);
      first(writer);
      writer.flush();
      const auto counts=writer.rows();
      if(counts!=std::vector<unsigned long long> {1}) {
        throw std::runtime_error("wrong committed row count");
      }
    } else if(mode=="automatic") {
      Writer writer(path,false);
      for(unsigned i=0; i<4096; ++i) {
        first(writer);
      }
      gate("automatic");
    } else {
      return 2;
    }
  } catch(const std::exception& error) {
    std::fprintf(stderr,"%s\n",error.what());
    return 1;
  }
  return 0;
}
