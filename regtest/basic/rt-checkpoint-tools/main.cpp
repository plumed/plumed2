#include "plumed/tools/Checkpoint.h"
#include "plumed/tools/Exception.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cstdio>
#include <functional>

using namespace PLMD;

// Every test prints one line to logfile, which is compared with logfile.reference

static std::ofstream ofs;

static void check(bool ok,const std::string& what) {
  ofs<<(ok?"OK   ":"FAIL ")<<what<<std::endl;
}

static std::string expectError(const std::string& what,const std::function<void()>& f) {
  try {
    f();
  } catch(const Exception& e) {
    std::string msg=e.what();
    // keep only the "checkpoint: ..." line, without the location added by plumed_merror
    const size_t pos=msg.find("checkpoint:");
    if(pos!=std::string::npos) {
      msg=msg.substr(pos);
      msg=msg.substr(0,msg.find('\n'));
    } else {
      msg="(unexpected message) "+msg;
    }
    ofs<<"OK   "<<what<<" -> "<<msg<<std::endl;
    return msg;
  }
  ofs<<"FAIL "<<what<<" (no error)"<<std::endl;
  return "";
}

static std::string buildSample() {
  CheckpointWriter w("2.x-test");
  w.write("step",(long long)1234567890123LL);
  w.write("time",0.125);
  w.write("float",1.5f);
  w.write("unsigned",(unsigned)7);
  w.write("negative",-3);
  w.write("flag",true);
  w.write("name","hello world");
  w.beginSection("a1");
  w.write("x",std::vector<double>{1.0,2.5,-3.0e-10});
  w.write("n",std::vector<int>{1,-2,3,4});
  w.write("m",std::vector<double>{1,2,3,4,5,6},{2,3});
  w.beginSection("inner");
  w.write("s",std::string("nested"));
  w.endSection();
  w.endSection();
  w.beginSection("empty");
  w.endSection();
  return w.finalize();
}

int main() {
  ofs.open("logfile");

  const std::string content=buildSample();
  ofs<<"size "<<content.size()<<std::endl;

  // round trip
  {
    CheckpointReader r(content);
    check(r.getFormatVersion()==1,"format version");
    check(r.getPlumedVersion()=="2.x-test","plumed version");
    check(r.get<long long>("step")==1234567890123LL,"int64 scalar");
    check(r.get<double>("time")==0.125,"double scalar");
    check(r.get<float>("float")==1.5f,"float scalar");
    check(r.get<unsigned>("unsigned")==7,"unsigned scalar");
    check(r.get<int>("negative")==-3,"negative int");
    check(r.get<bool>("flag")==true,"bool scalar");
    check(r.get<std::string>("name")=="hello world","string");
    check(r.has("step") && !r.has("missing"),"has");
    check(r.hasSection("a1") && !r.hasSection("a2"),"hasSection");
    r.beginSection("a1");
    std::vector<double> x;
    r.read("x",x);
    check(x.size()==3 && x[0]==1.0 && x[1]==2.5 && x[2]==-3.0e-10,"double array");
    std::vector<int> n;
    r.read("n",n);
    check(n.size()==4 && n[1]==-2 && n[3]==4,"int array");
    std::vector<unsigned long long> shape=r.getShape("m");
    check(shape.size()==2 && shape[0]==2 && shape[1]==3,"array shape");
    std::vector<float> m;
    r.read("m",m);
    check(m.size()==6 && m[5]==6.0f,"2d array as float");
    check(r.hasSection("inner"),"nested hasSection");
    r.beginSection("inner");
    check(r.get<std::string>("s")=="nested","nested string");
    check(!r.has("step"),"names are relative to the section");
    r.endSection();
    r.endSection();
    check(r.getSections().size()==3,"section list");
    check(r.getEntryNames().size()==11,"entry list");
    check(r.getEntryNames()[7]=="a1/x","entry full name");
  }

  // errors
  expectError("missing entry",[&]() {
    CheckpointReader r(content);
    r.get<double>("missing");
  });
  expectError("type mismatch",[&]() {
    CheckpointReader r(content);
    r.get<long long>("time");
  });
  expectError("array as scalar",[&]() {
    CheckpointReader r(content);
    r.beginSection("a1");
    r.get<double>("x");
  });
  expectError("missing section",[&]() {
    CheckpointReader r(content);
    r.beginSection("nope");
  });
  expectError("truncated file",[&]() {
    CheckpointReader r(content.substr(0,content.size()/2));
  });
  expectError("truncated file, one byte",[&]() {
    CheckpointReader r(content.substr(0,content.size()-1));
  });
  expectError("corrupted byte",[&]() {
    std::string bad=content;
    bad[bad.size()/2]^=0x55;
    CheckpointReader r(bad);
  });
  expectError("trailing garbage",[&]() {
    CheckpointReader r(content+"xx");
  });
  expectError("not a checkpoint",[&]() {
    CheckpointReader r("#! FIELDS time d1\n 0.0 1.0\n");
  });
  expectError("newer format version",[&]() {
    std::string bad=content;
    bad[8]=9; // format version byte
    CheckpointReader r(bad);
  });
  expectError("unclosed section on finalize",[&]() {
    CheckpointWriter w;
    w.beginSection("a");
    w.finalize();
  });
  expectError("invalid name",[&]() {
    CheckpointWriter w;
    w.write("a/b",1.0);
  });

  // entry counting (used to warn about actions with unported state)
  {
    CheckpointWriter w;
    w.beginSection("a");
    check(w.getNumberOfEntriesInSection()==0,"empty section count");
    w.write("x",1.0);
    w.beginSection("b");
    w.endSection();
    check(w.getNumberOfEntriesInSection()==2,"section count");
    w.endSection();
  }

  // file: atomic write with previous copy
  {
    std::remove("test.cpt");
    std::remove("test_prev.cpt");
    CheckpointWriter w("v1");
    w.write("step",1);
    w.writeFile("test.cpt","test_prev.cpt");
    CheckpointWriter w2("v2");
    w2.write("step",2);
    w2.writeFile("test.cpt","test_prev.cpt");
    CheckpointReader r(CheckpointReader::readFile("test.cpt"));
    CheckpointReader rp(CheckpointReader::readFile("test_prev.cpt"));
    check(r.get<int>("step")==2 && rp.get<int>("step")==1,"file and previous file");
    std::ifstream tmp("test.cpt.tmp");
    check(!tmp,"no temporary file left");
  }
  expectError("missing file",[&]() {
    CheckpointReader::readFile("nonexistent.cpt");
  });

  // empty arrays (e.g. no hills yet)
  {
    CheckpointWriter w;
    w.write("d",std::vector<double>());
    w.write("i",std::vector<int>());
    CheckpointReader r(w.finalize());
    std::vector<double> d(3,1.0);
    std::vector<int> i(2,1);
    r.read("d",d);
    r.read("i",i);
    check(d.empty() && i.empty() && r.getShape("d").size()==1 && r.getShape("d")[0]==0,"empty arrays");
  }

  return 0;
}
