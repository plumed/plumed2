#ifndef __PLUMED_PlumedFile_h
#define __PLUMED_PlumedFile_h

#include <cstdio>
#include <vector>
#include <string>
#include <sstream>

namespace PLMD{

class PlumedCommunicator;
class PlumedMain;
class Action;

class PlumedFileBase{
protected:
  class FieldBase{
  public:
    std::string name;
    double value;
    bool constant;
    FieldBase(): value(0.0),constant(false){}
  };

  static const int buflen=10000;
  FILE* fp;
  PlumedCommunicator* comm;
  PlumedMain* plumed;
  Action* action;
  bool cloned;
  PlumedFileBase();
public:
  PlumedFileBase& link(FILE*);
  PlumedFileBase& link(PlumedMain&);
  PlumedFileBase& link(PlumedCommunicator&);
  PlumedFileBase& link(Action&);
  PlumedFileBase& flush();
  PlumedFileBase& open(const std::string&name,const std::string& mode);
  void        close();
  virtual ~PlumedFileBase();
  static void test();
};

class PlumedOFile:
public virtual PlumedFileBase{
  PlumedOFile* linked;
  char* buffer;
  class Field:
  public FieldBase{
  public:
    bool set;
    std::string fmt;
    Field(): set(false) {}
  };
  size_t llwrite(const char*,size_t);
  bool fieldChanged;
  std::string fieldFmt;
  std::vector<Field> fields;
  std::string linePrefix;
  std::ostringstream oss;
public:
  PlumedOFile();
  ~PlumedOFile();
  using PlumedFileBase::link;
  PlumedOFile& link(PlumedOFile&);
  PlumedOFile& setLinePrefix(const std::string&);
  PlumedOFile& addField(const std::string&);
  PlumedOFile& addField(const std::string&,double);
  PlumedOFile& clearFields();
  PlumedOFile& fmtFields(const std::string&);
  PlumedOFile& fmtField(const std::string&,const std::string&);
  PlumedOFile& printField(const std::string&,double);
  PlumedOFile& printField();
  int printf(const char*fmt,...);
  template <class T>
  friend PlumedOFile& operator<<(PlumedOFile&,const T &);
};

class PlumedIFile:
public virtual PlumedFileBase{
  class Field:
  public FieldBase{
  public:
    bool read;
    Field(): read(false) {}
  };
  size_t llread(char*,size_t);
  std::vector<std::string> fields;
public:
  PlumedIFile();
  ~PlumedIFile();
  PlumedIFile& scanFieldList(std::vector<std::string>&);
  PlumedIFile& scanField(const std::string&,double&);
  PlumedIFile& scanField();
  PlumedIFile& getline(std::string&);
};

class PlumedFile:
public PlumedOFile,
public PlumedIFile
{
};

/// Write using << syntax
template <class T>
PlumedOFile& operator<<(PlumedOFile&of,const T &t){
  of.oss<<t;
  of.printf("%s",of.oss.str().c_str());
  of.oss.str("");
  return of;
}


}

#endif
