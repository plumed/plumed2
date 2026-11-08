/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2025 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "TrajectoryParser.h"

#include "FileBase.h"
#include "IFile.h"
#include "Tools.h"
#include "xdrfile/xdrfile.h"
#include "xdrfile/xdrfile_trr.h"
#include "xdrfile/xdrfile_xtc.h"

#include <algorithm>
#include <array>
#include <cstring>

// when using molfile plugin
#ifdef __PLUMED_HAS_MOLFILE_PLUGINS
#ifndef __PLUMED_HAS_EXTERNAL_MOLFILE_PLUGINS
/* Use the internal ones. Alternatively:
 *    ifeq (,$(findstring __PLUMED_HAS_EXTERNAL_MOLFILE_PLUGINS,$(CPPFLAGS)))
 *    CPPFLAGS+=-I../molfile
 */
#include "molfile/libmolfile_plugin.h"
#include "molfile/molfile_plugin.h"
using namespace PLMD::molfile;
#else
#include <libmolfile_plugin.h>
#include <molfile_plugin.h>
#endif

static int register_cb(void *v, vmdplugin_t *p);
struct molFilePlugins {
  std::vector<molfile_plugin_t *> plugins{};
  std::map <std::string, unsigned> pluginmap{};
  static molFilePlugins& get() {
    static molFilePlugins me;
    return me;
  }
  std::vector<std::string> getMolfilePluginsnames() {
    std::vector<std::string> toRet(plugins.size());
    for(unsigned i=0; i<plugins.size(); i++) {
      toRet[i] = std::string(plugins[i]->name);
    }
    return toRet;
  }
private:
  molFilePlugins() {};
};
static int register_cb(void *v, vmdplugin_t *p) {
  const auto ret = molFilePlugins::get().pluginmap.insert ( std::pair<std::string,unsigned>(
                     std::string(p->name),
                     molFilePlugins::get().plugins.size()) );

  if (ret.second) {
    molFilePlugins::get().plugins.push_back(reinterpret_cast<molfile_plugin_t *>(p));
#ifdef PLUMED_MOLFILE_PLUGIN_DEBUG
    std::cerr<<"MOLFILE: loading plugin "<<p->name<<" number "<<molFilePlugins::get().plugins.size()-1<<std::endl;
  } else {
    std::cerr<<"MOLFILE: found duplicate plugin for "<<p->name<<" : not inserted "<<std::endl;
#endif
  }
  return VMDPLUGIN_SUCCESS;
}
#endif


namespace PLMD {

/// A frame as read by a fileParser, in double precision and in nm and ps.
/// TrajectoryParser copies it to the buffers of the caller.
struct RawFrame {
  /// input of readAtoms(): see TrajectoryFrame::readBox
  bool readBox=true;
  /// output: the step stored in the file. Set by readHeader() (dlp4) or by readAtoms() (xtc, trr)
  std::optional<long long> step;
  /// output of readHeader(): the time step stored in the file (dlp4)
  std::optional<double> timestep;
  /// output of readAtoms(): the time of the frame stored in the file (xtc, trr)
  std::optional<double> time;
  /// output of readAtoms(): 3*natoms values, set to the size before the call
  std::vector<double> positions;
  /// output of readAtoms(): the cell vectors, one per row (zero if the file has no box)
  std::array<double,9> box{};
  /// output of readAtoms(): natoms values if the file stores masses and charges, empty otherwise
  std::vector<double> masses;
  /// output of readAtoms(): as masses
  std::vector<double> charges;
};

class fileParser {
protected:
  int natoms {-1};
public:
  virtual ~fileParser()=default;
  int nOfAtoms() const {
    return natoms;
  }
  virtual std::optional<std::string> init(std::string_view fmt,
                                          std::string_view fname,
                                          int command_line_natoms=-1)=0;
  virtual std::optional<std::string> init(FILE* /*fileHandle*/) {
    return "this trajectory format cannot be read from standard input";
  };
  /// Read the header of the next frame and set natoms
  virtual std::optional<std::string> readHeader(RawFrame& raw)=0;
  /// Read the atoms of the frame whose header was read last
  virtual std::optional<std::string> readAtoms(RawFrame& raw)=0;
  virtual std::optional<std::string> rewind() {
    return "rewind not (yet) implemented for this kind of trajectory file";
  }
};

namespace {

struct file_deleter {
  void operator()(FILE* fp) {
    std::fclose(fp);
  }
};

using safeFile=std::unique_ptr<std::FILE, file_deleter>;

/// The common part of the parsers of text files, which can also be read from standard input
class textParser : public fileParser {
protected:
  // not owning pointer for "-"mode
  std::FILE* fp{nullptr};
  // owning pointer otherwise
  safeFile fp_owned{nullptr,{}};
public:
  std::optional<std::string> init(std::string_view /*fmt*/,
                                  std::string_view fname,
                                  int /*command_line_natoms*/=-1) override {
    std::string trajectoryFile=std::string(fname);
    fp_owned.reset(std::fopen(trajectoryFile.c_str(),"r"));
    fp = fp_owned.get();
    if(!fp) {
      return "ERROR: Error opening trajectory file "+trajectoryFile;
    }
    return std::nullopt;
  }
  std::optional<std::string> init(FILE* fileHandle) override {
    fp = fileHandle;
    return std::nullopt;
  }
  std::optional<std::string> rewind() override {
    int error = std::fseek(fp,0,SEEK_SET);
    if (error) {
      return "ERROR: Error rewinding trajectory file";
    }
    return std::nullopt;
  }
};

class xyzParser final : public textParser {
public:
  std::optional<std::string> readHeader(RawFrame& /*raw*/) override {
    std::string line;
    if(!Tools::getline(fp,line)) {
      return "EOF";
    }
    std::sscanf(line.c_str(),"%100d",&natoms);
    if(natoms<0) {
      return "cannot read the number of atoms of the xyz file from the line "+line;
    }
    return std::nullopt;
  }

  std::optional<std::string> readAtoms(RawFrame& raw) override {
    std::string line;
    //cell
    if(!Tools::getline(fp,line)) {
      return "premature end of trajectory file";
    }
    if(raw.readBox) {
      std::vector<std::string> words=Tools::getWords(line);
      if(words.size()==3) {
        Tools::convert(words[0],raw.box[0]);
        Tools::convert(words[1],raw.box[4]);
        Tools::convert(words[2],raw.box[8]);
      } else if(words.size()==9) {
        for(unsigned i=0; i<9; i++) {
          Tools::convert(words[i],raw.box[i]);
        }
      } else {
        return "needed box in second line of xyz file";
      }
    }
    //atoms, finally
    for(int i=0; i<natoms; i++) {
      if(!Tools::getline(fp,line)) {
        return "premature end of trajectory file";
      }
      char dummy[1000];
      double* cc=&raw.positions[3*i];
      int ret=std::sscanf(line.c_str(),"%999s %100lf %100lf %100lf",
                          dummy,&cc[0],&cc[1],&cc[2]);
      if(ret!=4) {
        return "cannot read line"+line;
      }
    }
    return std::nullopt;
  }
};

class dlp4Parser final : public textParser {
  int lvl{0};
  std::optional<std::string> readFileHeader() {
    // dlp4 header
    std::string line;
    if(!Tools::getline(fp,line)) {
      return "error reading title";
    }
    if(!Tools::getline(fp,line)) {
      return "error reading atoms";
    }
    // in the original driver this pb is only read here, and then not used anymore
    int pb=1;
    std::sscanf(line.c_str(),"%d %d %d",&lvl,&pb,&natoms);
    return std::nullopt;
  }
public:
  std::optional<std::string> init(std::string_view fmt,
                                  std::string_view fname,
                                  int command_line_natoms=-1) override {
    auto error=textParser::init(fmt,fname,command_line_natoms);
    if(error) {
      return error;
    }
    return readFileHeader();
  }
  std::optional<std::string> init(FILE* fileHandle) override {
    // the header of the file was not read when reading from standard input
    textParser::init(fileHandle);
    return readFileHeader();
  }
  std::optional<std::string> rewind() override {
    // the header of the file is not a frame
    auto error=textParser::rewind();
    if(error) {
      return error;
    }
    return readFileHeader();
  }

  std::optional<std::string> readHeader(RawFrame& raw) override {
    std::string line;
    if(!Tools::getline(fp,line)) {
      return "EOF";
    }
    char xa[9];
    int xb,xc,xd;
    long long step;
    double t;
    auto res = std::sscanf(line.c_str(),"%8s %lld %d %d %d %lf",xa,&step,&xb,&xc,&xd,&t);
    if (res!=6) {
      return "error reading header";
    }
    raw.step=step;
    raw.timestep=t;
    return std::nullopt;
  }

  std::optional<std::string> readAtoms(RawFrame& raw) override {
    std::string line;
    //cell
    if(raw.readBox) {
      if(!Tools::getline(fp,line)) {
        return "error reading vector a of cell";
      }
      std::sscanf(line.c_str(),"%lf %lf %lf",&raw.box[0],&raw.box[1],&raw.box[2]);
      if(!Tools::getline(fp,line)) {
        return "error reading vector b of cell";
      }
      std::sscanf(line.c_str(),"%lf %lf %lf",&raw.box[3],&raw.box[4],&raw.box[5]);
      if(!Tools::getline(fp,line)) {
        return "error reading vector c of cell";
      }
      std::sscanf(line.c_str(),"%lf %lf %lf",&raw.box[6],&raw.box[7],&raw.box[8]);
      for(auto & b : raw.box) {
        //converting to nm
        b/=10.0;
      }
    }
    raw.masses.resize(natoms);
    raw.charges.resize(natoms);
    // reading the frame
    for(int i=0; i<natoms; i++) {
      double cc[3];
      bool ok=Tools::getline(fp,line);
      if(!ok) {
        return "premature end of trajectory file";
      }

      char dummy[9];
      int idummy;
      double m,c;
      auto res=std::sscanf(line.c_str(),"%8s %d %lf %lf",dummy,&idummy,&m,&c);
      if (res != 4) {
        return "cannot read line"+line;
      }
      raw.masses[i]=m;
      raw.charges[i]=c;
      if(!Tools::getline(fp,line)) {
        return "error reading coordinates";
      }
      std::sscanf(line.c_str(),"%lf %lf %lf",&cc[0],&cc[1],&cc[2]);
      raw.positions[3*i]=cc[0]*0.1;
      raw.positions[3*i+1]=cc[1]*0.1;
      raw.positions[3*i+2]=cc[2]*0.1;
      if(lvl>0) {
        if(!Tools::getline(fp,line)) {
          return "error skipping velocities";
        }
      }
      if(lvl>1) {
        if(!Tools::getline(fp,line)) {
          return "error skipping forces";
        }
      }
    }
    return std::nullopt;
  }
};

class groParser final : public textParser {
public:
  std::optional<std::string> readHeader(RawFrame& /*raw*/) override {
    std::string line;
    if(!Tools::getline(fp,line)) {
      return "EOF";
    }
    if(!Tools::getline(fp,line)) {
      return "premature end of trajectory file";
    }
    std::sscanf(line.c_str(),"%100d",&natoms);
    if(natoms<0) {
      return "cannot read the number of atoms of the gro file from the line "+line;
    }
    return std::nullopt;
  }

  std::optional<std::string> readAtoms(RawFrame& raw) override {
    std::string line;
    int ddist=0;
    for(int i=0; i<natoms; i++) {
      if(!Tools::getline(fp,line)) {
        return "premature end of trajectory file";
      }
      // do the gromacs way
      if(!i) {
        //
        // calculate the distance between dots (as in gromacs gmxlib/confio.c, routine get_w_conf )
        //
        const char      *p1, *p2, *p3;
        p1 = std::strchr(line.c_str(), '.');
        if (p1 == NULL) {
          return "seems there are no coordinates in the gro file";
        }
        p2 = std::strchr(&p1[1], '.');
        if (p2 == NULL) {
          return "seems there is only one coordinates in the gro file";
        }
        ddist = p2 - p1;
        p3 = std::strchr(&p2[1], '.');
        if (p3 == NULL) {
          return "seems there are only two coordinates in the gro file";
        }
        if (p3 - p2 != ddist) {
          return "not uniform spacing in fields in the gro file";
        }
      }
      Tools::convert(line.substr(20,ddist),raw.positions[3*i]);
      Tools::convert(line.substr(20+ddist,ddist),raw.positions[3*i+1]);
      Tools::convert(line.substr(20+ddist+ddist,ddist),raw.positions[3*i+2]);
    }
    if(!Tools::getline(fp,line)) {
      return "premature end of trajectory file";
    }
    std::vector<std::string> words=Tools::getWords(line);
    if(words.size()<3) {
      return "cannot understand box format";
    }
    // the order of the gro file: v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)
    const unsigned order[9]= {0,4,8,1,2,3,5,6,7};
    for(unsigned i=0; i<words.size() && i<9; i++) {
      Tools::convert(words[i],raw.box[order[i]]);
    }
    return std::nullopt;
  }
};

struct xd_deleter {
  void operator()(PLMD::xdrfile::XDRFILE* xd) {
    xdrfile::xdrfile_close(xd);
  }
};

enum class xdType {trr,xtc};
template <xdType is>
class xdParser final: public fileParser {
  // owning pointer
  std::unique_ptr<PLMD::xdrfile::XDRFILE, xd_deleter> xd{nullptr,{}};
  std::string trajectoryFile;
  // the positions of the last frame, sized at init
  // (a std::vector<rvec> does not compile, hence std::unique_ptr<rvec[]>)
  std::unique_ptr<xdrfile::rvec[]> pos;
  static constexpr const char* name = (is == xdType::xtc) ? "xtc" : "trr";
  std::string errorMessage(int ret) const {
    std::string msg = (ret>=0 && ret<xdrfile::exdrNR) ? xdrfile::exdr_message[ret] : "unknown error";
    return "cannot read a frame of the "+std::string(name)+" file "+trajectoryFile+" (xdrfile: "+msg+")";
  }
public:
  std::optional<std::string> init(std::string_view /*fmt*/,
                                  std::string_view fname,
                                  int /*command_line_natoms*/=-1) override {
    trajectoryFile=std::string(fname);
    xd.reset(xdrfile::xdrfile_open(trajectoryFile.c_str(),"r"));
    if(!xd) {
      return "ERROR: Error opening trajectory file "+trajectoryFile;
    }
    int ret=xdrfile::exdrOK;
    if constexpr(is == xdType::xtc) {
      ret=xdrfile::read_xtc_natoms(&trajectoryFile[0],&natoms);
    }
    if constexpr(is == xdType::trr) {
      ret=xdrfile::read_trr_natoms(&trajectoryFile[0],&natoms);
    }
    if(ret!=xdrfile::exdrOK || natoms<=0) {
      return "cannot read the number of atoms of the "+std::string(name)+" file "+trajectoryFile;
    }
    pos=Tools::make_unique<xdrfile::rvec[]>(natoms);
    return std::nullopt;
  }
  std::optional<std::string> readHeader(RawFrame& /*raw*/) override {
    // the header is read by readAtoms, together with the coordinates
    return std::nullopt;
  }
  std::optional<std::string> readAtoms(RawFrame& raw) override {
    int localstep;
    float time;
    xdrfile::matrix box;
    int ret=xdrfile::exdrOK;
    if constexpr(is == xdType::xtc) {
      float prec;
      ret=xdrfile::read_xtc(xd.get(),natoms,&localstep,&time,box,pos.get(),&prec);
    }
    if constexpr(is == xdType::trr) {
      float lambda;
      ret=xdrfile::read_trr(xd.get(),natoms,&localstep,&time,&lambda,box,pos.get(),NULL,NULL);
    }
    if(ret==xdrfile::exdrENDOFFILE) {
      return "EOF";
    }
    if(ret!=xdrfile::exdrOK) {
      return errorMessage(ret);
    }
    raw.step=localstep;
    raw.time=time;
    for(unsigned i=0; i<3; i++) {
      for(unsigned j=0; j<3; j++) {
        raw.box[3*i+j]=box[i][j];
      }
    }
    for(int i=0; i<natoms; i++) {
      for(unsigned j=0; j<3; j++) {
        raw.positions[3*i+j]=pos[i][j];
      }
    }
    return std::nullopt;
  }
};

#ifdef __PLUMED_HAS_MOLFILE_PLUGINS
class molfileParser final:  public fileParser {
  molfile::molfile_plugin_t *api=NULL;
  molfile_timestep_t ts_in;
  std::string trajectoryFile;
  std::string trajectory_fmt_str;
  //Initializing ts_in_coords should prevent that ts_in_coords.data() results in UB (and passes codecheck)
  std::vector<float> ts_in_coords{0.0f,0.0f,0.0f};
  struct molfile_deleter {
    molfile::molfile_plugin_t *api=NULL;
    void operator()(void* h) {
      if(h) {
        std::unique_ptr<std::lock_guard<std::mutex>> lck;
        if(api->is_reentrant==VMDPLUGIN_THREADUNSAFE) {
          lck=Tools::molfile_lock();
        }
        api->close_file_read(h);
      }
    }
  };
  std::unique_ptr<void, molfile_deleter> h_in{nullptr,{nullptr}};

public:
  molfileParser() {
    ts_in.coords=NULL;
    ts_in.velocities=NULL;
    ts_in.A=-1; // we use this to check whether cell is provided or not
  }
  std::optional<std::string> init(std::string_view fmt,
                                  std::string_view fname,
                                  int command_line_natoms=-1) override {
    trajectoryFile = std::string(fname);
    ts_in.coords=ts_in_coords.data();
    ts_in.velocities=NULL;
    ts_in.A=-1; // we use this to check whether cell is provided or not
    trajectory_fmt_str=std::string(fmt);
    bool found=false;
    for(unsigned i=0; i<molFilePlugins::get().plugins.size(); ++i) {
      if(fmt == molFilePlugins::get().plugins[i]->name) {
        found=true;
        api = molFilePlugins::get().plugins[i];
      }
    }
    if(!found) {
      return "trajectory format "+trajectory_fmt_str+" is not supported by the molfile plugin";
    }
    std::unique_ptr<std::lock_guard<std::mutex>> lck;
    if(api->is_reentrant==VMDPLUGIN_THREADUNSAFE) {
      lck=Tools::molfile_lock();
    }
    h_in = std::unique_ptr<void, molfile_deleter> {
      api->open_file_read(trajectoryFile.c_str(), trajectory_fmt_str.c_str(), &natoms),
      {api}
    };
    if(!h_in) {
      return "cannot read "+trajectoryFile+" with the molfile "+trajectory_fmt_str+" plugin";
    }
    if(natoms==MOLFILE_NUMATOMS_UNKNOWN) {
      if(command_line_natoms>=0) {
        natoms=command_line_natoms;
      } else {
        return "this file format does not provide number of atoms; use --natoms on the command line";
      }
    }
    ts_in_coords.resize(3*natoms);
    ts_in.coords = ts_in_coords.data();
    return std::nullopt;
  }

  std::optional<std::string> readHeader(RawFrame& /*raw*/) override {
    // the whole frame is read here; readAtoms converts it
    std::unique_ptr<std::lock_guard<std::mutex>> lck;
    if(api->is_reentrant==VMDPLUGIN_THREADUNSAFE) {
      lck=Tools::molfile_lock();
    }
    int rc;
    rc = api->read_next_timestep(h_in.get(), natoms, &ts_in);
    if(rc==MOLFILE_EOF) {
      return "EOF";
    }
    if(rc!=MOLFILE_SUCCESS) {
      return "error reading a frame of "+trajectoryFile+" with the molfile "+trajectory_fmt_str+" plugin";
    }
    return std::nullopt;
  }

  std::optional<std::string> readAtoms(RawFrame& raw) override {
    if(raw.readBox) {
      //ts_in has been read in the previous step
      if(ts_in.A>0.0) { // this is negative if molfile does not provide box
        // info on the cell: convert using pbcset.tcl from pbctools in vmd distribution
        constexpr double r180=180.;
        double cosBC=cos(double(ts_in.alpha)*pi/r180);
        double cosAC=std::cos(double(ts_in.beta)*pi/r180);
        double cosAB=std::cos(double(ts_in.gamma)*pi/r180);
        double sinAB=std::sin(double(ts_in.gamma)*pi/r180);
        double Ax=double(ts_in.A);
        double Bx=double(ts_in.B)*cosAB;
        double By=double(ts_in.B)*sinAB;
        double Cx=double(ts_in.C)*cosAC;
        double Cy=(double(ts_in.C)*double(ts_in.B)*cosBC-Cx*Bx)/By;
        double Cz=std::sqrt(double(ts_in.C)*double(ts_in.C)-Cx*Cx-Cy*Cy);
        raw.box[0]=Ax/10.;
        raw.box[1]=0.;
        raw.box[2]=0.;
        raw.box[3]=Bx/10.;
        raw.box[4]=By/10.;
        raw.box[5]=0.;
        raw.box[6]=Cx/10.;
        raw.box[7]=Cy/10.;
        raw.box[8]=Cz/10.;
      }
    }
    // info on coords
    // the order is xyzxyz...
    for(int i=0; i<3*natoms; i++) {
      raw.positions[i]=double(ts_in.coords[i])/double(10.); //convert to nm
    }
    return std::nullopt;
  }
};
#endif

std::unique_ptr<fileParser> parserFromFMT(TrajectoryParser::trajfmt fmt) {
  using trajfmt = TrajectoryParser::trajfmt;
  switch (fmt) {
  case trajfmt::molfile: {
#ifdef __PLUMED_HAS_MOLFILE_PLUGINS
    return std::make_unique<molfileParser>();
#else
    plumed_assert(false) <<"Plumed has not been compiled with molfile,\n"
                         "it is impossible to read a file with the molfile plugin";
#endif
  }
  case trajfmt::xdr_xtc: {
    return std::make_unique<xdParser<xdType::xtc>>();
  }
  case trajfmt::xdr_trr: {
    return std::make_unique<xdParser<xdType::trr>>();
  }
  case trajfmt::xyz: {
    return std::make_unique<xyzParser>();
  }
  case trajfmt::gro: {
    return std::make_unique<groParser>();
  }
  case trajfmt::dlp4: {
    return std::make_unique<dlp4Parser>();
  }
  default: {
    plumed_assert(false) <<"you should not be here" ;
  }
  }
}
}

void TrajectoryParser::registerKeywords(Keywords& keys) {
  keys.add("atoms","--ixyz","the trajectory in xyz format");
  keys.add("atoms","--igro","the trajectory in gro format");
  keys.add("atoms","--idlp4","the trajectory in DL_POLY_4 format");
  keys.add("atoms","--ixtc","the trajectory in xtc format (xdrfile implementation)");
  keys.add("atoms","--itrr","the trajectory in trr format (xdrfile implementation)");
#ifdef __PLUMED_HAS_MOLFILE_PLUGINS
  MOLFILE_INIT_ALL
  MOLFILE_REGISTER_ALL(NULL, register_cb)
  for(unsigned i=0; i<molFilePlugins::get().plugins.size(); i++) {
    std::string kk="--mf_"+std::string(molFilePlugins::get().plugins[i]->name);
    std::string mm=" molfile: the trajectory in "+std::string(molFilePlugins::get().plugins[i]->name)+" format " ;
    keys.add("atoms",kk,mm);
  }
#endif
}

TrajectoryParser::TrajectoryParser():
  raw(std::make_unique<RawFrame>()) {
}
TrajectoryParser::~TrajectoryParser()=default;
std::vector<std::string> TrajectoryParser::getMolfilePluginsnames() {
#ifdef __PLUMED_HAS_MOLFILE_PLUGINS
  return molFilePlugins::get().getMolfilePluginsnames();
#else
  return {};
#endif
}
std::vector<std::string> TrajectoryParser::trajectoryOptions() {
  return {"xyz","gro","dlp4","xtc","trr"};
}


TrajectoryParser::trajfmt TrajectoryParser::FMTfromString(std::string_view fmt) {
  if(fmt=="xyz") {
    return trajfmt::xyz;
  }
  if(fmt=="gro") {
    return trajfmt::gro;
  }
  if(fmt=="dlp4") {
    return trajfmt::dlp4;
  }
  if(fmt=="xtc") {
    return trajfmt::xdr_xtc;
  }
  if(fmt=="trr") {
    return trajfmt::xdr_trr;
  }
  return trajfmt::error;
}
std::string TrajectoryParser::toString(trajfmt fmt) {
  switch (fmt) {
  case trajfmt::molfile: {
    return "molfile";
  }
  case trajfmt::xdr_xtc: {
    return "xtc";
  }
  case trajfmt::xdr_trr: {
    return "trr";
  }
  case trajfmt::xyz: {
    return "xyz";
  }
  case trajfmt::gro: {
    return "gro";
  }
  case trajfmt::dlp4: {
    return "dlp4";
  }
  case trajfmt::error: {
    return "error";
  }
  }

  return "unknown";
}

std::optional<std::string> TrajectoryParser::init(std::string_view fmt, FILE* in) {
  auto trajectory_fmt =FMTfromString(fmt);
  if (trajectory_fmt==trajfmt::error) {
    return "Error: wrong input format specified";
  }
  parser = parserFromFMT(trajectory_fmt);
  //not owning, will not be deleted on destruction
  return parser->init(in);
}
std::optional<std::string> TrajectoryParser::init(std::string_view fmt,
    std::string_view fname,
    bool useMolfile,
    int command_line_natoms) {
  auto trajectory_fmt =FMTfromString(fmt);
  if (useMolfile) {
    trajectory_fmt=trajfmt::molfile;
  }
  if (trajectory_fmt==trajfmt::error) {
    return "Error: wrong input format specified";
  }

  parser = parserFromFMT(trajectory_fmt);
  return parser->init(fmt,fname,command_line_natoms);
}

template <typename real>
std::optional<std::string> TrajectoryParser::readHeader(TrajectoryFrame<real>& frame) {
  raw->step.reset();
  raw->timestep.reset();
  auto error=parser->readHeader(*raw);
  if(error) {
    return error;
  }
  if(raw->step) {
    frame.step=*raw->step;
  }
  frame.timestep=raw->timestep;
  return std::nullopt;
}

template <typename real>
std::optional<std::string> TrajectoryParser::readAtoms(TrajectoryFrame<real>& frame) {
  const int natoms=parser->nOfAtoms();
  raw->readBox=frame.readBox;
  raw->step.reset();
  raw->time.reset();
  raw->positions.resize(3*natoms);
  raw->box.fill(0.0);
  raw->masses.clear();
  raw->charges.clear();
  auto error=parser->readAtoms(*raw);
  if(error) {
    return error;
  }
  if(raw->step && frame.stride==0) {
    frame.step=*raw->step;
  }
  frame.time=raw->time;
  const int first=(frame.count<0) ? 0 : std::max(frame.first,0);
  const int last=(frame.count<0) ? natoms : std::min(natoms,frame.first+frame.count);
  for(int i=3*first; i<3*last; i++) {
    frame.positions[i]=real(raw->positions[i]);
  }
  if(frame.box) {
    for(unsigned i=0; i<9; i++) {
      frame.box[i]=real(raw->box[i]);
    }
  }
  frame.hasMassesAndCharges=!raw->masses.empty();
  if(frame.hasMassesAndCharges) {
    for(int i=0; i<natoms; i++) {
      frame.masses[i]=real(raw->masses[i]);
      frame.charges[i]=real(raw->charges[i]);
    }
  }
  return std::nullopt;
}

template <typename real>
std::optional<std::string> TrajectoryParser::readFrame(TrajectoryFrame<real>& frame) {
  auto error = readHeader(frame);
  if (error) {
    return error;
  }
  return readAtoms(frame);
}

template std::optional<std::string> TrajectoryParser::readHeader<double>(TrajectoryFrame<double>&);
template std::optional<std::string> TrajectoryParser::readHeader<float>(TrajectoryFrame<float>&);
template std::optional<std::string> TrajectoryParser::readAtoms<double>(TrajectoryFrame<double>&);
template std::optional<std::string> TrajectoryParser::readAtoms<float>(TrajectoryFrame<float>&);
template std::optional<std::string> TrajectoryParser::readFrame<double>(TrajectoryFrame<double>&);
template std::optional<std::string> TrajectoryParser::readFrame<float>(TrajectoryFrame<float>&);

//some getters
int TrajectoryParser::nOfAtoms() const {
  if (parser) {
    return parser->nOfAtoms();
  }
  return -1;
}

std::optional<std::string> TrajectoryParser::rewind() {
  return parser->rewind();
}

} // namespace PLMD
