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
  virtual std::optional<std::string> init(FILE* fileHandle) {
    return std::nullopt;
  };

  virtual std::optional<std::string> readHeader(
    long long int &step,
    double &timeStep
  )=0;

  virtual std::optional<std::string> readAtoms(    int stride,
      bool dont_read_pbc,
      bool debug_pd,
      int pd_start,
      int pd_nlocal,
      long long int &step,
      double* masses,
      double* charges,
      double* coordinates,
      double* cell )=0;

  virtual std::optional<std::string> readHeader(
    long long int &step,
    float &timeStep
  ) =0;
  virtual std::optional<std::string> readAtoms(int stride,
      bool dont_read_pbc,
      bool debug_pd,
      int pd_start,
      int pd_nlocal,
      long long int &step,
      float* masses,
      float* charges,
      float* coordinates,
      float* cell ) =0;

  virtual std::optional<std::string> rewind() {
    return "rewind not (yet) implemented for this kind of trajectory file";
  }

};
namespace {

#define READATOMS  \
  std::optional<std::string> readAtoms(int stride,bool dont_read_pbc, \
    bool debug_pd,int pd_start,int pd_nlocal,long long int &step, \
    double* masses,double* charges,double* coordinates,double* cell ) override { \
    return readAtoms_t<double>(stride,dont_read_pbc,debug_pd,pd_start,pd_nlocal, \
        step,masses,charges,coordinates,cell ); \
  } \
 \
  std::optional<std::string> readAtoms(int stride,bool dont_read_pbc, \
    bool debug_pd,int pd_start,int pd_nlocal,long long int &step, \
    float* masses,float* charges,float* coordinates,float* cell ) override { \
    return readAtoms_t<float>(stride,dont_read_pbc,debug_pd,pd_start,pd_nlocal, \
        step,masses,charges,coordinates,cell ); \
  }
struct file_deleter {
  void operator()(FILE* fp) {
    std::fclose(fp);
  }
};

using safeFile=std::unique_ptr<std::FILE, file_deleter>;
class xyzParser final :public fileParser  {
  // not owning pointer for "-"mode
  std::FILE* fp{nullptr};
  // owning pointer otherwise
  safeFile fp_owned{nullptr,{}};

  template <typename real>
  std::optional<std::string> readAtoms_t(int stride,
                                         bool dont_read_pbc,
                                         bool debug_pd,
                                         int pd_start,
                                         int pd_nlocal,
                                         long long int &step,
                                         real* masses,
                                         real* charges,
                                         real* coordinates,
                                         real* cell )  {
    std::string line;
    //cell
    if(!Tools::getline(fp,line)) {
      return "premature end of trajectory file";
    }
    std::vector<double> celld(9,0.0);
    if(!dont_read_pbc) {
      std::vector<std::string> words=Tools::getWords(line);
      if(words.size()==3) {
        Tools::convert(words[0],celld[0]);
        Tools::convert(words[1],celld[4]);
        Tools::convert(words[2],celld[8]);
      } else if(words.size()==9) {
        Tools::convert(words[0],celld[0]);
        Tools::convert(words[1],celld[1]);
        Tools::convert(words[2],celld[2]);
        Tools::convert(words[3],celld[3]);
        Tools::convert(words[4],celld[4]);
        Tools::convert(words[5],celld[5]);
        Tools::convert(words[6],celld[6]);
        Tools::convert(words[7],celld[7]);
        Tools::convert(words[8],celld[8]);
      } else {
        return "needed box in second line of xyz file";
      }

    }
    std::copy(celld.begin(),celld.end(),cell);
    //atoms, finally
    for(int i=0; i<natoms; i++) {
      double cc[3];
      if(!Tools::getline(fp,line)) {
        return "premature end of trajectory file";
      }
      char dummy[1000];
      int ret=std::sscanf(line.c_str(),"%999s %100lf %100lf %100lf",
                          dummy,&cc[0],&cc[1],&cc[2]);
      if(ret!=4) {
        return "cannot read line"+line;
      }
      if(!debug_pd || ( i>=pd_start && i<pd_start+pd_nlocal) ) {
        coordinates[3*i]=real(cc[0]);
        coordinates[3*i+1]=real(cc[1]);
        coordinates[3*i+2]=real(cc[2]);
      }
    }
    return std::nullopt;
  }
public:
  std::optional<std::string> init(std::string_view fmt,
                                  std::string_view fname,
                                  int command_line_natoms=-1) override {
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
  std::optional<std::string> readHeader(
    long long int &/*step*/,
    double &/*timeStep*/
  ) override {
    std::string line;
    if(!Tools::getline(fp,line)) {
      return "EOF";
    }
    std::sscanf(line.c_str(),"%100d",&natoms);
    return std::nullopt;
  }

  std::optional<std::string> readHeader(
    long long int &step,
    float &/*timeStep*/
  ) override {
    double ts;
    //ignore this -Wmaybe-uninitialized
    auto msg=readHeader(step,ts);
    return msg;
  }
  //see the template readAtoms_t for the implementation
  READATOMS;
  std::optional<std::string> rewind() override {
    int error = std::fseek(fp,0,SEEK_SET);
    if (error) {
      return "ERROR: Error rewinding trajectory file";
    }
    return std::nullopt;
  }
};

class dlp4Parser final:public fileParser {
  // not owning pointer for "-"mode
  std::FILE* fp{nullptr};
  // owning pointer otherwise
  safeFile fp_owned{nullptr,{}};
  int lvl{0};
  template <typename real>
  std::optional<std::string> readAtoms_t(int stride,
                                         bool dont_read_pbc,
                                         bool debug_pd,
                                         int pd_start,
                                         int pd_nlocal,
                                         long long int &step,
                                         real* masses,
                                         real* charges,
                                         real* coordinates,
                                         real* cell ) {
    std::string line;

    //cell
    std::vector<double> celld(9,0.0);
    if(!dont_read_pbc) {
      if(!Tools::getline(fp,line)) {
        return "error reading vector a of cell";
      }
      std::sscanf(line.c_str(),"%lf %lf %lf",&celld[0],&celld[1],&celld[2]);
      if(!Tools::getline(fp,line)) {
        return "error reading vector b of cell";
      }
      std::sscanf(line.c_str(),"%lf %lf %lf",&celld[3],&celld[4],&celld[5]);
      if(!Tools::getline(fp,line)) {
        return "error reading vector c of cell";
      }
      std::sscanf(line.c_str(),"%lf %lf %lf",&celld[6],&celld[7],&celld[8]);
      std::transform(celld.begin(),celld.end(),celld.begin(),
      [](double x) {
        //converting to nm
        return x/10.0;
      });
      std::copy(celld.begin(),celld.end(),cell);
    }
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
      masses[i]=real(m);
      charges[i]=real(c);
      if(!Tools::getline(fp,line)) {
        return "error reading coordinates";
      }
      std::sscanf(line.c_str(),"%lf %lf %lf",&cc[0],&cc[1],&cc[2]);
      cc[0]*=0.1;
      cc[1]*=0.1;
      cc[2]*=0.1;
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
      if(!debug_pd || ( i>=pd_start && i<pd_start+pd_nlocal) ) {
        coordinates[3*i]=real(cc[0]);
        coordinates[3*i+1]=real(cc[1]);
        coordinates[3*i+2]=real(cc[2]);
      }
    }

    return std::nullopt;
  }

public:
  std::optional<std::string> init(std::string_view fmt,
                                  std::string_view fname,
                                  int command_line_natoms=-1) override {
    std::string trajectoryFile=std::string(fname);
    fp_owned.reset(std::fopen(trajectoryFile.c_str(),"r"));
    fp = fp_owned.get();
    if(!fp) {
      return "ERROR: Error opening trajectory file "+trajectoryFile;
    }
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
  std::optional<std::string> init(FILE* fileHandle) override {
    fp = fileHandle;
    return std::nullopt;
  }

  std::optional<std::string> readHeader (
    long long int &step,
    double &timeStep
  ) override {
    std::string line;
    if(!Tools::getline(fp,line)) {
      return "EOF";
    }
    char xa[9];
    int xb,xc,xd;
    double t;
    auto res = std::sscanf(line.c_str(),"%8s %lld %d %d %d %lf",xa,&step,&xb,&xc,&xd,&t);
    if (res!=6) {
      return "error reading header";
    }
    timeStep=t;
    return std::nullopt;
  }

  std::optional<std::string> readHeader(
    long long int &step,
    float &timeStep
  ) override {
    double ts=timeStep;
    auto msg=readHeader(step,ts);
    timeStep = float(ts);
    return msg;
  }
  //see the template readAtoms_t for the implementation
  READATOMS;
  std::optional<std::string> rewind() override {
    int error = std::fseek(fp,0,SEEK_SET);
    if (error) {
      return "ERROR: Error rewinding trajectory file";
    }
    return std::nullopt;
  }
};

class groParser final:public fileParser {
  // not owning pointer for "-"mode
  std::FILE* fp{nullptr};
  // owning pointer otherwise
  safeFile fp_owned{nullptr,{}};
  template <typename real>
  std::optional<std::string> readAtoms_t(int stride,
                                         bool dont_read_pbc,
                                         bool debug_pd,
                                         int pd_start,
                                         int pd_nlocal,
                                         long long int &step,
                                         real* masses,
                                         real* charges,
                                         real* coordinates,
                                         real* cell ) {
    std::string line;
    int ddist=0;
    for(int i=0; i<natoms; i++) {
      if(!Tools::getline(fp,line)) {
        return "premature end of trajectory file";
      }
      double cc[3];
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
      Tools::convert(line.substr(20,ddist),cc[0]);
      Tools::convert(line.substr(20+ddist,ddist),cc[1]);
      Tools::convert(line.substr(20+ddist+ddist,ddist),cc[2]);
      if(!debug_pd || ( i>=pd_start && i<pd_start+pd_nlocal) ) {
        coordinates[3*i]=real(cc[0]);
        coordinates[3*i+1]=real(cc[1]);
        coordinates[3*i+2]=real(cc[2]);
      }

    }
    if(!Tools::getline(fp,line)) {
      return "premature end of trajectory file";
    }
    std::vector<std::string> words=Tools::getWords(line);
    if(words.size()<3) {
      return "cannot understand box format";
    }
    Tools::convert(words[0],cell[0]);
    Tools::convert(words[1],cell[4]);
    Tools::convert(words[2],cell[8]);
    if(words.size()>3) {
      Tools::convert(words[3],cell[1]);
    }
    if(words.size()>4) {
      Tools::convert(words[4],cell[2]);
    }
    if(words.size()>5) {
      Tools::convert(words[5],cell[3]);
    }
    if(words.size()>6) {
      Tools::convert(words[6],cell[5]);
    }
    if(words.size()>7) {
      Tools::convert(words[7],cell[6]);
    }
    if(words.size()>8) {
      Tools::convert(words[8],cell[7]);
    }
    return std::nullopt;
  }

public:
  std::optional<std::string> init(std::string_view fmt,
                                  std::string_view fname,
                                  int command_line_natoms=-1) override {
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

  std::optional<std::string> readHeader (
    long long int &/**/,
    double &/**/
  ) override {
    std::string line;
    if(!Tools::getline(fp,line)) {
      return "EOF";
    }
    if(!Tools::getline(fp,line)) {
      return "premature end of trajectory file";
    }
    std::sscanf(line.c_str(),"%100d",&natoms);

    return std::nullopt;
  }

  std::optional<std::string> readHeader(
    long long int &step,
    float &/**/
  ) override {
    double ts;
    auto msg=readHeader(step,ts);
    return msg;
  }
  //see the template readAtoms_t for the implementation
  READATOMS;
  std::optional<std::string> rewind() override {
    int error = std::fseek(fp,0,SEEK_SET);
    if (error) {
      return "ERROR: Error rewinding trajectory file";
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

  template <typename real>
  std::optional<std::string> readAtoms_t(int stride,
                                         bool dont_read_pbc,
                                         bool debug_pd,
                                         int pd_start,
                                         int pd_nlocal,
                                         long long int &step,
                                         real* masses,
                                         real* charges,
                                         real* coordinates,
                                         real* cell ) {
    int localstep;
    float time;
    xdrfile::matrix box;
    // here we cannot use a std::vector<rvec> since it does not compile.
    // we thus use a std::unique_ptr<rvec[]>
    auto pos=Tools::make_unique<xdrfile::rvec[]>(natoms);
    int ret=xdrfile::exdrOK;
    if constexpr(is == xdType::xtc) {
      float prec;
      ret=xdrfile::read_xtc(xd.get(),natoms,&localstep,&time,box,pos.get(),&prec);
    }
    if constexpr(is == xdType::trr) {
      float lambda;
      ret=xdrfile::read_trr(xd.get(),natoms,&localstep,&time,&lambda,box,pos.get(),NULL,NULL);
    }
    if(stride==0) {
      step=localstep;
    }
    if(ret==xdrfile::exdrENDOFFILE) {
      return "EOF";
    }
    if(ret!=xdrfile::exdrOK) {
      return "EOF";
    }
    for(unsigned i=0; i<3; i++) {
      for(unsigned j=0; j<3; j++) {
        cell[3*i+j]=box[i][j];
      }
    }
    //xdrfiles are not used with dd ?
    for(int i=0; i<natoms; i++) {
      for(unsigned j=0; j<3; j++) {
        coordinates[3*i+j]=real(pos[i][j]);
      }
    }
    return std::nullopt;
  }

public:
  std::optional<std::string> init(std::string_view fmt,
                                  std::string_view fname,
                                  int command_line_natoms=-1) override {
    using trajfmt=TrajectoryParser::trajfmt;
    std::string trajectoryFile=std::string(fname);
    auto trajectory_fmt =TrajectoryParser::FMTfromString(fmt);
    xd.reset(xdrfile::xdrfile_open(trajectoryFile.c_str(),"r"));
    if(!xd) {
      return "ERROR: Error opening trajectory file "+trajectoryFile;
    }
    if constexpr(is == xdType::xtc) {
      plumed_assert(trajectory_fmt==trajfmt::xdr_xtc) << "trajectory type should be xtc";
      xdrfile::read_xtc_natoms(&trajectoryFile[0],&natoms);
    }
    if constexpr(is == xdType::trr) {
      plumed_assert(trajectory_fmt==trajfmt::xdr_trr) << "trajectory type should be trr";
      xdrfile::read_trr_natoms(&trajectoryFile[0],&natoms);
    }
    return std::nullopt;
  }

  std::optional<std::string> init(FILE* /*fileHandle*/) override {
    plumed_assert(false) << "xd readers are not compatible with reading from stdin";
    return std::nullopt;
  }

  std::optional<std::string> readHeader (long long int &/*step*/,double &/*timeStep*/) override {
    return std::nullopt;
  }

  std::optional<std::string> readHeader(long long int &/*step*/,float &/*timeStep*/) override {
    return std::nullopt;
  }
  //see the template readAtoms_t for the implementation
  READATOMS;
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


  template <typename real>
  std::optional<std::string> readAtoms_t(int stride,
                                         bool dont_read_pbc,
                                         bool debug_pd,
                                         int pd_start,
                                         int pd_nlocal,
                                         long long int &step,
                                         real* masses,
                                         real* charges,
                                         real* coordinates,
                                         real* cell ) {
    if(!dont_read_pbc) {
      //ts_in has been read in the previous step
      if(ts_in.A>0.0) { // this is negative if molfile does not provide box
        // info on the cell: convert using pbcset.tcl from pbctools in vmd distribution
        constexpr real r180=real(180.);
        real cosBC=cos(real(ts_in.alpha)*pi/r180);
        real cosAC=std::cos(real(ts_in.beta)*pi/r180);
        real cosAB=std::cos(real(ts_in.gamma)*pi/r180);
        real sinAB=std::sin(real(ts_in.gamma)*pi/r180);
        real Ax=real(ts_in.A);
        real Bx=real(ts_in.B)*cosAB;
        real By=real(ts_in.B)*sinAB;
        real Cx=real(ts_in.C)*cosAC;
        real Cy=(real(ts_in.C)*real(ts_in.B)*cosBC-Cx*Bx)/By;
        real Cz=std::sqrt(real(ts_in.C)*real(ts_in.C)-Cx*Cx-Cy*Cy);
        cell[0]=Ax/10.;
        cell[1]=0.;
        cell[2]=0.;
        cell[3]=Bx/10.;
        cell[4]=By/10.;
        cell[5]=0.;
        cell[6]=Cx/10.;
        cell[7]=Cy/10.;
        cell[8]=Cz/10.;
      } else {
        cell[0]=0.0;
        cell[1]=0.0;
        cell[2]=0.0;
        cell[3]=0.0;
        cell[4]=0.0;
        cell[5]=0.0;
        cell[6]=0.0;
        cell[7]=0.0;
        cell[8]=0.0;
      }
    }
    // info on coords
    // the order is xyzxyz...
    for(int i=0; i<3*natoms; i++) {
      coordinates[i]=real(ts_in.coords[i])/real(10.); //convert to nm
      //cerr<<"COOR "<<coordinates[i]<<endl;
    }
    return std::nullopt;
  }

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

  std::optional<std::string> init(FILE* /*fileHandle*/) override {
    plumed_assert(false) << "molfile readers are not compatible with reading from stdin";
    return std::nullopt;
  }

  std::optional<std::string> readHeader (long long int &/*step*/,double &/*timeStep*/) override {
    std::unique_ptr<std::lock_guard<std::mutex>> lck;
    if(api->is_reentrant==VMDPLUGIN_THREADUNSAFE) {
      lck=Tools::molfile_lock();
    }
    int rc;
    rc = api->read_next_timestep(h_in.get(), natoms, &ts_in);
    if(rc==MOLFILE_EOF) {
      return "EOF";
    }
    return std::nullopt;
  }

  std::optional<std::string> readHeader(long long int &/*step*/,float &/*timeStep*/) override {
    return std::nullopt;
  }
//see the template readAtoms_t for the implementation
  READATOMS;
};
#endif

#undef READATOMS

std::unique_ptr<fileParser> parserFromFMT(TrajectoryParser::trajfmt fmt) {
  using trajfmt = TrajectoryParser::trajfmt;
  switch (fmt) {
  case trajfmt::molfile: {
#ifdef __PLUMED_HAS_MOLFILE_PLUGINS
    return std::make_unique<molfileParser>(molfileParser {});
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

TrajectoryParser::TrajectoryParser()=default;
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

std::optional<std::string> TrajectoryParser::readHeader(
  long long int &step,
  double &timeStep
) {
  return  parser->readHeader(step,timeStep);
}

std::optional<std::string> TrajectoryParser::readAtoms(
  int stride,
  bool dont_read_pbc,
  bool debug_pd,
  int pd_start,
  int pd_nlocal,
  long long int &step,
  double* masses,
  double* charges,
  double* coordinates,
  double* cell
) {
  return parser->readAtoms(
           stride,
           dont_read_pbc,
           debug_pd,
           pd_start,
           pd_nlocal,
           step,
           masses,
           charges,
           coordinates,
           cell
         );
}

std::optional<std::string> TrajectoryParser::readFrame(
  int stride,
  bool dont_read_pbc,
  bool debug_pd,
  int pd_start,
  int pd_nlocal,
  long long int &step,
  double &timeStep,
  double * masses,
  double * charges,
  double * coordinates,
  double * cell
) {
  auto error = readHeader(step,timeStep);
  if (error) {
    return error;
  }
  ////////////////////////////////////////////////////////////////////////////
  // here the first step settings
  ////////////////////////////////////////////////////////////////////////////
  //now the actual reading of the trajectory
  return readAtoms(
           stride,
           dont_read_pbc,
           debug_pd,
           pd_start,
           pd_nlocal,
           step,
           masses,
           charges,
           coordinates,
           cell
         );

}


std::optional<std::string> TrajectoryParser::readFrame(
  int stride,
  bool dont_read_pbc,
  bool debug_pd,
  int pd_start,
  int pd_nlocal,
  long long int &step,
  float &timeStep,
  float * masses,
  float * charges,
  float * coordinates,
  float * cell
) {
  auto error = readHeader(step,timeStep);
  if (error) {
    return error;
  }
  ////////////////////////////////////////////////////////////////////////////
  // here the first step settings
  ////////////////////////////////////////////////////////////////////////////
  //now the actual reading of the trajectory
  return readAtoms(
           stride,
           dont_read_pbc,
           debug_pd,
           pd_start,
           pd_nlocal,
           step,
           masses,
           charges,
           coordinates,
           cell
         );
}

std::optional<std::string> TrajectoryParser::readHeader(
  long long int &step,
  float &timeStep
) {
  return parser->readHeader(step,timeStep);
}

std::optional<std::string> TrajectoryParser::readAtoms(
  int stride,
  bool dont_read_pbc,
  bool debug_pd,
  int pd_start,
  int pd_nlocal,
  long long int &step,
  float* masses,
  float* charges,
  float* coordinates,
  float* cell
) {
  return parser->readAtoms(stride,
                           dont_read_pbc,
                           debug_pd,
                           pd_start,
                           pd_nlocal,
                           step,
                           masses,
                           charges,
                           coordinates,
                           cell);
}


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
