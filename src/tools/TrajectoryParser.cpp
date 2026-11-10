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
#include "FileTools.h"
#include "IFile.h"
#include "Tools.h"
#include "xdrfile/xdrfile.h"
#include "xdrfile/xdrfile_trr.h"
#include "xdrfile/xdrfile_xtc.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>

#ifdef __PLUMED_HAS_LIBCHEMFILES
#include <chemfiles.hpp>
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
  /// output of readAtoms(): set with "EOF" when the last frame is incomplete
  std::optional<std::string> warning;
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
                                          std::string_view fname)=0;
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
                                  std::string_view fname) override {
    std::string trajectoryFile=std::string(fname);
    fp_owned.reset(fopenUTF8(trajectoryFile,"r"));
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
                                  std::string_view fname) override {
    auto error=textParser::init(fmt,fname);
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
                                  std::string_view fname) override {
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
      // a frame that ends with the file is incomplete (the simulation stopped while writing it):
      // as gmx, it is the end of the trajectory, with a warning. Other errors are corrupted data
      if(xdrfile::xdrfile_eof(xd.get())) {
        raw.warning="the last frame of the "+std::string(name)+" file "+trajectoryFile+" is incomplete, it is ignored";
        return "EOF";
      }
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

/// The box of the cell with sides a, b and c and angles alpha (bc), beta (ac) and gamma (ab) (in degrees),
/// with a along x and b in the xy plane (as pbcset.tcl of the pbctools of VMD), in the units of a, b and c.
/// The box is zero (no box) if a is not positive
void boxFromLengthsAndAngles(double a, double b, double c,
                             double alpha, double beta, double gamma,
                             std::array<double,9>& box) {
  if(!(a>0.0)) {
    return;
  }
  constexpr double r180=180.;
  const double cosBC=std::cos(alpha*pi/r180);
  const double cosAC=std::cos(beta*pi/r180);
  const double cosAB=std::cos(gamma*pi/r180);
  const double sinAB=std::sin(gamma*pi/r180);
  const double Bx=b*cosAB;
  const double By=b*sinAB;
  const double Cx=c*cosAC;
  const double Cy=(c*b*cosBC-Cx*Bx)/By;
  const double Cz=std::sqrt(c*c-Cx*Cx-Cy*Cy);
  box= {a,0.,0.,
        Bx,By,0.,
        Cx,Cy,Cz
       };
}

/// pdb trajectories: frames of ATOM and HETATM records that end with an END or ENDMDL record.
/// As the molfile plugin did, only the coordinates (columns 31-54) and the CRYST1 record are read:
/// the other columns and records (serial numbers, occupancy, beta, REMARK ...) are not used, since the
/// structure is read elsewhere (MOLINFO, --pdb) with the PDB class, which checks them.
/// A frame without CRYST1 keeps the box of the previous frame. All the frames must have the same number of atoms
class pdbParser final : public textParser {
  std::string trajectoryFile{"standard input"};
  /// the first frame, read by init() for the number of atoms
  std::vector<double> pending;
  /// the positions of the frame whose header was read last, in nm
  std::vector<double> current;
  /// the box of the last CRYST1 record, in nm
  std::array<double,9> box{};
  long long frameNumber{0};
  long long lineNumber{0};
  std::string where() const {
    return "line "+std::to_string(lineNumber)+" of the pdb file "+trajectoryFile;
  }
  /// read the next frame (positions in nm): "EOF" if it has no atoms
  std::optional<std::string> readPDBFrame(std::vector<double>& pos) {
    pos.clear();
    ++frameNumber;
    std::string line;
    while(Tools::getline(fp,line)) {
      ++lineNumber;
      std::string record=line.substr(0,6);
      Tools::trim(record);
      if(record=="END" || record=="ENDMDL") {
        break;
      }
      if(record=="CRYST1") {
        line.resize(std::max<std::size_t>(line.size(),54),' ');
        double cell[6];
        const std::size_t start[6]= {6,15,24,33,40,47};
        const std::size_t width[6]= {9,9,9,7,7,7};
        for(unsigned i=0; i<6; i++) {
          if(!Tools::convertNoexcept(line.substr(start[i],width[i]),cell[i])) {
            return "cannot read the CRYST1 record (columns 7-54) at "+where()+":\n"+line;
          }
        }
        if(cell[0]>0.0) {
          boxFromLengthsAndAngles(cell[0],cell[1],cell[2],cell[3],cell[4],cell[5],box);
          for(auto & x : box) {
            x/=10.;
          }
        }
      }
      if(record=="ATOM" || record=="HETATM") {
        line.resize(std::max<std::size_t>(line.size(),54),' ');
        for(unsigned j=0; j<3; j++) {
          double x;
          if(!Tools::convertNoexcept(line.substr(30+8*j,8),x)) {
            return "cannot read the coordinates (columns 31-54) at "+where()+":\n"+line;
          }
          pos.push_back(x/10.);
        }
      }
    }
    if(pos.empty()) {
      return "EOF";
    }
    return std::nullopt;
  }
  std::optional<std::string> readFirstFrame() {
    frameNumber=0;
    lineNumber=0;
    box.fill(0.0);
    auto error=readPDBFrame(pending);
    if(error) {
      if(*error=="EOF") {
        return "the pdb file "+trajectoryFile+" contains no atoms";
      }
      return error;
    }
    natoms=pending.size()/3;
    return std::nullopt;
  }
public:
  std::optional<std::string> init(std::string_view fmt,
                                  std::string_view fname) override {
    auto error=textParser::init(fmt,fname);
    if(error) {
      return error;
    }
    trajectoryFile=std::string(fname);
    return readFirstFrame();
  }
  std::optional<std::string> init(FILE* fileHandle) override {
    textParser::init(fileHandle);
    return readFirstFrame();
  }
  std::optional<std::string> rewind() override {
    auto error=textParser::rewind();
    if(error) {
      return error;
    }
    return readFirstFrame();
  }
  std::optional<std::string> readHeader(RawFrame& /*raw*/) override {
    if(!pending.empty()) {
      current.swap(pending);
      pending.clear();
      return std::nullopt;
    }
    auto error=readPDBFrame(current);
    if(error) {
      return error;
    }
    if(current.size()!=3*std::size_t(natoms)) {
      return "frame "+std::to_string(frameNumber)+" of the pdb file "+trajectoryFile+" has "
             +std::to_string(current.size()/3)+" atoms, the first frame has "+std::to_string(natoms);
    }
    return std::nullopt;
  }
  std::optional<std::string> readAtoms(RawFrame& raw) override {
    std::copy(current.begin(),current.end(),raw.positions.begin());
    if(raw.readBox) {
      raw.box=box;
    }
    return std::nullopt;
  }
};

/// dcd trajectories as written by CHARMM, NAMD, X-PLOR and other codes: either endianness,
/// 32 or 64 bit (CHARMM -i8) record markers, fixed atoms, unit cell, 4th dimension (skipped).
/// Written after the dcd plugin of VMD (molfile), with every size of the file checked.
/// The coordinates and the cell are single precision values in angstrom (they give the same numbers as the plugin did).
class dcdParser final : public fileParser {
  safeFile fp{nullptr,{}};
  std::string trajectoryFile;
  /// true if the file has the other endianness
  bool swap{false};
  /// size of the record markers of the Fortran records (4, or 8 for CHARMM -i8)
  unsigned markerSize{4};
  bool charmm{false};
  bool hasUnitCell{false};
  bool hasFourDims{false};
  /// the indexes (from 0) of the atoms that are not fixed, if there are fixed atoms
  std::vector<std::int32_t> freeIndexes;
  int nfixed{0};
  /// x, y and z of all the atoms in the first frame (fixed atoms only)
  std::vector<float> fixedCoords;
  std::vector<float> buffer;
  long headerEnd{0};
  bool firstFrame{true};
  /// bytes read in the current frame
  std::size_t frameBytes{0};

  static void swap4(void* p, std::size_t n) {
    auto c=static_cast<unsigned char*>(p);
    for(std::size_t i=0; i<n; i++,c+=4) {
      std::swap(c[0],c[3]);
      std::swap(c[1],c[2]);
    }
  }
  std::int32_t int32At(const unsigned char* p) const {
    std::int32_t v;
    std::memcpy(&v,p,4);
    if(swap) {
      swap4(&v,1);
    }
    return v;
  }
  std::string corrupted(const std::string& what) const {
    return "the dcd file "+trajectoryFile+" is corrupted or is not a dcd file ("+what+")";
  }
  /// read n bytes of the header
  bool readHeaderBytes(void* p, std::size_t n) {
    return std::fread(p,1,n,fp.get())==n;
  }
  /// a record marker of the header
  bool readHeaderMarker(long long& v) {
    unsigned char m[8];
    if(!readHeaderBytes(m,markerSize)) {
      return false;
    }
    // as the plugin: with 64 bit markers, the sum of the two 32 bit halves
    v=int32At(m);
    if(markerSize==8) {
      v+=int32At(m+4);
    }
    return true;
  }
  /// read n bytes of a frame: at the end of the file, "EOF" if the frame was not started,
  /// otherwise "EOF" and the warning of an incomplete last frame
  std::optional<std::string> readFrameBytes(void* p, std::size_t n, RawFrame& raw) {
    const auto r=std::fread(p,1,n,fp.get());
    frameBytes+=r;
    if(r==n) {
      return std::nullopt;
    }
    if(std::ferror(fp.get())) {
      return "error reading the dcd file "+trajectoryFile;
    }
    if(frameBytes>0) {
      raw.warning="the last frame of the dcd file "+trajectoryFile+" is incomplete, it is ignored";
    }
    return "EOF";
  }
  std::optional<std::string> readFrameMarker(long long& v, RawFrame& raw) {
    unsigned char m[8];
    if(auto error=readFrameBytes(m,markerSize,raw)) {
      return error;
    }
    v=int32At(m);
    if(markerSize==8) {
      v+=int32At(m+4);
    }
    return std::nullopt;
  }
  /// a record of n 4-byte values
  std::optional<std::string> readFrameRecord(void* p, std::size_t n, const char* what, RawFrame& raw) {
    long long m;
    if(auto error=readFrameMarker(m,raw)) {
      return error;
    }
    if(m!=4*static_cast<long long>(n)) {
      return corrupted(std::string("size of the ")+what+" record");
    }
    if(auto error=readFrameBytes(p,4*n,raw)) {
      return error;
    }
    if(auto error=readFrameMarker(m,raw)) {
      return error;
    }
    if(m!=4*static_cast<long long>(n)) {
      return corrupted(std::string("size of the ")+what+" record");
    }
    if(swap) {
      swap4(p,n);
    }
    return std::nullopt;
  }
  std::optional<std::string> readFileHeader() {
    unsigned char first[8];
    if(!readHeaderBytes(first,8)) {
      return corrupted("the file is shorter than the header");
    }
    std::int32_t m0,m1;
    std::memcpy(&m0,first,4);
    std::memcpy(&m1,first+4,4);
    const bool cord=std::memcmp(first+4,"CORD",4)==0;
    auto swapped=[](std::int32_t v) {
      swap4(&v,1);
      return v;
    };
    if(m0==84 && cord) {
      swap=false;
      markerSize=4;
    } else if(swapped(m0)==84 && cord) {
      swap=true;
      markerSize=4;
    } else if(static_cast<long long>(m0)+m1==84) {
      swap=false;
      markerSize=8;
    } else if(static_cast<long long>(swapped(m0))+swapped(m1)==84) {
      swap=true;
      markerSize=8;
    } else {
      return corrupted("no CORD header");
    }
    if(markerSize==8) {
      char magic[4];
      if(!readHeaderBytes(magic,4) || std::memcmp(magic,"CORD",4)!=0) {
        return corrupted("no CORD header");
      }
    }
    unsigned char hdr[80];
    if(!readHeaderBytes(hdr,80)) {
      return corrupted("the file is shorter than the header");
    }
    // CHARMM writes its version in the last integer, X-PLOR writes zero
    charmm=int32At(hdr+76)!=0;
    hasUnitCell=charmm && int32At(hdr+40)!=0;
    hasFourDims=charmm && int32At(hdr+44)==1;
    nfixed=int32At(hdr+32);
    long long m;
    if(!readHeaderMarker(m) || m!=84) {
      return corrupted("end of the first record");
    }
    // title: the number of lines and 80 characters per line.
    // As the plugin, the lines are skipped with the number of lines, and the end of the record is not checked:
    // some writers give the record a size that does not match the number of lines
    long long titleSize;
    if(!readHeaderMarker(titleSize) || titleSize<4 || (titleSize-4)%80!=0) {
      return corrupted("size of the title record");
    }
    std::int32_t ntitle;
    if(!readHeaderBytes(&ntitle,4)) {
      return corrupted("the file is shorter than the header");
    }
    if(swap) {
      swap4(&ntitle,1);
    }
    if(ntitle<0) {
      return corrupted("number of title lines");
    }
    if(ntitle>1000) {
      // as the plugin: a broken file written by Vega ZZ 2.4.0 has 2 lines, other files none
      ntitle=(ntitle==1095062083) ? 2 : 0;
    }
    if(std::fseek(fp.get(),80L*ntitle,SEEK_CUR)!=0 || !readHeaderMarker(m)) {
      return corrupted("the file is shorter than the header");
    }
    std::int32_t n;
    if(!readHeaderMarker(m) || m!=4 || !readHeaderBytes(&n,4) || !readHeaderMarker(m) || m!=4) {
      return corrupted("number of atoms record");
    }
    if(swap) {
      swap4(&n,1);
    }
    // larger numbers of atoms are taken as a corrupted file (as the plugin did), the sizes of the records then overflow
    if(n<=0 || n>(1<<29)) {
      return corrupted("invalid number of atoms "+std::to_string(n));
    }
    natoms=n;
    if(nfixed<0 || nfixed>natoms) {
      return corrupted("invalid number of fixed atoms "+std::to_string(nfixed));
    }
    freeIndexes.clear();
    if(nfixed>0) {
      const long long nfree=natoms-nfixed;
      if(!readHeaderMarker(m) || m!=4*nfree) {
        return corrupted("size of the record of the free atoms");
      }
      freeIndexes.resize(nfree);
      if(!readHeaderBytes(freeIndexes.data(),4*nfree)) {
        return corrupted("the file is shorter than the header");
      }
      if(swap) {
        swap4(freeIndexes.data(),nfree);
      }
      for(auto & i : freeIndexes) {
        // 1-based in the file
        if(i<1 || i>natoms) {
          return corrupted("invalid index of a free atom "+std::to_string(i));
        }
        --i;
      }
      if(!readHeaderMarker(m) || m!=4*nfree) {
        return corrupted("end of the record of the free atoms");
      }
    }
    headerEnd=std::ftell(fp.get());
    // a file with a huge (corrupted) number of atoms is rejected before the buffers are allocated
    const long long rec=2LL*markerSize;
    const long long firstFrameSize=(hasUnitCell ? rec+48 : 0)+(hasFourDims ? 4 : 3)*(rec+4LL*natoms);
    if(headerEnd>=0 && std::fseek(fp.get(),0,SEEK_END)==0) {
      const long fileEnd=std::ftell(fp.get());
      if(std::fseek(fp.get(),headerEnd,SEEK_SET)!=0) {
        return "cannot read the dcd file "+trajectoryFile;
      }
      if(fileEnd>=0 && fileEnd-headerEnd<firstFrameSize) {
        return "the dcd file "+trajectoryFile+" does not contain a complete frame";
      }
    }
    buffer.assign(natoms,0.0f);
    fixedCoords.assign(nfixed>0 ? 3*natoms : 0,0.0f);
    firstFrame=true;
    return std::nullopt;
  }
public:
  std::optional<std::string> init(std::string_view /*fmt*/,
                                  std::string_view fname) override {
    trajectoryFile=std::string(fname);
    fp.reset(fopenUTF8(trajectoryFile,"rb"));
    if(!fp) {
      return "ERROR: Error opening trajectory file "+trajectoryFile;
    }
    return readFileHeader();
  }
  std::optional<std::string> rewind() override {
    if(std::fseek(fp.get(),headerEnd,SEEK_SET)!=0) {
      return "ERROR: Error rewinding trajectory file";
    }
    firstFrame=true;
    return std::nullopt;
  }
  std::optional<std::string> readHeader(RawFrame& /*raw*/) override {
    // the header is in the file header
    return std::nullopt;
  }
  std::optional<std::string> readAtoms(RawFrame& raw) override {
    frameBytes=0;
    // the unit cell, as the plugin read it: a, cos(gamma), b, cos(beta), cos(alpha), c
    // (NAMD 2.5 wrote the angles in degrees)
    float unitcell[6]= {0.0f,90.0f,0.0f,90.0f,90.0f,0.0f};
    if(hasUnitCell) {
      // six doubles, read as twelve 4-byte values
      unsigned char bytes[48];
      if(auto error=readFrameRecord(bytes,12,"unit cell",raw)) {
        return error;
      }
      for(unsigned i=0; i<6; i++) {
        unsigned char* b=bytes+8*i;
        if(swap) {
          // each 4-byte half was swapped by readFrameRecord: swap the two halves too
          std::swap_ranges(b,b+4,b+4);
        }
        double cell;
        std::memcpy(&cell,b,8);
        unitcell[i]=cell;
      }
    }
    const bool readAll=(nfixed==0 || firstFrame);
    const std::size_t nread=readAll ? natoms : freeIndexes.size();
    for(unsigned d=0; d<3; d++) {
      const char* what[3]= {"x","y","z"};
      if(auto error=readFrameRecord(buffer.data(),nread,what[d],raw)) {
        return error;
      }
      if(readAll) {
        for(int i=0; i<natoms; i++) {
          raw.positions[3*i+d]=double(buffer[i])/double(10.);
        }
        if(nfixed>0) {
          std::copy(buffer.begin(),buffer.end(),fixedCoords.begin()+std::size_t(d)*natoms);
        }
      } else {
        for(int i=0; i<natoms; i++) {
          raw.positions[3*i+d]=double(fixedCoords[std::size_t(d)*natoms+i])/double(10.);
        }
        for(std::size_t k=0; k<nread; k++) {
          raw.positions[3*freeIndexes[k]+d]=double(buffer[k])/double(10.);
        }
      }
    }
    if(hasFourDims) {
      // skipped: a record with one value per atom
      long long m;
      if(auto error=readFrameMarker(m,raw)) {
        return error;
      }
      if(m<0 || m>4LL*natoms) {
        return corrupted("size of the record of the 4th dimension");
      }
      if(auto error=readFrameBytes(buffer.data(),m,raw)) {
        return error;
      }
      long long mend;
      if(auto error=readFrameMarker(mend,raw)) {
        return error;
      }
      if(mend!=m) {
        return corrupted("size of the record of the 4th dimension");
      }
    }
    firstFrame=false;
    if(raw.readBox) {
      const float a=unitcell[0];
      const float b=unitcell[2];
      const float c=unitcell[5];
      float alpha,beta,gamma;
      auto isCosine=[](float x) {
        return x>=-1.0f && x<=1.0f;
      };
      if(isCosine(unitcell[1]) && isCosine(unitcell[3]) && isCosine(unitcell[4])) {
        // cosines of the angles (CHARMM, NAMD after 2.5)
        alpha=90.0-std::asin(unitcell[4])*90.0/(pi/2.0);
        beta=90.0-std::asin(unitcell[3])*90.0/(pi/2.0);
        gamma=90.0-std::asin(unitcell[1])*90.0/(pi/2.0);
      } else if(std::any_of(unitcell,unitcell+6,[](float x) {
      return x<0.0f;
    }) ||
    unitcell[1]>180.0f || unitcell[3]>180.0f || unitcell[4]>180.0f) {
        // not lengths and angles: the symmetric shape matrix h11, h21, h22, h31, h32, h33,
        // as written by CHARMM (at least since c36b2) for triclinic cells (the same guess as MDAnalysis).
        // The rows are the cell vectors, in the frame of the coordinates
        // (the plugin took these numbers for angles, which gave a box that is not valid)
        const unsigned order[9]= {0,1,3, 1,2,4, 3,4,5};
        for(unsigned i=0; i<9; i++) {
          raw.box[i]=double(unitcell[order[i]])/10.;
        }
        return std::nullopt;
      } else {
        // angles in degrees (NAMD 2.5)
        alpha=unitcell[4];
        beta=unitcell[3];
        gamma=unitcell[1];
      }
      boxFromLengthsAndAngles(a,b,c,alpha,beta,gamma,raw.box);
      for(auto & x : raw.box) {
        x/=10.;
      }
    }
    return std::nullopt;
  }
};

#ifdef __PLUMED_HAS_LIBCHEMFILES
/// Any trajectory format read by chemfiles (--itraj). chemfiles guesses the format from the extension
/// of the file (including the compression: .gz, .xz, .bz2), unless it is given (--traj-format).
/// The step, when the file stores it, and the time are passed as for xtc and trr files.
/// Masses and charges are passed only for the formats that store both of them (TPR, LAMMPS Data):
/// for the other formats chemfiles guesses the masses from the atom names and sets the charges to zero
class chemfilesParser final : public fileParser {
  std::string trajectoryFile;
  /// the format given by the user, empty to let chemfiles guess it
  std::string format;
  std::unique_ptr<chemfiles::Trajectory> trajectory;
  /// the frame whose header was read last
  chemfiles::Frame frame;
  long long frameNumber{0};
  bool massesAndCharges{false};
  std::string where() const {
    return "frame "+std::to_string(frameNumber)+" of the trajectory file "+trajectoryFile;
  }
  /// the name of the format (e.g. "XTC"), from the format given by the user or from the extension of the file
  std::string formatName() const {
    if(!format.empty()) {
      // the name is before the compression, as in "XYZ / GZ"
      std::string name=format.substr(0,format.find('/'));
      Tools::trim(name);
      return name;
    }
    std::string fname=trajectoryFile;
    const std::array<std::string,3> compressions{".gz",".xz",".bz2"};
    for(const auto & compression : compressions) {
      if(fname.size()>compression.size() && fname.compare(fname.size()-compression.size(),compression.size(),compression)==0) {
        fname.resize(fname.size()-compression.size());
        break;
      }
    }
    const auto dot=fname.rfind('.');
    if(dot==std::string::npos) {
      return "";
    }
    const std::string extension=fname.substr(dot);
    for(const auto & metadata : chemfiles::formats_list()) {
      if(metadata.get().extension && *metadata.get().extension==extension) {
        return metadata.get().name;
      }
    }
    return "";
  }
public:
  explicit chemfilesParser(std::string_view fmt):
    format(fmt) {
  }
  std::optional<std::string> init(std::string_view /*fmt*/,
                                  std::string_view fname) override {
    trajectoryFile=std::string(fname);
    frameNumber=0;
    natoms=-1;
    try {
      trajectory=std::make_unique<chemfiles::Trajectory>(trajectoryFile,'r',format);
    } catch(const std::exception & e) {
      return "cannot open the trajectory file "+trajectoryFile+" with chemfiles: "+e.what();
    }
    const std::string name=formatName();
    massesAndCharges=(name=="TPR" || name=="LAMMPS Data");
    return std::nullopt;
  }
  std::optional<std::string> rewind() override {
    const int nat=natoms;
    auto error=init("",trajectoryFile);
    natoms=nat;
    return error;
  }
  std::optional<std::string> readHeader(RawFrame& /*raw*/) override {
    // chemfiles reads the whole frame: it is kept until readAtoms()
    try {
      if(trajectory->done()) {
        return "EOF";
      }
      ++frameNumber;
      frame=trajectory->read();
    } catch(const std::exception & e) {
      return "cannot read "+where()+" with chemfiles: "+e.what();
    }
    const auto n=frame.size();
    if(n==0) {
      return where()+" has no atoms";
    }
    if(natoms>=0 && n!=std::size_t(natoms)) {
      return where()+" has "+std::to_string(n)+" atoms, the first frame has "+std::to_string(natoms);
    }
    if(n>std::size_t(std::numeric_limits<int>::max())) {
      return where()+" has too many atoms";
    }
    natoms=int(n);
    return std::nullopt;
  }
  std::optional<std::string> readAtoms(RawFrame& raw) override {
    // chemfiles works in A, ps and degrees
    const auto & pos=frame.positions();
    for(int i=0; i<natoms; i++) {
      for(unsigned j=0; j<3; j++) {
        raw.positions[3*i+j]=pos[i][j]/10.;
      }
    }
    if(raw.readBox) {
      const auto matrix=frame.cell().matrix();
      for(unsigned i=0; i<3; i++) {
        for(unsigned j=0; j<3; j++) {
#if CHEMFILES_VERSION_MAJOR==0 && CHEMFILES_VERSION_MINOR<11
          // chemfiles 0.10: the cell vectors are the columns of the matrix
          raw.box[3*i+j]=matrix[j][i]/10.;
#else
          raw.box[3*i+j]=matrix[i][j]/10.;
#endif
        }
      }
    }
#if CHEMFILES_VERSION_MAJOR==0 && CHEMFILES_VERSION_MINOR<11
    // chemfiles 0.10: the step of the simulation for the formats that store it, the index of the frame otherwise
    raw.step=static_cast<long long>(frame.step());
#else
    if(const auto step=frame.get("simulation_step")) {
      if(step->kind()==chemfiles::Property::DOUBLE) {
        raw.step=std::llround(step->as_double());
      }
    }
#endif
    if(const auto time=frame.get("time")) {
      if(time->kind()==chemfiles::Property::DOUBLE) {
        raw.time=time->as_double();
      }
    }
    if(massesAndCharges) {
      raw.masses.resize(natoms);
      raw.charges.resize(natoms);
      for(int i=0; i<natoms; i++) {
        raw.masses[i]=frame[i].mass();
        raw.charges[i]=frame[i].charge();
      }
    }
    return std::nullopt;
  }
};
#endif

std::unique_ptr<fileParser> parserFromFMT(TrajectoryParser::trajfmt fmt,[[maybe_unused]] std::string_view format) {
  using trajfmt = TrajectoryParser::trajfmt;
  switch (fmt) {
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
  case trajfmt::pdb: {
    return std::make_unique<pdbParser>();
  }
  case trajfmt::dcd: {
    return std::make_unique<dcdParser>();
  }
#ifdef __PLUMED_HAS_LIBCHEMFILES
  case trajfmt::chemfiles: {
    return std::make_unique<chemfilesParser>(format);
  }
#endif
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
  keys.add("atoms","--ipdb","the trajectory in pdb format (frames separated by END or ENDMDL records)");
  keys.add("atoms","--idcd","the trajectory in dcd format (CHARMM, NAMD, X-PLOR)");
  // without chemfiles the options are hidden (not in the help and in the manual), and --itraj stops with a message (see init())
  keys.add(hasChemfiles() ? "atoms" : "hidden","--itraj","the trajectory in any format read by chemfiles (needs PLUMED configured with --enable-libchemfiles): "
           "the format is guessed from the extension of the file, unless it is given with --traj-format");
  keys.add(hasChemfiles() ? "optional" : "hidden","--traj-format","with --itraj, the format of the trajectory as named by chemfiles (e.g. \"Amber NetCDF\", \"LAMMPS\", \"XYZ / GZ\"), "
           "when it cannot be guessed from the extension of the file");
  for(const auto & opt : deprecatedOptions()) {
    keys.addDeprecatedKeyword(opt.first,"--i"+opt.second);
  }
  for(const auto & opt : removedOptions()) {
    keys.add("hidden",opt,removedOptionMessage(opt));
  }
}

TrajectoryParser::TrajectoryParser():
  raw(std::make_unique<RawFrame>()) {
}
TrajectoryParser::~TrajectoryParser()=default;
std::vector<std::string> TrajectoryParser::trajectoryOptions() {
  return {"xyz","gro","dlp4","xtc","trr","pdb","dcd","traj"};
}

const std::vector<std::pair<std::string,std::string>> & TrajectoryParser::deprecatedOptions() {
  static const std::vector<std::pair<std::string,std::string>> options {
    {"--mf_xtc","xtc"},
    {"--mf_trr","trr"},
    {"--mf_gro","gro"},
    {"--mf_pdb","pdb"},
    {"--mf_dcd","dcd"}
  };
  return options;
}

const std::vector<std::string> & TrajectoryParser::removedOptions() {
  static const std::vector<std::string> options {"--mf_crd","--mf_crdbox","--mf_trj","--mf_g96"};
  return options;
}

std::string TrajectoryParser::removedOptionMessage(std::string_view option) {
  std::string message="the option "+std::string(option)+" was removed together with the molfile plugins (v2.11): ";
  if(option=="--mf_crd" || option=="--mf_crdbox") {
    message+="convert the AMBER trajectory, e.g. to xtc or dcd with cpptraj, and read it with --ixtc or --idcd";
  } else if(option=="--mf_g96") {
    message+="convert the trajectory, e.g. to xtc with gmx trjconv, and read it with --ixtc";
  } else {
    message+="convert the trajectory, e.g. to xtc, and read it with --ixtc";
  }
  return message;
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
  if(fmt=="pdb") {
    return trajfmt::pdb;
  }
  if(fmt=="dcd") {
    return trajfmt::dcd;
  }
  if(fmt=="traj") {
    return trajfmt::chemfiles;
  }
  return trajfmt::error;
}
std::string TrajectoryParser::toString(trajfmt fmt) {
  switch (fmt) {
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
  case trajfmt::pdb: {
    return "pdb";
  }
  case trajfmt::dcd: {
    return "dcd";
  }
  case trajfmt::chemfiles: {
    return "traj";
  }
  case trajfmt::error: {
    return "error";
  }
  }

  return "unknown";
}

bool TrajectoryParser::hasChemfiles() {
#ifdef __PLUMED_HAS_LIBCHEMFILES
  return true;
#else
  return false;
#endif
}

std::optional<std::string> TrajectoryParser::init(std::string_view fmt, FILE* in) {
  auto trajectory_fmt =FMTfromString(fmt);
  if (trajectory_fmt==trajfmt::error) {
    return "Error: wrong input format specified";
  }
  if (trajectory_fmt==trajfmt::chemfiles) {
    return "a trajectory read with chemfiles (--itraj) cannot be read from standard input";
  }
  parser = parserFromFMT(trajectory_fmt,"");
  //not owning, will not be deleted on destruction
  return parser->init(in);
}
std::optional<std::string> TrajectoryParser::init(std::string_view fmt,
    std::string_view fname,
    std::string_view format) {
  auto trajectory_fmt =FMTfromString(fmt);
  if (trajectory_fmt==trajfmt::error) {
    return "Error: wrong input format specified";
  }
  if (trajectory_fmt==trajfmt::chemfiles && !hasChemfiles()) {
    return "cannot read the trajectory "+std::string(fname)+" with --itraj: PLUMED was compiled without chemfiles (configure it with --enable-libchemfiles)";
  }
  if (trajectory_fmt!=trajfmt::chemfiles && !format.empty()) {
    return "the format of the trajectory (--traj-format) can only be given with --itraj";
  }

  parser = parserFromFMT(trajectory_fmt,format);
  return parser->init(fmt,fname);
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
  raw->warning.reset();
  auto error=parser->readAtoms(*raw);
  frame.warning=raw->warning;
  if(error) {
    return error;
  }
  // e.g. the cell of a dcd or pdb file with angles that do not make a cell (a box that is not finite made the Pbc class loop forever)
  if(!std::all_of(raw->box.begin(),raw->box.end(),[](double x) {
  return std::isfinite(x);
  })) {
    return "the box of a frame of the trajectory is not finite (the cell in the file is not valid)";
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
