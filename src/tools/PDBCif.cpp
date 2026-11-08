/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2026 The plumed team
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
#include "PDB.h"
#include "Tools.h"
#include "Exception.h"
#include "gemmi/cif.h"
#include <array>
#include <cctype>
#include <cmath>
#include <functional>
#include <limits>
#include <map>

// Reader of mmCIF files for the PDB class (the documentation is in PDB.cpp).
//
// The file is parsed with the grammar of gemmi (gemmi::cif::rules), with our own
// PEGTL actions: the values of the _atom_site loop are used as soon as a row is
// complete and nothing else is stored, so that the memory needed is that of the PDB
// object (gemmi's Document would hold every value of the file as a string).

namespace PLMD {

namespace {

namespace pegtl = tao::pegtl;
namespace rules = gemmi::cif::rules;

/// the columns of _atom_site that are read
enum Field : unsigned {
  ID, AUTH_ATOM, LABEL_ATOM, AUTH_COMP, LABEL_COMP, AUTH_ASYM, LABEL_ASYM, AUTH_SEQ, LABEL_SEQ,
  INS_CODE, ALT_ID, TYPE_SYMBOL, X, Y, Z, OCCUPANCY, B_ISO, MODEL, NFIELDS
};

const std::map<std::string,unsigned> & atomSiteFields() {
  static const std::map<std::string,unsigned> fields= {
    {"id",ID},
    {"auth_atom_id",AUTH_ATOM}, {"label_atom_id",LABEL_ATOM},
    {"auth_comp_id",AUTH_COMP}, {"label_comp_id",LABEL_COMP},
    {"auth_asym_id",AUTH_ASYM}, {"label_asym_id",LABEL_ASYM},
    {"auth_seq_id",AUTH_SEQ}, {"label_seq_id",LABEL_SEQ},
    {"pdbx_pdb_ins_code",INS_CODE}, {"label_alt_id",ALT_ID}, {"type_symbol",TYPE_SYMBOL},
    {"cartn_x",X}, {"cartn_y",Y}, {"cartn_z",Z},
    {"occupancy",OCCUPANCY}, {"b_iso_or_equiv",B_ISO}, {"pdbx_pdb_model_num",MODEL}
  };
  return fields;
}

std::string lower(std::string s) {
  for(auto & c : s) {
    c=static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }
  return s;
}

/// The values of a CIF file: ? and . (unquoted) are null, quotes and text-field delimiters are removed
bool isNull(const std::string& raw) {
  return raw.empty() || gemmi::cif::is_null(raw);
}

std::string value(const std::string& raw) {
  return gemmi::cif::as_string(raw);
}

/// What the actions fill while the file is parsed
struct CifReader {
  std::string file;
  /// called for each row of _atom_site, with the line where the row starts
  std::function<void(std::size_t)> addRow;

  unsigned nblocks=0;
  bool inFrame=false;
  /// the current loop
  bool inLoop=false;
  bool atomSiteLoop=false;
  bool atomSiteFound=false;
  std::vector<std::string> loopTags;
  std::vector<int> columnField;
  std::size_t nvalues=0;
  std::array<std::string,NFIELDS> row;
  std::array<bool,NFIELDS> present{};
  std::size_t rowLine=0;
  /// the current item (tag value pair)
  std::string itemTag;
  /// _cell and space group: tag -> value, from the first data block
  std::map<std::string,std::string> items;

  bool active() const {
    return nblocks==1 && !inFrame;
  }
  static bool interestingItem(const std::string& tag) {
    return Tools::startWith(tag,"_cell.") || Tools::startWith(tag,"_symmetry.") || Tools::startWith(tag,"_space_group.");
  }
  void startLoop() {
    inLoop=true;
    atomSiteLoop=false;
    loopTags.clear();
    columnField.clear();
    nvalues=0;
  }
  void loopTag(const std::string& tag) {
    loopTags.push_back(lower(tag));
    const std::string & t=loopTags.back();
    int f=-1;
    if(Tools::startWith(t,"_atom_site.")) {
      if(atomSiteFound && loopTags.size()==1) {
        plumed_merror("the mmCIF file " + file + " has more than one _atom_site table");
      }
      atomSiteLoop=true;
      const auto p=atomSiteFields().find(t.substr(11));
      if(p!=atomSiteFields().end()) {
        f=static_cast<int>(p->second);
      }
    }
    columnField.push_back(f);
  }
  void loopValue(const std::string& raw,std::size_t line) {
    // line is the line of the value (only used for the first value of a row)
    const std::size_t ncol=loopTags.size();
    const std::size_t col=nvalues%ncol;
    nvalues++;
    if(atomSiteLoop) {
      if(col==0) {
        rowLine=line;
        present.fill(false);
      }
      if(columnField[col]>=0) {
        const auto field=static_cast<std::size_t>(columnField[col]);
        row[field]=raw;
        present[field]=true;
      }
      if(col==ncol-1) {
        addRow(rowLine);
      }
    } else if(nvalues<=ncol) {
      // _cell and space group given as a table: the first row
      if(interestingItem(loopTags[col])) {
        items[loopTags[col]]=raw;
      }
    }
  }
  void endLoop(std::size_t line) {
    if(atomSiteLoop) {
      atomSiteFound=true;
      if(nvalues%loopTags.size()!=0) {
        plumed_merror("wrong number of values in the _atom_site table of the mmCIF file " + file
                      + " (" + std::to_string(nvalues) + " values for " + std::to_string(loopTags.size())
                      + " columns, table starting at line " + std::to_string(line) + ")");
      }
    }
    inLoop=false;
    atomSiteLoop=false;
  }
};

template<typename Rule> struct Action : pegtl::nothing<Rule> {};

template<> struct Action<rules::datablockname> {
  template<typename Input> static void apply(const Input&, CifReader& r) {
    r.nblocks++;
  }
};
template<> struct Action<rules::str_global> {
  template<typename Input> static void apply(const Input&, CifReader& r) {
    r.nblocks++;
  }
};
template<> struct Action<rules::framename> {
  template<typename Input> static void apply(const Input&, CifReader& r) {
    r.inFrame=true;
  }
};
template<> struct Action<rules::endframe> {
  template<typename Input> static void apply(const Input&, CifReader& r) {
    r.inFrame=false;
  }
};
template<> struct Action<rules::item_tag> {
  template<typename Input> static void apply(const Input& in, CifReader& r) {
    if(r.active()) {
      r.itemTag=lower(in.string());
      if(Tools::startWith(r.itemTag,"_atom_site.")) {
        plumed_merror("the _atom_site data of the mmCIF file " + r.file + " must be a table (loop_): line "
                      + std::to_string(in.iterator().line));
      }
    }
  }
};
template<> struct Action<rules::item_value> {
  template<typename Input> static void apply(const Input& in, CifReader& r) {
    if(r.active() && CifReader::interestingItem(r.itemTag)) {
      r.items[r.itemTag]=in.string();
    }
  }
};
template<> struct Action<rules::str_loop> {
  template<typename Input> static void apply(const Input&, CifReader& r) {
    if(r.active()) {
      r.startLoop();
    }
  }
};
template<> struct Action<rules::loop_tag> {
  template<typename Input> static void apply(const Input& in, CifReader& r) {
    if(r.active()) {
      r.loopTag(in.string());
    }
  }
};
template<> struct Action<rules::loop_value> {
  template<typename Input> static void apply(const Input& in, CifReader& r) {
    if(r.active()) {
      r.loopValue(in.string(),in.iterator().line);
    }
  }
};
template<> struct Action<rules::loop> {
  template<typename Input> static void apply(const Input& in, CifReader& r) {
    if(r.active() && r.inLoop) {
      r.endLoop(in.iterator().line);
    }
  }
};

}

void PDB::readCif(const std::string& file,double scale) {
  CifReader r;
  r.file=file;

  std::string firstModel;
  bool firstRow=true;
  unsigned nread=0;
  const auto & f=r.row;
  const auto & p=r.present;
  // the value of a column, empty if the column is missing or the value is null
  auto get=[&](unsigned field) {
    return (p[field] && !isNull(f[field])) ? value(f[field]) : std::string();
  };
  // the author value, or the label value if the author value is missing
  auto getAuthOrLabel=[&](unsigned auth,unsigned label) {
    return (p[auth] && !isNull(f[auth])) ? value(f[auth]) : get(label);
  };
  auto fail=[&](std::size_t line,const std::string& msg) {
    plumed_merror("in the _atom_site table of the mmCIF file " + file + " at line " + std::to_string(line) + ": " + msg);
  };
  auto getReal=[&](unsigned field,const std::string& name,double defaultValue,std::size_t line) {
    if(!p[field] || isNull(f[field])) {
      if(std::isnan(defaultValue)) {
        fail(line,"the " + name + " is missing");
      }
      return defaultValue;
    }
    double v=0.0;
    if(!Tools::convertNoexcept(value(f[field]),v)) {
      fail(line,"cannot read the " + name + " from " + f[field]);
    }
    return v;
  };

  r.addRow=[&](std::size_t line) {
    if(firstRow) {
      if(!p[X] || !p[Y] || !p[Z]) {
        plumed_merror("the _atom_site table of the mmCIF file " + file + " has no Cartn_x, Cartn_y or Cartn_z column");
      }
      if(numbering==AtomNumbering::serial && !p[ID]) {
        plumed_merror("the _atom_site table of the mmCIF file " + file + " has no id column: use ATOM_ID=order in MOLINFO");
      }
      firstModel=get(MODEL);
      firstRow=false;
    }
    // only the first model
    if(get(MODEL)!=firstModel) {
      return;
    }
    std::string chainid=getAuthOrLabel(AUTH_ASYM,LABEL_ASYM);
    if(chainid.empty()) {
      chainid=" "; // as a blank chain in a PDB file
    }
    int resnum=0;
    {
      const std::string s=getAuthOrLabel(AUTH_SEQ,LABEL_SEQ);
      if(!s.empty() && !Tools::convertNoexcept(s,resnum)) {
        fail(line,"cannot read the residue number from " + s);
      }
    }
    const std::string icode=get(INS_CODE);
    if( !keepAltLoc(chainid,resnum,icode,get(ALT_ID)) ) {
      skippedAltLocs++;
      return;
    }
    AtomNumber a;
    if(numbering==AtomNumbering::serial) {
      const std::string s=get(ID);
      long long id=0;
      if(!Tools::convertNoexcept(s,id) || id<1 || id>std::numeric_limits<int>::max()) {
        fail(line,"the atom id " + (s.empty() ? std::string("is missing") : s + " is not a positive integer")
             + ": use ATOM_ID=order in MOLINFO to number the atoms in the order of the file");
      }
      a.setSerial(static_cast<unsigned>(id));
      if(checkForAtom(a)) {
        fail(line,"the atom id " + s + " is used twice: use ATOM_ID=order in MOLINFO to number the atoms in the order of the file");
      }
    } else {
      a.setIndex(nread);
    }
    Vector pos;
    const double nan=std::numeric_limits<double>::quiet_NaN();
    pos[0]=getReal(X,"x coordinate (Cartn_x)",nan,line);
    pos[1]=getReal(Y,"y coordinate (Cartn_y)",nan,line);
    pos[2]=getReal(Z,"z coordinate (Cartn_z)",nan,line);
    const double o=getReal(OCCUPANCY,"occupancy",1.0,line);
    const double b=getReal(B_ISO,"beta factor (B_iso_or_equiv)",1.0,line);
    std::string resname=getAuthOrLabel(AUTH_COMP,LABEL_COMP);
    // as in the columns 18-20 of a PDB file, where short names are right-justified
    if(resname.length()<3) {
      resname=std::string(3-resname.length(),' ')+resname;
    }
    addAtom(a,getAuthOrLabel(AUTH_ATOM,LABEL_ATOM),resname,resnum,chainid,icode,get(TYPE_SYMBOL),pos*scale,o,b);
    nread++;
  };

  try {
    pegtl::file_input<> in(file);
    if(!pegtl::parse<rules::file,Action,gemmi::cif::Errors>(in,r)) {
      plumed_merror("the mmCIF file " + file + " is empty");
    }
  } catch(const pegtl::parse_error& e) {
    plumed_merror(std::string("cannot read the mmCIF file: ") + e.what());
  } catch(const std::system_error& e) {
    plumed_merror("cannot read the mmCIF file " + file + ": " + e.what());
  }
  if(!r.atomSiteFound) {
    plumed_merror("the mmCIF file " + file + " has no _atom_site table");
  }
  if(size()==0) {
    plumed_merror("the _atom_site table of the mmCIF file " + file + " has no atoms");
  }
  block_ends.push_back( size() );

  // unit cell, as the CRYST1 record of a PDB file
  const std::array<std::string,6> cellTags= {"_cell.length_a","_cell.length_b","_cell.length_c",
                                             "_cell.angle_alpha","_cell.angle_beta","_cell.angle_gamma"
                                            };
  unsigned ncell=0;
  for(const auto & t : cellTags) {
    const auto it=r.items.find(t);
    if(it!=r.items.end() && !isNull(it->second)) {
      ncell++;
    }
  }
  if(ncell>0) {
    std::array<double,6> cell;
    for(unsigned i=0; i<6; i++) {
      const auto it=r.items.find(cellTags[i]);
      if(it==r.items.end() || isNull(it->second) || !Tools::convertNoexcept(value(it->second),cell[i]) || !std::isfinite(cell[i])) {
        plumed_merror("the unit cell of the mmCIF file " + file + " cannot be read: " + cellTags[i] + " is "
                      + (it==r.items.end() ? std::string("missing") : it->second));
      }
    }
    for(unsigned i=0; i<3; i++) {
      if(!(cell[i]>0.0)) {
        plumed_merror("the unit cell of the mmCIF file " + file + " is not valid: " + cellTags[i] + " must be positive");
      }
      if(!(cell[3+i]>0.0 && cell[3+i]<180.0)) {
        plumed_merror("the unit cell of the mmCIF file " + file + " is not valid: " + cellTags[3+i] + " must be between 0 and 180 degrees");
      }
    }
    BoxXYZ=Vector(cell[0],cell[1],cell[2]);
    BoxXYZ*=scale;
    BoxABG=Vector(cell[3],cell[4],cell[5]);
    setBoxFromLengthsAndAngles();
    if(!(Box[2][2]>0.0)) {
      plumed_merror("the unit cell of the mmCIF file " + file + " is not valid: the angles do not define a cell");
    }
    const std::array<std::string,2> spaceGroupTags= {"_symmetry.space_group_name_h-m","_space_group.name_h-m_alt"};
    for(const auto & t : spaceGroupTags) {
      const auto it=r.items.find(t);
      if(it!=r.items.end() && !isNull(it->second)) {
        spaceGroup=value(it->second);
        Tools::trim(spaceGroup);
        break;
      }
    }
  }
}

}
