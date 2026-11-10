/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#include "FileTools.h"
#include "Log.h"
#include "h36.h"
#include <cctype>
#include <cstdio>
#include <set>
#include "Tensor.h"

//+PLUMEDOC INTERNAL pdbreader
/*
PLUMED can use the PDB format in several places

- To read molecular structure (\ref MOLINFO).
- To read reference conformations (\ref RMSD, but also many other methods in \ref dists, \ref FIT_TO_TEMPLATE, etc).

The implemented PDB reader expects a file formatted correctly according to the
[PDB standard](http://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html).
In particular, the following columns are read from ATOM records
\verbatim
columns | content
1-6     | record name (ATOM or HETATM)
7-11    | serial number of the atom (starting from 1)
13-16   | atom name
17      | alternate location
18-20   | residue name (and column 21, if not blank, for residue names of four characters, as written by GROMACS)
22      | chain id
23-26   | residue number
27      | insertion code
31-38   | x coordinate
39-46   | y coordinate
47-54   | z coordinate
55-60   | occupancy
61-66   | beta factor
77-78   | element symbol
\endverbatim

REMARK records are read only when they are not numbered: the numbered remarks of the wwPDB (`REMARK 465 ...`) are free text.

The serial numbers of a PDB file written by the wwPDB count the TER records, so that they skip a number at the end of each chain.
\ref MOLINFO has the keyword ATOM_ID: with ATOM_ID=order the atoms are numbered 1,2,3... in the order in which they are read,
and the serial numbers in the file are not used.

When atoms have alternate locations (column 17), only the first location found in each residue is read,
together with the atoms of the residue that have no alternate location; the other atoms are skipped
(\ref MOLINFO writes how many in the log). Residues are identified by chain, residue number and insertion code,
so that residues 52 and 52A of the same chain are different residues.

\par mmCIF files

Files with extension `.cif` or `.mmcif` are read as mmCIF (PDBx/mmCIF) files, with the CIF parser of
[gemmi](https://gemmi.readthedocs.io). The atoms are read from the `_atom_site` table of the first data block:

\verbatim
column                                      | content
_atom_site.id                               | serial number of the atom (see below)
_atom_site.auth_atom_id, label_atom_id      | atom name
_atom_site.auth_comp_id, label_comp_id      | residue name
_atom_site.auth_asym_id, label_asym_id      | chain id
_atom_site.auth_seq_id, label_seq_id        | residue number
_atom_site.pdbx_PDB_ins_code                | insertion code
_atom_site.label_alt_id                     | alternate location
_atom_site.type_symbol                      | element symbol
_atom_site.Cartn_x, Cartn_y, Cartn_z        | coordinates
_atom_site.occupancy                        | occupancy
_atom_site.B_iso_or_equiv                   | beta factor
_atom_site.pdbx_PDB_model_num               | model number: only the first model is read
\endverbatim

The author (`auth_`) columns are used, as in PDB files written from the same structure; the `label_` columns
are used for the atoms where the `auth_` value is missing (or when the column is not there). Chain ids and residue names
are read in full (chain ids with more than one character, residue names with four or five characters, which a PDB file
truncates to the three characters of columns 18-20); residue names shorter than three characters are right-justified
to three, as in the columns 18-20 of a PDB file.
Only the coordinates are compulsory: a missing occupancy or beta factor is 1, a missing residue number is 0.
The unit cell is read from `_cell.length_a`, `_cell.length_b`, `_cell.length_c`, `_cell.angle_alpha`,
`_cell.angle_beta` and `_cell.angle_gamma` (as the CRYST1 record of a PDB file), the space group from
`_symmetry.space_group_name_H-M` or `_space_group.name_H-M_alt`.
The PLUMED-specific REMARK lines of PDB files (arguments, TYPE, ...) cannot be given in an mmCIF file.

By default the serial number of an atom is the value of `_atom_site.id`, which must then be a positive integer.
With ATOM_ID=order in \ref MOLINFO the atoms are numbered 1,2,3... in the order in which they are read (after skipping
the other models and the alternate locations), for files where the ids do not match the numbering of the atoms in the simulation.
The same numbering is then obtained from the PDB and the mmCIF file of a wwPDB entry (whose ids do not count the TER records).

\par GROMACS TPR files

When PLUMED is configured with `--enable-libchemfiles`, files with extension `.tpr` are read as GROMACS TPR files,
with the [chemfiles](https://chemfiles.org) library (0.10.4 or later). The structure is that of the PDB file written by
`gmx editconf -f topol.tpr -o reference.pdb`:
- the atoms are numbered 1,2,3... in the order of the file, as in the simulation (ATOM_ID has no effect);
- the atom names, the residue names, numbers and insertion codes, the elements and the coordinates are those of the TPR file;
- a TPR file has no chain ids: the sequences of residues in which each residue is bonded to the next one (e.g. the protein and
  nucleic acid chains) are the chains A, B, C... (then AA, AB...) in the order of the file, and keep their residue numbers;
  the other residues (water, ions...) have no chain and are numbered after the largest residue number of the chains.
  This is what `gmx editconf` does with the molecules of more than one residue: the result is the same, unless a molecule
  of the topology has residues that are not bonded to each other (e.g. two chains in the same molecule, which are here two chains);
- the occupancy and the beta factor of all the atoms are 1;
- the box of the simulation is the unit cell (space group P 1).
PLUMED parser is slightly more permissive than the official PDB format
in the fact that the format of real numbers is not fixed. In other words,
any real number that can be parsed is OK and the dot can be placed anywhere. However,
__columns are interpret strictly__. A sample PDB should look like the following
\verbatim
ATOM      2  CH3 ACE     1      12.932 -14.718  -6.016  1.00  1.00
ATOM      5  C   ACE     1      21.312  -9.928  -5.946  1.00  1.00
ATOM      9  CA  ALA     2      19.462 -11.088  -8.986  1.00  1.00
\endverbatim

Notice that serial numbers need not to be consecutive. In the three-line example above,
only the coordinates of three atoms are provided. This is perfectly legal and indicates PLUMED
that information about these atoms only is available. This could be both for structural
information in \ref MOLINFO, where the other atoms would have no name assigned, and for
reference structures used in \ref RMSD, where only the provided atoms would be used to compute RMSD.

\par Occupancy and beta factors

PLUMED reads also occupancy and beta factors that however are given a very special meaning.
In cases where the PDB structure is used as a reference for an alignment (that's the case
for instance in \ref RMSD and in \ref FIT_TO_TEMPLATE), the occupancy column is used
to provide the weight of each atom in the alignment. In cases where, perhaps after alignment,
the displacement between running coordinates and the provided PDB is computed, the beta factors
are used as weight for the displacement.
Since setting the weights to zero is the same as __not__ including an atom in the alignment or
displacement calculation, the two following reference files would be equivalent when used in an \ref RMSD
calculation. First file:
\verbatim
ATOM      2  CH3 ACE     1      12.932 -14.718  -6.016  1.00  1.00
ATOM      5  C   ACE     1      21.312  -9.928  -5.946  1.00  1.00
ATOM      9  CA  ALA     2      19.462 -11.088  -8.986  0.00  0.00
\endverbatim
Second file:
\verbatim
ATOM      2  CH3 ACE     1      12.932 -14.718  -6.016  1.00  1.00
ATOM      5  C   ACE     1      21.312  -9.928  -5.946  1.00  1.00
\endverbatim
However notice that many extra atoms with zero weight might slow down the calculation, so
removing lines is better than setting their weights to zero.
In addition, weights for alignment need not to be equivalent to weights for displacement.
Starting with PLUMED 2.7, if all the weights are set to zero they will be normalized to be equal to the
inverse of the number of involved atoms. This means that it will be possible to use files with
the weight columns set to zero obtaining a meaningful result. In previous PLUMED versions,
setting all weights to zero was resulting in an error instead.


\par Systems with more than 100k atoms

Notice that it very likely does not make any sense to compute the \ref RMSD or any other structural
deviation __using__ so many atoms. However, if the protein for which you want to compute \ref RMSD
has atoms with large serial numbers (e.g. because it is located __after__ solvent in the sorted list of atoms)
you might end up with troubles with the limitations of the PDB format. Indeed, since there are 5
columns available for atom serial number, this number cannot be larger than 99999.
In addition, providing \ref MOLINFO with names associated to atoms with a serial larger than 99999 would be impossible.

Since PLUMED 2.4 we allow [hybrid 36](http://cci.lbl.gov/hybrid_36/) format to be used to specify atom numbers.
This format is not particularly widespread, but has the nice feature that it provides a one-to-one mapping
between numbers up to approximately 80 millions and strings with 5 characters, plus it is backward compatible
for numbers smaller than 100000. This is not true for notations like the hex notation exported by VMD.
Using the hybrid 36 format, the ATOM records for atom ranging from 99997 to 100002 would read like these:
\verbatim
ATOM  99997  Ar      X   1      45.349  38.631  15.116  1.00  1.00
ATOM  99998  Ar      X   1      46.189  38.631  15.956  1.00  1.00
ATOM  99999  Ar      X   1      46.189  39.471  15.116  1.00  1.00
ATOM  A0000  Ar      X   1      45.349  39.471  15.956  1.00  1.00
ATOM  A0000  Ar      X   1      45.349  38.631  16.796  1.00  1.00
ATOM  A0001  Ar      X   1      46.189  38.631  17.636  1.00  1.00
\endverbatim
There are tools that can be found to translate from integers to strings and back using hybrid 36 format
(a simple python script can be found [here](https://sourceforge.net/p/cctbx/code/HEAD/tree/trunk/iotbx/pdb/hybrid_36.py)).
In addition, as of PLUMED 2.5, we provide a \ref pdbrenumber "command line tool" that can be used to renumber atoms in a PDB file.

*/
//+ENDPLUMEDOC


namespace PLMD {

namespace {
/// the string without leading and trailing blanks (empty if it is all blanks)
std::string withoutBlanks(const std::string& s) {
  const auto first=s.find_first_not_of(" \t");
  if(first==std::string::npos) {
    return "";
  }
  return s.substr(first,s.find_last_not_of(" \t")-first+1);
}
}

bool PDB::getArgumentValue( const std::string& name, std::vector<double>& value ) const {
  std::map<std::string,std::vector<double> >::const_iterator it = arg_data.find(name);
  if( it!=arg_data.end() ) {
    if( value.size()!=it->second.size() ) {
      return false;
    }
    for(unsigned i=0; i<value.size(); ++i) {
      value[i] = it->second[i];
    }
    return true;
  }
  return false;
}

// bool PDB::hasRequiredProperties( const std::vector<std::string>& inproperties ){
//   bool hasprop=false;
//   for(unsigned i=0;i<remark.size();++i){
//       if( remark[i].find("PROPERTIES=")!=std::string::npos){ hasprop=true; break; }
//   }
//   if( !hasprop ){
//       std::string mypropstr="PROPERTIES=" + inproperties[0];
//       for(unsigned i=1;i<inproperties.size();++i) mypropstr += "," + inproperties[i];
//       remark.push_back( mypropstr );
//   }
//   // Now check that all required properties are there
//   for(unsigned i=0;i<inproperties.size();++i){
//       hasprop=false;
//       for(unsigned j=0;j<remark.size();++j){
//           if( remark[j].find(inproperties[i]+"=")!=std::string::npos){ hasprop=true; break; }
//       }
//       if( !hasprop ) return false;
//   }
//   return true;
// }

unsigned PDB::getNumberOfAtomBlocks()const {
  return block_ends.size();
}

const std::vector<unsigned> & PDB::getAtomBlockEnds()const {
  return block_ends;
}

const std::vector<Vector> & PDB::getPositions()const {
  return positions;
}

void PDB::setPositions(const std::vector<Vector> &v ) {
  plumed_assert( v.size()==positions.size() );
  positions=v;
}

const std::vector<double> & PDB::getOccupancy()const {
  return occupancy;
}

const std::vector<double> & PDB::getBeta()const {
  return beta;
}

void PDB::addRemark( std::vector<std::string>& v1 ) {
  Tools::parse(v1,"TYPE",mtype);
  Tools::parseVector(v1,"ARG",argnames);
  for(unsigned i=0; i<v1.size(); ++i) {
    if( v1[i].find("=")!=std::string::npos ) {
      std::size_t eq=v1[i].find_first_of('=');
      std::string name=v1[i].substr(0,eq);
      std::string sval=v1[i].substr(eq+1);
      std::vector<std::string> words=Tools::getWords(sval,"\t\n ,");
      std::vector<double> val( words.size() );
      for(unsigned ii=0; ii<val.size(); ++ii) {
        Tools::convert( words[ii], val[ii] );
      }
      arg_data.insert( std::pair<std::string,std::vector<double> >( name, val ) );
    } else {
      flags.push_back(v1[i]);
    }
  }
}


const std::vector<AtomNumber> & PDB::getAtomNumbers()const {
  return numbers;
}

const Vector & PDB::getBoxAxs()const {
  return BoxXYZ;
}

const Vector & PDB::getBoxAng()const {
  return BoxABG;
}

const Tensor & PDB::getBoxVec()const {
  return Box;
}

const std::string & PDB::getSpaceGroup()const {
  return spaceGroup;
}

std::string PDB::getAtomName(AtomNumber a)const {
  const auto p=number2index.find(a);
  if(p==number2index.end()) {
    std::string num;
    Tools::convert( a.serial(), num );
    plumed_merror("Name of atom " + num + " not found" );
  } else {
    return atomsymb[p->second];
  }
}

int PDB::getResidueNumber(AtomNumber a)const {
  const auto p=number2index.find(a);
  if(p==number2index.end()) {
    std::string num;
    Tools::convert( a.serial(), num );
    plumed_merror("Residue for atom " + num + " not found" );
  } else {
    return residue[p->second];
  }
}

std::string PDB::getResidueName(AtomNumber a) const {
  const auto p=number2index.find(a);
  if(p==number2index.end()) {
    std::string num;
    Tools::convert( a.serial(), num );
    plumed_merror("Residue for atom " + num + " not found" );
  } else {
    return residuenames[p->second];
  }
}

unsigned PDB::size()const {
  return positions.size();
}

bool PDB::readFromFilepointer(FILE *fp,bool naturalUnits,double scale) {
  bool file_is_alive=false;
  if(naturalUnits) {
    scale=1.0;
  }
  std::string line;
  fpos_t pos;
  bool between_ters=true;
  while(Tools::getline(fp,line)) {
    fgetpos (fp,&pos);
    while(line.length()<80) {
      line.push_back(' ');
    }
    std::string record=line.substr(0,6);
    std::string serial=line.substr(6,5);
    std::string atomname=line.substr(12,4);
    std::string altloc=line.substr(16,1);
    // residue names of four characters (as written by GROMACS) use column 21;
    // shorter names keep the three columns 18-20 as they are, with their blanks
    std::string residuename=line.substr(17,3) + (line[20]!=' ' ? line.substr(20,1) : "");
    std::string chainID=line.substr(21,1);
    std::string resnum=line.substr(22,4);
    std::string icode=line.substr(26,1);
    std::string elem=line.substr(76,2);
    std::string x=line.substr(30,8);
    std::string y=line.substr(38,8);
    std::string z=line.substr(46,8);
    std::string occ=line.substr(54,6);
    std::string bet=line.substr(60,6);
    std::string BoxX=line.substr(6,9);
    std::string BoxY=line.substr(15,9);
    std::string BoxZ=line.substr(24,9);
    std::string BoxA=line.substr(33,7);
    std::string BoxB=line.substr(40,7);
    std::string BoxG=line.substr(47,7);
    Tools::trim(record);
    if(record=="TER") {
      between_ters=false;
      block_ends.push_back( positions.size() );
    }
    if(record=="END") {
      file_is_alive=true;
      break;
    }
    if(record=="ENDMDL") {
      file_is_alive=true;
      break;
    }
    if(record=="REMARK") {
      std::vector<std::string> v1;
      v1=Tools::getWords(line.substr(6));
      // numbered remarks ("REMARK 465 ...") are those of the wwPDB, with free text: only the remarks of PLUMED are read
      int remarkNumber;
      if( v1.empty() || !Tools::convertNoexcept(v1[0],remarkNumber) ) {
        addRemark( v1 );
      }
    }
    if(record=="CRYST1") {
      Tools::convert(BoxX,BoxXYZ[0]);
      Tools::convert(BoxY,BoxXYZ[1]);
      Tools::convert(BoxZ,BoxXYZ[2]);
      Tools::convert(BoxA,BoxABG[0]);
      Tools::convert(BoxB,BoxABG[1]);
      Tools::convert(BoxG,BoxABG[2]);
      spaceGroup=line.substr(55,11);
      Tools::trim(spaceGroup);
      BoxXYZ*=scale;
      setBoxFromLengthsAndAngles();
    }
    if(record=="ATOM" || record=="HETATM") {
      between_ters=true;
      AtomNumber a;
      int resno=0; // GB: when resnum string is not present, we set res number to zero
      double o,b;
      Vector p;
      // with AtomNumbering::order the serial is not read (it is set below, when the atom is kept)
      if(numbering==AtomNumbering::serial) {
        int result;
        auto trimmed=serial;
        Tools::trim(trimmed);
        while(trimmed.length()<5) {
          trimmed = std::string(" ") + trimmed;
        }
        const char* errmsg = h36::hy36decode(5, trimmed.c_str(),trimmed.length(), &result);
        if(errmsg) {
          std::string msg(errmsg);
          plumed_merror(msg);
        }
        a.setSerial(result);
      }

      // allow skipping residue number
      {
        auto trimmed=resnum;
        Tools::trim(trimmed);
        if(trimmed.length()>0) {
          int result;
          while(trimmed.length()<4) {
            trimmed = std::string(" ") + trimmed;
          }
          const char* errmsg = h36::hy36decode(4, trimmed.c_str(),trimmed.length(), &result);
          if(errmsg) {
            std::string msg(errmsg);
            plumed_merror(msg);
          }
          resno=result;
        }
      }

      // a field that cannot be read is an error that shows the line (it was an assertion that showed only the text)
      auto readField=[&line](const std::string& field, double& val, const std::string& what) {
        if( !Tools::convertNoexcept(field,val) ) {
          plumed_merror("cannot read the " + what + " from the following line of a PDB file: in a PDB file each field must be in its own columns\n"
                        + line.substr(0,line.find_last_not_of(' ')+1));
        }
      };
      // occupancy and beta can be missing (a line that ends after the coordinates): they are then 1
      // (v2.10 accepted these lines and left the values undefined)
      auto blank=[](const std::string& field) {
        return field.find_first_not_of(" \t")==std::string::npos;
      };
      o=1.0;
      b=1.0;
      if(!blank(occ)) {
        readField(occ,o,"occupancy (columns 55-60)");
      }
      if(!blank(bet)) {
        readField(bet,b,"beta (columns 61-66)");
      }
      readField(x,p[0],"x coordinate (columns 31-38)");
      readField(y,p[1],"y coordinate (columns 39-46)");
      readField(z,p[2],"z coordinate (columns 47-54)");
      // scale into nm
      p*=scale;
      icode=withoutBlanks(icode);
      if( !keepAltLoc(chainID,resno,icode,withoutBlanks(altloc)) ) {
        skippedAltLocs++;
        continue;
      }
      if(numbering==AtomNumbering::order) {
        a.setIndex(size());
      }
      addAtom(a,withoutBlanks(atomname),residuename,resno,chainID,icode,withoutBlanks(elem),p,o,b);
    }
  }
  if( between_ters ) {
    block_ends.push_back( positions.size() );
  }
  return file_is_alive;
}

namespace {
/// the extension of a file name, in lower case (empty if there is none)
std::string lowerCaseExtension(const std::string& file) {
  const auto dot=file.find_last_of('.');
  if(dot==std::string::npos) {
    return "";
  }
  std::string ext=file.substr(dot+1);
  for(auto & c : ext) {
    c=static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }
  return ext;
}
}

bool PDB::isCifFile(const std::string& file) {
  const std::string ext=lowerCaseExtension(file);
  return ext=="cif" || ext=="mmcif";
}

bool PDB::isTprFile(const std::string& file) {
  return lowerCaseExtension(file)=="tpr";
}

bool PDB::read(const std::string&file,bool naturalUnits,double scale,AtomNumbering atomNumbering) {
  numbering=atomNumbering;
  FILE* fp=fopenUTF8(file,"r");
  if(!fp) {
    return false;
  }
  if(isCifFile(file)) {
    std::fclose(fp);
    readCif(file,(naturalUnits ? 1.0 : scale));
    return true;
  }
  if(isTprFile(file)) {
    std::fclose(fp);
    readTpr(file,(naturalUnits ? 1.0 : scale));
    return true;
  }
// call fclose when exiting this function
  auto deleter=[](auto f) {
    std::fclose(f);
  };
  std::unique_ptr<FILE,decltype(deleter)> fp_deleter(fp,deleter);
  readFromFilepointer(fp,naturalUnits,scale);
  return true;
}

void PDB::setBoxFromLengthsAndAngles() {
  double cosA=std::cos(BoxABG[0]*pi/180.);
  double cosB=std::cos(BoxABG[1]*pi/180.);
  double cosG=std::cos(BoxABG[2]*pi/180.);
  double sinG=std::sin(BoxABG[2]*pi/180.);
  for (unsigned i=0; i<3; i++) {
    Box[i][0]=0.;
    Box[i][1]=0.;
    Box[i][2]=0.;
  }
  Box[0][0]=BoxXYZ[0];
  Box[1][0]=BoxXYZ[1]*cosG;
  Box[1][1]=BoxXYZ[1]*sinG;
  Box[2][0]=BoxXYZ[2]*cosB;
  Box[2][1]=(BoxXYZ[2]*BoxXYZ[1]*cosA-Box[2][0]*Box[1][0])/Box[1][1];
  Box[2][2]=std::sqrt(BoxXYZ[2]*BoxXYZ[2]-Box[2][0]*Box[2][0]-Box[2][1]*Box[2][1]);
}

bool PDB::keepAltLoc(const std::string& chainid,int resnum,const std::string& icode,const std::string& altloc) {
  if(altloc.empty()) {
    return true;
  }
  // the key cannot be confused between residues: the chain id comes last and the other fields have no blanks
  const std::string key=std::to_string(resnum)+" "+icode+" "+chainid;
  const auto p=altLocKept.insert({key,altloc}).first;
  return p->second==altloc;
}

void PDB::addAtom(AtomNumber a,const std::string& atomname,const std::string& residuename,int resnum,
                  const std::string& chainid,const std::string& icode,const std::string& elem,
                  const Vector& p,double o,double b) {
  numbers.push_back(a);
  number2index[a]=positions.size();
  atomsymb.push_back(atomname);
  residue.push_back(resnum);
  chain.push_back(chainid);
  insertion.push_back(icode);
  element.push_back(elem);
  occupancy.push_back(o);
  beta.push_back(b);
  positions.push_back(p);
  residuenames.push_back(residuename);
}

std::string PDB::getInsertionCode(AtomNumber a) const {
  const auto p=number2index.find(a);
  if(p==number2index.end()) {
    plumed_merror("Insertion code for atom " + std::to_string(a.serial()) + " not found" );
  }
  return insertion[p->second];
}

std::string PDB::getElement(AtomNumber a) const {
  const auto p=number2index.find(a);
  if(p==number2index.end()) {
    plumed_merror("Element for atom " + std::to_string(a.serial()) + " not found" );
  }
  return element[p->second];
}

unsigned PDB::getNumberOfSkippedAlternateLocations() const {
  return skippedAltLocs;
}

std::vector<std::pair<int,std::string> > PDB::getResidueIdsInChain( const std::string& chainname ) const {
  std::vector<std::pair<int,std::string> > res;
  std::set<std::pair<int,std::string> > seen;
  for(unsigned i=0; i<size(); ++i) {
    if( chain[i]==chainname && seen.insert({residue[i],insertion[i]}).second ) {
      res.emplace_back(residue[i],insertion[i]);
    }
  }
  return res;
}

long PDB::findResidue(int resnum,const std::string& chainid,const std::string& icode) const {
  for(unsigned i=0; i<size(); ++i) {
    if( residue[i]==resnum && ( chainid=="*" || chain[i]==chainid) && ( icode=="*" || insertion[i]==icode ) ) {
      return i;
    }
  }
  return -1;
}

void PDB::getAdjacentResidue(int resnum,const std::string& chainid,const std::string& icode,int step,int& adjnum,std::string& adjicode) const {
  plumed_assert(step==1 || step==-1);
  const long i0=findResidue(resnum,chainid,icode);
  if(i0>=0) {
    const std::string & c=chain[i0];
    const std::string & ic=insertion[i0];
    // the residue before (after) it in the file, in the same chain
    for(long j=i0+step; j>=0 && j<static_cast<long>(size()); j+=step) {
      if( chain[j]!=c || ( residue[j]==resnum && insertion[j]==ic ) ) {
        continue;
      }
      if( residue[j]==resnum+step || ( residue[j]==resnum && insertion[j]!=ic ) ) {
        adjnum=residue[j];
        adjicode=insertion[j];
        return;
      }
      break;
    }
  }
  // as in previous versions
  adjnum=resnum+step;
  adjicode="*";
}

void PDB::getChainNames( std::vector<std::string>& chains ) const {
  chains.resize(0);
  if(chain.empty()) {
    return;
  }
  chains.push_back( chain[0] );
  for(unsigned i=1; i<size(); ++i) {
    if( chains[chains.size()-1]!=chain[i] ) {
      chains.push_back( chain[i] );
    }
  }
}

void PDB::getResidueRange( const std::string& chainname, int& res_start, int& res_end, std::string& errmsg ) const {
  bool inres=false, foundchain=false;
  for(unsigned i=0; i<size(); ++i) {
    if( chain[i]==chainname ) {
      if(!inres) {
        if(foundchain) {
          errmsg="found second start of chain named " + chainname;
        }
        res_start=residue[i];
      }
      inres=true;
      foundchain=true;
    } else if( inres && chain[i]!=chainname ) {
      inres=false;
      res_end=residue[i-1];
    }
  }
  if(inres) {
    res_end=residue[size()-1];
  }
}

std::vector<int> PDB::getResiduesInChain( const std::string& chainname ) const {
  std::vector<int> res;
  std::set<int> seen;
  for(unsigned i=0; i<size(); ++i) {
    if( chain[i]==chainname && seen.insert(residue[i]).second ) {
      res.push_back(residue[i]);
    }
  }
  return res;
}

void PDB::getAtomRange( const std::string& chainname, AtomNumber& a_start, AtomNumber& a_end, std::string& errmsg ) const {
  bool inres=false, foundchain=false;
  for(unsigned i=0; i<size(); ++i) {
    if( chain[i]==chainname ) {
      if(!inres) {
        if(foundchain) {
          errmsg="found second start of chain named " + chainname;
        }
        a_start=numbers[i];
      }
      inres=true;
      foundchain=true;
    } else if( inres && chain[i]!=chainname ) {
      inres=false;
      a_end=numbers[i-1];
    }
  }
  if(inres) {
    a_end=numbers[size()-1];
  }
}

std::string PDB::getResidueName(int resnum,const std::string& chainid ) const {
  return getResidueName(resnum,chainid,"*");
}

std::string PDB::getResidueName(int resnum,const std::string& chainid,const std::string& icode ) const {
  const long i=findResidue(resnum,chainid,icode);
  if(i>=0) {
    return residuenames[i];
  }
  std::string num;
  Tools::convert( resnum, num );
  if(icode!="*") {
    num+=icode;
  }
  plumed_merror("residue " + num + " not found in chain " + chainid );
}


AtomNumber PDB::getNamedAtomFromResidueAndChain( const std::string& aname, int resnum, const std::string& chainid ) const {
  return getNamedAtomFromResidueAndChain(aname,resnum,chainid,"*");
}

AtomNumber PDB::getNamedAtomFromResidueAndChain( const std::string& aname, int resnum, const std::string& chainid, const std::string& icode ) const {
  for(unsigned i=0; i<size(); ++i) {
    if( residue[i]==resnum && atomsymb[i]==aname && ( chainid=="*" || chain[i]==chainid) && ( icode=="*" || insertion[i]==icode ) ) {
      return numbers[i];
    }
  }
  std::string num;
  Tools::convert( resnum, num );
  if(icode!="*") {
    num+=icode;
  }
  plumed_merror("residue " + num + " from chain " + chainid + " does not contain an atom named " + aname );
}

std::vector<AtomNumber> PDB::getAtomsInResidue(int resnum,const std::string& chainid)const {
  return getAtomsInResidue(resnum,chainid,"*");
}

std::vector<AtomNumber> PDB::getAtomsInResidue(int resnum,const std::string& chainid,const std::string& icode)const {
  std::vector<AtomNumber> tmp;
  for(unsigned i=0; i<size(); ++i) {
    if( residue[i]==resnum && ( chainid=="*" || chain[i]==chainid) && ( icode=="*" || insertion[i]==icode ) ) {
      tmp.push_back(numbers[i]);
    }
  }
  if(tmp.size()==0) {
    std::string num;
    Tools::convert( resnum, num );
    if(icode!="*") {
      num+=icode;
    }
    plumed_merror("Cannot find residue " + num + " from chain " + chainid  );
  }
  return tmp;
}

std::string PDB::getChainID(AtomNumber a) const {
  const auto p=number2index.find(a);
  if(p==number2index.end()) {
    std::string num;
    Tools::convert( a.serial(), num );
    plumed_merror("Chain for atom " + num + " not found" );
  }
  return chain[p->second];
}

bool PDB::checkForAtom( AtomNumber a ) const {
  const auto p=number2index.find(a);
  return (p!=number2index.end());
}

Vector PDB::getPosition(AtomNumber a)const {
  const auto p=number2index.find(a);
  if(p==number2index.end()) {
    plumed_merror("atom not available");
  } else {
    return positions[p->second];
  }
}

std::vector<std::string> PDB::getArgumentNames()const {
  return argnames;
}

bool PDB::allowedResidue( const std::string& type, const std::string& residuename ) const {
  if( type=="protein" ) {
    if(residuename=="ALA") {
      return true;
    } else if(residuename=="ARG") {
      return true;
    } else if(residuename=="ASN") {
      return true;
    } else if(residuename=="ASP") {
      return true;
    } else if(residuename=="CYS") {
      return true;
    } else if(residuename=="GLN") {
      return true;
    } else if(residuename=="GLU") {
      return true;
    } else if(residuename=="GLY") {
      return true;
    } else if(residuename=="HIS") {
      return true;
    } else if(residuename=="ILE") {
      return true;
    } else if(residuename=="LEU") {
      return true;
    } else if(residuename=="LYS") {
      return true;
    } else if(residuename=="MET") {
      return true;
    } else if(residuename=="PHE") {
      return true;
    } else if(residuename=="PRO") {
      return true;
    } else if(residuename=="SER") {
      return true;
    } else if(residuename=="THR") {
      return true;
    } else if(residuename=="TRP") {
      return true;
    } else if(residuename=="TYR") {
      return true;
    } else if(residuename=="VAL") {
      return true;
    }
// Terminal groups
    else if(residuename=="ACE") {
      return true;
    } else if(residuename=="NME") {
      return true;
    } else if(residuename=="NH2") {
      return true;
    }
// Alternative residue names in common force fields
    else if(residuename=="GLH") {
      return true;  // neutral GLU
    } else if(residuename=="ASH") {
      return true;  // neutral ASP
    } else if(residuename=="HID") {
      return true;  // HIS-D amber
    } else if(residuename=="HSD") {
      return true;  // HIS-D charmm
    } else if(residuename=="HIE") {
      return true;  // HIS-E amber
    } else if(residuename=="HSE") {
      return true;  // HIS-E charmm
    } else if(residuename=="HIP") {
      return true;  // HIS-P amber
    } else if(residuename=="HSP") {
      return true;  // HIS-P charmm
    } else if(residuename=="CYX") {
      return true;  // disulfide bridge CYS
    }
// Weird amino acids
    else if(residuename=="NLE") {
      return true;
    } else if(residuename=="SFO") {
      return true;
    } else {
      return false;
    }
  } else if( type=="dna" ) {
    if(residuename=="A") {
      return true;
    } else if(residuename=="A5") {
      return true;
    } else if(residuename=="A3") {
      return true;
    } else if(residuename=="AN") {
      return true;
    } else if(residuename=="G") {
      return true;
    } else if(residuename=="G5") {
      return true;
    } else if(residuename=="G3") {
      return true;
    } else if(residuename=="GN") {
      return true;
    } else if(residuename=="T") {
      return true;
    } else if(residuename=="T5") {
      return true;
    } else if(residuename=="T3") {
      return true;
    } else if(residuename=="TN") {
      return true;
    } else if(residuename=="C") {
      return true;
    } else if(residuename=="C5") {
      return true;
    } else if(residuename=="C3") {
      return true;
    } else if(residuename=="CN") {
      return true;
    } else if(residuename=="DA") {
      return true;
    } else if(residuename=="DA5") {
      return true;
    } else if(residuename=="DA3") {
      return true;
    } else if(residuename=="DAN") {
      return true;
    } else if(residuename=="DG") {
      return true;
    } else if(residuename=="DG5") {
      return true;
    } else if(residuename=="DG3") {
      return true;
    } else if(residuename=="DGN") {
      return true;
    } else if(residuename=="DT") {
      return true;
    } else if(residuename=="DT5") {
      return true;
    } else if(residuename=="DT3") {
      return true;
    } else if(residuename=="DTN") {
      return true;
    } else if(residuename=="DC") {
      return true;
    } else if(residuename=="DC5") {
      return true;
    } else if(residuename=="DC3") {
      return true;
    } else if(residuename=="DCN") {
      return true;
    } else {
      return false;
    }
  } else if( type=="rna" ) {
    if(residuename=="A") {
      return true;
    } else if(residuename=="A5") {
      return true;
    } else if(residuename=="A3") {
      return true;
    } else if(residuename=="AN") {
      return true;
    } else if(residuename=="G") {
      return true;
    } else if(residuename=="G5") {
      return true;
    } else if(residuename=="G3") {
      return true;
    } else if(residuename=="GN") {
      return true;
    } else if(residuename=="U") {
      return true;
    } else if(residuename=="U5") {
      return true;
    } else if(residuename=="U3") {
      return true;
    } else if(residuename=="UN") {
      return true;
    } else if(residuename=="C") {
      return true;
    } else if(residuename=="C5") {
      return true;
    } else if(residuename=="C3") {
      return true;
    } else if(residuename=="CN") {
      return true;
    } else if(residuename=="RA") {
      return true;
    } else if(residuename=="RA5") {
      return true;
    } else if(residuename=="RA3") {
      return true;
    } else if(residuename=="RAN") {
      return true;
    } else if(residuename=="RG") {
      return true;
    } else if(residuename=="RG5") {
      return true;
    } else if(residuename=="RG3") {
      return true;
    } else if(residuename=="RGN") {
      return true;
    } else if(residuename=="RU") {
      return true;
    } else if(residuename=="RU5") {
      return true;
    } else if(residuename=="RU3") {
      return true;
    } else if(residuename=="RUN") {
      return true;
    } else if(residuename=="RC") {
      return true;
    } else if(residuename=="RC5") {
      return true;
    } else if(residuename=="RC3") {
      return true;
    } else if(residuename=="RCN") {
      return true;
    } else {
      return false;
    }
  } else if( type=="water" ) {
    if(residuename=="SOL") {
      return true;
    }
    if(residuename=="WAT") {
      return true;
    }
    return false;
  } else if( type=="ion" ) {
    if(residuename=="IB+") {
      return true;
    }
    if(residuename=="CA") {
      return true;
    }
    if(residuename=="CL") {
      return true;
    }
    if(residuename=="NA") {
      return true;
    }
    if(residuename=="MG") {
      return true;
    }
    if(residuename=="K") {
      return true;
    }
    if(residuename=="RB") {
      return true;
    }
    if(residuename=="CS") {
      return true;
    }
    if(residuename=="LI") {
      return true;
    }
    if(residuename=="ZN") {
      return true;
    }
    return false;
  }
  return false;
}

}

