/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2018,2019 The plumed team
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
#include "CLTool.h"
#include "CLToolRegister.h"
#include "tools/Tools.h"
#include "core/ActionRegister.h"
#include "tools/IFile.h"
#include "tools/OFile.h"
#include "tools/h36.h"
#include <cstdio>
#include <string>
#include <vector>
#include <array>
#include <limits>

using namespace std;

namespace PLMD {
namespace cltools {

//+PLUMEDOC TOOLS pdbrenumber
/*
Modify atom numbers in a PDB, possibly using hybrid-36 coding.

When reading a PDB files, PLUMED honors the serial number of each atom.
This command can be used to process a PDB file changing the atom serial numbers.
Notice that the resulting list might have gaps. It is however fundamental
that atom numbers correspond to those used within the MD code.
Importantly, if the serial number of an atom is greater than 99999, it is
written in hybrid-36 notation (see \ref pdbreader).
The main use of \ref pdbrenumber is thus that of producing files where atoms
are numbered using hybrid-36 convention.

The output PDB file is identical to the input PDB file, except for the atom number
field.
The rest of the line is written unchanged
to the output file, even if it is incorrectly formatted. Residue numbers are not touched,
and atom numbers in the input file are ignored.


\par Examples

By default, \ref pdbreader  just sets the numbers as progressive starting from 1.
For instance the following command:
\verbatim
> plumed pdbrenumber --ipdb input.pdb --opdb output.pdb
\endverbatim
will copy file `input.pdb` to `output.pdb` replacing all the serial atoms with
increasing numbers starting from one. Atoms that have an index that is greater than 99999 will be written
in the output PDB file in hybrid-36 code.

It is possible to set a different serial number for the first atom, letting the
following ones grow by one at each line. Here for instance the first atom
will be assigned serial 1000, the second serial 1001, etc:
\verbatim
> plumed pdbrenumber --ipdb input.pdb --opdb output.pdb --firstatomnumber 1000
\endverbatim
If the first atom number is >99999, it should be given as a decimal number (not in hybrid-36 code).
However, numbers >99999 in the output PDB file will be written in hybrid-36 code.

As an alternative, one can provide a list of atoms as one per line in an auxiliary file.
\verbatim
> plumed pdbrenumber --ipdb input.pdb --opdb output.pdb --atomnumbers list.txt
\endverbatim
The `list.txt` file might be something like this
\verbatim
120000
120001
120002
1
2
3
\endverbatim
Numbers >99999 in the list should be provided as decimal numbers (not in hybrid-36 code).
However, numbers >99999 in the output PDB file will be written in hybrid-36 code.
Notice that there should be at least enough lines in `list.txt` as many atoms in the PDB file.
Additional lines in `list.txt` will just be ignored.


*/
//+ENDPLUMEDOC

class PdbRenumber:
  public CLTool
{
public:
  static void registerKeywords( Keywords& keys );
  explicit PdbRenumber(const CLToolOptions& co );
  int main(FILE* in, FILE*out,Communicator& pc);
  string description()const {
    return "Modify atom numbers in a PDB, possibly using hybrid-36 coding";
  }
};

PLUMED_REGISTER_CLTOOL(PdbRenumber,"pdbrenumber")

void PdbRenumber::registerKeywords( Keywords& keys ) {
  CLTool::registerKeywords( keys );
  keys.add("compulsory","--ipdb","specify the name of the input PDB file");
  keys.add("compulsory","--opdb","specify the name of the output PDB file");
  keys.add("optional","--firstatomnumber","specify the desired serial number of the first atom of the output file");
  keys.add("optional","--atomnumbers","specify the desired serial numbers of the atoms of the output file using a separate list");
}

PdbRenumber::PdbRenumber(const CLToolOptions& co ):
  CLTool(co)
{
  inputdata=commandline;
}

int PdbRenumber::main(FILE* in, FILE*out,Communicator& pc) {

  std::string ipdb;
  parse("--ipdb",ipdb);
  std::string opdb;
  parse("--opdb",opdb);

  unsigned iat=0;

  parse("--firstatomnumber",iat);

  std::string atomnumbers;
  parse("--atomnumbers",atomnumbers);

  plumed_assert(ipdb.length()>0) << "please specify the input PDB with --ipdb";
  plumed_assert(opdb.length()>0) << "please specify the onput PDB with --opdb";
  fprintf(out,"  with input PDB: %s\n",ipdb.c_str());
  fprintf(out,"  with output PDB: %s\n",opdb.c_str());

  std::vector<unsigned> serials;

  if(atomnumbers.length()>0) {
    plumed_assert(iat==0) << "it is not possible to use both --atomnumbers and --firstatomnumber";
    fprintf(out,"  reading atom numbers from file %s\n",atomnumbers.c_str());
    IFile ifile;
    ifile.open(atomnumbers);
    std::string line;
    while(ifile.getline(line)) {
      int i;
      Tools::convert(line,i);
      serials.push_back(i);
    }
  } else {
    if(iat==0) iat=1;
    fprintf(out,"  with atoms starting from %u\n",iat);
  }

  IFile ifile;
  ifile.open(ipdb);

  OFile ofile;
  ofile.open(opdb);

  std::string line;
  while(ifile.getline(line)) {
    std::string record=line.substr(0,6);
    Tools::trim(record);

    if(record=="ATOM" || record=="HETATM") {
      std::array<char,6> at;
      unsigned ii=iat;
      if(serials.size()>0) {
        plumed_assert(iat<serials.size()) << "there are more atoms in the PDB than serials in the file";
        ii=serials[iat];
      }
      const char* msg = h36::hy36encode(5,ii,&at[0]);
      plumed_assert(msg==nullptr) << msg;
      at[5]=0;
      ofile << line.substr(0,6) << &at[0] << line.substr(11) << "\n";
      iat++;
    } else {
      if(record=="END" || record=="ENDMDL") iat=0;
      ofile << line << "\n";
    }
  }

  return 0;
}
}

} // End of namespace
