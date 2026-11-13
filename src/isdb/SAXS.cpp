/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2023 The plumed team
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
/*
 This class was originally written by Alexander Jussupow
 Arrayfire implementation (now removed) by Alexander Jussupow and CC
 Extension for the middleman algorithm (now removed) by Max Muehlbauer
 Refactoring for hySAXS Martini form factors for Nucleic Acids by Cristina Paissoni
 Refactoring for hySAS OneBead form factors with solvent correction by Federico Ballabio and Riccardo Capelli
 Resolution function by Henrique Musseli Cezar
*/

#include "SAXS.h"
#include "core/ActionRegister.h"
#include "core/ActionSet.h"
#include "core/GenericMolInfo.h"
#include "tools/MolDataClass.h"
#include "tools/Communicator.h"
#include "tools/Pbc.h"
#include "tools/Tools.h"
#include "tools/IFile.h"
#include "tools/OpenMP.h"

#include <map>
#include <iterator>
#include <iostream>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <memory>
#include <cstring>
#include <cstdint>

namespace PLMD {
namespace isdb {

//+PLUMEDOC ISDB_COLVAR SAXS
/*
Calculates SAXS intensity.

SAXS intensities are calculated for a set of scattering vectors using QVALUE keywords numbered from 1.
Form factors can be assigned either by polynomial expansion of any order by using the PARAMETERS keywords, or
automatically matched to atoms using the ATOMISTIC flag by reading a PDB file. Alternatively to the atomistic
representation, two types of coarse-grained mapping are available:
- MARTINI.
- ONEBEAD.

Whether for PARAMETERS, ATOMISTIC, and ONEBEAD the user must provide an all-atom PDB file via MOLINFO before the
SAXS instruction. MARTINI requires a mapping scheme consisting of a PDB file that contains both the all-atom
and MARTINI representations, and a bead position file (e.g., bead1: CENTER ATOMS=1,5,7,11,12 WEIGHTS=14,12,12,
12,16). The MARTINI_FF keyword selects the form-factor generation: MARTINI_FF=2 (default) uses the Martini 2.2
form factors for proteins and nucleic acids, while MARTINI_FF=3 uses Martini 3 form factors for proteins. The
Martini 3 protein form factors are described in the corresponding paper in the bibliography below.

ONEBEAD scheme consists in a single-bead per amino acid residue, three-bead for nucleic acid residue (one for
the phosphate group, one for the pentose sugar, one for the nucleobase), and a single bead per monosaccharide
for glycans. Every bead carries the sign of its net scattering amplitude, so a bead whose contrast is negative
with respect to the solvent enters the Debye sum with the correct sign. This matters for SANS close to a
contrast-match point and for nuclei with a negative coherent scattering length such as Mn; in SAXS every bead
has positive contrast and the sign is always positive. PLUMED creates a virtual bead on which
the SAXS calculations are performed, centred on the COM of all atoms belonging to the bead. It is possible to
account for the contribution of the solvation layer to the SAXS intensity by adding a correction term for the
solvent accessible beads only: the form factors of the amino acids / phosphate groups / pentose sugars /
nucleobases / monosaccharides with a SASA (computed via LCPO algorithm) greater than a threshold are corrected
according to an electron density term. Both the surface cut-off threshold and the electron density term can be
set by the user with the SASA_CUTOFF and SOLVATION_CORRECTION keywords. Moreover, SASA stride calculation can be
modified using SOLVATION_STRIDE, which is set to 10 steps by default.
ONEBEAD requires an additional PDB file to perform mapping conversion, which must be provided via TEMPLATE
keyword. This PDB file must contain exactly the atoms selected by ATOMS, listed in the same order, with the
atoms of each chain contiguous and the residues of each chain in ascending order. Any mismatch is reported as
an error. Note that the ONEBEAD mapping is taken entirely from the TEMPLATE: unlike ATOMISTIC and MARTINI, it
does not read residue or atom names from MOLINFO.
The AMBER OL3 (RNA) and OL15 (DNA) naming is required for nucleic acids.
Three additional bead types are available for DNA and RNA besides phosphate group, pentose sugar, and nucleobase:
- 5'-end pentose sugar capped with an hydroxyl moiety at C5' (the residue name in the PDB must be followed by
"5", e.g., DC5 or C5 for cytosine in DNA and RNA, respectively);
- 5'-phosphorylated end with an additional hydroxyl moiety at P (the residue name in the PDB must be followed by
"T", e.g., DCT or CT for cytosine in DNA and RNA, respectively);
- 3'-end pentose sugar capped with an hydroxyl moiety at C3' (the residue name in the PDB must be followed by
"3", e.g., DC3 or C3 for cytosine in DNA and RNA, respectively).
An additional bead type is available for proteins:
- Cysteine residue involved in disulfide bridge (the residue in the PDB must be named CYX).

Experimental reference intensities can be added using the EXPINT keywords. All these values must be normalised
to the SAXS intensity at q = 0. To facilitate this operation, the SCALE_EXPINT keyword can be used to provide
the intensity at q = 0. Each EXPINT is divided by SCALE_EXPINT.
The maximum QVALUE for ONEBEAD is set to 0.3 inverse angstroms.
The solvent density, that by default is set to 0.334 electrons per cubic angstrom (bulk water), can be modified
using the SOLVDENS keyword.

The ABSOLUTE flag can be used in order to calculate intensities in the absolute scale. It is only available for
the ATOMISTIC scheme and cannot be used with SCALE_EXPINT.

[METAINFERENCE](METAINFERENCE.md) can be activated using DOSCORE and the other relevant keywords.

##Examples

in the following example the SAXS intensities are calculated using single-bead per residue approximation, with a
SASA threshold of 1 square nanometer and a solvation term of 0.04. Each experimental intensity is divided by
1.4002, which is the corresponding theoretical intensity value at q = 0. The form factors are selected according
to the PDB file specified by TEMPLATE keyword.

```plumed
MOLINFO STRUCTURE=template_AA.pdb

SAXS ...
LABEL=SAXS
ATOMS=1-355
ONEBEAD
TEMPLATE=template_AA.pdb
SOLVDENS=0.334
SOLVATION_CORRECTION=0.04
SOLVATION_STRIDE=1
SASA_CUTOFF=1.0
SCALE_EXPINT=1.4002
QVALUE1=0.03 EXPINT1=1.0902
QVALUE2=0.06 EXPINT2=0.790632
QVALUE3=0.09 EXPINT3=0.453808
QVALUE4=0.12 EXPINT4=0.254737
QVALUE5=0.15 EXPINT5=0.154928
QVALUE6=0.18 EXPINT6=0.0921503
QVALUE7=0.21 EXPINT7=0.052633
QVALUE8=0.24 EXPINT8=0.0276557
QVALUE9=0.27 EXPINT9=0.0122775
QVALUE10=0.30 EXPINT10=0.00880634
... SAXS

PRINT ARG=(SAXS\.q-.*),(SAXS\.exp-.*) FILE=saxsdata STRIDE=1
```

The Debye sum can be computed on a GPU with the USEGPU flag. The GPU implementations are provided by
plugins, that must be loaded with [LOAD](LOAD.md) before SAXS/SANS (see plugins/openaccPTM in the
PLUMED sources). The GPU is the one assigned to PLUMED by the MD code, otherwise the one the MD code
is using, or the first one (see the documentation of the plugin). With MPI each rank computes a share
of the sum on its own GPU; with SERIAL only the first rank uses its GPU and the results are shared with the others.

````
LOAD FILE=/path/to/plumedOpenACC.so
MOLINFO STRUCTURE=template.pdb
SAXS ...
LABEL=SAXS
ATOMS=1-107
ATOMISTIC
USEGPU
QVALUE1=0.03
QVALUE2=0.06
QVALUE3=0.09
... SAXS
````

(The usual PLUMED syntax highlighting is not used above only because the plugin is not available when the manual is built.)

*/
//+ENDPLUMEDOC

//+PLUMEDOC ISDB_COLVAR SANS
/*
Calculates SANS intensity.

SANS intensities are calculated for a set of scattering vectors using QVALUE keywords numbered from 1.
Form factors are automatically assigned to atoms using the ATOMISTIC flag by reading a PDB file, by reading
the scattering lengths with the PARAMETERS keyword from input or with the PARAMETERSFILE keyword or, alternatively,
a ONEBEAD coarse-grained implementation is available.

Both for ATOMISTIC and ONEBEAD the user must provide an all-atom PDB file via MOLINFO before the SANS instruction.

ONEBEAD scheme consists in a single-bead per amino acid residue, three-bead for nucleic acid residue (one for
the phosphate group, one for the pentose sugar, one for the nucleobase), and a single bead per monosaccharide
for glycans. Every bead carries the sign of its net scattering amplitude, so a bead whose contrast is negative
with respect to the solvent enters the Debye sum with the correct sign. This matters for SANS close to a
contrast-match point and for nuclei with a negative coherent scattering length such as Mn; in SAXS every bead
has positive contrast and the sign is always positive. PLUMED creates a virtual bead on which
the SANS calculations are performed, centred on the COM of all atoms belonging to the bead. It is possible to
account for the contribution of the solvation layer to the SANS intensity by adding a correction term for the
solvent accessible beads only: the form factors of the amino acids / phosphate groups / pentose sugars /
nucleobases / monosaccharides with a SASA (computed via LCPO algorithm) greater than a threshold are corrected
according to an electron density term. Both the surface cut-off threshold and the electron density term can be set by the user
with the SASA_CUTOFF and SOLVATION_CORRECTION keywords. Moreover, SASA stride calculation can be modified using
SOLVATION_STRIDE, which is set to 10 steps by default. The deuteration of the solvent-exposed residues is chosen
with a probability equal to the deuterium concentration in the buffer. The deuterated residues are updated with a
stride equal to SOLVATION_STRIDE. The fraction of deuterated water can be set with DEUTER_CONC, the default value
is 0.
ONEBEAD requires an additional PDB file to perform mapping conversion, which must be provided via TEMPLATE
keyword. This PDB file must contain exactly the atoms selected by ATOMS, listed in the same order, with the
atoms of each chain contiguous and the residues of each chain in ascending order. Any mismatch is reported as
an error. Note that the ONEBEAD mapping is taken entirely from the TEMPLATE: unlike ATOMISTIC and MARTINI, it
does not read residue or atom names from MOLINFO.
The AMBER OL3 (RNA) and OL15 (DNA) naming is required for nucleic acids.
Three additional bead types are available for DNA and RNA besides phosphate group, pentose sugar, and nucleobase:
- 5'-end pentose sugar capped with an hydroxyl moiety at C5' (the residue name in the PDB must be followed by "5",
e.g., DC5 or C5 for cytosine in DNA and RNA, respectively);
- 5'-phosphorylated end with an additional hydroxyl moiety at P (the residue name in the PDB must be followed by
"T", e.g., DCT or CT for cytosine in DNA and RNA, respectively);
- 3'-end pentose sugar capped with an hydroxyl moiety at C3' (the residue name in the PDB must be followed by "3",
e.g., DC3 or C3 for cytosine in DNA and RNA, respectively).
An additional bead type is available for proteins:
- Cysteine residue involved in disulfide bridge (the residue in the PDB must be named CYX).

PLEASE NOTE: at the moment, we DO NOT explicitly take into account deuterated residues in the ATOMISTIC
representation, but we correct the solvent contribution via the DEUTER_CONC keyword.

Experimental reference intensities can be added using the EXPINT keywords. All these values must be normalised
to the SANS intensity at q = 0. To facilitate this operation, the SCALE_EXPINT keyword can be used to provide
the intensity at q = 0. Each EXPINT is divided by SCALE_EXPINT.
The maximum QVALUE for ONEBEAD is set to 0.3 inverse angstroms.
The solvent density, that by default is set to 0.334 electrons per cubic angstrom (bulk water), can be modified
using the SOLVDENS keyword.

The ABSOLUTE flag can be used in order to calculate intensities in the absolute scale. It is only available for
the ATOMISTIC scheme and cannot be used with SCALE_EXPINT.

[METAINFERENCE](METAINFERENCE.md) can be activated using DOSCORE and the other relevant keywords.

##Examples

in the following example the SANS intensities are calculated at atomistic resolution. The form factors are assigned
according to the PDB file specified in the MOLINFO. Each experimental intensity is divided by 1.4002, which is the
corresponding theoretical intensity value at q = 0. The deuterated water fraction is set to 48%.

```plumed
MOLINFO STRUCTURE=template_AA.pdb

SANS ...
LABEL=SANS
ATOMS=1-355
ATOMISTIC
SCALE_EXPINT=1.4002
DEUTER_CONC=0.48
QVALUE1=0.03 EXPINT1=1.0902
QVALUE2=0.06 EXPINT2=0.790632
QVALUE3=0.09 EXPINT3=0.453808
QVALUE4=0.12 EXPINT4=0.254737
QVALUE5=0.15 EXPINT5=0.154928
QVALUE6=0.18 EXPINT6=0.0921503
QVALUE7=0.21 EXPINT7=0.052633
QVALUE8=0.24 EXPINT8=0.0276557
QVALUE9=0.27 EXPINT9=0.0122775
QVALUE10=0.30 EXPINT10=0.00880634
... SANS

PRINT ARG=(SANS\.q-.*),(SANS\.exp-.*) FILE=sansdata STRIDE=1
```

As for [SAXS](SAXS.md), the Debye sum can be computed on a GPU with the USEGPU flag (and
GPU_BACKEND), using a backend provided by a plugin loaded with [LOAD](LOAD.md) (see plugins/openaccPTM in the
PLUMED sources). Everything else, including the deuteration and the solvation correction, is computed as usual.

*/
//+ENDPLUMEDOC


PLUMED_REGISTER_ACTION(SAXS,"SAXS")
PLUMED_REGISTER_ACTION(SAXS,"SANS")

void SAXS::registerKeywords(Keywords& keys) {
  MetainferenceBase::registerKeywords(keys);
  keys.addFlag("NOPBC",false,"Ignore the periodic boundary conditions when calculating distances");
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose. With USEGPU only the first rank uses its GPU");
  keys.addFlag("USEGPU",false,"Calculate the Debye sum on the GPU, using a backend provided by a plugin that must be loaded with LOAD");
  keys.add("optional","GPU_BACKEND","The GPU backend to be used with USEGPU, needed only if more than one is loaded");
  keys.addFlag("ABSOLUTE",false,"Absolute intensity: the intensities for each q-value are not normalised for the intensity at q=0.");
  keys.addFlag("ATOMISTIC",false,"Calculate SAXS for an atomistic model");
  keys.addFlag("MARTINI",false,"Calculate SAXS for a Martini model");
  keys.add("compulsory","MARTINI_FF","2","Martini form-factor generation to use with MARTINI: 2 (default, Martini 2.2 proteins and nucleic acids) or 3 (Martini 3 proteins)");
  keys.addFlag("ONEBEAD",false,"calculate SAXS for a single bead model");
  keys.add("compulsory","TEMPLATE","template.pdb","A PDB file is required for ONEBEAD mapping");
  keys.add("atoms","ATOMS","The atoms to be included in the calculation, e.g. the whole protein");
  keys.add("numbered","QVALUE","Selected scattering lengths in inverse angstroms are given as QVALUE1, QVALUE2, ...");
  keys.add("numbered","PARAMETERS","Used parameter Keywords like PARAMETERS1, PARAMETERS2. These are used to calculate the form factor for the $i$th atom/bead");
  keys.add("optional","PARAMETERSFILE","Read the PARAMETERS from a file");
  keys.add("compulsory","DEUTER_CONC","0.","Fraction of deuterated solvent");
  keys.add("compulsory","SOLVDENS","0.334","Density of the solvent to be used for the correction of atomistic form factors");
  keys.add("compulsory","SOLVATION_CORRECTION","0.0","Solvation layer electron density correction (ONEBEAD only)");
  keys.add("compulsory","SASA_CUTOFF","1.0","SASA value to consider a residue as exposed to the solvent (ONEBEAD only)");
  keys.add("numbered","EXPINT","Add an experimental value for each q value");
  keys.add("numbered","EXPERR","Add an experimental standard uncertainty for each q value. Values are divided by SCALE_EXPINT, exactly as EXPINT");
  keys.addFlag("CHI2",false,"Calculate the scale-optimized chi-squared between calculated and experimental intensities using EXPERR");
  keys.addFlag("CHI2ENSEMBLE",false,"Calculate chi-squared from arithmetic replica-averaged intensities for monitoring");
  keys.add("numbered","SIGMARES","Standard deviation of the resolution function for each q value, in inverse angstroms");
  keys.add("compulsory","N","10","Number of points in the resolution function integral");
  keys.add("compulsory","SOLVATION_STRIDE","10","Number of steps between every new residues solvation estimation via LCPO (ONEBEAD only)");
  keys.add("compulsory","SCALE_EXPINT","1.0","Scaling value for experimental data normalization");
  keys.addOutputComponent("q","default","scalar","The # SAXS of q");
  keys.addOutputComponent("exp","EXPINT","scalar","The # experimental intensity");
  keys.addOutputComponent("chi2","CHI2","scalar","Scale-optimized chi-squared against EXPINT using EXPERR");
  keys.addOutputComponent("chi2red","CHI2","scalar","Reduced scale-optimized chi-squared with N-1 degrees of freedom");
  keys.addOutputComponent("chi2scale","CHI2","scalar","Best-fit multiplicative scale applied to the calculated intensities");
  keys.addOutputComponent("chi2ensemble","CHI2ENSEMBLE","scalar","Chi-squared of replica-averaged intensities");
  keys.addOutputComponent("chi2ensemblered","CHI2ENSEMBLE","scalar","Reduced chi-squared of replica-averaged intensities");
  keys.addOutputComponent("chi2ensemblescale","CHI2ENSEMBLE","scalar","Best-fit scale for replica-averaged intensities");
}

SAXS::SAXS(const ActionOptions&ao):
  PLUMED_METAINF_INIT(ao),
  saxs(true),
  absolute(false),
  pbc(true),
  serial(false),
  first_row_(0),
  row_stride_(1),
  has_rows_(true),
  reduce_(false),
  usegpu_(false),
  onebead(false),
  chi2(false),
  chi2ensemble(false),
  chi2ensemble_replicas(1),
  isFirstStep(true) {
  if( getName().find("SAXS")!=std::string::npos) {
    saxs=true;
  } else if( getName().find("SANS")!=std::string::npos) {
    saxs=false;
  }

  std::vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  unsigned size = atoms.size();
  if(size == 0) {
    error("SAXS/SANS requires at least one atom");
  }
  {
    // an atom listed twice would be two particles at the same position (distance zero):
    // the Debye sum would give nan
    std::vector<AtomNumber> sorted(atoms);
    std::sort(sorted.begin(),sorted.end());
    const auto dup=std::adjacent_find(sorted.begin(),sorted.end());
    if(dup!=sorted.end()) {
      error("ATOMS contains duplicate atoms (atom " + std::to_string(dup->serial()) + " is listed more than once)");
    }
  }

  parseFlag("SERIAL",serial);

  parseFlag("USEGPU",usegpu_);
  const bool usegpu=usegpu_;
  std::string gpu_backend;
  parse("GPU_BACKEND",gpu_backend);
  // rows of the Debye sum of this rank: on the CPU, with SERIAL every rank computes
  // everything; on the GPU, with SERIAL only the first rank does
  if(!serial) {
    first_row_ = comm.Get_rank();
    row_stride_ = comm.Get_size();
  } else if(usegpu && comm.Get_rank()>0) {
    has_rows_ = false;
  }
  reduce_ = (!serial || usegpu) && comm.Get_size()>1;
  if(usegpu) {
    // a rank that fails to create the backend must not leave the others waiting
    std::string failure;
    if(has_rows_) {
      try {
        backend_ = saxsDebyeBackendRegister().create(gpu_backend,*this);
      } catch(const std::exception & e) {
        failure = e.what();
        if(failure.empty()) {
          failure = "unknown error";
        }
      }
    }
    int failed = failure.empty() ? 0 : 1;
    comm.Sum(failed);
    if(failed>0) {
      error("cannot use the GPU" + (failure.empty() ? std::string(" (failure on another rank)") : ": " + failure));
    }
    log.printf("  using GPU to calculate the Debye sum\n");
  } else if(!gpu_backend.empty()) {
    error("GPU_BACKEND can only be used together with USEGPU");
  }

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;
  if(pbc) {
    log.printf("  using periodic boundary conditions\n");
  } else {
    log.printf("  without periodic boundary conditions\n");
  }



  bool atomistic=false;
  parseFlag("ATOMISTIC",atomistic);
  if(atomistic) {
    log.printf("  using ATOMISTIC form factors\n");
  }
  bool martini=false;
  parseFlag("MARTINI",martini);
  unsigned martini_ff=2;
  parse("MARTINI_FF",martini_ff);
  if(martini) {
    if(martini_ff!=2&&martini_ff!=3) {
      error("MARTINI_FF must be 2 (Martini 2.2) or 3 (Martini 3)");
    }
    log.printf("  using MARTINI %u form factors\n",martini_ff);
  } else if(martini_ff!=2) {
    error("MARTINI_FF can only be used together with MARTINI");
  }
  onebead=false;
  parseFlag("ONEBEAD",onebead);
  if(onebead) {
    log.printf("  using ONEBEAD form factors\n");
  }
  bool fromfile=false;
  std::string parametersfile;
  parse("PARAMETERSFILE",parametersfile);
  if (parametersfile.length() != 0) {
    fromfile=true;
  }
  if(fromfile) {
    log.printf("  will read form factors from file\n");
  }
  parseFlag("ABSOLUTE",absolute);

  if(martini&&atomistic) {
    error("You cannot use MARTINI and ATOMISTIC at the same time");
  }
  if(martini&&onebead) {
    error("You cannot use MARTINI and ONEBEAD at the same time");
  }
  if(onebead&&atomistic) {
    error("You cannot use ONEBEAD and ATOMISTIC at the same time");
  }
  if((martini)&&(!saxs)) {
    error("MARTINI cannot be used with SANS");
  }
  if((fromfile)&&((atomistic)||(martini)||(onebead))) {
    error("You cannot read parameters from file and use ATOMISTIC/MARTINI/ONEBEAD");
  }

  unsigned ntarget=0;
  for(unsigned i=0;; ++i) {
    double t_list;
    if( !parseNumbered( "QVALUE", i+1, t_list) ) {
      break;
    }
    if(t_list<=0.) {
      error("QVALUE cannot be less or equal to zero!\n");
    }
    if(onebead&&t_list>0.3) {
      error("ONEBEAD mapping QVALUE must be smaller or equal to 0.3");
    }
    q_list.push_back(t_list);
    ntarget++;
  }
  const unsigned numq = ntarget;
  if(numq==0) {
    error("At least one QVALUE must be specified");
  }

  for(unsigned i=0; i<numq; ++i) {
    if(q_list[i]==0.) {
      error("it is not possible to set q=0\n");
    }
    if(i>0&&q_list[i]<q_list[i-1]) {
      error("QVALUE must be in ascending order");
    }
    log.printf("  my q: %lf \n",q_list[i]);
  }

  rho = 0.334;
  parse("SOLVDENS", rho);
  log.printf("  Solvent density: %lf\n", rho);

  double scale_expint=1.;
  parse("SCALE_EXPINT",scale_expint);
  if(scale_expint<=0.) {
    error("SCALE_EXPINT must be greater than zero");
  }
  parseFlag("CHI2",chi2);
  parseFlag("CHI2ENSEMBLE",chi2ensemble);
  if(chi2ensemble && !chi2) {
    error("CHI2ENSEMBLE requires CHI2");
  }
  if(chi2ensemble) {
    chi2ensemble_replicas = comm.Get_rank()==0 ? multi_sim_comm.Get_size() : 0;
    comm.Sum(&chi2ensemble_replicas,1);
  }

  if((!atomistic&&absolute)||(absolute&&scale_expint!=1)) {
    error("ABSOLUTE can be used only combined with ATOMISTIC without SCALE_EXPINT");
  }
  if(atomistic) {
    log.printf("  Scale for intensities: %s\n", absolute ? "absolute" : "normalised");
  }

  double correction = 0.00;
  parse("SOLVATION_CORRECTION", correction);
  rho_corr=rho-correction;
  if(onebead) {
    log.printf("  Solvation density contribution: %lf\n", correction);
  }
  if((atomistic||martini||fromfile)&&(rho_corr!=rho)) {
    log.printf("  Solvation density contribution is taken into account in ONEBEAD only\n");
  }

  solv_stride = 10;
  parse("SOLVATION_STRIDE", solv_stride);
  if(solv_stride < 1.) {
    error("SOLVATION_STRIDE must be greater than 0");
  }
  if(onebead&&(rho_corr!=rho)) {
    log.printf("  SASA calculation stride: %u\n", solv_stride);
  }

  sasa_cutoff = 1.0;
  parse("SASA_CUTOFF", sasa_cutoff);
  if(sasa_cutoff <= 0.) {
    error("SASA_CUTOFF must be greater than 0");
  }

  deuter_conc = 0.;
  parse("DEUTER_CONC", deuter_conc);
  if ((deuter_conc)&&(fromfile)) {
    error("DEUTER_CONC cannot be used with PARAMETERSFILE");
  }
  if(deuter_conc < 0. || deuter_conc > 1.) {
    error("DEUTER_CONC must be in 0-1 range");
  }
  if ((atomistic||onebead)&&(!saxs)) {
    log.printf("  Solvent deuterium fraction: %lf/1.000000\n", deuter_conc);
  }

  PDB pdb;
  if(onebead) {
    std::string template_name;
    parse("TEMPLATE",template_name);
    log.printf("  Template for ONEBEAD mapping conversion: %s\n", template_name.c_str());
    if( !pdb.read(template_name,usingNaturalUnits(),1.) ) {
      plumed_merror("missing input file " + template_name);
    }
  }

  // preliminary mapping for onebead representation
  if(onebead) {
    LCPOparam.resize(size);
    nres = getOnebeadMapping(pdb, atoms);
    bead_sign.assign(nres, 1.);
    if(saxs) {
      Iq0_vac.resize(nres);
      Iq0_solv.resize(nres);
      Iq0_mix.resize(nres);
    } else { // SANS
      Iq0_vac_H.resize(nres);
      Iq0_solv_H.resize(nres);
      Iq0_mix_H.resize(nres);
      Iq0_vac_D.resize(nres);
      Iq0_mix_D.resize(nres);
    }
    atoi.resize(nres);
  } else {
    atoi.resize(size);
  }

  Iq0=0;
  std::vector<std::vector<long double> > FF_tmp;
  std::vector<std::vector<long double> > FF_tmp_vac;
  std::vector<std::vector<long double> > FF_tmp_mix;
  std::vector<std::vector<long double> > FF_tmp_solv;
  // SANS
  std::vector<std::vector<long double> > FF_tmp_vac_H;
  std::vector<std::vector<long double> > FF_tmp_mix_H;
  std::vector<std::vector<long double> > FF_tmp_solv_H;
  std::vector<std::vector<long double> > FF_tmp_vac_D;
  std::vector<std::vector<long double> > FF_tmp_mix_D;
  std::vector<std::vector<long double> > parameter_H;
  std::vector<std::vector<long double> > parameter_D;

  if(!atomistic&&!martini&&!onebead&&!fromfile) { // read PARAMETERS from PLUMED file
    if (saxs) {
      // read in parameter std::vector
      std::vector<std::vector<long double> > parameter;
      parameter.resize(size);
      ntarget=0;
      for(unsigned i=0; i<size; ++i) {
        if( !parseNumberedVector( "PARAMETERS", i+1, parameter[i]) ) {
          break;
        }
        ntarget++;
      }
      if( ntarget!=size ) {
        error("found wrong number of parameter std::vectors");
      }
      FF_tmp.resize(numq,std::vector<long double>(size));
      for(unsigned i=0; i<size; ++i) {
        atoi[i]=i;
        for(unsigned k=0; k<numq; ++k) {
          for(unsigned j=0; j<parameter[i].size(); ++j) {
            FF_tmp[k][i]+= parameter[i][j]*std::pow(static_cast<long double>(q_list[k]),j);
          }
        }
      }
      for(unsigned i=0; i<size; ++i) {
        Iq0+=parameter[i][0];
      }
      Iq0 *= Iq0;
    } else { // SANS
      std::vector<long double> parameter;
      parameter.resize(size);
      ntarget=0;
      for(unsigned i=0; i<size; ++i) {
        if( !parseNumbered( "PARAMETERS", i+1, parameter[i]) ) {
          break;
        }
        ntarget++;
      }
      if( ntarget!=size ) {
        error("found wrong number of parameter std::vectors");
      }
      FF_tmp.resize(numq,std::vector<long double>(size));
      for(unsigned i=0; i<size; ++i) {
        atoi[i]=i;
        for(unsigned k=0; k<numq; ++k) {
          FF_tmp[k][i]+= parameter[i];
        }
      }
      for(unsigned i=0; i<size; ++i) {
        Iq0+=parameter[i];
      }
      Iq0 *= Iq0;
    }
  } else if (fromfile) { // read PARAMETERS from user-provided file
    log.printf("  Reading PARAMETERS from file: %s\n", parametersfile.c_str());
    if (saxs) {
      FF_tmp.resize(numq,std::vector<long double>(size));
      std::vector<std::vector<long double> > parameter;
      parameter.resize(size);

      IFile ifile;
      ifile.open(parametersfile);
      std::string line;

      ntarget=0;
      while(ifile.getline(line)) {
        Tools::ltrim(line);
        Tools::trimComments(line);
        if (line.empty()) {
          continue;
        }
        if (ntarget >= size) {
          error("PARAMETERSFILE has more PARAMETERS than there are scattering centers");
        }
        std::string num;
        Tools::convert(ntarget+1,num);
        std::vector<std::string> lineread{line};
        if (!Tools::parseVector(lineread, "PARAMETERS"+num, parameter[ntarget], -1)) {
          error("Missing PARAMETERS or PARAMETERS not sorted");
        }
        ntarget++;
      }
      if( ntarget!=size ) {
        error("found wrong number of PARAMETERS in file");
      }

      for(unsigned i=0; i<size; ++i) {
        atoi[i]=i;
        for(unsigned k=0; k<numq; ++k) {
          for(unsigned j=0; j<parameter[i].size(); ++j) {
            FF_tmp[k][i]+= parameter[i][j]*std::pow(static_cast<long double>(q_list[k]),j);
          }
        }
      }
      for(unsigned i=0; i<size; ++i) {
        Iq0+=parameter[i][0];
      }
      Iq0 *= Iq0;
    } else { // SANS
      FF_tmp.resize(numq,std::vector<long double>(size));

      IFile ifile;
      ifile.open(parametersfile);
      std::string line;

      ntarget=0;
      while(ifile.getline(line)) {
        Tools::ltrim(line);
        Tools::trimComments(line);
        if (line.empty()) {
          continue;
        }
        if (ntarget >= size) {
          error("PARAMETERSFILE has more PARAMETERS than there are scattering centers");
        }
        std::string num;
        Tools::convert(ntarget+1,num);
        std::vector<std::string> lineread{line};
        long double scatlen;
        atoi[ntarget]=ntarget;
        if (!Tools::parse(lineread, "PARAMETERS"+num, scatlen, -1)) {
          error("Missing PARAMETERS or PARAMETERS not sorted");
        }
        for(unsigned k=0; k<numq; ++k) {
          FF_tmp[k][ntarget] = scatlen;
        }
        ntarget++;
      }
      if( ntarget!=size ) {
        error("found wrong number of PARAMETERS in file");
      }
      for(unsigned i=0; i<size; ++i) {
        Iq0+=FF_tmp[0][i];
      }
      Iq0 *= Iq0;
    }
  } else if(onebead) {
    if(saxs) {
      // read built-in ONEBEAD parameters
      FF_tmp_vac.resize(numq,std::vector<long double>(NONEBEAD));
      FF_tmp_mix.resize(numq,std::vector<long double>(NONEBEAD));
      FF_tmp_solv.resize(numq,std::vector<long double>(NONEBEAD));
      std::vector<std::vector<long double> > parameter_vac(NONEBEAD);
      std::vector<std::vector<long double> > parameter_mix(NONEBEAD);
      std::vector<std::vector<long double> > parameter_solv(NONEBEAD);
      getOnebeadparam(pdb, atoms, parameter_vac, parameter_mix, parameter_solv, residue_atom);
      // ion bead types carry no polynomial coefficients; insert a q=0 placeholder so the
      // expansion and Iq0 loops below stay valid. Their ion entries are overwritten with
      // analytic values by getOnebeadIonparam_saxs just afterwards.
      for(unsigned s=OB_NA; s<NONEBEAD; ++s) {
        parameter_vac[s].assign(1, 0.0);
        parameter_mix[s].assign(1, 0.0);
        parameter_solv[s].assign(1, 0.0);
      }
      for(unsigned i=0; i<nres; ++i) {
        if(atoi[i] >= NONEBEAD || parameter_vac[atoi[i]].empty() || parameter_mix[atoi[i]].empty() || parameter_solv[atoi[i]].empty()) {
          error("ONEBEAD form-factor parameters are not available for bead type " + std::to_string(atoi[i]) + " in this build");
        }
      }
      for(unsigned i=0; i<NONEBEAD; ++i) {
        for(unsigned k=0; k<numq; ++k) {
          for(unsigned j=0; j<parameter_vac[i].size(); ++j) {
            FF_tmp_vac[k][i]+= parameter_vac[i][j]*std::pow(static_cast<long double>(q_list[k]),j);
          }
          for(unsigned j=0; j<parameter_mix[i].size(); ++j) {
            FF_tmp_mix[k][i]+= parameter_mix[i][j]*std::pow(static_cast<long double>(q_list[k]),j);
          }
          for(unsigned j=0; j<parameter_solv[i].size(); ++j) {
            FF_tmp_solv[k][i]+= parameter_solv[i][j]*std::pow(static_cast<long double>(q_list[k]),j);
          }
        }
      }
      for(unsigned i=0; i<nres; ++i) {
        Iq0_vac[i]=parameter_vac[atoi[i]][0];
        Iq0_mix[i]=parameter_mix[atoi[i]][0];
        Iq0_solv[i]=parameter_solv[atoi[i]][0];
      }
      // overwrite ion bead types with analytic Cromer-Mann + excluded-volume terms
      getOnebeadIonparam_saxs(FF_tmp_vac, FF_tmp_mix, FF_tmp_solv, Iq0_vac, Iq0_mix, Iq0_solv);
    } else { // SANS
      // read built-in ONEBEAD parameters
      FF_tmp_vac_H.resize(numq,std::vector<long double>(NONEBEAD));
      FF_tmp_mix_H.resize(numq,std::vector<long double>(NONEBEAD));
      FF_tmp_solv_H.resize(numq,std::vector<long double>(NONEBEAD));
      FF_tmp_vac_D.resize(numq,std::vector<long double>(NONEBEAD));
      FF_tmp_mix_D.resize(numq,std::vector<long double>(NONEBEAD));
      std::vector<std::vector<long double> > parameter_vac_H(NONEBEAD);
      std::vector<std::vector<long double> > parameter_mix_H(NONEBEAD);
      std::vector<std::vector<long double> > parameter_solv_H(NONEBEAD);
      std::vector<std::vector<long double> > parameter_vac_D(NONEBEAD);
      std::vector<std::vector<long double> > parameter_mix_D(NONEBEAD);
      getOnebeadparam_sansH(pdb, atoms, parameter_vac_H, parameter_mix_H, parameter_solv_H);
      getOnebeadparam_sansD(pdb, atoms, parameter_vac_D, parameter_mix_D);
      // q=0 placeholders for ion bead types (overwritten by getOnebeadIonparam_sans below)
      for(unsigned s=OB_NA; s<NONEBEAD; ++s) {
        parameter_vac_H[s].assign(1, 0.0);
        parameter_mix_H[s].assign(1, 0.0);
        parameter_solv_H[s].assign(1, 0.0);
        parameter_vac_D[s].assign(1, 0.0);
        parameter_mix_D[s].assign(1, 0.0);
      }
      for(unsigned i=0; i<nres; ++i) {
        if(atoi[i] >= NONEBEAD || parameter_vac_H[atoi[i]].empty() || parameter_mix_H[atoi[i]].empty() || parameter_solv_H[atoi[i]].empty() ||
            parameter_vac_D[atoi[i]].empty() || parameter_mix_D[atoi[i]].empty()) {
          error("ONEBEAD form-factor parameters are not available for bead type " + std::to_string(atoi[i]) + " in this build");
        }
      }
      for(unsigned i=0; i<NONEBEAD; ++i) {
        for(unsigned k=0; k<numq; ++k) {
          for(unsigned j=0; j<parameter_vac_H[i].size(); ++j) { // same number of parameters
            FF_tmp_vac_H[k][i]+= parameter_vac_H[i][j]*std::pow(static_cast<long double>(q_list[k]),j);
            FF_tmp_vac_D[k][i]+= parameter_vac_D[i][j]*std::pow(static_cast<long double>(q_list[k]),j);
          }
          for(unsigned j=0; j<parameter_mix_H[i].size(); ++j) {
            FF_tmp_mix_H[k][i]+= parameter_mix_H[i][j]*std::pow(static_cast<long double>(q_list[k]),j);
            FF_tmp_mix_D[k][i]+= parameter_mix_D[i][j]*std::pow(static_cast<long double>(q_list[k]),j);
          }
          for(unsigned j=0; j<parameter_solv_H[i].size(); ++j) {
            FF_tmp_solv_H[k][i]+= parameter_solv_H[i][j]*std::pow(static_cast<long double>(q_list[k]),j);
          }
        }
      }
      for(unsigned i=0; i<nres; ++i) {
        Iq0_vac_H[i]=parameter_vac_H[atoi[i]][0];
        Iq0_mix_H[i]=parameter_mix_H[atoi[i]][0];
        Iq0_solv_H[i]=parameter_solv_H[atoi[i]][0];
        Iq0_vac_D[i]=parameter_vac_D[atoi[i]][0];
        Iq0_mix_D[i]=parameter_mix_D[atoi[i]][0];
      }
      // overwrite ion bead types with analytic terms (neutron b + excluded volume)
      getOnebeadIonparam_sans(FF_tmp_vac_H, FF_tmp_mix_H, FF_tmp_solv_H, FF_tmp_vac_D, FF_tmp_mix_D,
                              Iq0_vac_H, Iq0_mix_H, Iq0_solv_H, Iq0_vac_D, Iq0_mix_D);
    }
  } else if(martini) {
    // read built-in MARTINI parameters
    FF_tmp.resize(numq,std::vector<long double>(NMARTINI));
    std::vector<std::vector<long double> > parameter;
    parameter.resize(NMARTINI);
    if(martini_ff==3) {
      getMartiniFFparam_v3(atoms, parameter);
    } else {
      getMartiniFFparam(atoms, parameter);
    }
    for(unsigned i=0; i<NMARTINI; ++i) {
      for(unsigned k=0; k<numq; ++k) {
        for(unsigned j=0; j<parameter[i].size(); ++j) {
          FF_tmp[k][i]+= parameter[i][j]*std::pow(static_cast<long double>(q_list[k]),j);
        }
      }
    }
    for(unsigned i=0; i<size; ++i) {
      Iq0+=parameter[atoi[i]][0];
    }
    Iq0 *= Iq0;
  } else if(atomistic) {
    FF_tmp.resize(numq,std::vector<long double>(NTT));
    if(saxs) {
      Iq0=calculateAFF(atoms, FF_tmp, rho);
    } else {
      Iq0=calculateAFFsans(atoms, FF_tmp, deuter_conc);
    }
    Iq0 *= Iq0;
  }


  expint.resize( numq );
  ntarget=0;
  for(unsigned i=0; i<numq; ++i) {
    if( !parseNumbered( "EXPINT", i+1, expint[i] ) ) {
      break;
    }
    ntarget++;
  }
  std::transform(expint.begin(), expint.begin() + ntarget, expint.begin(), [scale_expint](double x) {
    return x / scale_expint;
  });
  bool exp=false;
  if(ntarget!=numq && ntarget!=0) {
    error("found wrong number of EXPINT values");
  }
  if(ntarget==numq) {
    exp=true;
  }
  if(getDoScore()&&!exp) {
    error("with DOSCORE you need to set the EXPINT values");
  }

  experr.resize( numq );
  ntarget=0;
  for(unsigned i=0; i<numq; ++i) {
    if( !parseNumbered( "EXPERR", i+1, experr[i] ) ) {
      break;
    }
    ntarget++;
  }
  std::transform(experr.begin(), experr.begin() + ntarget, experr.begin(), [scale_expint](double x) {
    return x / scale_expint;
  });
  bool experr_complete=false;
  if(ntarget!=numq && ntarget!=0) {
    error("found wrong number of EXPERR values");
  }
  if(ntarget==numq) {
    experr_complete=true;
    for(unsigned i=0; i<numq; ++i) {
      if(!std::isfinite(experr[i]) || experr[i]<=0.) {
        error("EXPERR values must be finite and greater than zero");
      }
    }
  }
  if(chi2) {
    if(absolute) {
      error("CHI2 cannot be combined with ABSOLUTE in this first implementation because CHI2 fits a multiplicative intensity scale");
    }
    if(numq<2) {
      error("CHI2 requires at least two QVALUE values because one multiplicative scale is fitted");
    }
    if(!exp) {
      error("CHI2 requires one EXPINT value for every QVALUE");
    }
    if(!experr_complete) {
      error("CHI2 requires one EXPERR value for every QVALUE");
    }
    for(unsigned i=0; i<numq; ++i) {
      if(!std::isfinite(expint[i])) {
        error("EXPINT values must be finite when using CHI2");
      }
    }
  }

  sigma_res.resize( numq );
  resolution=false;
  ntarget=0;
  for(unsigned i=0; i<numq; ++i) {
    if( !parseNumbered( "SIGMARES", i+1, sigma_res[i] ) ) {
      break;
    }
    ntarget++;
  }
  if(ntarget!=numq && ntarget!=0) {
    error("found wrong number of SIGMARES values");
  }
  if(ntarget==numq) {
    resolution=true;
    if(numq==1) {
      error("Resolution function requires at least two QVALUE values");
    }
    for(unsigned i=0; i<numq; ++i) {
      if(sigma_res[i]<=0.) {
        error("SIGMARES values must be greater than zero");
      }
      if(i>0&&q_list[i]<=q_list[i-1]) {
        error("QVALUE values must be strictly increasing when using SIGMARES");
      }
    }
  }

  if(chi2 && resolution) {
    error("CHI2 with SIGMARES is not enabled in the first validation stage; resolution-smeared CHI2 will be tested separately");
  }

  Nj = 10;
  parse("N", Nj);
  if (Nj < 2) {
    error("N should be larger than 1");
  }
  if (resolution) {
    log.printf("  Resolution function with N: %d\n", Nj);
  }

  FF_rank.resize(numq);
  unsigned n_atom_types;
  if(onebead) {
    FF_value.resize(nres,std::vector<double>(numq));
    n_atom_types=NONEBEAD;
    if(saxs) {
      FF_value_vacuum.resize(n_atom_types,std::vector<double>(numq));
      FF_value_solv.resize(n_atom_types,std::vector<double>(numq));
      FF_value_mixed.resize(n_atom_types,std::vector<double>(numq));
    } else {
      FF_value_vacuum_H.resize(n_atom_types,std::vector<double>(numq));
      FF_value_solv_H.resize(n_atom_types,std::vector<double>(numq));
      FF_value_mixed_H.resize(n_atom_types,std::vector<double>(numq));
      FF_value_vacuum_D.resize(n_atom_types,std::vector<double>(numq));
      FF_value_mixed_D.resize(n_atom_types,std::vector<double>(numq));
    }
  } else {
    FF_value.resize(size,std::vector<double>(numq));
  }
  for(unsigned k=0; k<numq; ++k) {
    if(!onebead) {
      for(unsigned i=0; i<size; ++i) {
        FF_value[i][k] = static_cast<double>(FF_tmp[k][atoi[i]])/(std::sqrt(Iq0));
      }
      for(unsigned i=0; i<size; ++i) {
        FF_rank[k] += FF_value[i][k]*FF_value[i][k];
      }
    } else {
      if(saxs) {
        for(unsigned i=0; i<n_atom_types; ++i) {
          FF_value_vacuum[i][k] = static_cast<double>(FF_tmp_vac[k][i]);
          FF_value_mixed[i][k] = static_cast<double>(FF_tmp_mix[k][i]);
          FF_value_solv[i][k] = static_cast<double>(FF_tmp_solv[k][i]);
        }
      } else { // SANS
        for(unsigned i=0; i<n_atom_types; ++i) {
          FF_value_vacuum_H[i][k] = static_cast<double>(FF_tmp_vac_H[k][i]);
          FF_value_mixed_H[i][k] = static_cast<double>(FF_tmp_mix_H[k][i]);
          FF_value_solv_H[i][k] = static_cast<double>(FF_tmp_solv_H[k][i]);
          FF_value_vacuum_D[i][k] = static_cast<double>(FF_tmp_vac_D[k][i]);
          FF_value_mixed_D[i][k] = static_cast<double>(FF_tmp_mix_D[k][i]);
        }
      }
    }
  }


  if(!getDoScore()) {
    for(unsigned i=0; i<numq; ++i) {
      std::string num;
      Tools::convert(i,num);
      addComponentWithDerivatives("q-"+num);
      componentIsNotPeriodic("q-"+num);
    }
    if(exp) {
      for(unsigned i=0; i<numq; ++i) {
        std::string num;
        Tools::convert(i,num);
        addComponent("exp-"+num);
        componentIsNotPeriodic("exp-"+num);
        Value* comp=getPntrToComponent("exp-"+num);
        comp->set(expint[i]);
      }
    }
    if(chi2) {
      addComponentWithDerivatives("chi2");
      componentIsNotPeriodic("chi2");
      addComponentWithDerivatives("chi2red");
      componentIsNotPeriodic("chi2red");
      addComponent("chi2scale");
      componentIsNotPeriodic("chi2scale");
    }
  } else {
    for(unsigned i=0; i<numq; ++i) {
      std::string num;
      Tools::convert(i,num);
      addComponent("q-"+num);
      componentIsNotPeriodic("q-"+num);
    }
    for(unsigned i=0; i<numq; ++i) {
      std::string num;
      Tools::convert(i,num);
      addComponent("exp-"+num);
      componentIsNotPeriodic("exp-"+num);
      Value* comp=getPntrToComponent("exp-"+num);
      comp->set(expint[i]);
    }
    if(chi2) {
      addComponent("chi2");
      componentIsNotPeriodic("chi2");
      addComponent("chi2red");
      componentIsNotPeriodic("chi2red");
      addComponent("chi2scale");
      componentIsNotPeriodic("chi2scale");
    }
  }
  if(chi2ensemble) {
    addComponent("chi2ensemble");
    componentIsNotPeriodic("chi2ensemble");
    addComponent("chi2ensemblered");
    componentIsNotPeriodic("chi2ensemblered");
    addComponent("chi2ensemblescale");
    componentIsNotPeriodic("chi2ensemblescale");
  }

  // convert units to nm^-1
  for(unsigned i=0; i<numq; ++i) {
    q_list[i]=q_list[i]*10.0;    // factor 10 to convert from A^-1 to nm^-1
    if (resolution) {
      sigma_res[i]=sigma_res[i]*10.0;
    }
  }

  // q values as a uniform grid plus small deviations, used by sincos_q
  qref0_ = q_list[0];
  qref_dq_ = (numq>1) ? (q_list[numq-1]-q_list[0])/static_cast<double>(numq-1) : 0.0;
  qref_maxeps_ = 0.0;
  qref_eps_.resize(numq);
  for(unsigned k=0; k<numq; ++k) {
    qref_eps_[k] = q_list[k]-(qref0_+k*qref_dq_);
    qref_maxeps_ = std::max(qref_maxeps_, std::fabs(qref_eps_[k]));
  }

  // the deuteration of the exposed residues is drawn by the first rank and broadcast
  random_deut_.setSeed(-1);

  // compute resolution function after converting units
  if (resolution) {
    qj_list.resize(numq, std::vector<double>(Nj));
    Rij.resize(numq, std::vector<double>(Nj));
    // compute Rij and qj_list
    resolution_function();
    // matrix of the smearing (spline interpolation and integral with the resolution function)
    smear_.assign(static_cast<std::size_t>(numq)*numq,0.0);
    for(unsigned l=0; l<numq; ++l) {
      std::vector<double> unit(numq,0.0);
      unit[l] = 1.0;
      std::vector<SplineCoeffs> coeffs = spline_coeffs(q_list, unit);
      for(unsigned k=0; k<numq; ++k) {
        const double dq = qj_list[k][1] - qj_list[k][0];
        for(unsigned j=0; j<Nj; ++j) {
          smear_[k*numq+l] += Rij[k][j] * interpolation(coeffs, qj_list[k][j]) * dq;
        }
      }
    }
  }

  log<<"  Bibliography ";
  if(onebead) {
    log<<plumed.cite("Ballabio, Paissoni, Bollati, de Rosa, Capelli, Camilloni, J. Chem. Theory Comput., 19, 22, 8401-8413 (2023)");
  }
  if(martini) {
    log<<plumed.cite("Niebling, Björling, Westenhoff, J. Appl. Crystallogr., 47, 1190–1198 (2014)");
    log<<plumed.cite("Paissoni, Jussupow, Camilloni, J. Appl. Crystallogr., 52, 394-402 (2019)");
    if(martini_ff==3) {
      log<<plumed.cite("Vinterbladh, Vuillemot, Grudinin, Lund, ChemRxiv, doi:10.26434/chemrxiv.15006205/v1 (2026)");
    }
  }
  if(atomistic) {
    log<<plumed.cite("Fraser, MacRae, Suzuki, J. Appl. Crystallogr., 11, 693–694 (1978)");
    log<<plumed.cite("Brown, Fox, Maslen, O'Keefe, Willis, International Tables for Crystallography, C, 554–595 (International Union of Crystallography, 2006)");
  }
  if(resolution) {
    log<<plumed.cite("Pedersen, Posselt, Mortensen, J. Appl. Crystallogr., 23, 321–333 (1990)");
  }

  // ion parameters: cite the sources only if any ion is actually present in the system
  bool ions_present=false;
  if(onebead) {
    for(unsigned i=0; i<nres; ++i) {
      if(atoi[i]>=OB_NA && atoi[i]<NONEBEAD) {
        ions_present=true;
        break;
      }
    }
  } else if(atomistic) {
    for(unsigned i=0; i<atoi.size(); ++i) {
      if(atoi[i]>=AFF_NA && atoi[i]<NTT) {
        ions_present=true;
        break;
      }
    }
  }
  if(ions_present) {
    log<<plumed.cite("Cromer, Mann, Acta Crystallogr. A, 24, 321-324 (1968)");
    log<<plumed.cite("Brown, Fox, Maslen, O'Keefe, Willis, International Tables for Crystallography, C, 554–595 (International Union of Crystallography, 2006)");
    log<<plumed.cite("Shannon, Acta Crystallogr. A, 32, 751-767 (1976)");
    if(!saxs) {
      log<<plumed.cite("Sears, Neutron News, 3, 26-37 (1992)");
    }
  }

  log<< plumed.cite("Bonomi, Camilloni, Bioinformatics, 33, 3999 (2017)");
  log<<"\n";

  requestAtoms(atoms, false);

  if(getDoScore()) {
    setParameters(expint);
    Initialise(numq);
  }
  setDerivatives();
  checkRead();

  if(onebead) {
    // the atoms of each bead, in increasing order, with their weights in the centre
    // of mass of the bead (the masses are summed in the order of the atoms)
    bead_off_.assign(nres+1,0);
    for(unsigned a=0; a<size; ++a) {
      bead_off_[residue_atom[a]+1]++;
    }
    for(unsigned r=0; r<nres; ++r) {
      bead_off_[r+1] += bead_off_[r];
    }
    bead_atoms_.resize(size);
    std::vector<unsigned> fill(bead_off_.begin(), bead_off_.end()-1);
    std::vector<double> bead_mass(nres,0.);
    for(unsigned a=0; a<size; ++a) {
      bead_atoms_[fill[residue_atom[a]]++] = a;
      bead_mass[residue_atom[a]] += atoms_masses[a];
    }
    bead_inv_mass_.resize(nres);
    for(unsigned r=0; r<nres; ++r) {
      bead_inv_mass_[r] = 1.0/bead_mass[r];
    }
    atom_weight_.resize(size);
    for(unsigned a=0; a<size; ++a) {
      atom_weight_[a] = atoms_masses[a]*bead_inv_mass_[residue_atom[a]];
    }
  }
}

// LCPO SASA of the atoms (ONEBEAD only).
// Two atoms with LCPO parameters are neighbours if d_ij < r_i + r_j (Angstrom).
// A Verlet list of candidates (d_ij < r_i + r_j + skin) is kept between the calls
// and rebuilt only when an atom has moved by more than about skin/2 since it was
// built; at every call the exact test is applied to the candidates, so the
// neighbours are the same as with a new search. The rows are sorted, so all the
// sums are done in the same order as with an all-pairs search.
namespace {
// skin of the candidate list (Angstrom)
constexpr double lcpo_skin = 2.0;
}

void SAXS::lcpoBuildCandidates() {
  const unsigned natoms = getNumberOfAtoms();
  double maxr = 0.0;
  for(const auto i : lcpo_atoms_) {
    maxr = std::max(maxr, LCPOparam[i][0]);
  }
  // cell size (nm) not smaller than the largest cutoff
  const double rc = (2.0*maxr+lcpo_skin)/10.;
  Vector lo = getPosition(lcpo_atoms_[0]), hi = lo;
  for(const auto i : lcpo_atoms_) {
    for(unsigned k=0; k<3; ++k) {
      lo[k] = std::min(lo[k], getPosition(i)[k]);
      hi[k] = std::max(hi[k], getPosition(i)[k]);
    }
  }
  int ncell[3];
  double csize[3];
  for(unsigned k=0; k<3; ++k) {
    const double len = hi[k]-lo[k];
    ncell[k] = std::max(1, std::min(200, static_cast<int>(std::floor(len/rc))));
    csize[k] = std::max(len/static_cast<double>(ncell[k]), rc);
  }
  auto cell=[&](const unsigned i, int *c) {
    for(unsigned k=0; k<3; ++k) {
      c[k] = std::min(ncell[k]-1, std::max(0, static_cast<int>(std::floor((getPosition(i)[k]-lo[k])/csize[k]))));
    }
  };
  const unsigned ncells = ncell[0]*ncell[1]*ncell[2];
  const unsigned nlcpo = lcpo_atoms_.size();
  std::vector<unsigned> cell_start(ncells+1,0), cell_atoms(nlcpo), atom_cell(nlcpo);
  for(unsigned a=0; a<nlcpo; ++a) {
    int c[3];
    cell(lcpo_atoms_[a],c);
    atom_cell[a] = (c[0]*ncell[1]+c[1])*ncell[2]+c[2];
    cell_start[atom_cell[a]+1]++;
  }
  for(unsigned ic=0; ic<ncells; ++ic) {
    cell_start[ic+1] += cell_start[ic];
  }
  {
    std::vector<unsigned> fill(cell_start.begin(), cell_start.end()-1);
    for(unsigned a=0; a<nlcpo; ++a) {
      cell_atoms[fill[atom_cell[a]]++] = lcpo_atoms_[a];
    }
  }
  std::vector<std::vector<unsigned> > rows(natoms);
  #pragma omp parallel for num_threads(OpenMP::getNumThreads()) schedule(dynamic,64)
  for(unsigned a=0; a<nlcpo; ++a) {
    const unsigned i = lcpo_atoms_[a];
    int c[3];
    cell(i,c);
    for(int i0=std::max(0,c[0]-1); i0<=std::min(ncell[0]-1,c[0]+1); ++i0) {
      for(int i1=std::max(0,c[1]-1); i1<=std::min(ncell[1]-1,c[1]+1); ++i1) {
        for(int i2=std::max(0,c[2]-1); i2<=std::min(ncell[2]-1,c[2]+1); ++i2) {
          const unsigned ic = (i0*ncell[1]+i1)*ncell[2]+i2;
          for(unsigned jj=cell_start[ic]; jj<cell_start[ic+1]; ++jj) {
            const unsigned j = cell_atoms[jj];
            if(j==i) {
              continue;
            }
            const double cut = (LCPOparam[i][0]+LCPOparam[j][0]+lcpo_skin)/10.;
            if(modulo2(delta(getPosition(i), getPosition(j))) < cut*cut) {
              rows[i].push_back(j);
            }
          }
        }
      }
    }
    std::sort(rows[i].begin(), rows[i].end());
  }
  lcpo_cand_off_.assign(natoms+1,0);
  for(unsigned i=0; i<natoms; ++i) {
    lcpo_cand_off_[i+1] = lcpo_cand_off_[i]+rows[i].size();
  }
  lcpo_cand_idx_.resize(lcpo_cand_off_[natoms]);
  for(unsigned i=0; i<natoms; ++i) {
    std::copy(rows[i].begin(), rows[i].end(), lcpo_cand_idx_.begin()+lcpo_cand_off_[i]);
  }
  lcpo_ref_pos_.resize(nlcpo);
  for(unsigned a=0; a<nlcpo; ++a) {
    lcpo_ref_pos_[a] = getPosition(lcpo_atoms_[a]);
  }
}

void SAXS::lcpoUpdateList() {
  const unsigned natoms = getNumberOfAtoms();
  if(lcpo_ref_pos_.empty()) {
    for(unsigned i = 0; i < natoms; ++i) {
      if(LCPOparam[i].size()>0) {
        lcpo_atoms_.push_back(i);
        if(LCPOparam[i].size()>1 && LCPOparam[i][1]>0.0) {
          lcpo_area_atoms_.push_back(i);
        }
      }
    }
    if(lcpo_atoms_.empty()) {
      lcpo_off_.assign(natoms+1,0);
      return;
    }
    lcpoBuildCandidates();
  } else {
    // rebuild when an atom has moved by more than 0.45*skin (0.05*skin of margin
    // for the rounding of the candidate test)
    const double maxd = 0.45*lcpo_skin/10.;
    const unsigned nlcpo = lcpo_atoms_.size();
    bool moved = false;
    #pragma omp parallel for num_threads(OpenMP::getNumThreads()) reduction(||:moved)
    for(unsigned a=0; a<nlcpo; ++a) {
      if(modulo2(delta(lcpo_ref_pos_[a], getPosition(lcpo_atoms_[a]))) >= maxd*maxd) {
        moved = true;
      }
    }
    if(moved) {
      lcpoBuildCandidates();
    }
  }
  // the neighbours among the candidates, and the overlap term A(i,j) of each
  // directed pair, with the same expressions as in the original LCPO loop
  lcpo_cand_d_.resize(lcpo_cand_idx_.size());
  lcpo_off_.assign(natoms+1,0);
  #pragma omp parallel for num_threads(OpenMP::getNumThreads()) schedule(dynamic,64)
  for(unsigned i=0; i<natoms; ++i) {
    std::size_t n = 0;
    for(std::size_t e=lcpo_cand_off_[i]; e<lcpo_cand_off_[i+1]; ++e) {
      const unsigned j = lcpo_cand_idx_[e];
      const double d_ij = modulo(delta(getPosition(i), getPosition(j)))*10.;
      lcpo_cand_d_[e] = d_ij;
      if(d_ij < LCPOparam[i][0]+LCPOparam[j][0]) {
        n++;
      }
    }
    lcpo_off_[i+1] = n;
  }
  for(unsigned i=0; i<natoms; ++i) {
    lcpo_off_[i+1] += lcpo_off_[i];
  }
  lcpo_idx_.resize(lcpo_off_[natoms]);
  lcpo_term_.resize(lcpo_off_[natoms]);
  #pragma omp parallel for num_threads(OpenMP::getNumThreads()) schedule(dynamic,64)
  for(unsigned i=0; i<natoms; ++i) {
    std::size_t n = lcpo_off_[i];
    for(std::size_t e=lcpo_cand_off_[i]; e<lcpo_cand_off_[i+1]; ++e) {
      const unsigned j = lcpo_cand_idx_[e];
      const double d_ij = lcpo_cand_d_[e];
      if(d_ij < LCPOparam[i][0]+LCPOparam[j][0]) {
        const double ri = LCPOparam[i][0];
        const double rj = LCPOparam[j][0];
        lcpo_idx_[n] = j;
        lcpo_term_[n] = (2.*pi*ri*(ri-d_ij/2.-((ri*ri-rj*rj)/(2.*d_ij))));
        n++;
      }
    }
  }
}

// calculates SASA according to LCPO algorithm
void SAXS::sasa_calculate(std::vector<bool> &solv_res) {
  const unsigned natoms = getNumberOfAtoms();
  lcpoUpdateList();
  // the area of each atom; with MPI the atoms are split over the ranks
  std::vector<double> atom_area(natoms, 0.);
  const unsigned nt = OpenMP::getNumThreads();
  const bool split = !serial && comm.Get_size()>1;
  const unsigned first = split ? comm.Get_rank() : 0;
  const unsigned stride = split ? comm.Get_size() : 1;
  const unsigned narea = lcpo_area_atoms_.size();
  #pragma omp parallel num_threads(nt)
  {
    // all bits set for the neighbours of atom i, 0 otherwise
    std::vector<std::uint64_t> is_neighbor(natoms, 0);
    #pragma omp for schedule(dynamic,64)
    for (unsigned a = first; a < narea; a += stride) {
      const unsigned i = lcpo_area_atoms_[a];
      double Aij = 0.0;
      double Aijk = 0.0;
      double Ajk = 0.0;
      const double ri = LCPOparam[i][0];
      const double S1 = 4.*pi*ri*ri;
      for (std::size_t e = lcpo_off_[i]; e < lcpo_off_[i+1]; ++e) {
        is_neighbor[lcpo_idx_[e]] = ~std::uint64_t(0);
      }
      for (std::size_t e = lcpo_off_[i]; e < lcpo_off_[i+1]; ++e) {
        const unsigned j = lcpo_idx_[e];
        const double Aijt = lcpo_term_[e];
        // sum of A(j,k) over the neighbours k of j that are also neighbours of i,
        // in the order of the list of j. The mask selects A(j,k) or +0 without a
        // branch (about half of the k are neighbours of i, so a branch is
        // mispredicted often); adding +0 does not change the sum.
        double Ajkt = 0.0;
        for (std::size_t f = lcpo_off_[j]; f < lcpo_off_[j+1]; ++f) {
          std::uint64_t bits;
          std::memcpy(&bits, &lcpo_term_[f], sizeof(bits));
          bits &= is_neighbor[lcpo_idx_[f]];
          double term;
          std::memcpy(&term, &bits, sizeof(term));
          Ajkt += term;
        }
        Aijk += (Aijt * Ajkt);
        Aij += Aijt;
        Ajk += Ajkt;
      }
      for (std::size_t e = lcpo_off_[i]; e < lcpo_off_[i+1]; ++e) {
        is_neighbor[lcpo_idx_[e]] = 0;
      }
      atom_area[i] = (LCPOparam[i][1]*S1+LCPOparam[i][2]*Aij+LCPOparam[i][3]*Ajk+LCPOparam[i][4]*Aijk);
    }
  }
  if(split) {
    // the other ranks contribute zeros: the sum is exact
    comm.Sum(atom_area);
  }
  // residue areas: per-thread sums over a static partition of the atoms, combined
  // in thread order (reproducible results); the buffers of the threads that did
  // not run (e.g. without OpenMP, or when the runtime gives fewer threads than
  // requested) are empty and are skipped
  std::vector<double> sasares(nres, 0.);
  std::vector< std::vector<double> > sasares_t(nt);
  #pragma omp parallel num_threads(nt)
  {
    std::vector<double> & private_sasares = sasares_t[OpenMP::getThreadNum()];
    private_sasares.assign(nres, 0.);
    #pragma omp for schedule(static)
    for (unsigned i = 0; i < natoms; ++i) {
      if (LCPOparam[i].size() > 1 && LCPOparam[i][1] > 0.0) {
        const double sasai = atom_area[i];
        if (sasai > 0) {
          private_sasares[residue_atom[i]] += sasai / 100.0;
        }
      }
    }
  }
  for (unsigned t = 0; t < nt; ++t) {
    if(sasares_t[t].empty()) {
      continue;
    }
    for (unsigned i = 0; i < nres; ++i) {
      sasares[i] += sasares_t[t][i];
    }
  }
  for(unsigned i=0; i<nres; ++i) { // updating solv_res based on sasares
    if(sasares[i]>sasa_cutoff) {
      solv_res[i] = 1;
    } else {
      solv_res[i] = 0;
    }
  }
}


void SAXS::sincos_q(const double r, double *s, double *c) const {
  // sin(q_k r) and cos(q_k r) for all the q values. When the q values are close to a
  // uniform grid q0+k*dq (as usual for experimental data), sin and cos of the grid are
  // obtained by angle addition (exact restart every 16 values) and corrected for the
  // small deviations eps_k=q_k-q0-k*dq with a Taylor expansion (|eps_k r|<=0.1, error
  // smaller than 1e-16). Otherwise sin and cos are computed for every q.
  const unsigned numq = q_list.size();
  if(numq<4 || qref_maxeps_*r>0.1) {
    for(unsigned k=0; k<numq; ++k) {
      s[k] = std::sin(q_list[k]*r);
      c[k] = std::cos(q_list[k]*r);
    }
    return;
  }
  const double sd = std::sin(qref_dq_*r);
  const double cd = std::cos(qref_dq_*r);
  double sk = 0.0, ck = 1.0;
  for(unsigned k=0; k<numq; ++k) {
    if((k&15)==0) {
      const double a = (qref0_+k*qref_dq_)*r;
      sk = std::sin(a);
      ck = std::cos(a);
    }
    const double x = qref_eps_[k]*r;
    const double x2 = x*x;
    const double cx = 1.-x2*(1./2-x2*(1./24-x2*(1./720-x2*(1./40320))));
    const double sx = x*(1.-x2*(1./6-x2*(1./120-x2*(1./5040-x2*(1./362880)))));
    s[k] = sk*cx+ck*sx;
    c[k] = ck*cx-sk*sx;
    const double s2 = sk*cd+ck*sd;
    ck = ck*cd-sk*sd;
    sk = s2;
  }
}

void SAXS::debye_cpu(const std::vector<Vector> &pos, std::vector<double> &sum, std::vector<Vector> &deriv) {
  const unsigned size = pos.size();
  const unsigned numq = q_list.size();
  const unsigned rank = first_row_;
  const unsigned stride = row_stride_;
  const unsigned nt=OpenMP::getNumThreads();
  // partial sums of the threads other than the first one, added in thread order
  std::vector< std::vector<Vector> > deriv_t(nt);
  std::vector< std::vector<double> > sum_t(nt);
  #pragma omp parallel num_threads(nt)
  {
    const unsigned t = OpenMP::getThreadNum();
    if(t>0) {
      deriv_t[t].assign(deriv.size(),Vector(0,0,0));
      sum_t[t].assign(numq,0.0);
    }
    std::vector<Vector> & my_deriv = (t>0) ? deriv_t[t] : deriv;
    std::vector<double> & my_sum = (t>0) ? sum_t[t] : sum;
    std::vector<double> sq(numq), cq(numq);
    // the work of row i decreases with i: rows are dealt round robin
    #pragma omp for schedule(static,1)
    for (unsigned i=rank; i<size-1; i+=stride) {
      Vector posi = pos[i];
      for (unsigned j=i+1; j<size ; ++j) {
        Vector c_distances = delta(posi,pos[j]);
        double m_distances = c_distances.modulo();
        c_distances = c_distances/m_distances/m_distances;
        sincos_q(m_distances, sq.data(), cq.data());
        for (unsigned k=0; k<numq; ++k) {
          unsigned kdx=k*size;
          double qdist = q_list[k]*m_distances;
          double FFF = 2.*FF_value[i][k]*FF_value[j][k];
          double tsq = sq[k]/qdist;
          double tmp = FFF*(cq[k]-tsq);
          Vector dd  = c_distances*tmp;
          my_deriv[kdx+i] -= dd;
          my_deriv[kdx+j] += dd;
          my_sum[k] += FFF*tsq;
        }
      }
    }
  }
  // only the threads that actually ran have filled their buffers: without OpenMP
  // (or when the runtime gives fewer threads than requested) the others are empty
  std::vector<unsigned> filled;
  for(unsigned t=1; t<nt; ++t) {
    if(sum_t[t].empty()) {
      continue;
    }
    filled.push_back(t);
    for(unsigned k=0; k<numq; ++k) {
      sum[k] += sum_t[t][k];
    }
  }
  if(!filled.empty()) {
    #pragma omp parallel for num_threads(nt)
    for(unsigned i=0; i<deriv.size(); ++i) {
      for(const auto t : filled) {
        deriv[i] += deriv_t[t][i];
      }
    }
  }
}

void SAXS::calculate_cpu(std::vector<Vector> &pos, std::vector<Vector> &deriv, const bool reduce_deriv) {
  unsigned size;
  if(onebead) {
    size = nres;
  } else {
    size = getNumberOfAtoms();
  }
  const unsigned numq = q_list.size();
  if(size == 0) {
    error("SAXS/SANS requires at least one atom or bead");
  }

  std::vector<double> sum(numq,0);
  if(has_rows_) {
    if(usegpu_) {
      ff_flat_.resize(static_cast<std::size_t>(size)*numq);
      for(unsigned i=0; i<size; ++i) {
        std::copy(FF_value[i].begin(), FF_value[i].end(), ff_flat_.begin()+static_cast<std::size_t>(i)*numq);
      }
      // with DOSCORE the derivatives stay on the GPU and are contracted later
      backend_->compute(pos, ff_flat_, q_list, first_row_, row_stride_, sum, getDoScore() ? nullptr : &deriv);
    } else {
      debye_cpu(pos, sum, deriv);
    }
  }

  if(reduce_) {
    comm.Sum(&sum[0], numq);
  }

  if (resolution) {
    // the resolution function is applied to the whole intensity, self term included
    for (unsigned k=0; k<numq; ++k) {
      sum[k]+=FF_rank[k];
    }
    const std::vector<double> old_sum(sum);
    for (unsigned k=0; k<numq; ++k) {
      sum[k] = 0.;
      for (unsigned l=0; l<numq; ++l) {
        sum[k] += smear_[k*numq+l]*old_sum[l];
      }
    }
    // the smearing is linear in the data, so it is applied to the partial
    // derivatives of this rank (of all the particles) before summing them
    // (with USEGPU and DOSCORE the derivatives stay on the GPU: the weights of
    // the contraction are smeared instead, see calculate)
    const unsigned nsmear = (usegpu_ && getDoScore()) ? 0 : size;
    #pragma omp parallel for num_threads(OpenMP::getNumThreads()) schedule(static)
    for (unsigned i=0; i<nsmear; i++) {
      std::vector<Vector> old_deriv_i(numq);
      for (unsigned l=0; l<numq; l++) {
        old_deriv_i[l] = deriv[l*size+i];
      }
      for (unsigned k=0; k<numq; k++) {
        Vector d;
        for (unsigned l=0; l<numq; l++) {
          d += smear_[k*numq+l]*old_deriv_i[l];
        }
        deriv[k*size+i] = d;
      }
    }
  } else {
    for (unsigned k=0; k<numq; ++k) {
      sum[k]+=FF_rank[k];
    }
  }

  // with DOSCORE only the contraction with the score derivatives is summed (in calculate)
  if(reduce_ && reduce_deriv) {
    comm.Sum(&deriv[0][0], 3*deriv.size());
  }

  for (unsigned k=0; k<numq; ++k) {
    std::string num;
    Tools::convert(k,num);
    Value* val=getPntrToComponent("q-"+num);
    val->set(sum[k]);
    if(getDoScore()) {
      setCalcData(k, sum[k]);
    }
  }
}

void SAXS::calculate() {
  if(pbc) {
    makeWhole();
  }

  const size_t size = getNumberOfAtoms();
  const size_t numq = q_list.size();
  const size_t beads_size = onebead ? nres : size;

  // the derivatives of the intensities with respect to the beads, [k*beads_size+i]
  // (with USEGPU and DOSCORE they stay on the GPU)
  std::vector<Vector> bd_deriv((usegpu_ && getDoScore()) ? 0 : numq*beads_size);

  std::vector<Vector> beads_pos(beads_size);
  if(onebead) {
    // centres of mass of the beads; the derivative of a bead position with respect
    // to one of its atoms is atom_weight_ (mass of the atom / mass of the bead)
    #pragma omp parallel for num_threads(OpenMP::getNumThreads()) schedule(static)
    for(unsigned r=0; r<nres; ++r) {
      Vector sum_pos(0,0,0);
      for(unsigned e=bead_off_[r]; e<bead_off_[r+1]; ++e) {
        const unsigned a = bead_atoms_[e];
        sum_pos += atoms_masses[a] * getPosition(a);
      }
      beads_pos[r] = sum_pos*bead_inv_mass_[r];
    }
    // SASA
    std::vector<bool> solv_res(nres, 0);
    if(saxs) {
      if(getStep()%solv_stride == 0 || isFirstStep) {
        isFirstStep = 0;
        if(rho_corr!=rho) {
          sasa_calculate(solv_res);
        }
        Iq0=0.;
        for(unsigned i=0; i<nres; ++i) {
          // Iq0 accumulates the MAGNITUDE of each bead amplitude and is left as
          // sum |A| on purpose: it stays well conditioned near a match point,
          // where sum A passes through zero, and it equals |sum A| whenever
          // every bead shares a sign, which is always the case in SAXS.
          const double rho_i = (solv_res[i] == 1) ? rho_corr : rho;
          bead_sign[i] = beadSign(Iq0_solv[i], Iq0_mix[i], rho_i);
          Iq0 += std::sqrt(std::fabs(Iq0_vac[i] + rho_i*rho_i*Iq0_solv[i] - rho_i*Iq0_mix[i]));
        }
        // Form Factors
        for(unsigned k=0; k<numq; ++k) {
          for(unsigned i=0; i<nres; ++i) {
            if(solv_res[i] == 0) { // buried
              FF_value[i][k] = bead_sign[i]*std::sqrt(std::fabs(FF_value_vacuum[atoi[i]][k] + rho*rho*FF_value_solv[atoi[i]][k] - rho*FF_value_mixed[atoi[i]][k]))/Iq0;
            } else { // surface
              FF_value[i][k] = bead_sign[i]*std::sqrt(std::fabs(FF_value_vacuum[atoi[i]][k] + rho_corr*rho_corr*FF_value_solv[atoi[i]][k] - rho_corr*FF_value_mixed[atoi[i]][k]))/Iq0;
            }

          }
        }
        for(unsigned k=0; k<numq; ++k) {
          FF_rank[k]=0.;
          for(unsigned i=0; i<nres; ++i) {
            FF_rank[k]+=FF_value[i][k]*FF_value[i][k];
          }
        }

      }
    } else { // SANS
      std::vector<bool> deut_res(nres, 0);
      double solv_sc_length = 0.1*(0.580 + 2.*((1. - deuter_conc) * (-0.374) + deuter_conc * 0.667)); // per water electron (10 electrons)
      double rho_sans = rho * solv_sc_length;
      double rho_sans_corr = rho_corr * solv_sc_length;
      if(getStep()%solv_stride == 0 || isFirstStep) {
        isFirstStep = 0;
        if(deuter_conc!=0.||rho != rho_corr) {
          sasa_calculate(solv_res);
        }
        // the exposed residues are deuterated with probability deuter_conc: the random
        // numbers are drawn by the first rank (same sequence as before on every rank)
        std::vector<int> deut_draw(nres,0);
        if(comm.Get_rank()==0) {
          for(unsigned i=0; i<nres; ++i) {
            if(solv_res[i] == 1) {
              deut_draw[i] = (random_deut_.RandU01()<deuter_conc) ? 1 : 0;
            }
          }
        }
        if(comm.Get_size()>1) {
          comm.Bcast(deut_draw,0);
        }
        Iq0=0.;
        for(unsigned i=0; i<nres; ++i) {
          // The density and the H/D term set are selected once and used for both
          // the sign and the magnitude, so the two can never disagree. In SANS
          // the sign genuinely varies between beads near a contrast-match point.
          double vac_i, mix_i, rho_i;
          if(solv_res[i] == 1) {
            rho_i = rho_sans_corr;
            if(deut_draw[i]) {
              deut_res[i] = 1;
              vac_i = Iq0_vac_D[i];
              mix_i = Iq0_mix_D[i];
            } else {
              vac_i = Iq0_vac_H[i];
              mix_i = Iq0_mix_H[i];
            }
          } else {
            rho_i = rho_sans;
            vac_i = Iq0_vac_H[i];
            mix_i = Iq0_mix_H[i];
          }
          bead_sign[i] = beadSign(Iq0_solv_H[i], mix_i, rho_i);
          Iq0 += std::sqrt(std::fabs(vac_i + rho_i*rho_i*Iq0_solv_H[i] - rho_i*mix_i));
        }
        // Form Factors
        for(unsigned k=0; k<numq; ++k) {
          for(unsigned i=0; i<nres; ++i) {
            if(solv_res[i] == 0) { // hydrogen
              FF_value[i][k] = bead_sign[i]*std::sqrt(std::fabs(FF_value_vacuum_H[atoi[i]][k] + rho_sans*rho_sans*FF_value_solv_H[atoi[i]][k] - rho_sans*FF_value_mixed_H[atoi[i]][k]))/Iq0;
            } else {
              if(deut_res[i] == 0) {
                FF_value[i][k] = bead_sign[i]*std::sqrt(std::fabs(FF_value_vacuum_H[atoi[i]][k] + rho_sans_corr*rho_sans_corr*FF_value_solv_H[atoi[i]][k] - rho_sans_corr*FF_value_mixed_H[atoi[i]][k]))/Iq0;
              } else {
                FF_value[i][k] = bead_sign[i]*std::sqrt(std::fabs(FF_value_vacuum_D[atoi[i]][k] + rho_sans_corr*rho_sans_corr*FF_value_solv_H[atoi[i]][k] - rho_sans_corr*FF_value_mixed_D[atoi[i]][k]))/Iq0;
              }
            }

          }
        }
        for(unsigned k=0; k<numq; ++k) {
          FF_rank[k]=0.;
          for(unsigned i=0; i<nres; ++i) {
            FF_rank[k]+=FF_value[i][k]*FF_value[i][k];
          }
        }

      }
    }
    // not ONEBEAD
  } else {
    for(unsigned i=0; i<size; ++i) {
      beads_pos[i] = getPosition(i);
    }
  }

  // with DOSCORE the derivatives of the intensities are first contracted with the
  // derivatives of the score, and only then summed over the ranks
  calculate_cpu(beads_pos, bd_deriv, !getDoScore());

  std::vector<double> chi2_dq;
  if(chi2 && !getDoScore()) {
    chi2_dq.assign(numq,0.);
  }
  double chi2_raw=0.;
  double chi2_reduced=0.;
  double chi2_scale=1.;
  if(chi2) {
    double sum_xy=0.;
    double sum_xx=0.;
    double sum_yy=0.;
    std::vector<double> calcint(numq,0.);
    for(unsigned k=0; k<numq; ++k) {
      std::string num;
      Tools::convert(k,num);
      calcint[k]=getPntrToComponent("q-"+num)->get();
      if(!std::isfinite(calcint[k])) {
        error("Non-finite calculated intensity encountered while evaluating CHI2");
      }
      const double weight=1./(experr[k]*experr[k]);
      sum_xy += weight*calcint[k]*expint[k];
      sum_xx += weight*calcint[k]*calcint[k];
      sum_yy += weight*expint[k]*expint[k];
    }
    if(!std::isfinite(sum_xx) || sum_xx<=0.) {
      error("Cannot evaluate CHI2 because the weighted calculated-intensity norm is not positive and finite");
    }
    chi2_scale=sum_xy/sum_xx;
    chi2_raw=sum_yy-sum_xy*sum_xy/sum_xx;
    const double roundoff_tol=1.e-12*std::max(1.,std::fabs(sum_yy));
    if(chi2_raw<0. && chi2_raw>=-roundoff_tol) {
      chi2_raw=0.;
    }
    if(!std::isfinite(chi2_scale) || !std::isfinite(chi2_raw) || chi2_raw<0.) {
      error("Non-finite or negative CHI2 encountered after scale optimization");
    }
    chi2_reduced=chi2_raw/static_cast<double>(numq-1);
    if(!getDoScore()) {
      for(unsigned k=0; k<numq; ++k) {
        const double weight=1./(experr[k]*experr[k]);
        chi2_dq[k]=2.*weight*chi2_scale*(chi2_scale*calcint[k]-expint[k]);
      }
    }
    getPntrToComponent("chi2")->set(chi2_raw);
    getPntrToComponent("chi2red")->set(chi2_reduced);
    getPntrToComponent("chi2scale")->set(chi2_scale);
  }

  if(getDoScore()) {
    /* Metainference */
    double score = getScore();
    setScore(score);
  }

  if(chi2ensemble) {
    double ensemble_raw=chi2_raw;
    double ensemble_reduced=chi2_reduced;
    double ensemble_scale=chi2_scale;
    if(chi2ensemble_replicas>1) {
      std::vector<double> ensemble(numq,0.);
      if(comm.Get_rank()==0) {
        for(unsigned k=0; k<numq; ++k) {
          std::string num;
          Tools::convert(k,num);
          ensemble[k]=getPntrToComponent("q-"+num)->get();
        }
        multi_sim_comm.Sum(&ensemble[0],numq);
        for(unsigned k=0; k<numq; ++k) {
          ensemble[k]/=static_cast<double>(chi2ensemble_replicas);
        }
      }
      comm.Sum(&ensemble[0],numq);
      double sum_xy=0.;
      double sum_xx=0.;
      double sum_yy=0.;
      for(unsigned k=0; k<numq; ++k) {
        if(!std::isfinite(ensemble[k])) {
          error("Non-finite ensemble intensity encountered while evaluating CHI2ENSEMBLE");
        }
        const double weight=1./(experr[k]*experr[k]);
        sum_xy += weight*ensemble[k]*expint[k];
        sum_xx += weight*ensemble[k]*ensemble[k];
        sum_yy += weight*expint[k]*expint[k];
      }
      if(!std::isfinite(sum_xx) || sum_xx<=0.) {
        error("Cannot evaluate CHI2ENSEMBLE because the weighted calculated-intensity norm is not positive and finite");
      }
      ensemble_scale=sum_xy/sum_xx;
      ensemble_raw=sum_yy-sum_xy*sum_xy/sum_xx;
      const double roundoff_tol=1.e-12*std::max(1.,std::fabs(sum_yy));
      if(ensemble_raw<0. && ensemble_raw>=-roundoff_tol) {
        ensemble_raw=0.;
      }
      if(!std::isfinite(ensemble_scale) || !std::isfinite(ensemble_raw) || ensemble_raw<0.) {
        error("Non-finite or negative CHI2ENSEMBLE encountered after scale optimization");
      }
      ensemble_reduced=ensemble_raw/static_cast<double>(numq-1);
    }
    getPntrToComponent("chi2ensemble")->set(ensemble_raw);
    getPntrToComponent("chi2ensemblered")->set(ensemble_reduced);
    getPntrToComponent("chi2ensemblescale")->set(ensemble_scale);
  }

  std::vector<Vector> chi2_atom_deriv;
  if(chi2 && !getDoScore()) {
    chi2_atom_deriv.assign(size,Vector(0,0,0));
  }
  Tensor chi2_box_deriv;

  // the derivative with respect to an atom is the derivative with respect to its bead
  // times atom_weight_ (1 without ONEBEAD); the box derivatives are computed from the
  // beads, since sum_atoms x_a (x) w_a d = (sum_atoms w_a x_a) (x) d = x_bead (x) d
  auto weight=[&](const unsigned atom_id) {
    return onebead ? atom_weight_[atom_id] : 1.0;
  };
  auto bead_of=[&](const unsigned atom_id) {
    return onebead ? residue_atom[atom_id] : atom_id;
  };

  // derivatives of each intensity (without DOSCORE), the atoms split over the threads
  if(!getDoScore()) {
    std::vector<Value*> qval(numq);
    for(unsigned k=0; k<numq; ++k) {
      std::string num;
      Tools::convert(k,num);
      qval[k]=getPntrToComponent("q-"+num);
    }
    #pragma omp parallel for num_threads(OpenMP::getNumThreads()) schedule(static)
    for(unsigned atom_id=0; atom_id<size; ++atom_id) {
      const unsigned bead_id = bead_of(atom_id);
      const double w = weight(atom_id);
      Vector chi2_d;
      for(unsigned k=0; k<numq; ++k) {
        const Vector & bd = bd_deriv[k*beads_size+bead_id];
        const Vector qderiv(w*bd[0], w*bd[1], w*bd[2]);
        setAtomsDerivatives(qval[k], atom_id, qderiv);
        if(chi2) {
          chi2_d += chi2_dq[k]*qderiv;
        }
      }
      if(chi2) {
        chi2_atom_deriv[atom_id] = chi2_d;
      }
    }
    std::vector<Tensor> deriv_box(numq);
    #pragma omp parallel for num_threads(OpenMP::getNumThreads()) schedule(static)
    for(unsigned k=0; k<numq; ++k) {
      for(unsigned b=0; b<beads_size; ++b) {
        deriv_box[k] += Tensor(beads_pos[b], bd_deriv[k*beads_size+b]);
      }
    }
    for(unsigned k=0; k<numq; ++k) {
      if(chi2) {
        chi2_box_deriv += chi2_dq[k]*deriv_box[k];
      }
      setBoxDerivatives(qval[k], -deriv_box[k]);
    }
  }

  if(getDoScore()) {
    // derivative of the score with respect to the beads, from all the q values
    std::vector<Vector> score_deriv(beads_size, Vector(0,0,0));
    if(usegpu_) {
      // the derivatives of the intensities are on the GPU: the contraction is done there.
      // With SIGMARES the smeared intensities are I_k = sum_l smear_kl I_l, hence the weights
      // of the unsmeared derivatives are w_l = sum_k smear_kl dS/dI_k
      std::vector<double> w(numq);
      for(unsigned k=0; k<numq; ++k) {
        w[k] = getMetaDer(k);
      }
      if(resolution) {
        const std::vector<double> metader(w);
        for(unsigned l=0; l<numq; ++l) {
          w[l] = 0.;
          for(unsigned k=0; k<numq; ++k) {
            w[l] += smear_[k*numq+l]*metader[k];
          }
        }
      }
      if(has_rows_) {
        backend_->contract(w, score_deriv);
      }
    } else {
      #pragma omp parallel for num_threads(OpenMP::getNumThreads())
      for(unsigned i=0; i<beads_size; ++i) {
        for(unsigned k=0; k<numq; ++k) {
          score_deriv[i] += getMetaDer(k)*bd_deriv[k*beads_size+i];
        }
      }
    }
    // calculate_cpu does not sum the derivatives over the ranks with DOSCORE
    if(reduce_) {
      comm.Sum(&score_deriv[0][0], 3*score_deriv.size());
    }
    Value* val=getPntrToComponent("score");
    #pragma omp parallel for num_threads(OpenMP::getNumThreads()) schedule(static)
    for(unsigned atom_id=0; atom_id<size; ++atom_id) {
      const Vector & sd = score_deriv[bead_of(atom_id)];
      const double w = weight(atom_id);
      setAtomsDerivatives(val, atom_id, Vector(w*sd[0], w*sd[1], w*sd[2]));
    }
    Tensor deriv_box;
    for(unsigned b=0; b<beads_size; ++b) {
      deriv_box += Tensor(beads_pos[b], score_deriv[b]);
    }
    setBoxDerivatives(val, -deriv_box);
  }

  if(chi2 && !getDoScore()) {
    Value* chi2_val=getPntrToComponent("chi2");
    Value* chi2red_val=getPntrToComponent("chi2red");
    const double invdof=1./static_cast<double>(numq-1);
    for(unsigned atom_id=0; atom_id<size; ++atom_id) {
      setAtomsDerivatives(chi2_val, atom_id, chi2_atom_deriv[atom_id]);
      setAtomsDerivatives(chi2red_val, atom_id, invdof*chi2_atom_deriv[atom_id]);
    }
    setBoxDerivatives(chi2_val, -chi2_box_deriv);
    setBoxDerivatives(chi2red_val, -invdof*chi2_box_deriv);
  }
}



// ----------------------------------------------------------------------------
// Monoatomic-ion parameters, shared by the ATOMISTIC and ONEBEAD routes.
//
// X-ray Cromer-Mann coefficients (a1..a4, b1..b4, c) are taken from the
// charged-ion scattering-factor entries tabulated in the International Tables
// for Crystallography / International Tables for X-ray Crystallography, based
// on Cromer & Mann, Acta Cryst. A24, 321-324 (1968).
// For each ion, Sum(a)+c gives the q=0 electron count of the ion
// (atomic number minus charge), as expected for the forward X-ray factor.
//
// Neutron bound coherent scattering lengths are taken from the NIST NCNR table,
// based on Sears, Neutron News 3, No. 3, 1992. These values are nuclear/isotopic
// scattering lengths, so they do not depend on the ion oxidation state in this
// model; for example Fe2+ and Fe3+ use the same Fe neutron value.
//
// The displaced-solvent volume is estimated from Shannon effective ionic radii
// (Shannon, Acta Cryst. A32, 751-767 (1976); here using a fixed coordination
// choice, currently 6-fold) as the volume of a sphere with that radius.
// It is used instead of neutral-atom van der Waals radii because neutral vdW
// radii are too large for small cations and would over-subtract the solvent
// excluded-volume contribution.
// ----------------------------------------------------------------------------
void SAXS::getIonData(unsigned ion, std::vector<double> &cm_a, std::vector<double> &cm_b, double &cm_c, double &vol, double &neutron_b, double &mass) const {
  cm_a.resize(4);
  cm_b.resize(4);
  switch(ion) {
  case ION_NA: // Na1+
    cm_a = {3.2565, 3.9362, 1.3998, 1.0032};
    cm_b = {2.6671, 6.1153, 0.2001, 14.039};
    cm_c = 0.404;
    vol = 4.445;
    neutron_b = 0.363;
    mass = 22.990;
    break;
  case ION_K: // K1+
    cm_a = {7.9578, 7.4917, 6.359, 1.1915};
    cm_b = {12.6331, 0.7674, -0.0020, 31.9128};
    cm_c = -4.9978;
    vol = 11.008;
    neutron_b = 0.367;
    mass = 39.098;
    break;
  case ION_CL: // Cl1-
    cm_a = {18.2915, 7.2084, 6.5337, 2.3386};
    cm_b = {0.0066, 1.1717, 19.5424, 60.4486};
    cm_c = -16.378;
    vol = 24.838;
    neutron_b = 0.9577;
    mass = 35.453;
    break;
  case ION_CA: // Ca2+
    cm_a = {15.6348, 7.9518, 8.4372, 0.8537};
    cm_b = {-0.0074, 0.6089, 10.3116, 25.9905};
    cm_c = -14.875;
    vol = 4.189;
    neutron_b = 0.47;
    mass = 40.078;
    break;
  case ION_MG: // Mg2+
    cm_a = {3.4988, 3.8378, 1.3284, 0.8497};
    cm_b = {2.1676, 4.7542, 0.185, 10.1411};
    cm_c = 0.4853;
    vol = 1.563;
    neutron_b = 0.5375;
    mass = 24.305;
    break;
  case ION_ZN: // Zn2+
    cm_a = {11.9719, 7.3862, 6.4668, 1.394};
    cm_b = {2.9946, 0.2031, 7.0826, 18.0995};
    cm_c = 0.7807;
    vol = 1.697;
    neutron_b = 0.568;
    mass = 65.38;
    break;
  case ION_FE2: // Fe2+
    cm_a = {11.0424, 7.374, 4.1346, 0.4399};
    cm_b = {4.6538, 0.3053, 12.0546, 31.2809};
    cm_c = 1.0097;
    vol = 1.988;
    neutron_b = 0.945;
    mass = 55.845;
    break;
  case ION_FE3: // Fe3+
    cm_a = {11.1764, 7.3863, 3.3948, 0.0724};
    cm_b = {4.6147, 0.3005, 11.6729, 38.5566};
    cm_c = 0.9707;
    vol = 1.124;
    neutron_b = 0.945;
    mass = 55.845;
    break;
  case ION_MN: // Mn2+
    cm_a = {10.8061, 7.362, 3.5268, 0.2184};
    cm_b = {5.2796, 0.3435, 14.343, 41.3235};
    cm_c = 1.0874;
    vol = 2.395;
    neutron_b = -0.373;
    mass = 54.938;
    break;
  default:
    plumed_merror("getIonData called with an unknown ion id");
  }
}

// Detect a monoatomic ion from its residue name. Detection is by RESIDUE name only.
// The protein/nucleic-acid branches are matched before this is reached.
// Accepted aliases follow the RCSB PDB het codes plus the AMBER and CHARMM ion residue names.
int SAXS::ionTypeFromResname(const std::string &resname) const {
  std::string r = resname;
  r.erase(std::remove_if(r.begin(), r.end(), ::isspace), r.end());
  std::transform(r.begin(), r.end(), r.begin(), ::toupper);
  if(r=="NA"  || r=="SOD" || r=="NA+") {
    return ION_NA;
  }
  if(r=="K"   || r=="POT" || r=="K+") {
    return ION_K;
  }
  if(r=="CL"  || r=="CLA" || r=="CL-") {
    return ION_CL;
  }
  if(r=="CA"  || r=="CAL" || r=="CA2+") {
    return ION_CA;
  }
  if(r=="MG"  || r=="MG2+") {
    return ION_MG;
  }
  if(r=="ZN"  || r=="ZN2" || r=="ZN2+") {
    return ION_ZN;
  }
  if(r=="FE2") {
    return ION_FE2;
  }
  if(r=="FE3") {
    return ION_FE3;
  }
  if(r=="MN"  || r=="MN2+") {
    return ION_MN;
  }
  return -1;
}


double SAXS::calculateAFF(const std::vector<AtomNumber> &atoms, std::vector<std::vector<long double> > &FF_tmp, const double Rho) {
  std::map<std::string, unsigned> AA_map;
  AA_map["H"] = H;
  AA_map["C"] = C;
  AA_map["N"] = N;
  AA_map["O"] = O;
  AA_map["P"] = P;
  AA_map["S"] = S;

  std::vector<std::vector<double> > param_a;
  std::vector<std::vector<double> > param_b;
  std::vector<double> param_c;
  std::vector<double> param_v;

  param_a.resize(NTT, std::vector<double>(5));
  param_b.resize(NTT, std::vector<double>(5));
  param_c.resize(NTT);
  param_v.resize(NTT);

  param_a[H][0] = 0.493002;
  param_b[H][0] = 10.5109;
  param_c[H] = 0.003038;
  param_a[H][1] = 0.322912;
  param_b[H][1] = 26.1257;
  param_v[H] = 5.15;
  param_a[H][2] = 0.140191;
  param_b[H][2] = 3.14236;
  param_a[H][3] = 0.040810;
  param_b[H][3] = 57.7997;
  param_a[H][4] = 0.0;
  param_b[H][4] = 1.0;

  param_a[C][0] = 2.31000;
  param_b[C][0] = 20.8439;
  param_c[C] = 0.215600;
  param_a[C][1] = 1.02000;
  param_b[C][1] = 10.2075;
  param_v[C] = 16.44;
  param_a[C][2] = 1.58860;
  param_b[C][2] = 0.56870;
  param_a[C][3] = 0.86500;
  param_b[C][3] = 51.6512;
  param_a[C][4] = 0.0;
  param_b[C][4] = 1.0;

  param_a[N][0] = 12.2126;
  param_b[N][0] = 0.00570;
  param_c[N] = -11.529;
  param_a[N][1] = 3.13220;
  param_b[N][1] = 9.89330;
  param_v[N] = 2.49;
  param_a[N][2] = 2.01250;
  param_b[N][2] = 28.9975;
  param_a[N][3] = 1.16630;
  param_b[N][3] = 0.58260;
  param_a[N][4] = 0.0;
  param_b[N][4] = 1.0;

  param_a[O][0] = 3.04850;
  param_b[O][0] = 13.2771;
  param_c[O] = 0.250800 ;
  param_a[O][1] = 2.28680;
  param_b[O][1] = 5.70110;
  param_v[O] = 9.13;
  param_a[O][2] = 1.54630;
  param_b[O][2] = 0.32390;
  param_a[O][3] = 0.86700;
  param_b[O][3] = 32.9089;
  param_a[O][4] = 0.0;
  param_b[O][4] = 1.0;

  param_a[P][0] = 6.43450;
  param_b[P][0] = 1.90670;
  param_c[P] = 1.11490;
  param_a[P][1] = 4.17910;
  param_b[P][1] = 27.1570;
  param_v[P] = 5.73;
  param_a[P][2] = 1.78000;
  param_b[P][2] = 0.52600;
  param_a[P][3] = 1.49080;
  param_b[P][3] = 68.1645;
  param_a[P][4] = 0.0;
  param_b[P][4] = 1.0;

  param_a[S][0] = 6.90530;
  param_b[S][0] = 1.46790;
  param_c[S] = 0.866900;
  param_a[S][1] = 5.20340;
  param_b[S][1] = 22.2151;
  param_v[S] = 19.86;
  param_a[S][2] = 1.43790;
  param_b[S][2] = 0.25360;
  param_a[S][3] = 1.58630;
  param_b[S][3] = 56.1720;
  param_a[S][4] = 0.0;
  param_b[S][4] = 1.0;

  // monoatomic-ion X-ray Cromer-Mann parameters (charged-ion rows; see getIonData)
  for(unsigned t=0; t<NIONS; ++t) {
    std::vector<double> a, b;
    double c, v, bN, m;
    getIonData(t, a, b, c, v, bN, m);
    const unsigned idx = affSlot(t);
    for(unsigned j=0; j<4; ++j) {
      param_a[idx][j] = a[j];
      param_b[idx][j] = b[j];
    }
    param_c[idx] = c;
    param_v[idx] = v;
  }

  auto* moldat=plumed.getActionSet().selectLatest<GenericMolInfo*>(this);

  double Iq_0=0.;
  if( moldat ) {
    // cycle over the atom types
    for(unsigned i=0; i<NTT; ++i) {
      const double volr = std::pow(param_v[i], (2.0/3.0)) /(4. * pi);
      // cycle over q
      for(unsigned k=0; k<q_list.size(); ++k) {
        const double q = q_list[k];
        const double s = q / (4. * pi);
        FF_tmp[k][i] = param_c[i];
        // SUM [a_i * EXP( - b_i * (q/4pi)^2 )] Waasmaier and Kirfel (1995)
        for(unsigned j=0; j<4; ++j) {
          FF_tmp[k][i] += param_a[i][j]*std::exp(-param_b[i][j]*s*s);
        }
        // subtract solvation: rho * v_i * EXP( (- v_i^(2/3) / (4pi)) * q^2  ) // since  D in Fraser 1978 is 2*s
        FF_tmp[k][i] -= Rho*param_v[i]*std::exp(-volr*q*q);
      }
    }
    // cycle over the atoms to assign the atom type and calculate I0
    for(unsigned i=0; i<atoms.size(); ++i) {
      // ions are detected by residue name (never by the leading element letter,
      // which would misread CL/CA as carbon, NA as nitrogen, ...)
      int ion = ionTypeFromResname(moldat->getResidueName(atoms[i]));
      if(ion >= 0) {
        const unsigned index = affSlot(ion);
        atoi[i] = index;
        for(unsigned j=0; j<4; ++j) {
          Iq_0 += param_a[index][j];
        }
        Iq_0 -= Rho*param_v[index] - param_c[index];
        continue;
      }
      // get atom name
      std::string name = moldat->getAtomName(atoms[i]);
      char type;
      // get atom type
      char first = name.at(0);
      // GOLDEN RULE: type is first letter, if not a number
      if (!isdigit(first)) {
        type = first;
        // otherwise is the second
      } else {
        type = name.at(1);
      }
      std::string type_s = std::string(1,type);
      if(AA_map.find(type_s) != AA_map.end()) {
        const unsigned index=AA_map[type_s];
        atoi[i] = AA_map[type_s];
        for(unsigned j=0; j<4; ++j) {
          Iq_0 += param_a[index][j];
        }
        Iq_0 -= Rho*param_v[index] - param_c[index];
      } else {
        error("Wrong atom type "+type_s+" from atom name "+name+"\n");
      }
    }
  } else {
    error("MOLINFO DATA not found\n");
  }
  if(absolute) {
    Iq_0 = 1;
  }

  return Iq_0;
}

double SAXS::calculateAFFsans(const std::vector<AtomNumber> &atoms, std::vector<std::vector<long double> > &FF_tmp, const double deuterConc) {
  std::map<std::string, unsigned> AA_map;
  AA_map["H"] = H;
  AA_map["C"] = C;
  AA_map["N"] = N;
  AA_map["O"] = O;
  AA_map["P"] = P;
  AA_map["S"] = S;

  std::vector<double> param_b;
  std::vector<double> param_v;

  param_b.resize(NTT);
  param_v.resize(NTT);

  param_b[H] = -0.374;
  param_v[H] = 5.15;
  // param_b[D] = 0.667;
  param_b[C] =  0.665;
  param_v[C] = 16.44;
  param_b[N] =  0.94;
  param_v[N] = 2.49;
  param_b[O] =  0.580;
  param_v[O] = 9.13;
  param_b[P] =  0.51;
  param_v[P] = 5.73;
  param_b[S] =  0.28;
  param_v[S] = 19.86;

  // monoatomic-ion neutron scattering lengths (nuclear, charge-independent) and volumes
  for(unsigned t=0; t<NIONS; ++t) {
    std::vector<double> a, b;
    double c, v, bN, m;
    getIonData(t, a, b, c, v, bN, m);
    const unsigned idx = affSlot(t);
    param_b[idx] = bN;
    param_v[idx] = v;
  }

  double solv_sc_length = 0.1*(param_b[O] + 2.*((1. - deuterConc) * param_b[H] + deuterConc * 0.667)); // per water electron (10 electrons)

  auto* moldat=plumed.getActionSet().selectLatest<GenericMolInfo*>(this);

  double Iq_0=0.;
  if( moldat ) {
    // cycle over the atom types
    for(unsigned i=0; i<NTT; ++i) {
      double volr = std::pow(param_v[i], (2.0/3.0)) /(4. * pi);
      // cycle over q
      for(unsigned k=0; k<q_list.size(); ++k) {
        const double q = q_list[k];
        FF_tmp[k][i] = param_b[i];
        // subtract solvation: rho * v_i * EXP( (- v_i^(2/3) / (4pi)) * q^2  ) // since  D in Fraser 1978 is 2*s
        FF_tmp[k][i] -= solv_sc_length*rho*param_v[i]*std::exp(-volr*q*q);
      }
    }
    // cycle over the atoms to assign the atom type and calculate I0
    for(unsigned i=0; i<atoms.size(); ++i) {
      // ions are detected by residue name (never by the leading element letter)
      int ion = ionTypeFromResname(moldat->getResidueName(atoms[i]));
      if(ion >= 0) {
        const unsigned index = affSlot(ion);
        atoi[i] = index;
        Iq_0 += param_b[index]-solv_sc_length*rho*param_v[index];
        continue;
      }
      // get atom name
      std::string name = moldat->getAtomName(atoms[i]);
      char type;
      // get atom type
      char first = name.at(0);
      // GOLDEN RULE: type is first letter, if not a number
      if (!isdigit(first)) {
        type = first;
        // otherwise is the second
      } else {
        type = name.at(1);
      }
      std::string type_s = std::string(1,type);
      if(AA_map.find(type_s) != AA_map.end()) {
        const unsigned index=AA_map[type_s];
        atoi[i] = AA_map[type_s];
        Iq_0 += param_b[index]-solv_sc_length*rho*param_v[index];
      } else {
        error("Wrong atom type "+type_s+" from atom name "+name+"\n");
      }
    }
  } else {
    error("MOLINFO DATA not found\n");
  }
  if(absolute) {
    Iq_0 = 1;
  }

  return Iq_0;
}


void SAXS::resolution_function() {
  const unsigned numq = q_list.size();

  // only OpenMP because numq might be smaller than the number of ranks
  #pragma omp parallel for num_threads(OpenMP::getNumThreads())
  for (unsigned i=0; i<numq; i++) {
    double qi = q_list[i];
    double dq = 6*sigma_res[i]/(Nj-1);
    double sigma_sq = sigma_res[i]*sigma_res[i];
    double qstart = qi - 3*sigma_res[i];
    for (unsigned j=0; j<Nj; j++) {
      double qj = qstart + j*dq;
      double I0exp = i0e(qj*qi/sigma_sq);

      qj_list[i][j] = qj;
      Rij[i][j] = (qj/sigma_sq)*std::exp(-0.5*(qj - qi)*(qj - qi)/sigma_sq)*I0exp;
    }
  }
}

// i0e function from cephes
// compute I0(x) * exp (-x), with I0 being the modified Bessel function
// of first kind and zeroth order.
inline double SAXS::i0e(double x) {
  double y = 0.0;

  if (x < 0) {
    x = -x;
  }
  if (x <= 8.0) {
    y = (x/2.0) - 2.0;
    return chbevl(y, A);
  }

  return chbevl(32.0/x - 2.0, B) / sqrt(x);
}

double SAXS::chbevl(double x, const std::vector<double> &coeffs) {
  double b0, b1, b2;
  unsigned n = coeffs.size();

  b0 = coeffs[0];
  b1 = 0.0;
  b2 = 0.0;

  for (unsigned i = 1; i < n; i++) {
    b2 = b1;
    b1 = b0;
    b0 = x * b1 - b2 + coeffs[i];
  }
  return 0.5 * (b0 - b2);
}

inline double SAXS::interpolation(std::vector<SplineCoeffs> &coeffs, double x) {
  unsigned s = 0;
  while ((x >= q_list[s+1]) && (s+1 < q_list.size()-1)) {
    s++;
  }

  double dx = x - coeffs[s].x;
  return coeffs[s].a + coeffs[s].b*dx + coeffs[s].c*dx*dx + coeffs[s].d*dx*dx*dx;
}

// natural bc cubic spline implementation from the Wikipedia algorithm
// modified from https://stackoverflow.com/a/19216702/3254658
std::vector<SAXS::SplineCoeffs> SAXS::spline_coeffs(std::vector<double> &x, std::vector<double> &y) {
  unsigned n = x.size()-1;
  std::vector<double> a;
  a.insert(a.begin(), y.begin(), y.end());
  std::vector<double> b(n);
  std::vector<double> d(n);
  std::vector<double> h;

  for(unsigned i=0; i<n; i++) {
    h.push_back(x[i+1]-x[i]);
  }

  std::vector<double> alpha;
  alpha.push_back(0);
  for(unsigned i=1; i<n; i++) {
    alpha.push_back( 3*(a[i+1]-a[i])/h[i] - 3*(a[i]-a[i-1])/h[i-1]  );
  }

  std::vector<double> c(n+1);
  std::vector<double> l(n+1);
  std::vector<double> mu(n+1);
  std::vector<double> z(n+1);
  l[0] = 1;
  mu[0] = 0;
  z[0] = 0;

  for(unsigned i=1; i<n; i++) {
    l[i] = 2 *(x[i+1]-x[i-1])-h[i-1]*mu[i-1];
    mu[i] = h[i]/l[i];
    z[i] = (alpha[i]-h[i-1]*z[i-1])/l[i];
  }

  l[n] = 1;
  z[n] = 0;
  c[n] = 0;

  for(int j=n-1; j>=0; j--) {
    c[j] = z[j] - mu[j] * c[j+1];
    b[j] = (a[j+1]-a[j])/h[j]-h[j]*(c[j+1]+2*c[j])/3;
    d[j] = (c[j+1]-c[j])/3/h[j];
  }

  std::vector<SplineCoeffs> output_set(n);
  for(unsigned i=0; i<n; i++) {
    output_set[i].a = a[i];
    output_set[i].b = b[i];
    output_set[i].c = c[i];
    output_set[i].d = d[i];
    output_set[i].x = x[i];
  }

  return output_set;
}


}//namespace isdb
}//namespace PLMD
