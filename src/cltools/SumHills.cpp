/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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
#include "core/CLToolRegister.h"
#include "tools/Tools.h"
#include "core/Action.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "tools/Communicator.h"
#include "tools/Random.h"
#include <algorithm>
#include <cstdio>
#include <string>
#include <vector>
#include "tools/File.h"
#include "core/Value.h"
#include "tools/Matrix.h"

namespace PLMD {
namespace cltools {

//+PLUMEDOC TOOLS sum_hills
/*
sum_hills is a tool that sums the Gaussians in a hills file and writes the resulting free energy on a grid.

This tool is most frequently used after a [METAD](METAD.md) simulation to sum all the Gaussians in the hills file
and output the free energy surface:

```plumed
plumed sum_hills --hills PATHTOMYHILLSFILE
```

The default name for the output file is `fes.dat`. It contains, for each point of the grid, the free energy
and its derivatives with respect to the variables. The free energy is minus the sum of the Gaussians: with
well-tempered metadynamics the heights written in the hills file are already scaled by $\gamma/(\gamma-1)$
(where $\gamma$ is the bias factor), so that their sum is directly the estimate of minus the free energy.
As always, the free energy is defined up to an additive constant; use `--mintozero` to shift its minimum to zero.

sum_hills automatically detects the variables in the hills file and their periodicity. A column `X` of the file
is a variable if the file also contains its width (`sigma_X`, or `sigma_X_X` for multivariate Gaussians),
whatever the name of the variable. If you used flexible hills (multivariate Gaussian kernels), sum_hills
understands it from the hills file.

## Choosing the grid

By default the grid is chosen so as to contain all the Gaussians, with a spacing of about a third of the smallest width.
You can instead choose the number of bins:

```plumed
plumed sum_hills --bin 99,99 --hills PATHTOMYHILLSFILE
```

In the command above the minimum and maximum of the grid are still detected for you. If you want to set them yourself
you can do so as shown below:

```plumed
plumed sum_hills --bin 99,99 --min -pi,-pi --max pi,pi --hills PATHTOMYHILLSFILE
```

You can of course use numbers instead of -pi/pi. `--min` and `--max` must be given together. Instead of the number of bins
you can give the grid spacing with `--spacing`, with or without `--min` and `--max`. If you give both `--bin` and `--spacing`,
the larger number of bins is used for each variable.

## Multiple files and multiple walkers

sum_hills can also accept multiple files:

```plumed
plumed sum_hills --hills PATHTOMYHILLSFILE1,PATHTOMYHILLSFILE2,PATHTOMYHILLSFILE3
```

The files are integrated one after the other. This does not change the final free energy, but when
you use `--stride` (see below) with the hills files of multiple walkers the intermediate free energies would
not follow the simulation. With `--sort time` the hills of all the files are merged in order of their time
(hills with the same time are taken in the order of the files), so that each intermediate free energy
contains the hills deposited by all the walkers up to a given time:

```plumed
plumed sum_hills --hills HILLS.0,HILLS.1,HILLS.2 --sort time --stride 300
```

With `--sort clock` the hills are merged in order of the `clock` field (the wall-clock time at which they were
deposited, written with multiple walkers), which is more appropriate when the walkers run at different speeds.
Each file should be ordered by the chosen field.

## Integrating out variables

If you want to integrate out some variables from your free energy surface you can use `--idw`:

```plumed
plumed sum_hills --hills PATHTOMYHILLSFILE --idw t1 --kt 0.6
```

With `--idw` you list the variables that you want in the output free energy surface.
All the other variables in the hills file are integrated out, i.e.
$F(s_1) = -k_BT \log \sum_{s_2} e^{-F(s_1,s_2)/k_BT}$, where the sum runs over the points of the grid.
`--kt` is the temperature of the system in energy units: be consistent with the units of your hills file,
since plumed will not check this for you. If you need more variables you can list them separated by commas:

```plumed
plumed sum_hills --hills PATHTOMYHILLSFILE --idw t1,t2 --kt 0.6
```

## Following the convergence

With `--stride` you can write the free energy after each bunch of hills (the stride is a number of hills):

```plumed
plumed sum_hills --stride 300 --hills PATHTOMYHILLSFILE
```

The output files are then called `fes_0.dat`, `fes_1.dat`, etc. (see below), and the last one contains all the hills.

With `--nohistory` each output contains only the hills of its own block of hills, instead of all the hills up to that point:

```plumed
plumed sum_hills --stride 300 --hills PATHTOMYHILLSFILE --nohistory
```

In this way you can monitor the bias deposited in blocks of time as a metadynamics simulation goes on:
when the bias deposited in a block is constant over the whole domain, the simulation can be considered to be at convergence.

To check the convergence you can also follow the free-energy difference between basins as the
simulation goes on. With `--basins` you can define two or more basins as boxes, each given by the minimum and
the maximum of each variable of the output, with the basins separated by colons. Along a periodic variable a
basin can cross the boundary if its minimum is larger than its maximum. For example

```plumed
plumed sum_hills --hills PATHTOMYHILLSFILE --stride 500 --kt 2.5 --basins -pi,0,-pi,pi:0,pi,-pi,pi
```

writes, at each output, the time of the last hill read, the number of hills, and the free energy of each basin
relative to the first one, computed as $-k_BT \log \sum_{s\in\mathrm{basin}} e^{-F(s)/k_BT}$ over the points of
the grid, in the file `basins.dat` (the name can be set with `--outbasins`). With `--idw` the basins are defined on
the projected variables. With `--histo` (see below) the same is done for the free energy from the histogram, in the file `histo_basins.dat`.
A plateau of these differences is a necessary (not sufficient) indication of convergence.

## Time average of the bias

With standard (not well-tempered) metadynamics the bias does not converge but fluctuates around the free energy.
A better estimate is the time average of the bias after an initial transient, which you can obtain with `--average`:

```plumed
plumed sum_hills --hills PATHTOMYHILLSFILE --average 1000
```

computes the average of the bias between time 1000 (in the units of the time column of the hills file) and the time of
the last hill. This is done by scaling the height of each hill deposited at time $t$ after 1000 by $(t_1-t)/(t_1-1000)$, where
$t_1$ is the final time. You can also give the final time explicitly, e.g. `--average 1000,5000`: hills deposited after
5000 are ignored. The same average is useful for well-tempered metadynamics, where it reduces the fluctuations of the
estimate of the free energy. Combined with `--stride`, each output is the average from the initial time up to the last hill read.
This is the same average computed by the function `hills_time_average` of the plumed python module.

## Bias instead of free energy

If you have run well-tempered metadynamics you can also ask sum_hills to output the negative bias, $-V(s)=(\gamma-1)/\gamma\, F(s)$,
instead of the free energy by using `--negbias`:

```plumed
plumed sum_hills --negbias --hills PATHTOMYHILLSFILE
```

Here the default output file name is `negativebias.dat`. If the simulation was not well tempered the bias is minus the free energy,
and the output is the same as without `--negbias`.

## Histograms

From time to time you might need to use a hills or a colvar file as a set of points from which you want to build
a free energy as $-k_BT \log P(s)$. You can do this with `--histo`:

```plumed
plumed sum_hills --histo PATHTOMYCOLVARORHILLSFILE --sigma 0.2,0.2 --kt 0.6
```

The histogram is computed with Gaussian kernels whose widths are given with `--sigma`, one for each variable,
so it is a continuous histogram, and `--kt` is needed to convert it into a free energy. The default output file is `histo.dat`.
With a hills file the variables are those of the hills; with a colvar file all the columns but time are variables, so you
will usually select some of them with `--idw` (and give one width for each of them).
Multiple input files can be given, separated by commas.

All the points have the same weight: the result is the distribution sampled by the simulation, which is the
free energy only if the simulation was not biased. **sum_hills does not reweight the points.** To compute a free
energy from a biased simulation (e.g. from the COLVAR file of a metadynamics simulation) use [driver](driver.md) with
[READ](READ.md), [REWEIGHT_BIAS](REWEIGHT_BIAS.md) (or [REWEIGHT_METAD](REWEIGHT_METAD.md)), [HISTOGRAM](HISTOGRAM.md)
and [CONVERT_TO_FES](CONVERT_TO_FES.md), for example

```plumed
#SETTINGS INPUTFILES=regtest/basic/rt-fametad-2-restart/COLVAR
s: READ FILE=regtest/basic/rt-fametad-2-restart/COLVAR VALUES=sum_abs IGNORE_TIME IGNORE_FORCES
metad: READ FILE=regtest/basic/rt-fametad-2-restart/COLVAR VALUES=MetaD.bias IGNORE_TIME IGNORE_FORCES
w: REWEIGHT_BIAS TEMP=300 ARG=metad.bias
hh: HISTOGRAM ARG=s GRID_MIN=0 GRID_MAX=1 GRID_BIN=100 BANDWIDTH=0.02 LOGWEIGHTS=w NORMALIZATION=true
ff: CONVERT_TO_FES ARG=hh TEMP=300
DUMPGRID ARG=ff FILE=fes.dat
```

run with `plumed driver --plumed plumed.dat --noatoms`. With well-tempered metadynamics the bias changes with time:
use the `rbias` component of [METAD](METAD.md) (with `CALC_RCT`) instead of the bias, or [REWEIGHT_METAD](REWEIGHT_METAD.md).
Block averages, which give the statistical error of the free energy, are discussed in the documentation of [HISTOGRAM](HISTOGRAM.md)
and in the tutorial linked there.

You can compute the free energy from the hills and the histogram at the same time, e.g. from the same file:

```plumed
plumed sum_hills --hills PATHTOMYHILLSFILE --histo PATHTOMYHILLSFILE --sigma 0.2,0.2 --kt 0.6
```

`--nohistory` and `--stride` can also be used with `--histo`, to compute the histogram in blocks of points:

```plumed
plumed sum_hills --histo PATHTOMYCOLVARORHILLSFILE --sigma 0.2,0.2 --kt 0.6 --stride 1000 --nohistory
```

## Output

The output format can be controlled by using the `--fmt` option as shown below

```plumed
plumed sum_hills --hills PATHTOMYHILLSFILE --fmt %8.3f
```

where here we chose a float with length of 8 and 3 digits.

The output file can be named in an arbitrary way by using the `--outfile` option as shown below:

```plumed
plumed sum_hills --hills PATHTOMYHILLSFILE --outfile myfes.dat
```

This command produces a file `myfes.dat` which contains the free energy surface.

If you use `--stride`, the name specified with `--outfile` is the prefix of the files, so this command:

```plumed
plumed sum_hills --hills PATHTOMYHILLSFILE --outfile myfes_ --stride 100
```

produces output files called `myfes_0.dat`, `myfes_1.dat`, `myfes_2.dat` etc.

The same is true for the output of the histogram, which is named with `--outhisto`. So this example

```plumed
plumed sum_hills --histo HILLS --kt 2.5 --sigma 0.01 --outhisto myhisto.dat
```

produces a file `myhisto.dat`, while this input

```plumed
plumed sum_hills --histo HILLS --kt 2.5 --sigma 0.01 --outhisto myhisto_ --stride 100
```

produces output files called `myhisto_0.dat`, `myhisto_1.dat`, `myhisto_2.dat` etc.

*/
//+ENDPLUMEDOC

class CLToolSumHills : public CLTool {
public:
  static void registerKeywords( Keywords& keys );
  explicit CLToolSumHills(const CLToolOptions& co );
  int main(FILE* in,FILE*out,Communicator& pc) override;
  std::string description()const override;
/// find a list of variables present, if they are periodic and which is the period
/// return false if the file does not exist
  static bool findCvsAndPeriodic(const std::string & filename, std::vector< std::vector <std::string> > &cvs,std::vector<std::string> &pmin,std::vector<std::string> &pmax, bool &multivariate, std::string &lowI_, std::string &uppI_);
};

void CLToolSumHills::registerKeywords( Keywords& keys ) {
  CLTool::registerKeywords( keys );
  keys.addFlag("--help-debug",false,"print special options that can be used to create regtests");
  keys.add("optional","--hills","the hills files to sum, separated by commas");
  keys.add("optional","--histo","the files (hills or colvar) from which a histogram is computed, separated by commas; the points are not reweighted");
  keys.add("optional","--stride","write the free energy (or histogram) every this number of hills (or points); by default only at the end");
  keys.add("optional","--min","the lower bounds of the grid, one for each variable; by default they are detected from the hills");
  keys.add("optional","--max","the upper bounds of the grid, one for each variable; by default they are detected from the hills");
  keys.add("optional","--bin","the number of bins of the grid, one for each variable");
  keys.add("optional","--spacing","the grid spacing, one for each variable, alternative to the number of bins");
  keys.add("optional","--idw","the variables of the output (default all). With --hills the other variables are integrated out, with --histo the other variables are not considered");
  keys.add("optional","--outfile","the output file for the free energy (fes.dat, or negativebias.dat with --negbias); with --stride the prefix of the files");
  keys.add("optional","--outhisto","the output file for the histogram (histo.dat); with --stride the prefix of the files");
  keys.add("optional","--kt","the temperature in energy units, needed with --idw, --histo and --basins");
  keys.add("optional","--sigma","the widths of the Gaussian kernels of the histogram, one for each variable (needed with --histo)");
  keys.addFlag("--negbias",false,"write minus the bias instead of the free energy (they differ only with well-tempered metadynamics)");
  keys.addFlag("--nohistory",false,"with --stride, each output contains only the hills (or points) of its own block");
  keys.addFlag("--mintozero",false,"shift the free energy (or histogram) so that its minimum is zero");
  keys.add("optional","--fmt","the format of the numbers in the output, e.g. %14.9f");
  keys.add("optional","--sort","time or clock: with multiple files (e.g. one per walker), merge the hills in order of time or clock instead of reading the files one after the other");
  keys.add("optional","--average","time average of the bias from time t0 (t0) or between times t0 and t1 (t0,t1); without t1 the average ends at the last hill read");
  keys.add("optional","--basins","basins separated by colons, each given by the minimum and maximum of each output variable (e.g. -3,-1:1,3); the free energy of each basin relative to the first one is written at each output");
  keys.add("optional","--outbasins","the output file for the free energies of the basins (basins.dat)");
}

CLToolSumHills::CLToolSumHills(const CLToolOptions& co ):
  CLTool(co) {
  inputdata=inputType::commandline;
}

std::string CLToolSumHills::description()const {
  return "sum the hills of metadynamics to compute the free energy";
}

int CLToolSumHills::main(FILE* in,FILE*out,Communicator& pc) {

// Read the hills input file name
  std::vector<std::string> hillsFiles;
  bool dohills;
  dohills=parseVector("--hills",hillsFiles);
// Read the histogram file
  std::vector<std::string> histoFiles;
  bool dohisto;
  dohisto=parseVector("--histo",histoFiles);

  plumed_massert(dohisto || dohills,"you should use --histo or/and --hills command");

  std::vector< std::vector<std::string> > vcvs;
  std::vector<std::string> vpmin;
  std::vector<std::string> vpmax;
  std::string lowI_, uppI_;
  if(dohills) {
    // parse it as it was a restart
    bool vmultivariate;
    if(!findCvsAndPeriodic(hillsFiles[0], vcvs, vpmin, vpmax, vmultivariate, lowI_, uppI_)) {
      plumed_merror("cannot find file "+hillsFiles[0]);
    }
    if(vcvs.empty()) {
      plumed_merror("no collective variable found in file "+hillsFiles[0]);
    }
  }

  std::vector< std::vector<std::string> > hcvs;
  std::vector<std::string> hpmin;
  std::vector<std::string> hpmax;

  std::vector<std::string> sigma;
  if(dohisto) {
    bool hmultivariate;
    if(!findCvsAndPeriodic(histoFiles[0], hcvs, hpmin, hpmax, hmultivariate, lowI_, uppI_)) {
      plumed_merror("cannot find file "+histoFiles[0]);
    }
    if(hcvs.empty()) {
      plumed_merror("no collective variable found in file "+histoFiles[0]);
    }
    // here need also the vector of sigmas
    parseVector("--sigma",sigma);
    if(sigma.size()==0) {
      plumed_merror("you should define --sigma vector when using histogram");
    }
    lowI_=uppI_="-1.";  // Interval is not use for histograms
  }

  if(dohisto && dohills) {
    plumed_massert(vcvs==hcvs,"variables for histogram and bias should have the same labels");
    plumed_massert(hpmin==vpmin,"variables for histogram and bias should have the same min for periodicity");
    plumed_massert(hpmax==vpmax,"variables for histogram and bias should have the same max for periodicity");
  }

  // now put into a neutral vector

  std::vector< std::vector<std::string> > cvs;
  std::vector<std::string> pmin;
  std::vector<std::string> pmax;

  if(dohills) {
    cvs=vcvs;
    pmin=vpmin;
    pmax=vpmax;
  }
  if(dohisto) {
    cvs=hcvs;
    pmin=hpmin;
    pmax=hpmax;
  }


  // setup grids
  unsigned grid_check=0;
  std::vector<std::string> gmin(cvs.size());
  if(parseVector("--min",gmin)) {
    if(gmin.size()!=cvs.size() && gmin.size()!=0) {
      plumed_merror("not enough values for --min");
    }
    grid_check++;
  }
  std::vector<std::string> gmax(cvs.size() );
  if(parseVector("--max",gmax)) {
    if(gmax.size()!=cvs.size() && gmax.size()!=0) {
      plumed_merror("not enough values for --max");
    }
    grid_check++;
  }
  std::vector<std::string> gbin(cvs.size());
  bool grid_has_bin;
  grid_has_bin=false;
  if(parseVector("--bin",gbin)) {
    if(gbin.size()!=cvs.size() && gbin.size()!=0) {
      plumed_merror("not enough values for --bin");
    }
    grid_has_bin=true;
  }
  std::vector<std::string> gspacing(cvs.size());
  bool grid_has_spacing;
  grid_has_spacing=false;
  if(parseVector("--spacing",gspacing)) {
    if(gspacing.size()!=cvs.size() && gspacing.size()!=0) {
      plumed_merror("not enough values for --spacing");
    }
    grid_has_spacing=true;
  }
  // allowed: no grids only bin
  // not allowed: partial grid definition
  plumed_massert( gmin.size()==gmax.size() && (gmin.size()==0 ||  gmin.size()==cvs.size() ),"you should specify --min and --max together with same number of components");



  PlumedMain plumed;
  std::string ss;
  unsigned nn=1;
  ss="setNatoms";
  plumed.cmd(ss,&nn);
  if(Communicator::initialized()) {
    plumed.cmd("setMPIComm",&pc.Get_comm());
  }
  plumed.cmd("init",&nn);
  std::vector <bool> isdone(cvs.size(),false);
  for(unsigned i=0; i<cvs.size(); i++) {
    if(!isdone[i]) {
      isdone[i]=true;
      std::vector<std::string> actioninput;
      std::vector <unsigned> inds;
      actioninput.push_back("FAKE");
      actioninput.push_back("ATOMS=1");
      actioninput.push_back("LABEL="+cvs[i][0]);
      std::vector<std::string> comps, periods;
      if(cvs[i].size()>1) {
        comps.push_back(cvs[i][1]);
        inds.push_back(i);
      }
      periods.push_back(pmin[i]);
      periods.push_back(pmax[i]);
      for(unsigned j=i+1; j<cvs.size(); j++) {
        if(cvs[i][0]==cvs[j][0] && !isdone[j]) {
          if(cvs[i].size()==1 || cvs[j].size()==1  ) {
            plumed_merror("you cannot have twice the same label and no components ");
          }
          if(cvs[j].size()>1) {
            comps.push_back(cvs[j][1]);
            periods.push_back(pmin[j]);
            periods.push_back(pmax[j]);
            isdone[j]=true;
            inds.push_back(j);
          }
        }

      }
      // drain all the components
      std::string addme;
      if(comps.size()>0) {
        addme="COMPONENTS=";
        for(unsigned ii=0; ii<comps.size()-1; ii++) {
          addme+=comps[ii]+",";
        }
        addme+=comps.back();
        actioninput.push_back(addme);
      }
      // periodicity (always explicit here)
      addme="PERIODIC=";
      for(unsigned j=0; j<periods.size()-1; j++) {
        addme+=periods[j]+",";
      }
      addme+=periods.back();
      actioninput.push_back(addme);
      for(unsigned j=0; j<inds.size(); j++) {
        unsigned jj;
        jj=inds[j];
        if(grid_check==2) {
          double gm;
          double pm;
          if(pmin[jj]!="none") {
            Tools::convert(gmin[jj],gm);
            Tools::convert(pmin[jj],pm);
            if(  gm<pm  ) {
              plumed_merror("Periodicity issue : GRID_MIN value ( "+gmin[jj]+" ) is less than periodicity in HILLS file in "+cvs[jj][0]+ " ( "+pmin[jj]+" ) ");
            }
          }
          if(pmax[jj]!="none") {
            Tools::convert(gmax[jj],gm);
            Tools::convert(pmax[jj],pm);
            if(  gm>pm ) {
              plumed_merror("Periodicity issue : GRID_MAX value ( "+gmax[jj]+" ) is more than periodicity in HILLS file in "+cvs[jj][0]+ " ( "+pmax[jj]+" ) ");
            }
          }
        }
      }

      plumed.readInputWords(actioninput,false);
    }

  }
  unsigned ncv=cvs.size();
  std::vector<std::string> actioninput;
  std::vector<std::string> idw;
  // check if the variables to be used are correct
  if(parseVector("--idw",idw)) {
    for(unsigned i=0; i<idw.size(); i++) {
      bool found=false;
      for(unsigned j=0; j<cvs.size(); j++) {
        if(cvs[j].size()>1) {
          if(idw[i]==cvs[j][0]+"."+cvs[j][1]) {
            found=true;
          }
        } else {
          if(idw[i]==cvs[j][0]) {
            found=true;
          }
        }
      }
      if(!found) {
        plumed_merror("variable "+idw[i]+" is not found in the bunch of cvs: revise your --idw option" );
      }
    }
    plumed_massert( idw.size()<=cvs.size(),"the number of variables to be integrated should be at most equal to the total number of cvs  ");
    // in this case you need a beta factor!
  }

  std::string basins;
  parse("--basins",basins);

  std::string kt;
  kt=std::string("1.");// assign an arbitrary value just in case that idw.size()==cvs.size()
  if ( dohisto || idw.size()!=0 || basins.length()>0 ) {
    plumed_massert(parse("--kt",kt),"if you make a dimensionality reduction (--idw), a histogram (--histo) or use --basins then you need to define --kt ");
  }

  std::string addme;

  actioninput.push_back("FUNCSUMHILLS");
  actioninput.push_back("ISCLTOOL");

  // set names
  std::string outfile;
  if(parse("--outfile",outfile)) {
    actioninput.push_back("OUTHILLS="+outfile);
  }
  std::string outhisto;
  if(parse("--outhisto",outhisto)) {
    actioninput.push_back("OUTHISTO="+outhisto);
  }


  addme="ARG=";
  for(unsigned i=0; i<(ncv-1); i++) {
    if(cvs[i].size()==1) {
      addme+=std::string(cvs[i][0])+",";
    } else {
      addme+=std::string(cvs[i][0])+"."+std::string(cvs[i][1])+",";
    }
  }
  if(cvs[ncv-1].size()==1) {
    addme+=std::string(cvs[ncv-1][0]);
  } else {
    addme+=std::string(cvs[ncv-1][0])+"."+std::string(cvs[ncv-1][1]);
  }
  actioninput.push_back(addme);
  if(dohills) {
    addme="HILLSFILES=";
    for(unsigned i=0; i<hillsFiles.size()-1; i++) {
      addme+=hillsFiles[i]+",";
    }
    addme+=hillsFiles[hillsFiles.size()-1];
    actioninput.push_back(addme);
    // set the grid
  }
  if(grid_check==2) {
    addme="GRID_MAX=";
    for(unsigned i=0; i<(ncv-1); i++) {
      addme+=gmax[i]+",";
    }
    addme+=gmax[ncv-1];
    actioninput.push_back(addme);
    addme="GRID_MIN=";
    for(unsigned i=0; i<(ncv-1); i++) {
      addme+=gmin[i]+",";
    }
    addme+=gmin[ncv-1];
    actioninput.push_back(addme);
  }
  if(grid_has_bin) {
    addme="GRID_BIN=";
    for(unsigned i=0; i<(ncv-1); i++) {
      addme+=gbin[i]+",";
    }
    addme+=gbin[ncv-1];
    actioninput.push_back(addme);
  }
  if(grid_has_spacing) {
    addme="GRID_SPACING=";
    for(unsigned i=0; i<(ncv-1); i++) {
      addme+=gspacing[i]+",";
    }
    addme+=gspacing[ncv-1];
    actioninput.push_back(addme);
  }
  std::string  stride;
  stride="";
  if(parse("--stride",stride)) {
    actioninput.push_back("INITSTRIDE="+stride);
    bool  nohistory;
    parseFlag("--nohistory",nohistory);
    if(nohistory) {
      actioninput.push_back("NOHISTORY");
    }
  }
  bool  mintozero;
  parseFlag("--mintozero",mintozero);
  if(mintozero) {
    actioninput.push_back("MINTOZERO");
  }
  if(idw.size()!=0) {
    addme="PROJ=";
    for(unsigned i=0; i<idw.size()-1; i++) {
      addme+=idw[i]+",";
    }
    addme+=idw.back();
    actioninput.push_back(addme);
  }

  if(dohisto) {
    if(idw.size()==0) {
      if(sigma.size()!=hcvs.size()) {
        plumed_merror("you should define as many --sigma vector as the number of collective variable used for the histogram ");
      }
    } else {
      if(idw.size()!=sigma.size()) {
        plumed_merror("you should define as many --sigma vector as the number of collective variable used for the histogram ");
      }
    }
  }

  if(idw.size()!=0 || dohisto || basins.length()>0) {
    actioninput.push_back("KT="+kt);
  }
  if(basins.length()>0) {
    actioninput.push_back("BASINS="+basins);
    std::string outbasins;
    if(parse("--outbasins",outbasins)) {
      actioninput.push_back("OUTBASINS="+outbasins);
    }
  }
  std::vector<std::string> average;
  if(parseVector("--average",average)) {
    if(average.size()!=1 && average.size()!=2) {
      plumed_merror("--average needs the initial time and optionally the final time");
    }
    actioninput.push_back("AVERAGE="+average[0]+(average.size()==2 ? ","+average[1] : ""));
  }
  if(dohisto) {
    addme="HISTOFILES=";
    for(unsigned i=0; i<histoFiles.size()-1; i++) {
      addme+=histoFiles[i]+",";
    }
    addme+=histoFiles[histoFiles.size()-1];
    actioninput.push_back(addme);

    addme="HISTOSIGMA=";
    for(unsigned i=0; i<sigma.size()-1; i++) {
      addme+=sigma[i]+",";
    }
    addme+=sigma.back();
    actioninput.push_back(addme);
  }

  bool negbias;
  parseFlag("--negbias",negbias);
  if(negbias) {
    actioninput.push_back("NEGBIAS");
  }

  if(lowI_!=uppI_) {
    addme="INTERVAL=";
    addme+=lowI_+",";
    addme+=uppI_;
    actioninput.push_back(addme);
  }

  std::string sortby;
  if(parse("--sort",sortby)) {
    if(sortby!="time" && sortby!="clock") {
      plumed_merror("--sort can only be time or clock");
    }
    actioninput.push_back("SORT="+sortby);
  }

  std::string fmt;
  fmt="";
  parse("--fmt",fmt);
  if(fmt!="") {
    actioninput.push_back("FMT="+fmt);
  }


  plumed.readInputWords(actioninput,false);
  // if not a grid, then set it up automatically
  return 0;
}

bool CLToolSumHills::findCvsAndPeriodic(const std::string & filename,
                                        std::vector< std::vector<std::string>  > &cvs,
                                        std::vector<std::string> &pmin,
                                        std::vector<std::string> &pmax,
                                        bool &multivariate,
                                        std::string &lowI_,
                                        std::string &uppI_) {
  IFile ifile;
  ifile.allowIgnoredFields();
  if(!ifile.FileExist(filename)) {
    return false;
  }
  cvs.clear();
  pmin.clear();
  pmax.clear();
  ifile.open(filename);
  std::vector<std::string> fields;
  ifile.scanFieldList(fields);
  auto hasField=[&](const std::string & name) {
    return std::find(fields.begin(),fields.end(),name)!=fields.end();
  };
  // In a hills file a column X is a collective variable if its width is also
  // there (sigma_X, or sigma_X_X for multivariate kernels). Names are compared
  // exactly, so that variables such as sigma_1 or lifetime are accepted (see #1226).
  // In a file without widths (e.g. a COLVAR file used with --histo) all the
  // columns but time are variables.
  bool hasSigma=false;
  for(const auto & f : fields) {
    if(f!="time" && (hasField("sigma_"+f) || hasField("sigma_"+f+"_"+f))) {
      hasSigma=true;
      break;
    }
  }
  for(const auto & f : fields) {
    if(f=="time") {
      continue;
    }
    if(hasSigma && !hasField("sigma_"+f) && !hasField("sigma_"+f+"_"+f)) {
      continue;
    }
    // separate the components
    auto dot=f.find_first_of('.');
    if(dot!=std::string::npos) {
      cvs.emplace_back(std::vector<std::string> {f.substr(0,dot),f.substr(dot+1)});
    } else {
      cvs.emplace_back(std::vector<std::string> {f});
    }
    // get periodicity
    pmin.push_back("none");
    pmax.push_back("none");
    if(ifile.FieldExist("min_"+f)) {
      ifile.scanField("min_"+f,pmin.back());
    }
    if(ifile.FieldExist("max_"+f)) {
      ifile.scanField("max_"+f,pmax.back());
    }
  }
  // is multivariate ???
  multivariate=false;
  if(ifile.FieldExist("multivariate")) {
    std::string sss;
    ifile.scanField("multivariate",sss);
    multivariate=(sss=="true");
  }
  // do interval?
  if(ifile.FieldExist("lower_int")) {
    ifile.scanField("lower_int",lowI_);
    ifile.scanField("upper_int",uppI_);
  } else {
    lowI_="-1.";
    uppI_="-1.";
  }
  ifile.scanField();
  return true;
}


PLUMED_REGISTER_CLTOOL(CLToolSumHills,"sum_hills")



}
}
