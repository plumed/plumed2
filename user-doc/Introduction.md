\mainpage Introduction

PLUMED is a plugin that works with a large number of molecular dynamics codes (\ref codes ). 
It can be used to analyze features of the dynamics on-the-fly or to perform a wide variety of free energy methods.
PLUMED can also work as a \ref tools to perform analysis on trajectories saved in most of the
existing formats. If PLUMED is useful for your work please read and cite \cite plumed2, if you are interested in 
the PLUMED 1 original publication please read and cite \cite plumed1 .

To follow the development of PLUMED 2, you can look at the detailed \ref ChangeLog .

To install PLUMED, see this page: \ref Installation , while in \ref Syntax you can find a brief introduction on how to write your first PLUMED input file.

\ref tutorials are available to introduce basic as well as more advanced features of PLUMED.
 
\section AboutManual About this manual

@VERSION@

This is the user manual -  if you want to modify PLUMED or to understand how it works internally, have a look at the 
<a href="../../developer-doc/html/index.html"> developer manual </a>.

@PDFMANUAL@

\section codes Codes interfaced with PLUMED 

PLUMED can be incorporated into an MD code and used to analyze or bias a molecular dynamics run on the fly.
Some MD code could already include calls to the PLUMED library
and be PLUMED-ready in its original distribution.
As far as we know, the following MD codes can be used with PLUMED out of the box:
- [Amber](http://ambermd.org/), pmemd module, since version 20.
- [AmberTools](http://ambermd.org/), sander module, since version 15.
- [CP2K](http://www.cp2k.org), since Feb 2015.
- [ESPResSo](http://espressomd.org), in a version that has been patched with PLUMED can be found
  [here](http://davidebr.github.io/espresso/).
- [PINY-MD](http://github.com/TuckermanGroup/PINY), in its plumed branch.
- [IPHIGENIE](http://sourceforge.net/projects/iphigenie/).
- [AceMD](http://www.multiscalelab.org/acemd/), see [this link](https://github.com/tonigi/ACEMD-PLUMED).
- [OpenMM](http://openmm.org), using the [openmm-plumed plugin](http://github.com/peastman/openmm-plumed).
- [DL_POLY4](https://www.scd.stfc.ac.uk/Pages/DL_POLY.aspx).
- [VNL-ATK](https://quantumwise.com), see [this link](https://docs.quantumwise.com/tutorials/metadynamics_with_plumed/metadynamics_with_plumed.html).
- [ABIN](https://github.com/PHOTOX/ABIN).
- [i-pi](https://github.com/i-pi/i-pi).
- [LAMMPS](https://lammps.sandia.gov/) since Nov 2018.
- [Yaff](https://github.com/molmod/yaff), since Jul 2019.
- [DFTB+](https://www.dftbplus.org/), since release 20.1.
- [Metalwalls](https://gitlab.com/ampere2/metalwalls)

Please refer to the documentation of the MD code to know how to use it with the latest PLUMED release.
If you maintain another MD code that is PLUMED-ready let us know and we will add it to this list.

Additionally, we provide patching procedures for the following codes:

@CODESL@

Alternatively, one
can use PLUMED as a \ref tools for post processing the results from molecular dynamics 
or enhanced sampling calculations.  Notice that PLUMED can be used as an analysis tool
also from the following packages:
- [PLUMED-GUI](http://github.com/tonigi/vmd_plumed) is a [VMD](http://www.ks.uiuc.edu/Research/vmd/) plugin that computes PLUMED collective variables.
- [HTMD](http://www.htmd.org/) can use PLUMED collective variables for analysis.
- [OpenPathSampling](http://openpathsampling.org/), using the [PLUMED Wrapper for OpenPathSampling](https://e-cam.readthedocs.io/en/latest/Classical-MD-Modules/modules/OpenPathSampling/ops_plumed_wrapper/readme.html).

