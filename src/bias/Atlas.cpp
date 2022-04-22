/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2018 The plumed team
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
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionWithValue.h"
#include "tools/IFile.h"
#include "gridtools/Histogram.h"
#include "MetadShortcut.h"
#include "core/ReweightBase.h"

namespace PLMD {
namespace bias {
//+PLUMEDOC BIAS ATLAS
/*
Used to performed the Adaptive Topography of Landscape for Accelerated Sampling(ATLAS) on one or more collective variables.

In an ATLAS simulations a history dependent bias composed of
intermittently added Gaussian functions is added to the potential \cite atlas.

The bias is constructed following a *divide-et-impera* idea, where a high-dimensional CVs space
is firstly partitioned into clusters. Then, a dimensionality reduction technique is applied to the
CVs to create a low-dimensional locally-valid projection that describe in a simple and effective
way each cluster. The global repulsive potential is constructed as a smooth sum of the local
potential acting on each individual cluster.

At the moment, ATLAS requires a Gaussian Mixture Model (GMM) to work. Starting from the
definition of GMM:

\f[
P({\bf s}) = \pi_0 + \sum_{k=1}^M \pi_k \ G({\bf s}|{\bf \mu}_k,{\bf \Sigma}_k),
\f]

where the sum run on clusters, \f$G({\bf s}|{\bf \mu}_k,{\bf \Sigma}_k)\f$ is a normalized gaussian, and \f$pi_0\f$ is a
baseline probability, we can introduce the following quantities

\f[
\theta_k({\bf s}) =  \frac{\pi_k \ G({\bf s}|{\bf \mu}_k,{\bf \Sigma}_k)}{\pi_0+\sum_{l=1}^M \pi_l \ G({\bf s}|{\bf \mu}_k,{\bf \Sigma}_k) }
\f]

called Probabilistic Motif Identifies, introduced by Gasparotto *et al.* \cite pamm. these functions
are used to identify in which cluster of the GMM the system is actually located. It approaches
1 if the system is in minimum \f$k\f$, and approaches 0 if otherwise. A similar  expression
can be also defined for points outside the GMM

\f[
\theta_0({\bf s}) = \frac{\pi_0}{\pi_0+\sum_{l=1}^M \pi_l \ G({\bf s}|{\bf \mu}_k,{\bf \Sigma}_k)  }.
\f]

These functions are used to identify in which cluster the system is, and allow a smooth and
continuous sum of the local potentials that create the global one.

The term \f$pi_0 \f$ is important because it avoids the PMIs to be undetermined when numerator and denominator
are close to zero, and decide the probability threshold below which we do not consider
points in the GMM (this is true also for cluster, so cluster smaller than the threshold will not be considere)

To select the threshold, we can safely used

\f[
\pi_0(f_0) = \frac{\pi_k}{\sqrt{(2 \pi)^{n_s} |{\bf \Sigma}_k|}}\ e^{-z_0^2/2}
\f]

where \f$z_0=ICDF_{\chi^2(n_s)}(1-f_0)\f$ is the value of the Gaussian
exponent that corresponds to the isocontour that discards a fraction \f$f_0\f$ of the probability.

To reduce the dimensionality of the CVs space there is ample freedom. At the moment, ATLAS
support the formulation of local coordinates \f$\vec{c}\f$ using a projection on the major
Principal Components of the Covariance matrix \f$\vec{Sigma}_k\f$. Once that a number of PC has been
decided (1, 2, 3 etc...) the local coordinates \f$\vec{c}\f$ are obtained  with the following
projection

\f[
c_k^l = \frac{{\bf s}^{T} {\bf U}_k^l}{\sqrt{\lambda_k^l}},
\f]

where \f$\lambda_k^l\f$ is the \f$l\f$ eigenvalue  of the \f$k\f$ cluster, and \f${\bf U}_k^l\f$ is the
corresponding eigenvector. For periodic variables the projection is not-normalized.

Once that the CVs space is partitioned, the PMIs function are defined and the low-dimensional space`
are constructed, the global potential at a given time \f$t\f$ can be expressed as

\f[
V({\bf s},t) = \sum_{k=1}^M v_k({\bf s},t)\ \theta_k({\bf s}) + v_0({\bf s},t) \theta_0({\bf s}).
\f]

The local potentials \f$v_k({\bf s},t)\f$ are formulated as

\f[
v_k({\bf s},t) = h   \sum_{t'=0}^t e^{-V({\bf s}(t'),t')\Delta T^{-1}} g({\bf c}_k({\bf s})-{\bf c}_k(t')) \frac{ \theta_k({\bf s}(t'))}{\sum_{l=0}^M \theta_l({\bf s}(t'))^2}.
\f]

In this expression, \f$g({\bf c}_k-{\bf c}_k(t'))\f$ is a non-normalized Gaussian function
computed relative to the local, low-dimensional variables. The indicator functions act so
that bias is only added to the basin the system is in at a given time.

A similar expression, but without repuslive gaussian, can be defined for \f$v_0({\bf s},t)\f$
 which is a potential acting if the system escape the GMM, and is responsible of bringing it back.

 \f[
v_0({\bf s},t) =\ h   \sum_{t'=0}^T e^{-V({\bf s}(t'),t')\Delta T^{-1}} \frac{ \theta_0({\bf s}(t'))}{\sum_{l=0}^M \theta_l({\bf s}(t'))^2}.
 \f]

To calculate unbiased quantities after the ATLAS calculation is necessary to weight the sampled microstates
according to the potential deposited during the simulation. To his extent, we suggest to use the Iterative Trajectory Reweigthing
method, as it does not depends on the dimensionality of the space sampled during the calculation.

As the simulation proceed, ATLAS automatically produce two files, THETA and LOWD_CVS.

The THETA file contains the values of the PMIs during the calculation, including the one
identifying the points outside the GMM.

The LOWD_CVS file contains the valued of the local low-dimensional CVs. The last
column of the file contains the height of the repulsive gaussian h, required by ITRE
to evaluate the repulsive potential.

\par Examples

As first example, let's see a simple example for a langevin particle in a 2D potential
that can be simulated using the pesmd command. To mimic the one particle in 2D, plumed
require the definition of 2 atoms, and the CVs that are going to be used are the
first two components of the distance between the two atoms

\plumedfile
UNITS NATURAL
d1: DISTANCE ATOMS=1,2 COMPONENTS

ff: MATHEVAL ...
ARG=d1.x,d1.y VAR=x0,x1 PERIODIC=NO
FUNC=30.0*exp(-4.0*(x0-1)^2-4.0*(x1-1)^2)+30.0*exp(-4.0*(x0-1)^2-4.0*(x1+1)^2)+30.0*exp(-4.0*(x0+1)^2-4.0*(x1+1)^2)+30.0/(1.0/((5.0*x0+5.0)^2+(5.0*x1+5.0)^2+1)+1.0/((5.0*x0-5.0)^2+(5.0*x1+5.0)^2+1)+1.0/((5.0*x0-5.0)^2+(5.0*x1-5.0)^2+1)+1.0/((-3.5355339059*x0+3.5355339059*x1)^2+(0.500000000834386*x0+0.500000000834386*x1)^8)+1.0/(1.0*x1^8+(5.0*x0-5.0)^2)+1.0/(1.0*x0^8+(5.0*x1+5.0)^2))
...

bb: BIASVALUE ARG=ff

at: ATLAS ...
ARG=d1.x,d1.y REFERENCE=cluster.dat PACE=500
SIGMA=0.20 BIASFACTOR=10 HEIGHT=2.0
GRID_MAX=6.0 GRID_BIN=600 TEMP=1 TRUNCATE_GRIDS
REGULARISE=1E-8
STATIC_WALL=0.0
ADAPTIVE_WALL=1.0
...

PRINT ARG=at,at_wtfact FILE=rr.gbias STRIDE=500
PRINT ARG=at_adaptive_wall FILE=rr.wall STRIDE=500
PRINT ARG=d1.x,d1.y,bb.bias FILE=TRAJ STRIDE=500
\endplumedfile

A few explanation are required. Each local potential is estimated on a grid so that
the method does not scale with \f$t^2\f$. To do so, we need to be sure that the
grid will cover values of the local CVs for which the PMIs are not zeros. All the points
that are far from the cluster and for which the PMI is zero can be ignored and not
saved, hende we discard them with TRUNCATE_GRIDS. REGULARISE is the keyword that set
the probability threshold \f$\pi_0\f$. There is only one SIGMA parameter and only
one GRID declaration since we assume the local CVs to span the same amount of space
since they are normalized. The STATIC_WALL and ADAPTIVE_WALL control the absence of presence
of the wall.

The GMM is declared in a file, cluster.dat, which reads as

\auxfile{cluster_plumed.dat}
#! FIELDS d1.x d1.y sigma_d1.x_d1.x sigma_d1.x_d1.y sigma_d1.y_d1.x sigma_d1.y_d1.y height
#! SET kerneltype gaussian
#! SET neigv 2
#! SET residual false
0.01671000 0.05344100 0.14059788 0.13116870 0.13116870 0.14015185 0.48990991
0.07364000 -1.01967600 0.13797482 0.00230512 0.00230512 0.01167570 0.26019881
1.01244000 -0.10820500 0.01063001 -0.00137699 -0.00137699 0.13882161 0.24989128
\endauxfile

The first line contains the declaration for PLUMED of the fields. The subsequent three
lines declare the type of kernel (gaussian or von-misses), the number of eigenvectors
to use in the dimensionality reduction technique (1, 2, 3 etc...). The last one, if set to true, instruct
PLUMED to use the first n eigenvectors and the residual projection.

To illustrate how this works in real atomistic systems, here we report two other cases,
Alanine dipeptide as well as a cluster of 38 LJ atoms. The parameters are the same as those
just explained for the 2D Langevin particle.

The following input is for a ATLAS calculation using as
collective variables the number of atoms coordinated with 4, 5, 6, 7, 8, 9, 10, 11, for
a cluster of 38 atoms in vacuum interacting through a LennardJones potential.

\plumedfile
UNITS NATURAL
nsa: COORDINATIONNUMBER SPECIES=1-38 SWITCH={CUBIC D_0=1.25 D_MAX=1.5}
n4: MATHEVAL ARG=nsa FUNC=exp(-(x-4)*(x-4)/(2*0.5*0.5)) PERIODIC=NO
n5: MATHEVAL ARG=nsa FUNC=exp(-(x-5)*(x-5)/(2*0.5*0.5)) PERIODIC=NO
n6: MATHEVAL ARG=nsa FUNC=exp(-(x-6)*(x-6)/(2*0.5*0.5)) PERIODIC=NO
n7: MATHEVAL ARG=nsa FUNC=exp(-(x-7)*(x-7)/(2*0.5*0.5)) PERIODIC=NO
n8: MATHEVAL ARG=nsa FUNC=exp(-(x-8)*(x-8)/(2*0.5*0.5)) PERIODIC=NO
n9: MATHEVAL ARG=nsa FUNC=exp(-(x-9)*(x-9)/(2*0.5*0.5)) PERIODIC=NO
n10: MATHEVAL ARG=nsa FUNC=exp(-(x-10)*(x-10)/(2*0.5*0.5)) PERIODIC=NO
n11: MATHEVAL ARG=nsa FUNC=exp(-(x-11)*(x-11)/(2*0.5*0.5)) PERIODIC=NO

at: ATLAS ...
ARG=n4,n5,n6,n7,n8,n9,n10,n11
REFERENCE=cluster_plumed.dat
PACE=500
SIGMA=0.1
BIASFACTOR=10
HEIGHT=0.5
TEMP=0.12
ADAPTIVE_WALL=1.0 STATIC_WALL=0.0
REGULARISE=1E-12
GRID_MAX=7.5 GRID_BIN=750 TRUNCATE_GRIDS
...

PRINT ARG=n4,n5,n6,n7,n8,n9,n10,n11 FMT=%8.4f FILE=colvar STRIDE=500
PRINT ARG=at_adaptive_wall FMT=%8.4f FILE=wall STRIDE=500
PRINT ARG=at FMT=%8.4f FILE=bias STRIDE=500
\endplumedfile

To work, ATLAS requires a file containing the clusters of the GMM, cluster_plumed.dat.
The file would read as:

\auxfile{cluster_plumed.dat}
#! FIELDS n4 n5 n6 n7 n8 n9 n10 n11 sigma_n4_n4 sigma_n4_n5 sigma_n4_n6 sigma_n4_n7 sigma_n4_n8 sigma_n4_n9 sigma_n4_n10 sigma_n4_n11 sigma_n5_n4 sigma_n5_n5 sigma_n5_n6 sigma_n5_n7 sigma_n5_n8 sigma_n5_n9 sigma_n5_n10 sigma_n5_n11 sigma_n6_n4 sigma_n6_n5 sigma_n6_n6 sigma_n6_n7 sigma_n6_n8 sigma_n6_n9 sigma_n6_n10 sigma_n6_n11 sigma_n7_n4 sigma_n7_n5 sigma_n7_n6 sigma_n7_n7 sigma_n7_n8 sigma_n7_n9 sigma_n7_n10 sigma_n7_n11 sigma_n8_n4 sigma_n8_n5 sigma_n8_n6 sigma_n8_n7 sigma_n8_n8 sigma_n8_n9 sigma_n8_n10 sigma_n8_n11 sigma_n9_n4 sigma_n9_n5 sigma_n9_n6 sigma_n9_n7 sigma_n9_n8 sigma_n9_n9 sigma_n9_n10 sigma_n9_n11 sigma_n10_n4 sigma_n10_n5 sigma_n10_n6 sigma_n10_n7 sigma_n10_n8 sigma_n10_n9 sigma_n10_n10 sigma_n10_n11 sigma_n11_n4 sigma_n11_n5 sigma_n11_n6 sigma_n11_n7 sigma_n11_n8 sigma_n11_n9 sigma_n11_n10 sigma_n11_n11 height
#! SET kerneltype gaussian
#! SET neigv 2
#! SET residual false
2.380300 4.518200 8.394300 5.288800 4.755600 6.038500 2.343800 1.876900 1.259338 0.130634 -0.869997 0.016250 -0.438566 -0.375105 0.127170 0.318129 0.130634 3.423514 -1.742815 -0.284503 -0.862195 -0.628138 0.210197 0.765419 -0.869997 -1.742815 3.586251 0.304950 -0.252102 -0.098865 -0.176326 0.174368 0.016250 -0.284503 0.304950 2.967988 -0.967389 -1.782510 0.232075 0.677206 -0.438566 -0.862195 -0.252102 -0.967389 3.245160 1.193654 -0.917513 -1.716616 -0.375105 -0.628138 -0.098865 -1.782510 1.193654 3.492971 -0.738052 -1.995773 0.127170 0.210197 -0.176326 0.232075 -0.917513 -0.738052 1.116574 0.469500 0.318129 0.765419 0.174368 0.677206 -1.716616 -1.995773 0.469500 2.419072 0.016174
0.007800 3.196600 23.995600 3.305800 1.089700 8.000000 1.085800 0.771400 0.035839 0.140432 -0.217355 0.037738 0.082495 -0.097664 0.012069 0.039429 0.140432 1.467490 -0.240125 -1.238690 0.349054 -0.571226 0.109911 0.228116 -0.217355 -0.240125 8.377959 -7.203111 -1.654707 1.548503 -0.184798 -0.386602 0.037738 -1.238690 -7.203111 7.673440 1.112053 -0.842305 0.051252 0.116922 0.082495 0.349054 -1.654707 1.112053 0.628199 -0.531707 -0.057895 0.124313 -0.097664 -0.571226 1.548503 -0.842305 -0.531707 0.803523 -0.203294 -0.174698 0.012069 0.109911 -0.184798 0.051252 -0.057895 -0.203294 0.262751 0.042041 0.039429 0.228116 -0.386602 0.116922 0.124313 -0.174698 0.042041 0.114536 0.012073
0.417400 5.011500 13.144200 6.298700 12.873800 1.546100 0.140800 1.734900 0.221036 -0.010111 -0.233996 0.149060 -0.292659 0.091049 0.087019 -0.058476 -0.010111 2.441686 -1.692420 -0.737146 0.545939 -0.486565 -0.213106 0.103008 -0.233996 -1.692420 2.677948 -0.256248 -0.315485 -0.020455 -0.015110 -0.036793 0.149060 -0.737146 -0.256248 2.830081 -2.118312 0.048805 0.215480 0.048522 -0.292659 0.545939 -0.315485 -2.118312 3.728157 -0.974893 -0.498718 0.101869 0.091049 -0.486565 -0.020455 0.048805 -0.974893 1.071804 0.234893 -0.159643 0.087019 -0.213106 -0.015110 0.215480 -0.498718 0.234893 0.177750 -0.045712 -0.058476 0.103008 -0.036793 0.048522 0.101869 -0.159643 -0.045712 0.222320 0.966387
0.836600 9.477500 9.255700 3.956400 4.583800 9.945900 1.105400 0.555700 1.029512 -0.366703 -1.129590 -0.528825 0.565086 -0.375293 0.060877 0.109134 -0.366703 2.876908 -1.717325 0.162965 -0.858175 0.016411 0.170222 0.015101 -1.129590 -1.717325 4.261461 0.874028 -0.649926 0.266438 -0.329186 -0.183487 -0.528825 0.162965 0.874028 1.599660 -1.144236 -0.277194 0.053264 -0.003173 0.565086 -0.858175 -0.649926 -1.144236 2.784413 -1.044301 -0.367285 0.050473 -0.375293 0.016411 0.266438 -0.277194 -1.044301 2.367485 -0.140604 -0.386127 0.060877 0.170222 -0.329186 0.053264 -0.367285 -0.140604 0.366430 0.122679 0.109134 0.015101 -0.183487 -0.003173 0.050473 -0.386127 0.122679 0.215191 0.005272
1.979700 7.842800 13.770300 1.825200 3.084100 5.325800 3.046100 3.826900 0.526268 -0.278447 -0.071992 -0.214195 -0.256734 0.108387 0.115535 -0.064924 -0.278447 1.825603 -1.628482 0.140279 0.494932 -0.541937 -0.099556 0.266131 -0.071992 -1.628482 2.646550 -0.116345 -0.756682 0.094646 0.089681 0.218240 -0.214195 0.140279 -0.116345 0.291698 0.296833 -0.218239 -0.034190 -0.008709 -0.256734 0.494932 -0.756682 0.296833 0.655983 -0.159005 -0.155757 -0.163024 0.108387 -0.541937 0.094646 -0.218239 -0.159005 0.959316 -0.191525 -0.597566 0.115535 -0.099556 0.089681 -0.034190 -0.155757 -0.191525 0.202282 0.146922 -0.064924 0.266131 0.218240 -0.008709 -0.163024 -0.597566 0.146922 0.614569 0.000095
\endauxfile

To finish with the examples, here is the PLUMED input and GMM.dat files for
Alanine Dipeptide in vacuum.

\plumedfile
UNITS ENERGY=kj/mol
t1: TORSION ATOMS=5,7,9,15
t2: TORSION ATOMS=7,9,15,17

at: ATLAS ...
REFERENCE=cluster_plumed.dat PACE=500
ARG=t1,t2
SIGMA=0.10 BIASFACTOR=10 HEIGHT=2.0
GRID_MAX=5.0 GRID_BIN=500 TEMP=300 TRUNCATE_GRIDS
REGULARISE=1E-6
STATIC_WALL=0.0
ADAPTIVE_WALL=1.0
...

PRINT ARG=t1,t2 FILE=colvar STRIDE=500
PRINT ARG=at,at_wtfact FILE=at.gbias STRIDE=500
PRINT ARG=at_adaptive_wall FILE=wall STRIDE=500
\endplumedfile

With the GMM file readings:

\auxfile{cluster_plumed.dat}
#! FIELDS t1 t2 sigma_t1_t1 sigma_t1_t2 sigma_t2_t1 sigma_t2_t2 height
#! SET neigv 2
#! SET kerneltype von-misses
#! SET residual false
-1.482918 1.047509 0.301925 0.000000 0.000000 0.646060 0.254361
1.054132 -0.726058 0.267969 0.000000 0.000000 1.327605 0.179060
-2.474009 2.817684 0.355761 0.000000 0.000000 0.339309 0.236480
-2.870169 -0.701767 0.061753 0.000000 0.000000 0.244917 0.330099
\endauxfile


*/
//+ENDPLUMEDOC


class Atlas : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit Atlas(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(Atlas,"ATLAS")

void Atlas::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","ARG","the arguments that should be used as input to this method");
  keys.add("compulsory","REFERENCE","the input file containing the definitions of the clusters");
  keys.add("compulsory","PACE","the frequency with which hills should be added");
  keys.add("compulsory","SIGMA","the widths of the Gaussians that we are using");
  keys.add("compulsory","HEIGHT","the heights of the hills that should be added");
  keys.add("compulsory","BIASFACTOR","the bias factor for the well tempered metadynamics");
  keys.add("optional","GRID_MAX","the maximum value to use for all the grids");
  keys.add("optional","GRID_BIN","the number of bins to use for all the grids");
  keys.add("compulsory","REGULARISE","0.001","don't allow the denominator to be smaller then this value");
  keys.add("compulsory","STATIC_WALL","the force constant of the wall applied outside the GMM");
  keys.add("compulsory","ADAPTIVE_WALL","the force constant of the wall applied outside the GMM");
  keys.add("optional","TEMP","the system temperature - this is only needed if you are doing well-tempered metadynamics");
  keys.addFlag("TRUNCATE_GRIDS",false,"set all histograms equal to zero outside specified range");
  keys.add("compulsory","THETA_FILE","THETA","print a file containing the kernel values with this name");
  keys.add("compulsory","LOWD_CVS_FILE","LOWD_CVS","print a file containing the values of the low dimensional CVs");
}

Atlas::Atlas(const ActionOptions& ao):
Action(ao),
ActionShortcut(ao)
{
  // Read the reference file and determine how many clusters we have
  bool truncate=false; parseFlag("TRUNCATE_GRIDS",truncate);
  std::string ktype, argstr; parse("ARG",argstr); std::vector<unsigned> neigv; std::vector<bool> resid;
  std::string fname; parse("REFERENCE",fname); std::vector<double> weights;
  IFile ifile; ifile.open(fname); ifile.allowIgnoredFields(); double h;
  for(unsigned k=0;; ++k) {
     if( !ifile.scanField("height",h) ) break;
     int meig; ifile.scanField("neigv",meig); neigv.push_back( meig );
     if( meig>0 ) {
         std::string ires; ifile.scanField("residual",ires);
         if( ires=="true" ) resid.push_back( true );
         else if( ires=="false" ) resid.push_back( false );
         else error("residual flag should be set to true/false");
     } else resid.push_back( false );
     // Create a Kernel for this cluster
     std::string num, wstr; Tools::convert( k+1, num ); Tools::convert( h, wstr ); ifile.scanField("kerneltype",ktype);
     readInputLine( getShortcutLabel() + "_kernel-" + num + ": KERNEL NORMALIZED ARG=" + argstr + " NUMBER=" + num + " REFERENCE=" + fname + " WEIGHT=" + wstr + " TYPE=" + ktype );
     // Compute eigenvalues and eigenvectors for the input covariance matrix if required
     if( meig>0 ) {
         std::string seig="1"; for(int j=1;j<meig;++j) { std::string eignum; Tools::convert( j+1, eignum ); seig += "," + eignum; }
         readInputLine( getShortcutLabel() + "_eigv" + num + ": DIAGONALIZE ARG=" + getShortcutLabel() + "_kernel-" + num + "_cov VECTORS=" + seig );
     }
     // Store the weights as we will use these when constructing the bias later in the input
     weights.push_back(h); ifile.scanField();
  }
  ifile.close();

  // Now build the basins
  for(unsigned k=0;k<weights.size();++k) {
      std::string num; Tools::convert( k+1, num );
      // Compute the distance between the center of the basin and the current configuration
      readInputLine( getShortcutLabel() + "_dist-" + num + ": MATHEVAL ARG=" + getShortcutLabel() + "_kernel-" + num + "_dist_2 FUNC=sqrt(x) PERIODIC=NO");
      // Get the negative of the distance from the center of the basin
      if( neigv[k]==0 ) {
        // And the reflection of the distance
        readInputLine( getShortcutLabel() + "_pdist-" + num + ": MATHEVAL ARG1=" + getShortcutLabel() + "_dist-" + num + " FUNC=0-x PERIODIC=NO");
      } else {
        // This computes the projections of the difference between the current point and the origin on the various eigenvectors
        std::string argstr = "ARG1=" + getShortcutLabel() + "_dist-" + num; std::string coeffstr="COEFFICIENTS=1"; std::string powstr="POWERS=2";
        for(unsigned i=0;i<neigv[k];++i) {
            coeffstr +=",-1"; powstr +=",2"; std::string anum, eignum; Tools::convert( i+1, eignum );
            Tools::convert( i+2, anum ); argstr += " ARG" + anum + "=" + getShortcutLabel() + "_proj" + eignum + "-" + num + " ";
            // Multiply difference in CVs by eigenvector - returns a vector
	    ActionWithValue* av=plumed.getActionSet().selectWithLabel<ActionWithValue*>(getShortcutLabel() + "_kernel-" + num + "_dist_2_diff" ); plumed_assert( av ); //////
	    // By default, we set the low dimensional CVs to be non-periodic. As at this stage periodic CVs has a diagonal covariance matrix, this affect in a
	    // minimum way the projection of periodic variable
            readInputLine( getShortcutLabel() + "_udproj" + eignum + "-" + num + ": DOT DIAGONAL_ELEMENTS_ONLY ARG1=" + getShortcutLabel() + "_kernel-" + num + "_dist_2_diff"
               + " ARG2=" + getShortcutLabel() + "_eigv" + num + ".vecs-" + eignum );
            // Divide the projection on the eigenvector by the eigenvalue so that gaussian widths are in units of covariance
	    // However, since it seems quite complex to normalize the periodic boundary too, we do not normalize the non-periodic boundaries for the sqrt(eigval)
	    // As a matter of fact, for periodic CVs this procedure is basically a selection of the most important modes in the basins
	    if( av->copyOutput(0)->isPeriodic()) {
		    // Periodic CVs -> not normalized
	            std::string min, max; av->copyOutput(0)->getDomain(min,max);
	            readInputLine( getShortcutLabel() + "_proj" + eignum + "-" + num + ": MATHEVAL ARG1="+  getShortcutLabel() + "_udproj" + eignum + "-" + num
	        		    +  " FUNC=x PERIODIC="+min+","+max );
	    } else {
		    // Non periodic CVs -> normalized
	            readInputLine( getShortcutLabel() + "_proj" + eignum + "-" + num + ": MATHEVAL ARG1="+  getShortcutLabel() + "_udproj" + eignum + "-" + num
	        		    + " ARG2=" + getShortcutLabel() + "_eigv" + num + ".vals-" + eignum + " FUNC=x/sqrt(y) PERIODIC=NO");
	    } ////
        }
        // Add this command to compute the residual distance
        if( resid[k] ) {
            readInputLine( getShortcutLabel() + "_resid2-" + num + ": COMBINE PERIODIC=NO " + argstr + coeffstr + " " + powstr );
            readInputLine( getShortcutLabel() + "_resid-" + num + ": MATHEVAL ARG1=" + getShortcutLabel() + "_resid2-" + num + " FUNC=sqrt(x) PERIODIC=NO");
            readInputLine( getShortcutLabel() + "_presid-" + num + ": MATHEVAL ARG1=" + getShortcutLabel() + "_resid-" + num + " FUNC=0-x PERIODIC=NO");
        }
      }
  }

  // Create the well-tempered weight
  std::string biasfactor, height; parse("HEIGHT",height); parse("BIASFACTOR",biasfactor);
  double temp=0.0; parse("TEMP",temp); std::string tempstr="";
  if( temp>0.0 ) { std::string tstr; Tools::convert( temp, tstr); tempstr = " TEMP=" + tstr; }
  readInputLine( getShortcutLabel() + "_wtfact: REWEIGHT_WELLTEMPERED HEIGHT=" + height + " BIASFACTOR=" + biasfactor + tempstr);

  // And sum the kernels
  std::string cinput = getShortcutLabel() + "_ksum: COMBINE PERIODIC=NO ARG=" + getShortcutLabel() + "_kernel-1"; 
  for(unsigned k=1;k<weights.size();++k) { std::string num; Tools::convert( k+1, num ); cinput += "," + getShortcutLabel() + "_kernel-" + num; } 
  readInputLine( cinput ) ; 

  // Add a small number to regularize the sum
  std::string regparam; parse("REGULARISE",regparam);
  readInputLine( getShortcutLabel() + "_rksum: MATHEVAL ARG1=" + getShortcutLabel() + "_ksum FUNC=x+" + regparam + " PERIODIC=NO");

  // Normalize the weights for each of the kernels
  for(unsigned k=0;k<weights.size();++k) {
      std::string num; Tools::convert( k+1, num );
      // And now compute the final weights of the basins
      readInputLine( getShortcutLabel() + "_wkernel-" + num + ": MATHEVAL ARG1=" + getShortcutLabel() + "_kernel-" + num + " ARG2=" + getShortcutLabel() +
                     "_rksum FUNC=x/y PERIODIC=NO");
  }

  // And compute the wkernel outside the GMM
  readInputLine( getShortcutLabel() + "_ext_wkernel: MATHEVAL ARG1=" + getShortcutLabel() + "_ksum FUNC=" + regparam + "/(x+" + regparam + ") PERIODIC=NO");

  // And sum the wkernels for the renormalization
  std::string ccinput = getShortcutLabel() + "_wksum: COMBINE PERIODIC=NO ARG=" + getShortcutLabel() + "_wkernel-1", ppwrs=" POWERS=2";
  for(unsigned k=1;k<weights.size();++k) { std::string num; Tools::convert( k+1, num ); ccinput += "," + getShortcutLabel() + "_wkernel-" + num; ppwrs += ",2"; }
  ccinput +="," + getShortcutLabel() + "_ext_wkernel"; ppwrs += ",2";
  readInputLine( ccinput + ppwrs );

  // Setup the histograms that will store the bias potential for each basin and compute the instantaneous bias from each basin
  std::string truncflag1="", truncflag2=""; if( truncate ) { truncflag1="IGNORE_IF_OUT_OF_RANGE"; truncflag2="ZERO_OUTSIDE_GRID_RANGE"; }
  std::string gmax, grid_nbins, pacestr, hstring; std::vector<std::string> sigma(1); std::vector<std::string> kargs,eargs,tgmin,tgmax,tgbins;
  parse("GRID_MAX",gmax); parse("GRID_BIN",grid_nbins); parse("SIGMA",sigma[0]); parse("PACE",pacestr);
  // Build the histograms for the bias potential
  readInputLine( getShortcutLabel() + "_height: CONSTANT VALUE=1.0");
  for(unsigned k=0;k<weights.size();++k) {
      std::string num; Tools::convert( k+1, num ); kargs.resize(0); eargs.resize(0); tgmin.resize(0); tgmax.resize(0); tgbins.resize(0);
      readInputLine(getShortcutLabel() + "_logwkernel-" + num + ": MATHEVAL ARG1=" + getShortcutLabel() + "_wkernel-" + num +
                                                                          " ARG2=" + getShortcutLabel() + "_wksum FUNC=log(x/y) PERIODIC=NO");
      readInputLine(getShortcutLabel() + "-" + num + "_wtfact: MATHEVAL ARG1=" + getShortcutLabel() + "_wtfact ARG2=" + getShortcutLabel() + "_logwkernel-" +
                    num + " FUNC=x+y PERIODIC=NO"); hstring = getShortcutLabel() + "-" + num + "_wtfact";
      if( neigv[k]==0 ) {
          readInputLine( getShortcutLabel() + "_cdist-" + num + ": CONCATENATE " + 
                            "ARG=" + getShortcutLabel() + "_dist-" + num + "," + getShortcutLabel() + "_pdist-" + num );
          kargs.push_back( getShortcutLabel() + "_cdist-" + num ); eargs.push_back( getShortcutLabel() + "_dist-" + num ); 
          // Convert the bandwidth to something constant actions
          gridtools::HistogramTools::convertBandwiths( getShortcutLabel() + "-" + num, sigma, this );
          if( grid_nbins.size()>0 ) {
              if( gmax.size()==0 ) error("you must set GRID_MAX if you set GRID_BIN");
              tgmin.push_back("0"); tgmax.push_back(gmax); tgbins.push_back( grid_nbins );
          } else {
              readInputLine( getShortcutLabel() + "-" + num + "_nwtfact: MATHEVAL ARG1=" + getShortcutLabel() + "-" + num + "_wtfact FUNC=x-log(2) PERIODIC=NO");
              readInputLine( getShortcutLabel() + "-" + num + "_hnwtfact: CONCATENATE ARG=" + getShortcutLabel() + "-" + num + "_nwtfact," + getShortcutLabel() + "-" + num + "_nwtfact");
              hstring = getShortcutLabel() + "-" + num + "_hnwtfact";
          }
      } else {
          std::vector<std::string> bw_str( neigv[k], sigma[0] ); if( resid[k] ) bw_str.push_back( sigma[0] );
          // Convert the bandwidth to something constant actions
          gridtools::HistogramTools::convertBandwiths( getShortcutLabel() + "-" + num, bw_str, this ); kargs.resize(0); eargs.resize(0);
          for(unsigned i=0;i<neigv[k];++i) {
              std::string eignum; Tools::convert( i+1, eignum );
              if( resid[k] ) { 
                   readInputLine( getShortcutLabel() + "_cproj" + eignum + "-" + num + ": CONCATENATE " +
                                     "ARG=" + getShortcutLabel() + "_proj" + eignum + "-" + num + "," + getShortcutLabel() + "_proj" + eignum + "-" + num );
                   kargs.push_back( getShortcutLabel() + "_cproj" + eignum + "-" + num ); eargs.push_back( getShortcutLabel() + "_proj" + eignum + "-" + num );
              } else { kargs.push_back( getShortcutLabel() + "_proj" + eignum + "-" + num ); eargs.push_back( getShortcutLabel() + "_proj" + eignum + "-" + num ); }
              if( grid_nbins.size()>0 ) { 
                  ActionWithValue* av=plumed.getActionSet().selectWithLabel<ActionWithValue*>(getShortcutLabel() + "_kernel-" + num + "_dist_2_diff" );
                  if( av->copyOutput(0)->isPeriodic() ) {
                      std::string min, max; av->copyOutput(0)->getDomain(min,max);
                      if( gmax.size()==0 ) warning("ignoring specified GRID_MAX and using domain of periodic variable instead");
                      tgmin.push_back( min ); tgmax.push_back( max );
                  } else { 
                      if( gmax.size()==0 ) error("you must set GRID_MAX if you set GRID_BIN");
                      tgmin.push_back( "-" + gmax ); tgmax.push_back( gmax ); 
                  }
                  tgbins.push_back( grid_nbins ); 
              }
          }
          if( resid[k] ) {
              readInputLine( getShortcutLabel() + "_cresid-" + num + ": CONCATENATE ARG=" + getShortcutLabel() + "_resid-" + num + "," + getShortcutLabel() + "_presid-" + num );
              kargs.push_back( getShortcutLabel() + "_cresid-" + num ); eargs.push_back( getShortcutLabel() + "_resid-" + num );
              if( grid_nbins.size()>0 ) {
                  if( gmax.size()==0 ) error("you must set GRID_MAX if you set GRID_BIN");
                  tgmin.push_back( "-" + gmax ); tgmax.push_back( gmax ); tgbins.push_back( grid_nbins );
              } else {
                  readInputLine( getShortcutLabel() + "-" + num + "_nwtfact: MATHEVAL ARG1=" + getShortcutLabel() + "-" + num + "_wtfact FUNC=x-log(2) PERIODIC=NO");
                  readInputLine( getShortcutLabel() + "-" + num + "_hnwtfact: CONCATENATE ARG=" + getShortcutLabel() + "-" + num + "_nwtfact," + getShortcutLabel() + "-" + num + "_nwtfact");
                  hstring = getShortcutLabel() + "-" + num + "_hnwtfact"; 
              }
          }
      }
      MetadShortcut::createMetadBias( getShortcutLabel() + "-" + num, pacestr, kargs, eargs, tgmin, tgmax, tgbins, hstring, truncflag1, truncflag2, this );
  }

  // Normalize the weights for each of the kernels and compute the final bias
  for(unsigned k=0;k<weights.size();++k) {
      std::string num; Tools::convert( k+1, num );
      // And the bias due to each basin (product of bias due to basin and kernel weight)
      readInputLine( getShortcutLabel() + "_wbias-" + num + ": MATHEVAL MIX_HISTORY_DEPENDENCE ARG1=" + getShortcutLabel() + "-" + num + "_bias ARG2=" +
                     getShortcutLabel() + "_wkernel-" + num + " FUNC=x*y PERIODIC=NO");
  }

  // And calculate the static_wall
  std::string static_wall; parse("STATIC_WALL",static_wall);
  readInputLine( getShortcutLabel() + "_static_wall: MATHEVAL ARG1=" + getShortcutLabel() + "_ext_wkernel FUNC=" + static_wall + "*x/(1-x) PERIODIC=NO");

  // And calculate the adaptive_wall
  std::string adaptive_wall; parse("ADAPTIVE_WALL",adaptive_wall);
  readInputLine( getShortcutLabel() + "_height_adaptive_wall: MATHEVAL ARG1=" + getShortcutLabel() + "_wksum "+
                                                                      "ARG2=" + getShortcutLabel() + "_ext_wkernel "+
                                                                      "ARG3=" + getShortcutLabel() + "_wtfact FUNC="+adaptive_wall+"*y/x*exp(z) PERIODIC=NO");
  readInputLine( getShortcutLabel() + "_cum_adaptive_wall: AVERAGE NORMALIZATION=false CLEAR=0 STRIDE="+pacestr+" ARG=" + getShortcutLabel() + "_height_adaptive_wall");
  readInputLine( getShortcutLabel() + "_adaptive_wall: MATHEVAL MIX_HISTORY_DEPENDENCE ARG1=" + getShortcutLabel() + "_cum_adaptive_wall ARG2=" + getShortcutLabel() + "_ext_wkernel FUNC=x*y PERIODIC=NO");

  // This is for the sum of these quantities
  std::string combstr = getShortcutLabel() + ": COMBINE MIX_HISTORY_DEPENDENCE PERIODIC=NO ARG=" + getShortcutLabel() + "_static_wall," + getShortcutLabel() + "_adaptive_wall," + getShortcutLabel() + "_wbias-1";
  for(unsigned k=1;k<weights.size();++k) { std::string num; Tools::convert( k+1, num ); combstr += "," + getShortcutLabel() + "_wbias-" + num; }

  // And the final bias
  readInputLine( combstr ); readInputLine("BIASVALUE ARG=" + getShortcutLabel() );
  readInputLine( getShortcutLabel() + "_wtheight: MATHEVAL PERIODIC=NO ARG=" + getShortcutLabel() + "_wtfact" + " FUNC=exp(x)");

  // Print the theta values to the THETA file
  std::string theta_file; parse("THETA_FILE",theta_file);
  std::string theta_str = "PRINT FILE=" + theta_file + " STRIDE="+pacestr+" ARG=" + getShortcutLabel() + "_wkernel-1" ;
  for(unsigned k=1;k<weights.size();++k) {
    std::string num; Tools::convert( k+1, num );
    theta_str += "," + getShortcutLabel() + "_wkernel-" + num;
  }
  theta_str += "," + getShortcutLabel() + "_ext_wkernel";
  readInputLine( theta_str );

  // Print the reduced CVs to a file
  std::string lowd_cvs_f; parse("LOWD_CVS_FILE",lowd_cvs_f);
  std::string cvs_str = "PRINT FILE=" + lowd_cvs_f + " STRIDE="+pacestr+" ARG=";

  for(unsigned k=0;k<weights.size();++k) {
    std::string num; Tools::convert( k+1, num );
    if( neigv[k]==0 ) {
      cvs_str += getShortcutLabel() + "_pdist-" + num + ",";
    } else {
      for(unsigned i=0;i<neigv[k];++i) {
        std::string eignum; Tools::convert( i+1, eignum );
        cvs_str +=  getShortcutLabel() + "_proj" + eignum + "-" + num + ",";
      }
      if( resid[k] ) {
        cvs_str +=  getShortcutLabel() + "_resid-" + num + ",";
      }
    }
  }
  cvs_str += getShortcutLabel() + "_wtheight";
  readInputLine( cvs_str ); 

  // Complete setup of the well tempered weights
  std::vector<std::string> args(1); args[0] = getShortcutLabel();
  ReweightBase* rwb = plumed.getActionSet().selectWithLabel<ReweightBase*>( getShortcutLabel() + "_wtfact" );
  plumed_assert( rwb ); rwb->setArguments( args );
}

}
}
