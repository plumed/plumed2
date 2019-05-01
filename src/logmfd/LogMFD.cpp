/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2015-2018
National Institute of Advanced Industrial Science and Technology (AIST), Japan.
This file contains module for LogMFD method proposed by Tetsuya Morishita(AIST).
then written by Mizuho information and Research Institute, Inc.

The LogMFD module is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The LogMFD module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

//+PLUMEDOC LOGMFDMOD_BIAS LOGMFD
/*
Used to perform LogMFD, LogPD, and TAMD/d-AFED.

\section LogMFD LogMFD

Consider a physical system of \f$N_q\f$ particles, for which the Hamiltonian is given as

\f[
  {H_{\rm MD}}\left( {\bf{\Gamma}} \right) = \sum\limits_{j = 1}^{{N_q}} {\frac{{{\bf{p}}_j^2}}{{2{m_j}}}}  + \Phi \left( {\bf{q}} \right)
\f]

where \f${\bf q}_j\f$, \f${\bf p}_j\f$ (\f$\bf{\Gamma}={\bf q},{\bf p}\f$), and \f$m_j\f$ are the position, momentum, and mass of particle \f$j\f$ respectively,
and \f$\Phi\f$ is the potential energy function for \f${\bf q}\f$.
The free energy \f$F({\bf X})\f$ as a function of a set of \f$N\f$ collective variables (CVs) is given as

\f{eqnarray*}{
  F\left( {{\bf X}} \right) &=&  - {k_B}T\log \int {\exp \left[ { - \beta {H_{\rm MD}}} \right]\prod\limits_{i = 1}^N {\delta \left( {{s_i}\left( {{\bf q}} \right) - {X_i}} \right)} d{\bf{\Gamma }}} \\
  &\simeq&  - {k_B}T\log \int {\exp \left[ { - \beta \left\{ {{H_{\rm MD}} + \sum\limits_i^N {\frac{{{K_i}}}{2}{{\left( {{s_i}\left( {{\bf q}} \right) - {X_i}} \right)}^2}} } \right\}} \right]} d{\bf{\Gamma }}
\f}

where \f$\bf{s}\f$ are the CVs, \f$k_B\f$ is Boltzmann constant, \f$\beta=k_BT\f$,
and \f$K_i\f$ is the spring constant which is large enough to invoke

\f[
 \delta \left( x \right) = \lim_{k \to \infty } \sqrt {\beta k/2\pi} \exp \left( -\beta kx^2/2 \right)
\f]

In mean-force dynamics, \f${\bf X}\f$ are treated as fictitious dynamical variables, which are associated with the following Hamiltonian,

\f[
 {H_{\rm log}} = \sum\limits_{i = 1}^N {\frac{{P_{{X_i}}^2}}{{2{M_i}}} + \Psi \left( {{\bf X}} \right)}
\f]

where \f${P_{{X_i}}}\f$  and \f$M_i\f$ are the momentum and mass of \f$X_i\f$, respectively, and \f$\Psi\f$ is the potential function for \f${\bf X}\f$.
We assume that \f$\Psi\f$ is a functional of \f$F({\bf X})\f$; \f$Ψ[F({\bf X})]\f$. The simplest form of \f$\Psi\f$ is \f$\Psi = F({\bf X})\f$,
which corresponds to TAMD/d-AFED \cite AbramsJ2008, \cite Maragliano2006 (or the extended Lagrangian dynamics in the limit of the adiabatic decoupling between \f$\bf{q}\f$ and \f${\bf X}\f$).
 In the case of LogMFD, the following form of \f$\Psi\f$ is introduced \cite MorishitaLogMFD;


\f[
  {\Psi _{\rm log}}\left( {{\bf X}} \right) = \gamma \log \left[ {\alpha F\left( {{\bf X}} \right) + 1} \right]
\f]

where \f$\alpha\f$ (ALPHA) and \f$\gamma\f$ (GAMMA) are positive parameters. The logarithmic form of \f$\Psi_{\rm log}\f$ ensures the dynamics of \f${\bf X}\f$ on a much smoother energy surface [i.e., \f$\Psi_{\rm log}({\bf X})\f$] than \f$F({\bf X})\f$, thus enhancing the sampling in the \f${\bf X}\f$-space. The parameters \f$\alpha\f$ and \f$\gamma\f$ determine the degree of flatness of \f$\Psi_{\rm log}\f$, but adjusting only \f$\alpha\f$ is normally sufficient to have a relatively flat surface (with keeping the relation \f$\gamma=1/\alpha\f$).

The equation of motion for \f$X_i\f$ in LogMFD (no thermostat) is

\f[
 {M_i}{\ddot X_i} =  - \left( {\frac{{\alpha \gamma }}{{\alpha F + 1}}} \right)\frac{{\partial F}}{{\partial {X_i}}}
\f]

where \f$-\partial F/\partial X_i\f$  is evaluated as a canonical average under the condition that \f${\bf X}\f$ is fixed;

\f{eqnarray*}{
 - \frac{{\partial F}}{{\partial {X_i}}} &\simeq& \frac{1}{Z}\int {{K_i}\left( {{s_i}\left( {{\bf q}} \right) - {X_i}} \right)\exp \left[ { - \beta \left\{ {{H_{\rm MD}} + \sum\limits_i^N {\frac{{{K_i}}}{2}{{\left( {{s_i}\left( {{\bf q}} \right) - {X_i}} \right)}^2}} } \right\}} \right]} d{\bf{\Gamma }}\\
 &\equiv& {\left\langle {{K_i}\left( {{s_i}\left( {{\bf q}} \right) - {X_i}} \right)} \right\rangle _{{\bf X}}}
\f}

where

\f[
 Z = \int {\exp \left[ { - \beta \left\{ {{H_{\rm MD}} + \sum\limits_i^N {\frac{{{K_i}}}{2}{{\left( {{s_i}\left( {{\bf q}} \right) - {X_i}} \right)}^2}} } \right\}} \right]} d{\bf{\Gamma }}
\f]

The mean-force (MF) is practically evaluated by performing a shot-time canonical MD run each time \f${\bf X}\f$ is updated according to the equation of motion for \f${\bf X}\f$.

If the canonical average for the MF is effectively converged, the dynamical variables \f${\bf q}\f$ and \f${\bf X}\f$ are decoupled and they evolve adiabatically, which can be exploited for the on-the-fly evaluation of \f$F({\bf X})\f$. I.e., \f$H_{\rm log}\f$ should be a constant of motion in this case, thus \f$F({\bf X})\f$ can be evaluated each time \f${\bf X}\f$ is updated as


\f[
 F\left( {{{\bf X}}\left( t \right)} \right) = \frac{1}{\alpha} \left[
  \exp \frac{1}{\gamma} \left\{ \left( H_{\rm log} - \sum_i \frac{P_{X_i}^2}{2M_i} \right) \right\} - 1 \right]
\f]


This means that \f$F({\bf X})\f$ can be constructed without post processing (on-the-fly free energy reconstruction). Note that the on-the-fly free energy reconstruction is also possible in TAMD/d-AFED if the Hamiltonian-like conserved quantity is available (e.g., the Nose-Hoover type dynamics).



\section LogPD LogPD


The accuracy in the MF is critical to the on-the-fly free energy reconstruction. To improve the evaluation of the MF, parallel-dynamics (PD) is incorporated into LogMFD, leading to logarithmic parallel-dynamics (LogPD) \cite MorishitaLogPD.


In PD, the MF is evaluated by a non-equilibrium path-ensemble based on the Crooks-Jarzynski non-equilibrium work relation. To this end, multiple replicas of the MD system which run in parallel are introduced. The CVs [\f${\bf s}({\bf q})\f$] in each replica is restrained to the same value of \f${\bf X}(t)\f$. A canonical MD run with \f$N_m\f$ steps is performed in each replica, then the MF on \f$X_i\f$ is evaluated using the MD trajectories from all replicas.
The MF is practically calculated as


\f[
 - \frac{{\partial F}}{{\partial {X_i}}} = \sum\limits_{k = 1}^{{N_r}} {{W_k}} \sum\limits_{n = 1}^{{N_m}} {\frac{1}{{{N_m}}}{K_i}\left[ {{s_i}\left( {{{\bf q}}_n^k} \right) - {X_i}} \right]}
\f]



where \f${\bf q}^k_n\f$  indicates the \f${\bf q}\f$-configuration at time step \f$n\f$ in the \f$N_m\f$-step shot-time MD run in the \f$k\f$th replica among \f$N_r\f$ replicas. \f$W_k\f$ is given as

\f[
 {W_k} = \frac{{{e^{ - \beta {w_k}\left( t \right)}}}}{{\sum\limits_{k=1}^{{N_r}} {{e^{ - \beta {w_k}\left( t \right)}}} }}
\f]


where


\f[
 {w_k}\left( t \right) = \int\limits_{{X_0}}^{X\left( t \right)} {\sum\limits_{i=1}^N {\frac{{\partial H_{\rm MD}^k}}{{\partial {X_i}}}d{X_i}} }
\f]

\f[
 H_{\rm MD}^k\left( {{\bf{\Gamma }},{{\bf X}}} \right) = {H_{\rm MD}}\left( {{{\bf{\Gamma }}^k}} \right) + \sum\limits_{i = 1}^N {\frac{{{K_i}}}{2}{{\left( {s_i^k - {X_i}} \right)}^2}}
\f]

and \f$s^k_i\f$ is the \f$i\f$th CV in the \f$k\f$th replica.

\f$W_k\f$ comes from the Crooks-Jarzynski non-equilibrium work relation by which we can evaluate an equilibrium ensemble average from a set of non-equilibrium trajectories. Note that, to avoid possible numerical errors in the exponential function, the following form of \f$W_k\f$ is instead used in PLUMED,

\f[
 {W_k}\left( t \right) = \frac{{\exp \left[ { - \beta \left\{ {{w_k}\left( t \right) - {w_{\min }}\left( t \right)} \right\}} \right]}}{{\sum\nolimits_k {\exp \left[ { - \beta \left\{ {{w_k}\left( t \right) - {w_{\min }}\left( t \right)} \right\}} \right]} }}
\f]

where

\f[
 {w_{\min }}\left( t \right) = {\rm Min}_k\left[ {{w_k}\left( t \right)} \right]
\f]


With the MF evaluated using the PD approach, reconstructing free energy profiles can be performed more efficiently (requiring less elapsed computing time) in LogPD than with a single MD system in LogMFD. In the case that there exists more than one stable state separated by high energy barriers in the hidden subspace orthogonal to the CV-subspace, LogPD is particularly of use to incorporate all the contributions from such hidden states with appropriate weights (in the limit \f$N_r\to\infty\f$ ).

Note that LogPD calculations should always be initiated with an equilibrium \f${\bf q}\f$-configuration in each replica, because the Crooks-Jarzynski non-equilibrium work relation is invoked. Also note that LogPD is currently available only with Gromacs, while LogMFD can be performed with LAMMPS, Gromacs, and NAMD.

\section Thermostat Using LogMFD with a thermostat /PD

Introducing a thermostat on \f${\bf X}\f$ is often recommended in LogMFD/PD to maintain the adiabatic decoupling between \f${\bf q}\f$ and \f${\bf X}\f$. In the framework of the LogMFD approach, the Nose-Hoover type thermostat and the Gaussian isokinetic (velocity scaling) thermostat can be used to control the kinetic energy of \f${\bf X}\f$.

\subsection Nose-Hoover Nose-Hoover LogMFD/PD

The equation of motion for \f$X_i\f$ coupled to a Nose-Hoover thermostat variable \f$\eta\f$ (single heat bath) is

\f[
 {M_i}{\ddot X_i} =  - \left( {\frac{{\alpha \gamma }}{{\alpha F + 1}}} \right)\frac{{\partial F}}{{\partial {X_i}}} - {M_i}{\dot X_i}\dot \eta
\f]

The equation of motion for \f$\eta\f$ is

\f[
 Q\ddot \eta  = \sum\limits_{i = 1}^N {\frac{{P_{{X_i}}^2}}{{{M_i}}} - N{k_B}T}
\f]

where \f$N\f$ is the number of the CVs. Since the following pseudo-Hamiltonian is a constant of motion in Nose-Hoover LogMFD/PD,

\f[
 H_{\rm log}^{\rm NH} = \sum\limits_{i = 1}^N {\frac{{P_{{X_i}}^2}}{{2{M_i}}} + \gamma \log \left[ {\alpha F\left( {{\bf X}} \right) + 1} \right] + \frac{1}{2}Q{{\dot \eta }^2} + N{k_B}T\eta}
\f]

\f$F({\bf X}(t))\f$ is obtained at each MFD step as

\f[
 F\left( {{{\bf X}}\left( t \right)} \right) = \frac{1}{\alpha }\left[ {\exp \left\{ {{{ \frac{1}{\gamma} \left( {H_{\rm log}^{\rm NH} - \sum_i {\frac{{P_{{X_i}}^2}}{{2{M_i}}}}  - \frac{1}{2}Q{{\dot \eta }^2} - N{k_B}T\eta} \right)}  }} \right\} - 1} \right]
\f]



\subsection VS Velocity scaling LogMFD/PD

The velocity scaling algorithm (which is equivalent to the Gaussian isokinetic dynamics in the limit \f$\Delta t\to 0\f$) can also be employed to control the velocity of \f${\bf X}\f$, \f$\bf{V}_x\f$.

The following algorithm is introduced to perform isokinetic LogMFD calculations \cite MorishitaVsLogMFD;

\f{eqnarray*}{
{V_{{X_i}}}\left( {{t_{n + 1}}} \right) &=& V_{X_i}^\prime \left( {{t_n}} \right) + \Delta t\left[ { - \left( {\frac{{\alpha \gamma }}{{\alpha F\left( {{t_n}} \right) + 1}}} \right)\frac{{\partial F\left( {{t_n}} \right)}}{{\partial {X_i}}}} \right]\\
S\left( {{t_{n + 1}}} \right) &=& \sqrt {\frac{{N{k_B}T}}{{\sum\limits_i {{M_i}V_{{X_i}}^2\left( {{t_{n + 1}}} \right)} }}} \\
{V_{{X_i}}}^\prime \left( {{t_{n + 1}}} \right) &=& S\left( {{t_{n + 1}}} \right){V_{{X_i}}}\left( {{t_{n + 1}}} \right)\\
{X_i}\left( {{t_{n + 1}}} \right) &=& {X_i}\left( {{t_n}} \right) + \Delta t V_{X_i}^\prime \left( {{t_{n + 1}}} \right)\\
{\Psi_{\rm log}}\left( {{t_{n + 1}}} \right) &=& N{k_B}T\log S\left( {{t_{n + 1}}} \right) + {\Psi_{\rm log}}\left( {{t_n}} \right)\\
F\left( {{t_{n + 1}}} \right) &=& \frac{1}{\alpha} \left[
    \exp \left\{ \Psi_{\rm log} \left( t_{n+1} \right) / \gamma \right\} - 1 \right]
\f}

Note that \f$V_{X_i}^\prime\left( {{t_0}} \right)\f$ is assumed to be initially given, which meets the following relation,

\f[
  \sum\limits_{i = 1}^N M_i V_{X_i}^{\prime 2} \left( t_0 \right)  = N{k_B}{T}
\f]

It should be stressed that, in the same way as in the NVE and Nose-Hoover LogMFD/PD, \f$F({\bf X}(t))\f$ can be evaluated at each MFD step (on-the-fly free energy reconstruction) in Velocity Scaling LogMFD/PD.


\par Examples
\section Examples Examples

\subsection Example-LoGMFD Example of LogMFD

The following input file tells plumed to restrain collective variables
to the fictitious dynamical variables in LogMFD/PD.

plumed.dat
\plumedfile
UNITS TIME=fs LENGTH=1.0 ENERGY=kcal/mol Mass=1.0 CHARGE=1.0
phi: TORSION ATOMS=5,7,9,15
psi: TORSION ATOMS=7,9,15,17

# LogMFD
LOGMFD ...
LABEL=logmfd
ARG=phi,psi
KAPPA=100.0,100.0
DELTA_T=0.5
INTERVAL=500
TEMP=300.0
FLOG=5.0
MFICT=5000000.0,5000000.0
# FICT=0.8,0.8
VFICT=3.5e-4,3.5e-4
ALPHA=4.0
THERMOSTAT=NVE
VETA=0.0
META=20000.0
FICT_MAX=3.1,3.1
FICT_MIN=-3.1,-3.1
... LOGMFD
\endplumedfile

To submit this simulation with Gromacs, use the following command line
to execute a LogMFD run with Gromacs-MD.
Here TOPO/topol0.tpr is an input file
which contains atomic coordinates and Gromacs parameters.

\plumedfile
gmx_mpi mdrun -s TOPO/topol0.tpr -plumed
\endplumedfile

This command will output files named logmfd.out and replica.out.

The output file logmfd.out records free energy and all fictitious dynamical variables at every MFD step.

logmfd.out

\plumedfile
# LogMFD
# CVs : phi psi
# Mass for CV particles : 5000000.000000000 5000000.000000000
# Mass for thermostat   :   20000.000000000
# 1:iter_md, 2:Flog, 3:2*Ekin/gkb[K], 4:eta, 5:Veta,
# 6:phi_fict(t), 7:phi_vfict(t), 8:phi_force(t),
# 9:psi_fict(t), 10:psi_vfict(t), 11:psi_force(t),
       1       4.99918574     308.24149708       0.00000000       0.00000000      -2.85605938       0.00035002       5.19074544       2.79216364       0.00035000      -0.53762989
       2       4.99836196     308.26124159       0.00000000       0.00000000      -2.85588436       0.00035005       4.71247605       2.79233863       0.00035000      -0.00532474
       3       4.99743572     308.28344595       0.00000000       0.00000000      -2.85570932       0.00035007       5.34358230       2.79251363       0.00035000      -0.05119816
...
\endplumedfile

The output file replica.out records all collective variables at every MFD step.

replica.out

\plumedfile
# Replica No. 0 of 1.
# 1:iter_md, 2:work, 3:weight,
# 4:phi(q)
# 5:psi(q)
       1   -8.142952e-04     1.000000e+00       -2.80432694       2.78661234
       2   -1.638105e-03     1.000000e+00       -2.80893462       2.79211039
       3   -2.564398e-03     1.000000e+00       -2.80244854       2.79182665
...
\endplumedfile

\subsection Example-LogPD Example of LogPD

Use the following command line to execute a LogPD run using two MD replicas (note that only Gromacs is currently available for LogPD).
Here TOPO/topol0.tpr and TOPO/topol1.tpr are input files
which contain atomic coordinates of each replica and Gromacs parameters.

\plumedfile
mpirun -np 2 gmx_mpi mdrun -s TOPO/topol -plumed -multi 2
\endplumedfile

This command will output files named logmfd.out, replica.out.0 and replica.out.1.

The output file logmfd.out records free energy and all fictitious dynamical variables at every MFD step.

logmfd.out

\plumedfile
# LogPD, replica parallel of LogMFD
# number of replica : 2
# CVs : phi psi
# Mass for CV particles : 5000000.000000000 5000000.000000000
# Mass for thermostat   :   20000.000000000
# 1:iter_md, 2:Flog, 3:2*Ekin/gkb[K], 4:eta, 5:Veta,
# 6:phi_fict(t), 7:phi_vfict(t), 8:phi_force(t),
# 9:psi_fict(t), 10:psi_vfict(t), 11:psi_force(t),
       1       5.00224715     308.16814691       0.00000000       0.00000000      -0.95937173       0.00034994     -12.91277494       0.78923967       0.00035000       0.07353010
       2       5.00476934     308.10774854       0.00000000       0.00000000      -0.95919679       0.00034989     -11.20093553       0.78941467       0.00034999      -3.21098229
       3       5.00702463     308.05376594       0.00000000       0.00000000      -0.95902187       0.00034983     -10.81712171       0.78958965       0.00034998      -2.07196718
...
\endplumedfile


The output file replica.out.0 records all collective variables of replica No.0 at every MFD step.

replica.out.0

\plumedfile
# Replica No. 0 of 2.
# 1:iter_md, 2:work, 3:weight,
# 4:phi(q)
# 5:psi(q)
       1    1.843110e-03     5.003389e-01       -1.10929125       0.83348865
       2    3.466179e-03     5.010942e-01       -1.05020764       0.78731283
       3    4.927870e-03     5.017619e-01       -1.04968867       0.79635198
...
\endplumedfile

The output file replica.out.1 records all collective variables of replica No.1 at every MFD step.

replica.out.1

\plumedfile
# Replica No. 1 of 2.
# 1:iter_md, 2:work, 3:weight,
# 4:phi(q)
# 5:psi(q)
       1    2.651173e-03     4.996611e-01       -1.06802968       0.74605205
       2    6.075530e-03     4.989058e-01       -1.09264741       0.72681448
       3    9.129358e-03     4.982381e-01       -1.08517238       0.74084241
...
\endplumedfile

*/
//+ENDPLUMEDOC

#include "bias/Bias.h"
#include "core/ActionRegister.h"
#include "core/Atoms.h"
#include "core/PlumedMain.h"

#include <iostream>

using namespace std;
using namespace PLMD;
using namespace bias;

namespace PLMD {
namespace logmfd {
/**
   \brief class for LogMFD parameters, variables and subroutines.
 */
class LogMFD : public Bias {
  bool firsttime;               ///< flag that indicates first MFD step or not.
  int    interval;              ///< input parameter, period of MD steps when fictitious dynamical variables are evolved.
  double delta_t;               ///< input parameter, one time step of MFD when fictitious dynamical variables are evolved.
  string thermostat;            ///< input parameter, type of thermostat for canonical dyanamics.
  double kbt;                   ///< k_B*temerature

  int    TAMD;                  ///< input parameter, perform TAMD instead of LogMFD.
  double alpha;                 ///< input parameter, alpha parameter for LogMFD.
  double gamma;                 ///< input parameter, gamma parameter for LogMFD.
  std::vector<double> kappa;    ///< input parameter, strength of the harmonic restraining potential.

  std::vector<double> fict_max; ///< input parameter, maximum of each fictitous dynamical variable.
  std::vector<double> fict_min; ///< input parameter, minimum of each fictitous dynamical variable.

  std::vector<double>  fict;    ///< current values of each fictitous dynamical variable.
  std::vector<double> vfict;    ///< current velocity of each fictitous dynamical variable.
  std::vector<double> mfict;    ///< mass of each fictitous dynamical variable.

  double xeta;                  ///< current eta variable of thermostat.
  double veta;                  ///< current velocity of eta variable of thermostat.
  double meta;                  ///< mass of eta variable of thermostat.

  double flog;                  ///< current free energy

  double hlog;                  ///< value invariant
  double phivs;                 ///< potential used in VS method
  double work;                  ///< current works done by fictitious dynamical variables in this replica.
  double weight;                ///< current weight of this replica.

  std::vector<double> ffict;    ///< current force of each fictitous dynamical variable.
  std::vector<double> fict_ave; ///< averaged values of each collective variable.

  std::vector<Value*>  fictValue; ///< pointers to fictitious dynamical variables
  std::vector<Value*> vfictValue; ///< pointers to velocity of fictitious dynamical variables

public:
  static void registerKeywords(Keywords& keys);

  explicit LogMFD(const ActionOptions&);
  void calculate();
  void update();
  void updateNVE();
  void updateNVT();
  void updateVS();
  void calcMeanForce();
  double calcEkin();
  double calcFlog();
  double calcClog();

private:
  double sgn( double x ) {
    return x>0.0 ? 1.0 : x<0.0 ? -1.0 : 0.0;
  }
};

PLUMED_REGISTER_ACTION(LogMFD,"LOGMFD")

/**
   \brief instruction of parameters for Plumed manual.
*/
void LogMFD::registerKeywords(Keywords& keys) {
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","INTERVAL",
           "Period of MD steps (\\f$N_m\\f$) to update fictitious dynamical variables." );
  keys.add("compulsory","DELTA_T",
           "Time step for the fictitious dynamical variables (MFD step)." );
  keys.add("compulsory","THERMOSTAT",
           "Type of thermostat for the fictitious dynamical variables. NVE, NVT, VS are available." );
  keys.add("optional","TEMP",
           "Temperature of the fictitious dynamical variables in LogMFD/PD thermostat. "
           "If not provided or provided as 0, it will be taken from the temperature of the MD system." );

  keys.add("optional","TAMD",
           "When TAMD=1, TAMD/d-AFED calculations can be performed instead of LogMFD. Otherwise, the LogMFD protocol is switched on (default)." );

  keys.add("optional","ALPHA",
           "Alpha parameter for LogMFD. "
           "If not provided or provided as 0, it will be taken as 1/gamma. "
           "If gamma is also not provided, Alpha is set as 4, which is a sensible value when the unit of kcal/mol is used." );
  keys.add("optional","GAMMA",
           "Gamma parameter for LogMFD. "
           "If not provided or provided as 0, it will be taken as 1/alpha. "
           "If alpha is also not provided, Gamma is set as 0.25, which is a sensible value when the unit of kcal/mol is used." );
  keys.add("compulsory","KAPPA",
           "Spring constant of the harmonic restraining potential for the fictitious dynamical variables." );

  keys.add("compulsory","FICT_MAX",
           "Maximum values reachable for the fictitious dynamical variables. The variables will elastically bounce back at the boundary (mirror boundary)." );
  keys.add("compulsory","FICT_MIN",
           "Minimum values reachable for the fictitious dynamical variables. The variables will elastically bounce back at the boundary (mirror boundary)." );

  keys.add("optional","FICT",
           "The initial values of the fictitious dynamical variables. "
           "If not provided, they are set equal to their corresponding CVs for the initial atomic configuration." );
  keys.add("optional","VFICT",
           "The initial velocities of the fictitious dynamical variables. "
           "If not provided, they will be taken as 0." );
  keys.add("optional","MFICT",
           "Masses of each fictitious dynamical variable. "
           "If not provided, they will be taken as 10000." );

  keys.add("optional","XETA",
           "The initial eta variable of the Nose-Hoover thermostat "
           "for the fictitious dynamical variables. "
           "If not provided, it will be taken as 0." );
  keys.add("optional","VETA",
           "The initial velocity of eta variable. "
           "If not provided, it will be taken as 0." );
  keys.add("optional","META",
           "Mass of eta variable. "
           "If not provided, it will be taken as \\f$N*kb*T*100*100\\f$." );

  keys.add("compulsory","FLOG",
           "The initial free energy value in the LogMFD/PD run."
           "The origin of the free energy profile is adjusted by FLOG to "
           "realize \\f$F({\\bf X}(t)) > 0\\f$ at any \\f${\\bf X}(t)\\f$, "
           "resulting in enhanced barrier-crossing. "
           "(The value of \\f$H_{\\rm log}\\f$ is automatically "
           "set according to FLOG).");

  keys.add("optional","WORK",
           "The initial value of the work done by fictitious dynamical "
           "variables in each replica. "
           "If not provided, it will be taken as 0.");

  componentsAreNotOptional(keys);
  keys.addOutputComponent("_fict","default",
                          "For example, the fictitious collective variable for LogMFD is specified as "
                          "ARG=dist12 and LABEL=logmfd in LOGMFD section in Plumed input file, "
                          "the associated fictitious dynamical variable can be specified as "
                          "PRINT ARG=dist12,logmfd.dist12_fict FILE=COLVAR");
  keys.addOutputComponent("_vfict","default",
                          "For example, the fictitious collective variable for LogMFD is specified as "
                          "ARG=dist12 and LABEL=logmfd in LOGMFD section in Plumed input file, the "
                          "velocity of the associated fictitious dynamical variable can be specified as "
                          "PRINT ARG=dist12,logmfd.dist12_vfict FILE=COLVAR");
}


/**
   \brief constructor of LogMFD class
   \details This constructor initializes all parameters and variables,
   reads input file and set value of parameters and initial value of variables,
   and writes messages as Plumed log.
*/
LogMFD::LogMFD( const ActionOptions& ao ):
  PLUMED_BIAS_INIT(ao),
  firsttime(true),
  interval(10),
  delta_t(1.0),
  thermostat("NVE"),
  kbt(-1.0),
  TAMD(0),
  alpha(0.0),
  gamma(0.0),
  kappa(getNumberOfArguments(),0.0),
  fict_max(getNumberOfArguments(),0.0),
  fict_min(getNumberOfArguments(),0.0),
  fict (getNumberOfArguments(),0.0),
  vfict(getNumberOfArguments(),0.0),
  mfict(getNumberOfArguments(),10000.0),
  xeta(0.0),
  veta(0.0),
  meta(0.0),
  flog(0.0),
  hlog(0.0),
  phivs(0.0),
  work(0.0),
  weight(0.0),
  ffict(getNumberOfArguments(),0.0),
  fict_ave(getNumberOfArguments(),0.0),
  fictValue(getNumberOfArguments(),NULL),
  vfictValue(getNumberOfArguments(),NULL)
{
  parse("INTERVAL",interval);
  parse("DELTA_T",delta_t);
  parse("THERMOSTAT",thermostat);
  parse("TEMP",kbt); // read as temperature

  parse("TAMD",TAMD);
  parse("ALPHA",alpha);
  parse("GAMMA",gamma);
  parseVector("KAPPA",kappa);

  parseVector("FICT_MAX",fict_max);
  parseVector("FICT_MIN",fict_min);

  parseVector("FICT",fict);
  parseVector("VFICT",vfict);
  parseVector("MFICT",mfict);

  parse("XETA",xeta);
  parse("VETA",veta);
  parse("META",meta);

  parse("FLOG",flog);

  // read initial value of work for each replica of LogPD
  if( multi_sim_comm.Get_size()>1 ) {
    vector<double> vwork(multi_sim_comm.Get_size(),0.0);
    parseVector("WORK",vwork);
    // initial work of this replica
    work = vwork[multi_sim_comm.Get_rank()];
  }
  else {
    work = 0.0;
  }

  if( kbt>=0.0 ) {
    kbt *= plumed.getAtoms().getKBoltzmann();
  }
  else {
    kbt = plumed.getAtoms().getKbT();
  }

  if( meta == 0.0 ) {
    const double nkt = getNumberOfArguments()*kbt;
    meta = nkt*100.0*100.0;
  }

  if(alpha == 0.0 && gamma == 0.0) {
    alpha = 4.0;
    gamma = 1 / alpha;
  }
  else if(alpha != 0.0 && gamma == 0.0) {
    gamma = 1 / alpha;
  }
  else if(alpha == 0.0 && gamma != 0.0) {
    alpha = 1 / gamma;
  }

  checkRead();

  // output messaages to Plumed's log file
  if( multi_sim_comm.Get_size()>1 ) {
    if( TAMD ) {
      log.printf("TAMD-PD, replica parallel of TAMD, no logarithmic flattening.\n");
    }
    else {
      log.printf("LogPD, replica parallel of LogMFD.\n");
    }
    log.printf("number of replica : %d.\n", multi_sim_comm.Get_size() );
  }
  else {
    if( TAMD ) {
      log.printf("TAMD, no logarithmic flattening.\n");
    }
    else {
      log.printf("LogMFD, logarithmic flattening.\n");
    }
  }

  log.printf("  with harmonic force constant      ");
  for(unsigned i=0; i<kappa.size(); i++) log.printf(" %f",kappa[i]);
  log.printf("\n");

  log.printf("  with interval of cv(ideal) update ");
  log.printf(" %d", interval);
  log.printf("\n");

  log.printf("  with time step of cv(ideal) update ");
  log.printf(" %f", delta_t);
  log.printf("\n");

  if( !TAMD ) {
    log.printf("  with alpha, gamma                 ");
    log.printf(" %f %f", alpha, gamma);
    log.printf("\n");
  }

  log.printf("  with Thermostat for cv(ideal)     ");
  log.printf(" %s", thermostat.c_str());
  log.printf("\n");

  log.printf("  with initial free energy          ");
  log.printf(" %f", flog);
  log.printf("\n");

  log.printf("  with mass of cv(ideal)");
  for(unsigned i=0; i<mfict.size(); i++) log.printf(" %f", mfict[i]);
  log.printf("\n");

  log.printf("  with initial value of cv(ideal)");
  for(unsigned i=0; i<fict.size(); i++) log.printf(" %f", fict[i]);
  log.printf("\n");

  log.printf("  with initial velocity of cv(ideal)");
  for(unsigned i=0; i<vfict.size(); i++) log.printf(" %f", vfict[i]);
  log.printf("\n");

  log.printf("  with maximum value of cv(ideal)    ");
  for(unsigned i=0; i<fict_max.size(); i++) log.printf(" %f",fict_max[i]);
  log.printf("\n");

  log.printf("  with minimum value of cv(ideal)    ");
  for(unsigned i=0; i<fict_min.size(); i++) log.printf(" %f",fict_min[i]);
  log.printf("\n");

  log.printf("  and kbt                           ");
  log.printf(" %f",kbt);
  log.printf("\n");

  // setup Value* variables
  for(unsigned i=0; i<getNumberOfArguments(); i++) {
    std::string comp = getPntrToArgument(i)->getName()+"_fict";
    addComponentWithDerivatives(comp);

    if(getPntrToArgument(i)->isPeriodic()) {
      std::string a,b;
      getPntrToArgument(i)->getDomain(a,b);
      componentIsPeriodic(comp,a,b);
    }
    else {
      componentIsNotPeriodic(comp);
    }
    fictValue[i] = getPntrToComponent(comp);

    comp = getPntrToArgument(i)->getName()+"_vfict";
    addComponent(comp);

    componentIsNotPeriodic(comp);
    vfictValue[i] = getPntrToComponent(comp);
  }
}

/**
   \brief calculate forces for fictitious variables at every MD steps.
   \details This function calculates initial values of fictitious variables
   and write header messages to LogMFD log files at the first MFD step,
   calculates restraining fources comes from difference between the fictious variable
   and collective variable at every MD steps.
*/
void LogMFD::calculate() {
  if( firsttime ) {
    firsttime = false;

    // set initial values of fictitious variables if they were not specified.
    for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      if( fict[i] != 0.0 ) continue;

      // use the collective variables as the initial of the fictitious variable.
      fict[i] = getArgument(i);

      // average values of fictitious variables by all replica.
      if( multi_sim_comm.Get_size()>1 ) {
        multi_sim_comm.Sum(fict[i]);
        fict[i] /= multi_sim_comm.Get_size();
      }
    }

    // initialize accumulation value to zero
    for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      fict_ave[i] = 0.0;
    }

    // calculate invariant for NVE
    if(thermostat == "NVE") {
      // kinetic energy
      const double ekin = calcEkin();
      // potential energy
      const double pot = TAMD ? flog : sgn(flog)*gamma * std::log1p( alpha*fabs(flog) );
      // invariant
      hlog = pot + ekin;
    }
    else if(thermostat == "NVT") {
      const double nkt = getNumberOfArguments()*kbt;
      // kinetic energy
      const double ekin = calcEkin();
      // bath energy
      const double ekin_bath = 0.5*veta*veta*meta + xeta*nkt;
      // potential energy
      const double pot = TAMD ? flog : sgn(flog)*gamma * std::log1p( alpha*fabs(flog) );
      // invariant
      hlog = pot + ekin + ekin_bath;
    }
    else if(thermostat == "VS") {
      // initial velocities
      for(unsigned i=0; i<getNumberOfArguments(); ++i) {
        vfict[i] = sqrt(kbt/mfict[i]);
      }
      // initial VS potential
      phivs = TAMD ? flog : sgn(flog)* gamma*std::log1p( alpha*fabs(flog) );

      // invariant
      hlog = 0.0;
    }

    weight = 1.0; // for replica parallel

    // open LogMFD's log file
    if( multi_sim_comm.Get_rank()==0 && comm.Get_rank()==0 ) {
      FILE *outlog = std::fopen("logmfd.out", "w");

      // output messages to LogMFD's log file
      if( multi_sim_comm.Get_size()>1 ) {
        fprintf(outlog, "# LogPD, replica parallel of LogMFD\n");
        fprintf(outlog, "# number of replica : %d\n", multi_sim_comm.Get_size() );
      }
      else {
        fprintf(outlog, "# LogMFD\n");
      }

      fprintf(outlog, "# CVs :");
      for(unsigned i=0; i<getNumberOfArguments(); ++i) {
        fprintf(outlog, " %s",  getPntrToArgument(i)->getName().c_str() );
      }
      fprintf(outlog, "\n");

      fprintf(outlog, "# Mass for CV particles :");
      for(unsigned i=0; i<getNumberOfArguments(); ++i) {
        fprintf(outlog, "%18.9f", mfict[i]);
      }
      fprintf(outlog, "\n");

      fprintf(outlog, "# Mass for thermostat   :");
      fprintf(outlog, "%18.9f", meta);
      fprintf(outlog, "\n");
      fprintf(outlog, "# 1:iter_md, 2:Flog, 3:2*Ekin/gkb[K], 4:eta, 5:Veta,\n");

      for(unsigned i=0; i<getNumberOfArguments(); ++i) {
        fprintf(outlog, "# %u:%s_fict(t), %u:%s_vfict(t), %u:%s_force(t),\n",
                6+i*3, getPntrToArgument(i)->getName().c_str(),
                7+i*3, getPntrToArgument(i)->getName().c_str(),
                8+i*3, getPntrToArgument(i)->getName().c_str() );
      }

      fclose(outlog);
    }

    if( comm.Get_rank()==0 ) {
      // the number of replica is added to file name to distingwish replica.
      FILE *outlog2 = fopen("replica.out", "w");
      fprintf(outlog2, "# Replica No. %d of %d.\n",
              multi_sim_comm.Get_rank(), multi_sim_comm.Get_size() );

      fprintf(outlog2, "# 1:iter_md, 2:work, 3:weight,\n");
      for(unsigned i=0; i<getNumberOfArguments(); ++i) {
        fprintf(outlog2, "# %u:%s(q)\n",
                4+i, getPntrToArgument(i)->getName().c_str() );
      }
      fclose(outlog2);
    }

    // output messages to Plumed's log file
    //    log.printf("LOGMFD thermostat parameters Xeta Veta Meta");
    //    log.printf(" %f %f %f", xeta, veta, meta);
    //    log.printf("\n");
    //    log.printf("# 1:iter_md, 2:Flog, 3:2*Ekin/gkb[K], 4:eta, 5:Veta,");
    //    log.printf("# 6:X1(t), 7:V1(t), 8:F1(t), 9:X2(t), 10:V2(t), 11:F2(t), ...");

  } // firsttime

  // calculate force for fictitious variable
  double ene=0.0;
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    // difference between fictitious variable and collective variable.
    const double diff = difference(i,fict[i],getArgument(i));
    // restraining force.
    const double f = -kappa[i]*diff;
    setOutputForce(i,f);

    // restraining energy.
    ene += 0.5*kappa[i]*diff*diff;

    // accumulate force, later it will be averaged.
    ffict[i] += -f;

    // accumulate varience of collective variable, later it will be averaged.
    fict_ave[i] += diff;
  }

  setBias(ene);
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    // correct fict so that it is inside [min:max].
    fict[i] = fictValue[i]->bringBackInPbc(fict[i]);
    fictValue[i]->set(fict[i]);
    vfictValue[i]->set(vfict[i]);
  }
} // calculate

/**
   \brief update fictitious variables.
   \details This function manages evolution of fictitious variables.
   This function calculates mean force, updates fictitious variables by one MFD step,
   bounces back variables, updates free energy, and record logs.
*/
void LogMFD::update() {
  if(getStep() == 0 || getStep()%interval != 0 ) return;

  // calc mean force for fictitious variables
  calcMeanForce();

  // update fictitious variables
  if(thermostat == "NVE") {
    updateNVE();
  }
  else if(thermostat == "NVT") {
    updateNVT();
  }
  else if(thermostat == "VS") {
    updateVS();
  }

  // bounce back variables if they are beyond their min and max.
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    if(fict[i] > fict_max[i]) {
      fict[i] = fict_max[i] - (fict[i] - fict_max[i]);
      vfict[i] *= -1.0;
    }
    if(fict[i] < fict_min[i]) {
      fict[i] = fict_min[i] + (fict_min[i] - fict[i]);
      vfict[i] *= -1.0;
    }
  }

  // update free energy
  flog = calcFlog();

  // record log for fictitious variables
  if( multi_sim_comm.Get_rank()==0 && comm.Get_rank()==0 ) {
    FILE *outlog = std::fopen("logmfd.out", "a");

    const double ekin = calcEkin();
    const double temp = 2.0*ekin/getNumberOfArguments()/plumed.getAtoms().getKBoltzmann();

    fprintf(outlog, "%*d", 8, (int)getStep()/interval);
    fprintf(outlog, "%17.8f", flog);
    fprintf(outlog, "%17.8f", temp);
    fprintf(outlog, "%17.8f", xeta);
    fprintf(outlog, "%17.8f", veta);
    for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      fprintf(outlog, "%17.8f", fict[i]);
      fprintf(outlog, "%17.8f", vfict[i]);
      fprintf(outlog, "%17.8f", ffict[i]);
    }
    fprintf(outlog," \n");
    fclose(outlog);
  }

  // record log for collective variables
  if( comm.Get_rank()==0 ) {
    // the number of replica is added to file name to distingwish replica.
    FILE *outlog2 = fopen("replica.out", "a");
    fprintf(outlog2, "%*d", 8, (int)getStep()/interval);
    fprintf(outlog2, "%16.6e ", work);
    fprintf(outlog2, "%16.6e ", weight);
    for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      fprintf(outlog2, "%17.8f", fict_ave[i]);
    }
    fprintf(outlog2," \n");
    fclose(outlog2);
  }

  // reset mean force
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    ffict[i] = 0.0;
    fict_ave[i] = 0.0;
  }

} // update

/**
   \brief update fictitious variables by NVE mechanics.
   \details This function updates ficitious variables by one NVE-MFD step using mean forces
   and flattening coefficient and free energy.
 */
void LogMFD::updateNVE() {
  const double dt = delta_t;

  // get latest free energy and flattening coefficient
  flog = calcFlog();
  const double clog = calcClog();

  // update all ficitious variables by one MFD step
  for( unsigned i=0; i<getNumberOfArguments(); ++i ) {
    // update velocity (full step)
    vfict[i]+=clog*ffict[i]*dt/mfict[i];
    // update position (full step)
    fict[i]+=vfict[i]*dt;
  }
} // updateNVE

/**
   \brief update fictitious variables by NVT mechanics.
   \details This function updates ficitious variables by one NVT-MFD step using mean forces
   and flattening coefficient and free energy.
 */
void LogMFD::updateNVT() {
  const double dt = delta_t;
  const double nkt = getNumberOfArguments()*kbt;

  // backup vfict
  std::vector<double> vfict_backup(getNumberOfArguments());
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    vfict_backup[i] = vfict[i];
  }

  const int niter=5;
  for(unsigned j=0; j<niter; ++j) {
    // get latest free energy and flattening coefficient
    flog = calcFlog();
    const double clog = calcClog();

    // restore vfict from vfict_backup
    for(unsigned l=0; l<getNumberOfArguments(); ++l) {
      vfict[l] = vfict_backup[l];
    }

    // evolve vfict from vfict_backup by dt/2
    for(unsigned m=0; m<getNumberOfArguments(); ++m) {
      vfict[m] *= exp(-0.25*dt*veta);
      vfict[m] += 0.5*dt*clog*ffict[m]/mfict[m];
      vfict[m] *= exp(-0.25*dt*veta);
    }
  }

  // get latest free energy and flattening coefficient
  flog = calcFlog();
  const double clog = calcClog();

  // evolve vfict by dt/2, and evolve fict by dt
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    vfict[i] *= exp(-0.25*dt*veta);
    vfict[i] += 0.5*dt*clog*ffict[i]/mfict[i];
    vfict[i] *= exp(-0.25*dt*veta);
    fict[i]  += dt*vfict[i];
  }

  // evolve xeta and veta by dt
  xeta += 0.5*dt*veta;
  const double ekin = calcEkin();
  veta += dt*(2.0*ekin-nkt)/meta;
  xeta += 0.5*dt*veta;
} // updateNVT

/**
   \brief update fictitious variables by VS mechanics.
   \details This function updates ficitious variables by one VS-MFD step using mean forces
   and flattening coefficient and free energy.
 */
void LogMFD::updateVS() {
  const double dt = delta_t;
  const double nkt = getNumberOfArguments()*kbt;

  // get latest free energy and flattening coefficient
  flog = calcFlog();
  const double clog = calcClog();

  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    // update velocity (full step)
    vfict[i] += clog*ffict[i]*dt/mfict[i];
  }

  const double ekin = calcEkin();
  const double svs = sqrt(nkt/ekin/2);

  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    // update position (full step)
    vfict[i] *= svs;
    fict[i] += vfict[i]*dt;
  }

  // evolve VS potential
  phivs += nkt*std::log(svs);
} // updateVS

/**
   \brief calculate mean force for fictitious variables.
   \details This function calculates mean forces by averaging forces accumulated during one MFD step,
   update work variables done by fictitious variables by one MFD step,
   calculate weight variable of this replica for LogPD.
*/
void LogMFD::calcMeanForce() {
  // cale temporal mean force for each CV
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    ffict[i] /= interval;
    // average of diff (getArgument(i)-fict[i])
    fict_ave[i] /= interval;
    // average of getArgument(i)
    fict_ave[i] += fict[i];

    // correct fict_ave so that it is inside [min:max].
    fict_ave[i] = fictValue[i]->bringBackInPbc(fict_ave[i]);
  }

  // accumulate work, it was initialized as 0.0
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    work -= ffict[i] * vfict[i] * delta_t; // modified sign
  }

  // for replica parallel
  if( multi_sim_comm.Get_size()>1 ) {
    // find the minimum work among all replicas
    double work_min = work;
    multi_sim_comm.Min(work_min);

    // weight of this replica.
    // here, work is reduced by work_min to avoid all exp(-work/kbt)s disconverge
    if( kbt == 0.0 ) {
      weight = work==work_min ? 1.0 : 0.0;
    }
    else {
      weight = exp(-(work-work_min)/kbt);
    }

    // normalize the weight
    double sum_weight = weight;
    multi_sim_comm.Sum(sum_weight);
    weight /= sum_weight;

    // weighting force of this replica
    for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      ffict[i] *= weight;
    }

    // mean forces of all replica.
    for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      multi_sim_comm.Sum(ffict[i]);
    }
    // now, mean force is obtained.
  }
} // calcMeanForce

/**
   \brief calculate kinetic energy of fictitious variables.
   \retval kinetic energy.
   \details This function calculates sum of kinetic energy of all fictitious variables.
 */
double LogMFD::calcEkin() {
  double ekin=0.0;
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    ekin += mfict[i]*vfict[i]*vfict[i]*0.5;
  }
  return ekin;
} // calcEkin

/**
   \brief calculate free energy of fictitious variables.
   \retval free energy.
   \details This function calculates free energy by using invariant of canonical mechanics.
 */
double LogMFD::calcFlog() {
  const double nkt = getNumberOfArguments()*kbt;
  const double ekin = calcEkin();
  double pot;

  if( false ) {
  }
  else if (thermostat == "NVE") {
    pot = hlog - ekin;
  }
  else if (thermostat == "NVT") {
    const double ekin_bath = 0.5*veta*veta*meta+xeta*nkt;
    pot = hlog - ekin - ekin_bath;
  }
  else if (thermostat == "VS") {
    pot = phivs;
  }
  else {
    pot = 0.0; // never occurs
  }

  return TAMD ? pot : sgn(pot)*expm1(fabs(pot)/gamma)/alpha;
} // calcFlog

/**
   \brief calculate coefficient for flattening.
   \retval flattering coefficient.
   \details This function returns 1.0 for TAMD, flattening coefficient for LogMFD.
 */
double LogMFD::calcClog() {
  return TAMD ? 1.0 : alpha*gamma/(alpha*fabs(flog)+1.0);
} // calcClog

} // logmfd
} // PLMD
