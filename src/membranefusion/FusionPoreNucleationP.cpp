/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2022.

CVs originally developed by the Jochen Hub group from the University of Saarland (Germany)
and adapted and implemented in PLUMED by Ary Lautaro Di Bartolo and Diego Masone from the
National University of Cuyo (Argentina).

The membranefusion module is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The membranefusion module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "colvar/Colvar.h"
#include "core/ActionRegister.h"
#include <cmath>
#ifdef _OPENMP
#if _OPENMP >= 201307
#include <omp.h>
#endif
#endif

namespace PLMD
{
namespace membranefusion
{
//+PLUMEDOC MEMBRANEFUSIONMOD_COLVAR FUSIONPORENUCLEATIONP
/*
A CV for inducing the nucleation of the fusion pore from a hemifusion stalk.

Calculate the collective variable designed by Hub and collaborators \cite Hub2017 and
implemented into PLUMED by Masone and collaborators.
This CV is capable of inducing the nucleation of the fusion pore from a hemifusion stalk.

\f[
\xi_n = \frac{1}{N_{sn}} \sum_{s=0}^{N_{sn}-1} \delta_{sn} (N_{sn}^{(p)})
\f]

Where \f$\xi_n\f$ is the CV, \f$N_{sn}\f$ is the number of slices of the cylinder that make up the CV,
\f$\delta_{sn}\f$ is a continuos function in the interval [0 1] (\f$\delta_{sf} = 0\f$ for no beads in the slice s, and
\f$\delta_{sf} = 1\f$ for 1 or more beads in the slice s) and \f$N_{sf}^{(p)}\f$ accounts for the number of water and
phosphateoxygens beads within the slice s.

\par Examples

This example induces the nucleation of the fusion pore (\f$\xi_n = 1.0\f$) from a hemifusion stalk (\f$\xi_n = 0.2\f$).

\plumedfile

lMem: GROUP ATOMS=1-10752,21505-22728,23953-24420 #All the lower membrane beads.
uMem: GROUP ATOMS=10753-21504,22729-23952,24421-24888 #All the upper membrane beads.
tails: GROUP ATOMS=8-23948:12,12-23952:12,23966-24884:18,23970-24888:18 #All the lipid tails beads (from the lower and upper membrane).
waters: GROUP ATOMS=24889-56490 #All the water beads.
po4: GROUP ATOMS=2-23942:12,23957-24875:18 #All the lipid phosphateoxygens beads.

fusionPoreNucleation: FUSIONPORENUCLEATIONP UMEMBRANE=uMem LMEMBRANE=lMem TAILS=tails WATERS=waters PHOSPHATEOXYGENS=po4 NSMEM=85 NS=45

MOVINGRESTRAINT ...
    ARG=fusionPoreNucleation
    STEP0=0 AT0=0.2 KAPPA0=10000.0
    STEP1=500000 AT1=1.0 KAPPA1=10000.0
...

PRINT ARG=fusionPoreNucleation FILE=COLVAR STRIDE=1

\endplumedfile

*/
//+ENDPLUMEDOC
class fusionPoreNucleationP : public Colvar
{
  std::vector<AtomNumber> UMEM, LMEM, TAILS, WATERS, POXYGENS;
  std::vector<double> NSMEM, DSMEM, HMEM, NS, DS, HCH, RCYL, ZETA, ONEOVERS2C2CUTOFF, XCYL, YCYL;

public:
  explicit fusionPoreNucleationP(const ActionOptions &);
  void calculate() override;
  static void registerKeywords(Keywords &keys);
};

PLUMED_REGISTER_ACTION(fusionPoreNucleationP, "FUSIONPORENUCLEATIONP")

void fusionPoreNucleationP::registerKeywords(Keywords &keys)
{
  Colvar::registerKeywords(keys);
  keys.add("atoms", "UMEMBRANE", "all the beads of the upper membrane.");
  keys.add("atoms", "LMEMBRANE", "all the beads of the lower membrane.");
  keys.add("atoms", "TAILS", "all the tail beads of the system.");
  keys.add("atoms", "WATERS", "all the water beads of the system.");
  keys.add("atoms", "PHOSPHATEOXYGENS", "all the lipid phosphateoxygens beads of the system.");
  keys.add("compulsory", "NSMEM", "the number of slices of the membrane fusion cylinder.");
  keys.add("optional", "DSMEM", "( default=0.1 ) thickness of the slices of the membrane fusion cylinder.");
  keys.add("optional", "HMEM", "( default=0.25 ) parameter of the step function θ(x,h) for the membrane fusion.");
  keys.add("compulsory", "NS", "the number of slices of the membrane-spanning cylinder in such a way that when the bilayers are flat and parallel the CV is equal to 0.2.");
  keys.add("optional", "DS", "( default=0.25 ) thickness of the slices of the membrane-spanning cylinder.");
  keys.add("optional", "HCH", "( default=0.25 ) parameter of the step function θ(x,h) for the CV.");
  keys.add("optional", "RCYL", "( default=0.8 ) the radius of the membrane-spanning cylinder.");
  keys.add("optional", "ZETA", "( default=0.75 ) parameter of the switch function ψ(x,ζ).");
  keys.add("optional", "ONEOVERS2C2CUTOFF", "( default=500 ) cut off large values for the derivative of the atan2 function to avoid violate energy.");
  keys.add("optional", "XCYL", "X coordinate of the fixed cylinder, if not present this will be calculated.");
  keys.add("optional", "YCYL", "X coordinate of the fixed cylinder, if not present this will be calculated.");
}

fusionPoreNucleationP::fusionPoreNucleationP(const ActionOptions &ao) : PLUMED_COLVAR_INIT(ao)
{
  parseAtomList("UMEMBRANE", UMEM);
  if (UMEM.size() == 0)
    error("UMEMBRANE has not any atom specified.");

  parseAtomList("LMEMBRANE", LMEM);
  if (LMEM.size() == 0)
    error("LMEMBRANE has not any atom specified.");

  parseAtomList("TAILS", TAILS);
  if (TAILS.size() == 0)
    error("TAILS has not any atom specified.");

  parseAtomList("WATERS", WATERS);
  if (WATERS.size() == 0)
    error("WATERS has not any atom specified.");

  parseAtomList("PHOSPHATEOXYGENS", POXYGENS);
  if (POXYGENS.size() == 0)
    error("PHOSPHATEOXYGENS has not any atom specified.");

  parseVector("NSMEM", NSMEM);
  if (NSMEM.size() > 1)
    error("NSMEM cannot take more than one value.");

  parseVector("DSMEM", DSMEM);
  if (DSMEM.size() > 1)
    error("DSMEM cannot take more than one value.");
  if (DSMEM.size() == 0)
    DSMEM.push_back(0.1);

  parseVector("HMEM", HMEM);
  if (HMEM.size() > 1)
    error("HMEM cannot take more than one value.");
  if (HMEM.size() == 0)
    HMEM.push_back(0.25);

  parseVector("NS", NS);
  if (NS.size() > 1)
    error("NS cannot take more than one value.");

  parseVector("DS", DS);
  if (DS.size() > 1)
    error("DS cannot take more than one value.");
  if (DS.size() == 0)
    DS.push_back(0.25);

  parseVector("HCH", HCH);
  if (HCH.size() > 1)
    error("H cannot take more than one value.");
  if (HCH.size() == 0)
    HCH.push_back(0.25);

  parseVector("RCYL", RCYL);
  if (RCYL.size() > 1)
    error("RCYL cannot take more than one value.");
  if (RCYL.size() == 0)
    RCYL.push_back(0.8);

  parseVector("ZETA", ZETA);
  if (ZETA.size() > 1)
    error("ZETA cannot take more than one value.");
  if (ZETA.size() == 0)
    ZETA.push_back(0.75);

  parseVector("ONEOVERS2C2CUTOFF", ONEOVERS2C2CUTOFF);
  if (ONEOVERS2C2CUTOFF.size() > 1)
    error("ONEOVERS2C2CUTOFF cannot take more than one value.");
  if (ONEOVERS2C2CUTOFF.size() == 0)
    ONEOVERS2C2CUTOFF.push_back(500);

  parseVector("XCYL", XCYL);
  if (XCYL.size() > 1)
    error("XCYL cannot take more than one value.");
  if (XCYL.size() == 0)
    XCYL.push_back(-1.0);

  parseVector("YCYL", YCYL);
  if (YCYL.size() > 1)
    error("YCYL cannot take more than one value.");
  if (YCYL.size() == 0)
    YCYL.push_back(-1.0);

  checkRead();

  std::vector<AtomNumber> atoms;
  for (unsigned i = 0; i < UMEM.size(); i++)
  {
    atoms.push_back(UMEM[i]);
  }
  for (unsigned i = 0; i < LMEM.size(); i++)
  {
    atoms.push_back(LMEM[i]);
  }
  for (unsigned i = 0; i < TAILS.size(); i++)
  {
    atoms.push_back(TAILS[i]);
  }
  for (unsigned i = 0; i < WATERS.size(); i++)
  {
    atoms.push_back(WATERS[i]);
  }
  for (unsigned i = 0; i < POXYGENS.size(); i++)
  {
    atoms.push_back(POXYGENS[i]);
  }

  addValueWithDerivatives();
  setNotPeriodic();
  requestAtoms(atoms);
}

void fusionPoreNucleationP::calculate()
{
  /*************************
  *                        *
  *         System         *
  *                        *
  **************************/

  // Box dimensions.
  double Lx = getBox()[0][0], Ly = getBox()[1][1], Lz = getBox()[2][2];

  // Z center of the upper membrane (uMem) and lower membrane (lMem) for systems with PBC: https://en.wikipedia.org/wiki/Center_of_mass#Systems_with_periodic_boundary_conditions .
  double ZuMem, ZuMemcos = 0.0, ZuMemsin = 0.0, uMemAngle, ZlMem, ZlMemcos = 0.0, ZlMemsin = 0.0, lMemAngle;

#ifdef _OPENMP
#if _OPENMP >= 201307
  #pragma omp parallel for private(uMemAngle, lMemAngle) reduction(+:ZuMemcos, ZuMemsin, ZlMemcos, ZlMemsin)
#endif
#endif
  for (unsigned i = 0; i < UMEM.size(); i++)
  {
    uMemAngle = 2.0 * M_PI * getPbc().realToScaled(pbcDistance(Vector(0.0, 0.0, 0.0), getPosition(i)))[2];
    lMemAngle = 2.0 * M_PI * getPbc().realToScaled(pbcDistance(Vector(0.0, 0.0, 0.0), getPosition(i + UMEM.size())))[2];
    ZuMemcos += cos(uMemAngle);
    ZuMemsin += sin(uMemAngle);
    ZlMemcos += cos(lMemAngle);
    ZlMemsin += sin(lMemAngle);
  }
  ZuMemcos = ZuMemcos / UMEM.size();
  ZuMemsin = ZuMemsin / UMEM.size();
  ZuMem = Lz * (atan2(-ZuMemsin, -ZuMemcos) + M_PI) / (2.0 * M_PI);
  ZlMemcos = ZlMemcos / UMEM.size();
  ZlMemsin = ZlMemsin / UMEM.size();
  ZlMem = Lz * (atan2(-ZlMemsin, -ZlMemcos) + M_PI) / (2.0 * M_PI);

  // Z center of the boths membranes (upper and lower).
  double ZMems = (ZuMem + ZlMem) / 2.0;

  /**************************
  *                         *
  *   Xcyl_Mem & Ycyl_Mem   *
  *                         *
  ***************************/

  // Quantity of beads of the membranes.
  unsigned membraneBeads = UMEM.size() + LMEM.size();

  // Z distance from the lipid tail to the geometric center of both membranes.
  double ZTailDistance;

  // Z position of the first slice.
  double firstSliceZ_Mem = ZMems + (0.0 + 0.5 - NSMEM[0] / 2.0) * DSMEM[0];

  // Z distance between the first slice and the Z center of the membrane.
  double firstSliceZDist_Mem = pbcDistance(Vector(0.0, 0.0, firstSliceZ_Mem), Vector(0.0, 0.0, ZMems))[2];

  // Position in the cylinder.
  double PositionS_Mem;

  // Slices to analyze per particle.
  unsigned s1_Mem, s2_Mem;

  // Eq. 7 Hub & Awasthi JCTC 2017.
  std::vector<double> faxial_Mem(TAILS.size() * NSMEM[0]);

  // Eq. 10 Hub & Awasthi JCTC 2017.
  std::vector<double> Fs_Mem(NSMEM[0]);

  // Eq. 11 Hub & Awasthi JCTC 2017.
  std::vector<double> ws_Mem(NSMEM[0]);

  // Eq. 10 Hub & Awasthi JCTC 2017.
  double W_Mem = 0.0;

  // Eq. 21 and 22 Hub & Awasthi JCTC 2017.
  std::vector<double> sx_Mem(NSMEM[0]), sy_Mem(NSMEM[0]), cx_Mem(NSMEM[0]), cy_Mem(NSMEM[0]);

  // Eq. 10 Hub & Awasthi JCTC 2017.
  double Xsc_Mem = 0.0, Xcc_Mem = 0.0, Ysc_Mem = 0.0, Ycc_Mem = 0.0;

  // Aux.
  double x, aux;

  // Scaled position of the lipid tail respect the origin of coordinates.
  Vector TailPosition;

#ifdef _OPENMP
#if _OPENMP >= 201307
  #pragma omp declare reduction(vec_double_plus : std::vector<double> : \
  std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
  initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))
#endif
#endif

#ifdef _OPENMP
#if _OPENMP >= 201307
  #pragma omp parallel for private(ZTailDistance, PositionS_Mem, s1_Mem, s2_Mem, TailPosition, x, aux) reduction(vec_double_plus:Fs_Mem, sx_Mem, sy_Mem, cx_Mem, cy_Mem)
#endif
#endif
  for (unsigned i = 0; i < TAILS.size(); i++)
  {
    ZTailDistance = pbcDistance(Vector(0.0, 0.0, ZMems), getPosition(i + membraneBeads))[2];
    PositionS_Mem = (ZTailDistance + firstSliceZDist_Mem) / DSMEM[0];
    // If the following condition is met the particle is in the Z space of the cylinder.
    if ((PositionS_Mem >= (-0.5 - HMEM[0])) && (PositionS_Mem <= (NSMEM[0] + 0.5 - 1.0 + HMEM[0])))
    {
      //Defining the slices to analyze each particle.
      if (PositionS_Mem < 1)
      {
        s1_Mem = 0;
        s2_Mem = 2;
      }
      else if (PositionS_Mem <= (NSMEM[0] - 2.0))
      {
        s1_Mem = floor(PositionS_Mem) - 1;
        s2_Mem = floor(PositionS_Mem) + 1;
      }
      else
      {
        s1_Mem = NSMEM[0] - 3;
        s2_Mem = NSMEM[0] - 1;
      }

      TailPosition = getPbc().realToScaled(pbcDistance(Vector(0.0, 0.0, 0.0), getPosition(i + membraneBeads)));

      for (unsigned s = s1_Mem; s <= s2_Mem; s++)
      {
        x = (ZTailDistance - (s + 0.5 - NSMEM[0] / 2.0) * DSMEM[0]) * 2.0 / DSMEM[0];
        if (!((x <= -1.0 - HMEM[0]) || (x >= 1.0 + HMEM[0])))
        {
          if (((-1.0 + HMEM[0]) <= x) && (x <= (1.0 - HMEM[0])))
          {
            faxial_Mem[i + TAILS.size() * s] = 1.0;
            Fs_Mem[s] += 1.0;
            sx_Mem[s] += sin(2.0 * M_PI * TailPosition[0]);
            sy_Mem[s] += sin(2.0 * M_PI * TailPosition[1]);
            cx_Mem[s] += cos(2.0 * M_PI * TailPosition[0]);
            cy_Mem[s] += cos(2.0 * M_PI * TailPosition[1]);
          }
          else if (((1.0 - HMEM[0]) < x) && (x < (1.0 + HMEM[0])))
          {
            aux = 0.5 - ((3.0 * x - 3.0) / (4.0 * HMEM[0])) + (pow((x - 1.0), 3) / (4.0 * pow(HMEM[0], 3)));
            faxial_Mem[i + TAILS.size() * s] = aux;
            Fs_Mem[s] += aux;
            sx_Mem[s] += aux * sin(2.0 * M_PI * TailPosition[0]);
            sy_Mem[s] += aux * sin(2.0 * M_PI * TailPosition[1]);
            cx_Mem[s] += aux * cos(2.0 * M_PI * TailPosition[0]);
            cy_Mem[s] += aux * cos(2.0 * M_PI * TailPosition[1]);
          }
          else if (((-1.0 - HMEM[0]) < x) && (x < (-1.0 + HMEM[0])))
          {
            aux = 0.5 + ((3.0 * x + 3.0) / (4.0 * HMEM[0])) - (pow((x + 1.0), 3) / (4.0 * pow(HMEM[0], 3)));
            faxial_Mem[i + TAILS.size() * s] = aux;
            Fs_Mem[s] += aux;
            sx_Mem[s] += (aux * sin(2.0 * M_PI * TailPosition[0]));
            sy_Mem[s] += (aux * sin(2.0 * M_PI * TailPosition[1]));
            cx_Mem[s] += (aux * cos(2.0 * M_PI * TailPosition[0]));
            cy_Mem[s] += (aux * cos(2.0 * M_PI * TailPosition[1]));
          }
        }
      }
    }
  }

  for (unsigned s = 0; s < NSMEM[0]; s++)
  {
    if (Fs_Mem[s] != 0.0)
    {
      ws_Mem[s] = tanh(Fs_Mem[s]);
      W_Mem += ws_Mem[s];
      sx_Mem[s] = sx_Mem[s] / Fs_Mem[s];
      sy_Mem[s] = sy_Mem[s] / Fs_Mem[s];
      cx_Mem[s] = cx_Mem[s] / Fs_Mem[s];
      cy_Mem[s] = cy_Mem[s] / Fs_Mem[s];
      Xsc_Mem += sx_Mem[s] * ws_Mem[s];
      Ysc_Mem += sy_Mem[s] * ws_Mem[s];
      Xcc_Mem += cx_Mem[s] * ws_Mem[s];
      Ycc_Mem += cy_Mem[s] * ws_Mem[s];
    }
  }

  Xsc_Mem = Xsc_Mem / W_Mem;
  Ysc_Mem = Ysc_Mem / W_Mem;
  Xcc_Mem = Xcc_Mem / W_Mem;
  Ycc_Mem = Ycc_Mem / W_Mem;

  // Eq. 12 Hub & Awasthi JCTC 2017.
  double Xcyl_Mem, Ycyl_Mem;

  if ((XCYL[0] > 0.0) && (YCYL[0] > 0.0))
  {
    Xcyl_Mem = XCYL[0];
    Ycyl_Mem = YCYL[0];
  }
  else
  {
    Xcyl_Mem = (atan2(-Xsc_Mem, -Xcc_Mem) + M_PI) * Lx / (2 * M_PI);
    Ycyl_Mem = (atan2(-Ysc_Mem, -Ycc_Mem) + M_PI) * Ly / (2 * M_PI);
  }

  /*************************
  *                        *
  *        Xi_n            *
  *                        *
  **************************/

  // Eq. 1 Hub & Awasthi JCTC 2017. This is the CV that describes de Pore Nucleation.
  double Xi_n = 0.0;

  // Quantity of beads that could participate in the calculation of the Xi_Chain
  unsigned chainBeads = WATERS.size() + POXYGENS.size();

  // Quantity of beads that don't participate in the calculation of the Xi_Chain
  unsigned noChainBeads = (UMEM.size() * 2) + TAILS.size();

  // Z Distances from the oxygens to the geometric center of the membranes.
  double ZMemDistance;

  // Scaled positions of the oxygens to respect of the origin of coordinates.
  Vector Position;

  // Distance from the water/phosphate group to the defect cylinder.
  Vector distCylinder;

  // Center of the cylinder. XY components are calculated (or defined), Z is the Z geometric center of the membranes of the system.
  Vector xyzCyl_Mem = pbcDistance(Vector(0.0, 0.0, 0.0), Vector(Xcyl_Mem, Ycyl_Mem, ZMems));

  // Average of the radius of the water and lipid cylinder.
  double RCYLAVERAGE = RCYL[0] * (1 + HCH[0]);

  // Conditions.
  bool condition1, condition2, condition3;

  // Z position of the first slice.
  double firstSliceZ = ZMems + (0.0 + 0.5 - NS[0] / 2.0) * DS[0];

  // Z distance between the first slice and the Z center of the membrane.
  double firstSliceZDist = pbcDistance(Vector(0.0, 0.0, firstSliceZ), Vector(0.0, 0.0, ZMems))[2];

  // Position in the cylinder.
  double PositionS;

  // Mark the particles to analyze.
  std::vector<double> analyzeThisParticle(chainBeads);

  // Slices to analyze per particle.
  std::vector<unsigned> s1(chainBeads), s2(chainBeads);

  // Eq. 7 Hub & Awasthi JCTC 2017.
  std::vector<double> faxial(chainBeads * NS[0]);

  // Eq. 16 Hub & Awasthi JCTC 2017.
  std::vector<double> d_faxial_dz(chainBeads * NS[0]);

  // Eq. 10 Hub & Awasthi JCTC 2017.
  std::vector<double> Fs(NS[0]);

  // Eq. 11 Hub & Awasthi JCTC 2017.
  std::vector<double> ws(NS[0]);

  // Eq. 10 Hub & Awasthi JCTC 2017.
  double W = 0.0;

  // Eq. 21 and 22 Hub & Awasthi JCTC 2017.
  std::vector<double> sx(NS[0]), sy(NS[0]), cx(NS[0]), cy(NS[0]);

  // Eq. 10 Hub & Awasthi JCTC 2017.
  double Xsc = 0.0, Xcc = 0.0, Ysc = 0.0, Ycc = 0.0;

#ifdef _OPENMP
#if _OPENMP >= 201307
  #pragma omp parallel for private(distCylinder, aux, condition1, condition2, condition3, ZMemDistance, PositionS, Position, x) reduction(vec_double_plus:Fs, sx, sy, cx, cy)
#endif
#endif
  for (unsigned i = 0; i < chainBeads; i++)
{
  distCylinder = pbcDistance(xyzCyl_Mem, pbcDistance(Vector(0.0, 0.0, 0.0), getPosition(i + noChainBeads)));
    aux = sqrt(pow(distCylinder[0], 2) + pow(distCylinder[1], 2));
    condition1 = ((aux / RCYLAVERAGE) < 1.0);
    condition2 = ((pbcDistance(Vector(0.0, 0.0, ZuMem), getPosition(i + noChainBeads))[2] > 0) && (aux / RCYLAVERAGE) < 2.0);
    condition3 = ((pbcDistance(getPosition(i + noChainBeads), Vector(0.0, 0.0, ZlMem))[2] > 0) && (aux / RCYLAVERAGE) < 2.0);
    if (condition1 || condition2 || condition3)
    {
      ZMemDistance = pbcDistance(Vector(0.0, 0.0, ZMems), getPosition(i + noChainBeads))[2];
      PositionS = (ZMemDistance + firstSliceZDist) / DS[0];
      // If the following condition is met the particle is in the Z space of the cylinder.
      if ((PositionS >= (-0.5 - HCH[0])) && (PositionS <= (NS[0] + 0.5 - 1.0 + HCH[0])))
      {
        analyzeThisParticle[i] = 1.0;

        //Defining the slices to analyze each particle.
        if (PositionS < 1)
        {
          s1[i] = 0;
          s2[i] = 2;
        }
        else if (PositionS <= (NS[0] - 2.0))
        {
          s1[i] = floor(PositionS) - 1;
          s2[i] = floor(PositionS) + 1;
        }
        else
        {
          s1[i] = NS[0] - 3;
          s2[i] = NS[0] - 1;
        }

        Position = getPbc().realToScaled(pbcDistance(Vector(0.0, 0.0, 0.0), getPosition(i + noChainBeads)));

        for (unsigned s = s1[i]; s <= s2[i]; s++)
        {
          x = (ZMemDistance - (s + 0.5 - NS[0] / 2.0) * DS[0]) * 2.0 / DS[0];
          if (!((x <= -1.0 - HCH[0]) || (x >= 1.0 + HCH[0])))
          {
            if (((-1.0 + HCH[0]) <= x) && (x <= (1.0 - HCH[0])))
            {
              faxial[i + chainBeads * s] = 1.0;
              Fs[s] += 1.0;
              sx[s] += sin(2.0 * M_PI * Position[0]);
              sy[s] += sin(2.0 * M_PI * Position[1]);
              cx[s] += cos(2.0 * M_PI * Position[0]);
              cy[s] += cos(2.0 * M_PI * Position[1]);
            }
            else if (((1.0 - HCH[0]) < x) && (x < (1.0 + HCH[0])))
            {
              aux = 0.5 - ((3.0 * x - 3.0) / (4.0 * HCH[0])) + (pow((x - 1.0), 3) / (4.0 * pow(HCH[0], 3)));
              faxial[i + chainBeads * s] = aux;
              d_faxial_dz[i + chainBeads * s] = ((-3.0 / (4.0 * HCH[0])) + ((3.0 * pow((x - 1), 2)) / (4.0 * pow(HCH[0], 3)))) * 2.0 / DS[0];
              Fs[s] += aux;
              sx[s] += aux * sin(2.0 * M_PI * Position[0]);
              sy[s] += aux * sin(2.0 * M_PI * Position[1]);
              cx[s] += aux * cos(2.0 * M_PI * Position[0]);
              cy[s] += aux * cos(2.0 * M_PI * Position[1]);
            }
            else if (((-1.0 - HCH[0]) < x) && (x < (-1.0 + HCH[0])))
            {
              aux = 0.5 + ((3.0 * x + 3.0) / (4.0 * HCH[0])) - (pow((x + 1.0), 3) / (4.0 * pow(HCH[0], 3)));
              faxial[i + chainBeads * s] = aux;
              d_faxial_dz[i + chainBeads * s] = ((3.0 / (4.0 * HCH[0])) - ((3.0 * pow((x + 1), 2)) / (4.0 * pow(HCH[0], 3)))) * 2.0 / DS[0];
              Fs[s] += aux;
              sx[s] += (aux * sin(2.0 * M_PI * Position[0]));
              sy[s] += (aux * sin(2.0 * M_PI * Position[1]));
              cx[s] += (aux * cos(2.0 * M_PI * Position[0]));
              cy[s] += (aux * cos(2.0 * M_PI * Position[1]));
            }
          }
        }
      }
    }
  }

  for (unsigned s = 0; s < NS[0]; s++)
  {
    if (Fs[s] != 0.0)
    {
      ws[s] = tanh(Fs[s]);
      W += ws[s];
      sx[s] = sx[s] / Fs[s];
      sy[s] = sy[s] / Fs[s];
      cx[s] = cx[s] / Fs[s];
      cy[s] = cy[s] / Fs[s];
      Xsc += sx[s] * ws[s];
      Ysc += sy[s] * ws[s];
      Xcc += cx[s] * ws[s];
      Ycc += cy[s] * ws[s];
    }
  }

  Xsc = Xsc / W;
  Ysc = Ysc / W;
  Xcc = Xcc / W;
  Ycc = Ycc / W;

  // Eq. 12 Hub & Awasthi JCTC 2017.
  double Xcyl, Ycyl;

  Xcyl = Xcyl_Mem;
  Ycyl = Ycyl_Mem;

  // Eq. 25, 26 and 27 Hub & Awasthi JCTC 2017.
  double d_sx_dx, d_sx_dz, d_sy_dy, d_sy_dz, d_cx_dx, d_cx_dz, d_cy_dy, d_cy_dz;

  // Eq. 29 Hub & Awasthi JCTC 2017.
  double d_ws_dz;

  // Eq. 31, 32 and 33 Hub & Awasthi JCTC 2017
  double d_Xsc_dx, d_Xsc_dz, d_Xcc_dx, d_Xcc_dz, d_Ysc_dy, d_Ysc_dz, d_Ycc_dy, d_Ycc_dz;

  // Center of the cylinder. X and Y are calculated (or defined), Z is the Z component of the geometric center of the membranes.
  Vector xyzCyl = pbcDistance(Vector(0.0, 0.0, 0.0), Vector(Xcyl, Ycyl, ZMems));

  // Distances from the oxygens to center of the cylinder.
  std::vector<Vector> CylDistances(chainBeads);

  // Modulo of the XY distances from the oxygens to the center of the cylinder.
  double ri;

  // Eq. 8 Hub & Awasthi JCTC 2017.
  double fradial;

  // Eq. 15 Hub & Awasthi JCTC 2017.
  std::vector<double> d_fradial_dx(chainBeads), d_fradial_dy(chainBeads);

  // Eq. 35, 36, 37 and 38 Hub & Awasthi JCTC 2017.
  std::vector<double> d_Xcyl_dx(chainBeads), d_Xcyl_dz(chainBeads), d_Ycyl_dy(chainBeads), d_Ycyl_dz(chainBeads);

  // To avoid rare instabilities auxX and auxY are truncated at a configurable value (default 500).
  double auxX = (1 / (pow(Xsc, 2) + pow(Xcc, 2))), auxY = (1 / (pow(Ysc, 2) + pow(Ycc, 2)));

  if (auxX > ONEOVERS2C2CUTOFF[0])
  {
    auxX = Lx * ONEOVERS2C2CUTOFF[0] / (2 * M_PI);
  }
  else
  {
    auxX = Lx * auxX / (2 * M_PI);
  }

  if (auxY > ONEOVERS2C2CUTOFF[0])
  {
    auxY = Ly * ONEOVERS2C2CUTOFF[0] / (2 * M_PI);
  }
  else
  {
    auxY = Ly * auxY / (2 * M_PI);
  }

  //Number of oxygens within the slice s of the membrane-spanning cylinder.
  std::vector<double> Nsp(NS[0]), psi(NS[0]), d_psi(NS[0]);

  // Eq. 3 Hub & Awasthi JCTC 2017.
  double b = (ZETA[0] / (1.0 - ZETA[0])), c = ((1.0 - ZETA[0]) * exp(b));

  // Eq. 19 Hub & Awasthi JCTC 2017.
  std::vector<double> fradial_d_faxial_dz(chainBeads * NS[0]);

  // Eq. 20 Hub & Awasthi JCTC 2017.
  std::vector<double> Axs(NS[0]), Ays(NS[0]);

#ifdef _OPENMP
#if _OPENMP >= 201307
  #pragma omp parallel for private(d_Xsc_dx,d_Xcc_dx,d_Ysc_dy,d_Ycc_dy,d_Xsc_dz,d_Xcc_dz,d_Ysc_dz,d_Ycc_dz,d_sx_dx,d_sy_dy,d_cx_dx,d_cy_dy,d_sx_dz,d_sy_dz,d_cx_dz,d_cy_dz,d_ws_dz,ri,x,fradial) reduction(vec_double_plus: Nsp, Axs, Ays)
#endif
#endif
  for (unsigned i = 0; i < chainBeads; i++)
{
  CylDistances[i] = pbcDistance(xyzCyl, pbcDistance(Vector(0.0, 0.0, 0.0), getPosition(i + noChainBeads)));
    if (analyzeThisParticle[i])
    {
      Position = getPbc().realToScaled(pbcDistance(Vector(0.0, 0.0, 0.0), getPosition(i + noChainBeads)));
      d_Xsc_dx = 0.0;
      d_Xcc_dx = 0.0;
      d_Ysc_dy = 0.0;
      d_Ycc_dy = 0.0;
      d_Xsc_dz = 0.0;
      d_Xcc_dz = 0.0;
      d_Ysc_dz = 0.0;
      d_Ycc_dz = 0.0;
      for (unsigned s = s1[i]; s <= s2[i]; s++)
      {
        if (Fs[s] != 0.0)
        {
          d_sx_dx = faxial[i + chainBeads * s] * 2.0 * M_PI * cos(2.0 * M_PI * Position[0]) / (Lx * Fs[s]);
          d_sy_dy = faxial[i + chainBeads * s] * 2.0 * M_PI * cos(2.0 * M_PI * Position[1]) / (Ly * Fs[s]);
          d_cx_dx = -faxial[i + chainBeads * s] * 2.0 * M_PI * sin(2.0 * M_PI * Position[0]) / (Lx * Fs[s]);
          d_cy_dy = -faxial[i + chainBeads * s] * 2.0 * M_PI * sin(2.0 * M_PI * Position[1]) / (Ly * Fs[s]);
          d_Xsc_dx += ws[s] * d_sx_dx / W;
          d_Xcc_dx += ws[s] * d_cx_dx / W;
          d_Ysc_dy += ws[s] * d_sy_dy / W;
          d_Ycc_dy += ws[s] * d_cy_dy / W;

          d_sx_dz = d_faxial_dz[i + chainBeads * s] * (sin(2.0 * M_PI * Position[0]) - sx[s]) / Fs[s];
          d_sy_dz = d_faxial_dz[i + chainBeads * s] * (sin(2.0 * M_PI * Position[1]) - sy[s]) / Fs[s];
          d_cx_dz = d_faxial_dz[i + chainBeads * s] * (cos(2.0 * M_PI * Position[0]) - cx[s]) / Fs[s];
          d_cy_dz = d_faxial_dz[i + chainBeads * s] * (cos(2.0 * M_PI * Position[1]) - cy[s]) / Fs[s];
          d_ws_dz = (1 - pow(ws[s], 2)) * d_faxial_dz[i + chainBeads * s];
          d_Xsc_dz += (ws[s] * d_sx_dz + d_ws_dz * (sx[s] - Xsc)) / W;
          d_Xcc_dz += (ws[s] * d_cx_dz + d_ws_dz * (cx[s] - Xcc)) / W;
          d_Ysc_dz += (ws[s] * d_sy_dz + d_ws_dz * (sy[s] - Ysc)) / W;
          d_Ycc_dz += (ws[s] * d_cy_dz + d_ws_dz * (cy[s] - Ycc)) / W;
        }
      }
      d_Xcyl_dx[i] = auxX * (-Xsc * d_Xcc_dx + Xcc * d_Xsc_dx);
      d_Xcyl_dz[i] = auxX * (-Xsc * d_Xcc_dz + Xcc * d_Xsc_dz);
      d_Ycyl_dy[i] = auxY * (-Ysc * d_Ycc_dy + Ycc * d_Ysc_dy);
      d_Ycyl_dz[i] = auxY * (-Ysc * d_Ycc_dz + Ycc * d_Ysc_dz);

      ri = sqrt(pow(CylDistances[i][0], 2) + pow(CylDistances[i][1], 2));
      x = ri / RCYL[0];
      if (!((x <= -1.0 - HCH[0]) || (x >= 1.0 + HCH[0])))
      {
        if (((-1.0 + HCH[0]) <= x) && (x <= (1.0 - HCH[0])))
        {
          fradial = 1.0;
        }
        else if (((1.0 - HCH[0]) < x) && (x < (1.0 + HCH[0])))
        {
          fradial = 0.5 - ((3.0 * x - 3.0) / (4.0 * HCH[0])) + (pow((x - 1.0), 3) / (4.0 * pow(HCH[0], 3)));
          d_fradial_dx[i] = ((-3.0 / (4.0 * HCH[0])) + ((3.0 * pow((x - 1), 2)) / (4.0 * pow(HCH[0], 3)))) * CylDistances[i][0] / (RCYL[0] * ri);
          d_fradial_dy[i] = ((-3.0 / (4.0 * HCH[0])) + ((3.0 * pow((x - 1), 2)) / (4.0 * pow(HCH[0], 3)))) * CylDistances[i][1] / (RCYL[0] * ri);
        }
        else if (((-1.0 - HCH[0]) < x) && (x < (-1.0 + HCH[0])))
        {
          fradial = 0.5 + ((3.0 * x + 3.0) / (4.0 * HCH[0])) - (pow((x + 1.0), 3) / (4.0 * pow(HCH[0], 3)));
          d_fradial_dx[i] = ((3.0 / (4.0 * HCH[0])) - ((3.0 * pow((x + 1), 2)) / (4.0 * pow(HCH[0], 3)))) * CylDistances[i][0] / (RCYL[0] * ri);
          d_fradial_dy[i] = ((3.0 / (4.0 * HCH[0])) - ((3.0 * pow((x + 1), 2)) / (4.0 * pow(HCH[0], 3)))) * CylDistances[i][1] / (RCYL[0] * ri);
        }

        for (unsigned s = s1[i]; s <= s2[i]; s++)
        {
          Nsp[s] += fradial * faxial[i + chainBeads * s];
          Axs[s] += faxial[i + chainBeads * s] * d_fradial_dx[i];
          Ays[s] += faxial[i + chainBeads * s] * d_fradial_dy[i];
          fradial_d_faxial_dz[i + chainBeads * s] = fradial * d_faxial_dz[i + chainBeads * s];
        }
      }
    }
  }

  for (unsigned s = 0; s < NS[0]; s++)
  {
    if (Nsp[s] <= 1.0)
    {
      psi[s] = ZETA[0] * Nsp[s];
      d_psi[s] = ZETA[0];
      Xi_n += psi[s];
    }
    else
    {
      psi[s] = 1.0 - c * exp(-b * Nsp[s]);
      d_psi[s] = b * c * exp(-b * Nsp[s]);
      Xi_n += psi[s];
    }
  }

  Xi_n = Xi_n / NS[0];

  // Eq. 18 Hub & Awasthi JCTC 2017.
  std::vector<double> faxial_d_fradial_dx(chainBeads * NS[0]), faxial_d_fradial_dy(chainBeads * NS[0]), faxial_d_fradial_dz(chainBeads * NS[0]);

  // Eq. 13 Hub & Awasthi JCTC 2017 modified to considere the Heaviside_Chain step function (this only affect during the transition).
  std::vector<Vector> derivatives_Chain(chainBeads);

#ifdef _OPENMP
#if _OPENMP >= 201307
  #pragma omp parallel for private(aux)
#endif
#endif
  for (unsigned i = 0; i < chainBeads; i++)
{
  if (analyzeThisParticle[i])
    {
      for (unsigned s = s1[i]; s <= s2[i]; s++)
      {
        if (faxial[i + chainBeads * s])
        {
          faxial_d_fradial_dx[i + chainBeads * s] = faxial[i + chainBeads * s] * d_fradial_dx[i] - d_Xcyl_dx[i] * Axs[s];
          faxial_d_fradial_dy[i + chainBeads * s] = faxial[i + chainBeads * s] * d_fradial_dy[i] - d_Ycyl_dy[i] * Ays[s];
          faxial_d_fradial_dz[i + chainBeads * s] = -d_Xcyl_dz[i] * Axs[s] - d_Ycyl_dz[i] * Ays[s];
        }
      }

      for (unsigned s = s1[i]; s <= s2[i]; s++)
      {
        aux = d_psi[s] / NS[0];
        derivatives_Chain[i][0] += aux * faxial_d_fradial_dx[i + chainBeads * s];
        derivatives_Chain[i][1] += aux * faxial_d_fradial_dy[i + chainBeads * s];
        derivatives_Chain[i][2] += aux * (faxial_d_fradial_dz[i + chainBeads * s] + fradial_d_faxial_dz[i + chainBeads * s]);
      }
    }
  }

  Tensor virial;
  for (unsigned i = 0; i < chainBeads; i++)
{
  setAtomsDerivatives((i + noChainBeads), derivatives_Chain[i]);
    virial -= Tensor(CylDistances[i], derivatives_Chain[i]);
  }

  setValue(Xi_n);
  setBoxDerivatives(virial);
}
}
}
