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
//+PLUMEDOC MEMBRANEFUSIONMOD_COLVAR FUSIONPOREEXPANSIONP
/*
A CV for inducing the expansion of a fusion pore from a nucleated fusion pore.

Calculate the collective variable designed by Hub \cite Hub2021  and implemented into PLUMED by Masone and collaborators.
This CV is capable of inducing the expansion of the fusion pore from a nucleated fusion pore.

\f[
\xi_e = \frac{R(r) - R_0}{R_0}
\f]

Where \f$\xi_e\f$ is the CV, \f$R_0\f$ is a normalization constant that makes zero the initial value of \f$\xi_e\f$, and
\f$R(r)\f$ is the approximate radius of the fusion pore, which is defined by the number of waters and phosphateoxygens
beads within a horizontal layer in the center of both membranes.

\par Examples

This example induces the expansion of a nucleated fusion pore (\f$\xi_e = 0.75\f$) from a just nucleated fusion pore (\f$\xi_e = 0.00\f$).

\plumedfile
lMem: GROUP ATOMS=1-10752,21505-22728,23953-24420 #All the lower membrane beads.
uMem: GROUP ATOMS=10753-21504,22729-23952,24421-24888 #All the upper membrane beads.
tails: GROUP ATOMS=8-23948:12,12-23952:12,23966-24884:18,23970-24888:18 #All the lipid tails beads (from the lower and upper membrane).
waters: GROUP ATOMS=24889-56589  #All the water beads.
po4: GROUP ATOMS=2-23942:12,23957-24875:18 #All the lipid phosphateoxygens beads.

fusionPoreExpansion: FUSIONPOREEXPANSIONP UMEMBRANE=uMem LMEMBRANE=lMem TAILS=tails WATERS=waters PHOSPHATEOXYGENS=po4 NSMEM=85 D=7.0 R0=0.57

MOVINGRESTRAINT ...
    ARG=fusionPoreExpansion
    STEP0=0 AT0=0.0 KAPPA0=10000.0
    STEP1=500000 AT1=0.75 KAPPA1=10000.0
...

PRINT ARG=fusionPoreExpansion FILE=COLVAR STRIDE=1

\endplumedfile

*/
//+ENDPLUMEDOC
class fusionPoreExpansionP : public Colvar
{
  std::vector<AtomNumber> UMEM, LMEM, TAILS, WATERS, POXYGENS;
  std::vector<double> NSMEM, DSMEM, HMEM, VO, D, H, RMAX, R0, XCYL, YCYL;

public:
  explicit fusionPoreExpansionP(const ActionOptions &);
  void calculate() override;
  static void registerKeywords(Keywords &keys);
};

PLUMED_REGISTER_ACTION(fusionPoreExpansionP, "FUSIONPOREEXPANSIONP")

void fusionPoreExpansionP::registerKeywords(Keywords &keys)
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
  keys.add("optional", "VO", "( default=0.076879 ) beads' molecular volume.");
  keys.add("compulsory", "D", "horizontal layer thickness, it depends on the Z separation of the membranes.");
  keys.add("optional", "H", "( default=0.1 ) parameter of the step function θ(x,h) for the fusion pore expansion.");
  keys.add("optional", "RMAX", "( default=2.5 ) to avoid effects of membrane undulations in large membranes (more than 256 lipids).");
  keys.add("compulsory", "R0", "normalization constant that makes 0 the initial value of the CV.");
  keys.add("optional", "XCYL", "X coordinate of the fixed cylinder, if not present this will be calculated.");
  keys.add("optional", "YCYL", "X coordinate of the fixed cylinder, if not present this will be calculated.");
}

fusionPoreExpansionP::fusionPoreExpansionP(const ActionOptions &ao) : PLUMED_COLVAR_INIT(ao)
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

  parseVector("VO", VO);
  if (VO.size() > 1)
    error("VO cannot take more than one value.");
  if (VO.size() == 0)
    VO.push_back(0.076879);

  parseVector("D", D);
  if (D.size() > 1)
    error("D cannot take more than one value.");

  parseVector("H", H);
  if (H.size() > 1)
    error("H cannot take more than one value.");
  if (H.size() == 0)
    H.push_back(0.1);

  parseVector("RMAX", RMAX);
  if (RMAX.size() > 1)
    error("RMAX cannot take more than one value.");
  if (RMAX.size() == 0)
    RMAX.push_back(2.5);

  parseVector("R0", R0);
  if (R0.size() > 1)
    error("R0 cannot take more than one value.");

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
void fusionPoreExpansionP::calculate()
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
  #pragma omp parallel for private(ZTailDistance, PositionS_Mem, TailPosition, x, aux, s1_Mem, s2_Mem) reduction(vec_double_plus:Fs_Mem, sx_Mem, sy_Mem, cx_Mem, cy_Mem)
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
  *         Xi_Exp         *
  *                        *
  **************************/

  // Quantity of beads that could participate in the calculation of the Xi_Chain
  unsigned chainBeads = WATERS.size() + POXYGENS.size();

  // Quantity of beads that don't participate in the calculation of the Xi_Chain
  unsigned noChainBeads = (UMEM.size() * 2) + TAILS.size();

  // Center of the cylinder. X and Y are calculated (or defined), Z is the Z component of the geometric center of the membranes.
  Vector xyzCyl = pbcDistance(Vector(0.0, 0.0, 0.0), Vector(Xcyl_Mem, Ycyl_Mem, ZMems));

  // Estimation of RO with the Hub 2021 JCTC method. Only needed for the expansion.
  double RO = R0[0];

  // Number of polar atoms inside the horizontal layer. Eq. 3 Hub 2021 JCTC.
  double np = 0.0, fz, fr, fz_prime, fr_prime;

  // Derivative of np. Eq. 8 Hub 2021 JCTC.
  std::vector<double> d_np_dx(chainBeads), d_np_dy(chainBeads), d_np_dz(chainBeads);

  // Pore radius of the defect. Eq. 2 Hub 2021 JCTC.
  double poreR = 1.0;

  // Z center of the Membrane in the RMAX radius.
  double ZMemRMAX, ZMemRMAXcos = 0.0, ZMemRMAXsin = 0.0, countAux = 0.0, auxcos, auxsin;

  ZMemRMAX = ZMems;

  // The curvature of large membranes (1024 lipids) makes the Z-center of the membranes not to be representative
  // in some sectors, particularly in the region of ​​the defect.
  //
  // To solve this, the center Z of the membranes in the defect sector is calculated and used to calculate
  // the number of polar atoms within the horizontal layer AND in the radious of the defect.
  //
  // ________       | |       ________
  // ________ \_____| |______/ _______<-- Top membrane.
  //         \______|P|_______/
  //                |O|
  //                | |               <-- Z-center of the membranes in the region of the defect.
  //          ______|R|_______        <-- Z-center of the membranes
  //         / _____|E|______ \ 
  //        / /     | |      \ \ 
  // ______/ /      | |       \ \______
  // _______/                  \_______<-- Bottom membrane.

  // Center of mass for systems with PBC: https://en.wikipedia.org/wiki/Center_of_mass#Systems_with_periodic_boundary_conditions
  Vector MemCylDistances, distCylinder;
  double angle, ri;

#ifdef _OPENMP
#if _OPENMP >= 201307
  #pragma omp parallel for private(MemCylDistances, x, angle, auxcos, auxsin) reduction(+:ZMemRMAXcos, ZMemRMAXsin, countAux)
#endif
#endif
  for (unsigned i = 0; i < membraneBeads; i++)
{
  MemCylDistances = pbcDistance(xyzCyl, pbcDistance(Vector(0.0, 0.0, 0.0), getPosition(i)));
    x = sqrt(pow(MemCylDistances[0], 2) + pow(MemCylDistances[1], 2)) / RMAX[0];
    if (!((x <= -1.0 - H[0]) || (x >= 1.0 + H[0])))
    {
      angle = 2.0 * M_PI * getPbc().realToScaled(pbcDistance(Vector(0.0, 0.0, 0.0), getPosition(i)))[2];
      auxcos = cos(angle);
      auxsin = sin(angle);
      if (((-1.0 + H[0]) <= x) && (x <= (1.0 - H[0])))
      {
        ZMemRMAXcos += 1.0 * auxcos;
        ZMemRMAXsin += 1.0 * auxsin;
        countAux += 1.0;
      }
      else if (((1.0 - H[0]) < x) && (x < (1.0 + H[0])))
      {
        ZMemRMAXcos += (0.5 - 0.75 * (x - 1.0) / H[0] + 0.25 * pow((x - 1.0), 3) / pow(H[0], 3)) * auxcos;
        ZMemRMAXsin += (0.5 - 0.75 * (x - 1.0) / H[0] + 0.25 * pow((x - 1.0), 3) / pow(H[0], 3)) * auxsin;
        countAux += (0.5 - 0.75 * (x - 1.0) / H[0] + 0.25 * pow((x - 1.0), 3) / pow(H[0], 3));
      }
      else if (((-1.0 - H[0]) < x) && (x < (-1.0 + H[0])))
      {
        ZMemRMAXcos += (0.5 + 0.75 * (x + 1.0) / H[0] - 0.25 * pow((x + 1.0), 3) / pow(H[0], 3)) * auxcos;
        ZMemRMAXsin += (0.5 + 0.75 * (x + 1.0) / H[0] - 0.25 * pow((x + 1.0), 3) / pow(H[0], 3)) * auxsin;
        countAux += (0.5 + 0.75 * (x + 1.0) / H[0] - 0.25 * pow((x + 1.0), 3) / pow(H[0], 3));
      }
    }
  }

  ZMemRMAXcos = ZMemRMAXcos / countAux;
  ZMemRMAXsin = ZMemRMAXsin / countAux;
  ZMemRMAX = Lz * (atan2(-ZMemRMAXsin, -ZMemRMAXcos) + M_PI) / (2.0 * M_PI);

  xyzCyl = pbcDistance(Vector(0.0, 0.0, 0.0), Vector(Xcyl_Mem, Ycyl_Mem, ZMemRMAX));

#ifdef _OPENMP
#if _OPENMP >= 201307
  #pragma omp parallel for private(distCylinder, fz, fz_prime, fr, fr_prime, ri, x) reduction(+:np)
#endif
#endif
  for (unsigned i = 0; i < chainBeads; i++)
{
  distCylinder = pbcDistance(xyzCyl, pbcDistance(Vector(0.0, 0.0, 0.0), getPosition(i + noChainBeads)));
    fz = 0.0;
    fz_prime = 0.0;
    fr = 0.0;
    fr_prime = 0.0;

    ri = sqrt(pow(distCylinder[0], 2) + pow(distCylinder[1], 2));
    x = ri / RMAX[0];
    if (!((x <= -1.0 - H[0]) || (x >= 1.0 + H[0])))
    {
      if (((-1.0 + H[0]) <= x) && (x <= (1.0 - H[0])))
      {
        fr = 1.0;
      }
      else if (((1.0 - H[0]) < x) && (x < (1.0 + H[0])))
      {
        fr = 0.5 - 0.75 * (x - 1.0) / H[0] + 0.25 * pow((x - 1.0), 3) / pow(H[0], 3);
        fr_prime = (-0.75 / H[0] + 0.75 * pow((x - 1.0), 2) / pow(H[0], 3)) / (RMAX[0] * ri);
      }
      else if (((-1.0 - H[0]) < x) && (x < (-1.0 + H[0])))
      {
        fr = 0.5 + 0.75 * (x + 1.0) / H[0] - 0.25 * pow((x + 1.0), 3) / pow(H[0], 3);
        fr_prime = (0.75 / H[0] - 0.75 * pow((x + 1), 2) / pow(H[0], 3)) / (RMAX[0] * ri);
      }

      x = distCylinder[2] * 2.0 / D[0];
      if (!((x <= -1.0 - H[0]) || (x >= 1.0 + H[0])))
      {
        if (((-1.0 + H[0]) <= x) && (x <= (1.0 - H[0])))
        {
          fz = 1.0;
        }
        else if (((1.0 - H[0]) < x) && (x < (1.0 + H[0])))
        {
          fz = 0.5 - 0.75 * (x - 1.0) / H[0] + 0.25 * pow((x - 1.0), 3) / pow(H[0], 3);
          fz_prime = (-0.75 / H[0] + 0.75 * pow((x - 1.0), 2) / pow(H[0], 3)) * 2.0 / D[0];
        }
        else if (((-1.0 - H[0]) < x) && (x < (-1.0 + H[0])))
        {
          fz = 0.5 + 0.75 * (x + 1.0) / H[0] - 0.25 * pow((x + 1.0), 3) / pow(H[0], 3);
          fz_prime = (0.75 / H[0] - 0.75 * pow((x + 1), 2) / pow(H[0], 3)) * 2.0 / D[0];
        }

        np += fz * fr;
        d_np_dx[i] = fz * fr_prime * distCylinder[0];
        d_np_dy[i] = fz * fr_prime * distCylinder[1];
        d_np_dz[i] = fz_prime * fr;
      }
    }
  }
  poreR = sqrt(np * VO[0] / (M_PI * D[0]));

  // This is the CV that describes the Pore Expansion.
  double Xi_Exp = (poreR - RO) / RO;

  // Derivatives vector.
  std::vector<Vector> derivatives(chainBeads);

  // Aux for the derivatives calculations. Eq. 7 Hub 2021 JCTC.
  double fact2 = 0.0;

  if (poreR != 0.0)
{
  fact2 = VO[0] / (2.0 * M_PI * RO * D[0] * poreR);
  }

  // Distances from the oxygens to center of the cylinder.
  std::vector<Vector> CylDistances(chainBeads);

#ifdef _OPENMP
#if _OPENMP >= 201307
  #pragma omp parallel for
#endif
#endif
  for (unsigned i = 0; i < chainBeads; i++)
{
  derivatives[i][0] = fact2 * d_np_dx[i];
    derivatives[i][1] = fact2 * d_np_dy[i];
    derivatives[i][2] = fact2 * d_np_dz[i];
    CylDistances[i] = pbcDistance(xyzCyl, pbcDistance(Vector(0.0, 0.0, 0.0), getPosition(i + noChainBeads)));
  }

  Tensor virial;
  for (unsigned i = 0; i < chainBeads; i++)
{
  setAtomsDerivatives((i + noChainBeads), derivatives[i]);
    virial -= Tensor(CylDistances[i], derivatives[i]);
  }

  setValue(Xi_Exp);
  setBoxDerivatives(virial);
}
}
}
