/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2016 The plumed team
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
#ifndef __PLUMED_tools_Kearsley_h
#define __PLUMED_tools_Kearsley_h

#include "Vector.h"
#include "Tensor.h"
#include <vector>

namespace PLMD{

class Log;

/*! A class that implements Kearsley's calculation
 which is optimal alignment via quaternion and
 analytical derivatives via perturbation theory

In Kearsley algorithm (see, S. K. Kearsley, Acta Crystallogr., Sect. A: Found. Crystallogr. 45, 208 1989 ),
here adapted to included the COM reset,
one first calculates the COM reset position respect to frame  \f$ \{x_0,y_0,z_0\} \f$ (running frame),
being \f$ w_{tot}^{al}=\sum_i w_{al}^i \f$

\f{eqnarray}
\tilde x_0^i=x_0^i-\sum_i \frac{w_{al}^i}{w_{tot}^{al}} x_0^i\\
\tilde y_0^i=y_0^i-\sum_i \frac{w_{al}^i}{w_{tot}^{al}} y_0^i\\
\tilde z_0^i=z_0^i-\sum_i \frac{w_{al}^i}{w_{tot}^{al}} z_0^i
\f}

and the same is done with the reference one \f$ \{x_1,y_1,z_1\} \f$
\f{eqnarray}
\tilde x_1^i=x_1^i-\sum_i \frac{w_{al}^i}{w_{tot}^{al}} x_1^i\\
\tilde y_1^i=y_1^i-\sum_i \frac{w_{al}^i}{w_{tot}^{al}} y_1^i\\
\tilde z_1^i=z_1^i-\sum_i \frac{w_{al}^i}{w_{tot}^{al}} z_1^i
\f}

Then one can build the \f$ \{x_p,y_p,z_p\} \f$ and \f$ \{x_m,y_m,z_m\} \f$
vectors of weighted summation and differences:
\f{eqnarray}
x_m^i=w_{al}^i(\tilde x_0^i-\tilde x_1^i)\\
y_m^i=w_{al}^i(\tilde y_0^i-\tilde y_1^i)\\
z_m^i=w_{al}^i(\tilde z_0^i-\tilde z_1^i)
\f}

\f{eqnarray}
x_p^i=w_{al}^i(x_0^i+x_1^i)\\
y_p^i=w_{al}^i(y_0^i+y_1^i)\\
z_p^i=w_{al}^i(z_0^i+z_1^i)
\f}


Then one build the COM-resetted matrix

\f{displaymath}
\mathbf{M}=\left[
\begin{array}{cccc}
\sum ( {x}_{m}^{2}+{y}_{m}^{2}+{z}_{m}^{2})      &
 \sum (y_{p}z_{m} -y_{m}z_{p})   &
 \sum ( x_{m}z_{p} -x_{p}z_{m}) &
 \sum (x_{p}y_{m}-x_{m}y_{p} ) \\
  \sum ( y_{p}z_{m} -y_{m}z_{p}) &
\sum (  {x}_{m}^{2}+{y}_{p}^{2}+{z}_{p}^{2}) &
\sum ( x_{m}y_{m} -x_{p}y_{p} ) &
 \sum (x_{m}z_{m}-x_{p}z_{p} ) \\
  \sum (x_m z_p - x_p z_m ) &
\sum ( x_m y_m -x_p y_p) &
\sum (  {x}_{p}^{2}+{y}_{m}^{2}+{z}_{p}^{2}) &
\sum ( y_m z_m -y_p z_p) \\
  \sum (x_p y_m -x_m y_p ) &
\sum (x_m z_m - x_p z_p ) &
\sum (y_m z_m- y_p z_p ) &
\sum (  {x}_{p}^{2}+{y}_{p}^{2}+{z}_{m}^{2}) \\
\end{array}
\right]
\f}

by diagonalizing one obtains the mean square deviation by using the lowest eigenvalue \f$ \lambda_0 \f$
\f{equation}
MSD=  \frac{\lambda_{0}}{w_{tot}^{al}}
\f}

The rotation matrix is obtained from the eigenvector corresponding to \f$ \lambda_0 \f$ eigenvalue
having components \f$ q_1, q_2, q_3, q_4 \f$

\f{displaymath}
\mathbf{R}=\left[
\begin{array}{ccc}
q_1 ^2 + q_2 ^2 - q_3 ^2 - q_4^2 &
2(q_2 q_3 + q_1 q_4) &
2(q_1 q_4 -q_1 q_3 )\\
2(q_2 q_3 - q_1 q_4) &
q_1 ^2 +q_3 ^2 -q_2 ^2 -q_4^2  &
2(q_3 q_4 - q_1 q_2)\\
2( q_2 q_4 + q_1 q_3) &
2( q_3 q_4 - q_1 q_2) &
q_1^2 +q_4 ^2 - q_2^2 - q_3 ^2 \\
\end{array}
\right]
\f}

by using the perturbation theory one can retrieve the various derivatives:

In derivative calculation we exploited the classical Perturbation Theory up to the first order.
In extensive manner, we introduce a perturbation over \f$\lambda_{0}\f$ correlated with
a pertubation of the states \f$\vert q_{0}\rangle \f$ (in bra-ket notation):
\f{displaymath}
[\mathbf{M}+d\mathbf{M}][\vert q_{0}\rangle + \vert dq_{0}\rangle ]=
[\lambda_{0}+d\lambda_{0}][\vert q_{0}\rangle +\vert dq_{0}\rangle ]
\f}
Grouping the zero order we recollect the unperturbed equation(see before).
To the first order:
\f{displaymath}
d\mathbf{M}q_{0}+\mathbf{M}\vert dq_{0}\rangle =d\lambda_{0}\vert q_{0}\rangle +\lambda_{0} \vert dq_{0}\rangle
\f}
Now we express \f$dq_{0}\f$ as linear combination of the other ortogonal eigenvectors:
\f{displaymath}
\vert dq_{0}\rangle =\sum_{j\neq0}c_{j}\vert q_{j}\rangle
\f}
thus we have
\f{displaymath}
d\mathbf{M}\vert q_{0}\rangle +\sum_{j\neq0}c_{j}\mathbf{M}\vert q_{j}\rangle=
d\lambda_{0}\vert q_{0}\rangle+\lambda_{0}\sum_{j\neq0}c_{j}\vert q_{j}\rangle
\f}
projecting onto the \f$q_{0}\f$ state and deleting the projection onto \f$\vert dq_{0}\rangle\f$ beacuse
of ortogonality:
\f{displaymath}
\langle q_{0}\vert d\mathbf{M}\vert q_{0}\rangle +\sum_{j\neq0}c_{j}\lambda_{j}\langle q_{0} \vert q_{j}\rangle=
d\lambda_{0}\langle q_{0}\vert q_{0}\rangle+\lambda_{0}\sum_{j\neq0}c_{j}\langle q_{0}\vert q_{j}\rangle
\f}
we get
\f{displaymath}
\langle q_{0}\vert d\mathbf{M}\vert q_{0}\rangle=d\lambda_{0}
\f}
So, using simple chain rules:
\f{displaymath}
\langle q_{0}\vert \frac{d\mathbf{M}}{dr_{k}^{\gamma}}\vert q_{0}\rangle
dr_{k}^{\gamma}=d\lambda_{0}
\f}
where here we used the notation \f$r_{k}^{\gamma}\f$ to denote an arbitrary position which can be
\f$\tilde x_0 ,\tilde y_0,\tilde z_0\f$ or   \f$\tilde x_1 ,\tilde y_1,\tilde z_1\f$
we get
\f{displaymath}
\langle q_{0}\vert \frac{d\mathbf{M}}{dr_{k}^{\gamma}}\vert q_{0}\rangle
=\frac{d\lambda_{0}}{dr_{k}^{\gamma}}
\f}

The derivatives of the matrix \f$\frac{d\mathbf{M}}{dr_{k}^{\gamma}} \f$ can be readily obtained via the
chain rule
\f{displaymath}
\frac{d\mathbf{M}}{dr_{k}^{\gamma}}=\sum_{\l}^{nat}\sum_{\alpha}^{x,y,z} \frac{d\mathbf{M}}{dP_{l}^{\alpha}}\frac{dP_{l}^{\alpha}}{dr_{k}^{\gamma}} +\\
\frac{d\mathbf{M}}{dM_{l}^{\alpha}}\frac{dM_{l}^{\alpha}}{dr_{k}^{\gamma}}
\f}

where \f$ M_{l}^{\alpha} \f$ corresponds to \f$ x_m^{l},y_m^{l},z_m^{l} \f$   and
\f$ P_{l}^{\alpha} \f$ corresponds to \f$ x_p^{l},y_p^{l},z_p^{l} \f$ according to the \f$ \alpha \f$ component.
*/


class Kearsley
{
  /// general log reference that needs to be initialized when constructed
  Log* log;
  /// position of atoms (first frame. In md is the running frame)
  std::vector<Vector> p0;
  /// position of atoms (second frame. In md is the  reference frame)
  std::vector<Vector> p1;
  /// alignment weight: the rmsd/msd that it provides is only based on this scalar
  std::vector<double> align;

  bool com0_is_removed;
  bool com1_is_removed;

public:
  /// error: the distance between two frames (might be rmsd/msd. See below)
  double err;
  /// displacement: the vector that goes from the p0 onto p1
  std::vector<Vector> diff0on1;
  /// displacement: the vector that goes from the p1 onto p0 (via inverse rotation)
  std::vector<Vector> diff1on0;

  /// center of mass of p0
  Vector com0;
  /// center of mass of p1
  Vector com1;
  /// position resetted wrt coms p0
  std::vector<Vector> p0reset;
  /// position resetted wrt coms p1
  std::vector<Vector> p1reset;
  /// position rotated: p0
  std::vector<Vector> p0rotated;
  /// position rotated: p1
  std::vector<Vector> p1rotated;
  /// rotation matrices p0 on p1 and reverse (p1 over p0)
  Tensor rotmat0on1,rotmat1on0;
  /// derivatives: derivative of the error respect p0
  std::vector<Vector> derrdp0;
  /// derivatives: derivative of the error respect p1
  std::vector<Vector> derrdp1;
  /// derivative of the rotation matrix
  /// note the dimension 3x3 x 3 x N
  std::vector<double> dmatdp0;
  std::vector<double> dmatdp1;

  /// constructor: need the two structure, the alignment vector and  the log reference
  Kearsley(  const std::vector<Vector> &p0, const std::vector<Vector> &p1,  const std::vector<double> &align , Log* &log);
  /// switch the assignment of the structure p0 (e.g. at each md step)
  void assignP0(const std::vector<Vector> & p0);
  /// derivatives: derivative of the error respect p1
  void assignP1(const std::vector<Vector> & p1);
  /// transfer the alignment vector
  void assignAlign(const std::vector<double> & align);
  /// finite differences of all the relevant quantities: takes a bool which decides if giving back rmsd or not (msd in this case)
  void finiteDifferenceInterface(bool rmsd);
  // this makes the real calculation: the rmsd bool decides wether doing rmsd or msd
  double calculate( bool rmsd );
};

}

#endif

