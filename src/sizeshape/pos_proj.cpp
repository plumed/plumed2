/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2024 by Glen Hocky, New York University on behalf of authors

The sizeshape module is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The sizeshape module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "colvar/Colvar.h"
#include "core/ActionRegister.h"
#include "tools/Pbc.h"
#include "tools/File.h"           // Input and output from files 
#include "tools/Matrix.h"         // Linear Algebra operations
#include <sstream>
#include <cmath>
#include "tools/Communicator.h"   // All the MPI related stuffs

namespace PLMD {
namespace sizeshape {

//+PLUMEDOC COLVAR SIZESHAPE_POSITION_LINEAR_PROJ
/*
Calculates a linear projection in the space of a given reference configurational distribution in size-and-shape space.

The linear projection is given by:

$$
    l(\mathbf{x}) = \mathbf{v}\cdot(\mathbf{R}\cdot(\mathbf{x}(t) - \vec{{\zeta}}(t)) - \mathbf{\mu}),
$$

Where $\mathbf{v}$ is a vector of linear coefficients, $\mathbf{x}(t)$ is the configuration at time t, $\vec{\zeta}(t)$ is the difference in the geometric mean of the current configuration and that of the reference configuration $\mathbf{\mu}$. $\vec{\zeta}(t) = \frac{1}{N} \sum_{i=1}^N \vec{x_{i}}(t) - \frac{1}{N} \sum_{i=1}^N \vec{\mu_{i}}(t)$, for N atoms.

$\mathbf{R}$ is an optimal rotation matrix that minimizes the Mahalanobis distance between the current configuration and reference. $\mathbf{R}$ is obtained by using Kabsch algorithm within the code. The Mahalanobis distance is given as:

$$
d(\mathbf{x}, \mathbf{\mu}, \mathbf{\Sigma}) = \sqrt{(\mathbf{x}-\mathbf{\mu})^T \mathbf{\Sigma}^{-1} (\mathbf{x}-\mathbf{\mu})}
$$

where, $\mathbf{\Sigma}^{-1}$ is the $N\times N$ precision matrix. See also [SIZESHAPE_POSITION_MAHA_DIST](SIZESHAPE_POSITION_MAHA_DIST.md) for information about calculating Mahalanobis distance in size-and-shape space.

Size-and-shape Gaussian Mixture Model (shapeGMM) \cite Heidi-shapeGMM-2022 is a probabilistic clustering technique that is used to perform structural clusteing on ensemble of molecular configurations and to obtain reference $(\mathbf{\mu})$ and precision $(\mathbf{\Sigma}^{-1})$ corresponding to each of the cluster centers.
Please chcek out <a href="https://github.com/mccullaghlab/shapeGMMTorch">shapeGMMTorch-GitHub</a> and <a href="https://pypi.org/project/shapeGMMTorch/"> shapeGMMTorch-PyPI</a> for examples and informations on preforming shapeGMM clustering.

##Examples

In the following example, a group is defined with atom indices of selected atoms and then linear projection is calculated for the given reference, precision and coefficients. Each file is a space separated list of 3N floating point numbers.

```plumed
#SETTINGS INPUTFILES=regtest/sizeshape/rt-sizeshape/global_avg.txt
#SETTINGS INPUTFILES=regtest/sizeshape/rt-sizeshape/global_precision.txt
#SETTINGS INPUTFILES=regtest/sizeshape/rt-sizeshape/ld1_scalings.txt

UNITS LENGTH=A TIME=ps ENERGY=kcal/mol
GROUP ATOMS=18,20,22,31,33,35,44,46,48,57,59,61,70,72,74,83,85,87,96,98,100,109,111 LABEL=ga_list
proj: SIZESHAPE_POSITION_LINEAR_PROJ ...
   REFERENCE=regtest/sizeshape/rt-sizeshape/global_avg.txt
   PRECISION=regtest/sizeshape/rt-sizeshape/global_precision.txt
   COEFFS=regtest/sizeshape/rt-sizeshape/ld1_scalings.txt
   GROUP=ga_list
...

PRINT ARG=proj STRIDE=1 FILE=COLVAR FMT=%8.8f
```

*/
//+ENDPLUMEDOC


class position_linear_proj : public Colvar {

private:
  bool pbc, serial;
  std::string prec_f_name;      		// precision file name
  std::string ref_f_name;       		// reference file name
  std::string coeffs_f_name; 			// file containing linear coeffs
  IFile in_;             			// create an object of class IFile
  Log out_;
  Matrix<double> ref_str;       	        // coords of reference
  Matrix<double> mobile_str;    		// coords of mobile
  Matrix<double> prec;        			// precision data
  Matrix<double> rotation;
  std::vector<double> linear_coeffs;            // Linear Coefficients
  Matrix<double> derv_numeric;
  void readinputs();            		// reads the input data
  double proj;					// projection value
  std::vector<AtomNumber> atom_list;            // list of atoms
  const double SMALL = 1.0E-30;
  const double delta = 0.00001;
public:
  static void registerKeywords( Keywords& keys );
  explicit position_linear_proj(const ActionOptions&);
  double determinant(int n, const std::vector<std::vector<double>>* B);
  void kabsch_rot_mat();   		// gives rotation matrix
  double cal_position_linear_proj();    // calculates the linear projection
  void numeric_grad();        		// calculates the numeric gradient
  // active methods:
  void calculate() override;
};

PLUMED_REGISTER_ACTION(position_linear_proj, "SIZESHAPE_POSITION_LINEAR_PROJ")

// register keywords function
void position_linear_proj::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords( keys );
  keys.add("compulsory", "PRECISION", "Precision Matrix (inverse of covariance)." );
  keys.add("compulsory", "REFERENCE", "Coordinates of the reference structure.");
  keys.add("atoms","GROUP","Group of atoms being used");
  keys.add("compulsory", "COEFFS", "Vector of linear coefficients.");
  keys.addFlag("SERIAL",false,"Perform the calculation in serial, for debug purposes only.");
  keys.setValueDescription("scalar","the linear projection");
}

// constructor function
position_linear_proj::position_linear_proj(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc(true),
  serial(false),
  proj(0),
  prec_f_name(""),
  ref_f_name(""),
  coeffs_f_name("") { // Note! no comma here in the last line.
  parseFlag("SERIAL",serial);
  parseAtomList("GROUP",atom_list);
  parse("REFERENCE", ref_f_name);
  parse("PRECISION", prec_f_name);
  parse("COEFFS", coeffs_f_name);
  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  checkRead();

  log.printf("  of %u atoms\n",static_cast<unsigned>(atom_list.size()));
  for(unsigned int i=0; i<atom_list.size(); ++i) {
    log.printf("  %d", atom_list[i].serial());
  }

  if(pbc) {
    log.printf("\n using periodic boundary conditions\n");
  } else {
    log.printf("\n without periodic boundary conditions\n");
  }

  addValueWithDerivatives();
  setNotPeriodic();

  requestAtoms(atom_list);

  // call the readinputs() function here
  readinputs();

}

// read inputs function
void position_linear_proj::readinputs() {
  unsigned N=getNumberOfAtoms();
  // read ref coords
  in_.open(ref_f_name);

  ref_str.resize(N,3);
  prec.resize(N,N);
  derv_numeric.resize(N,3);

  std::string line_, val_;
  unsigned c_=0;

  while (c_ < N) {
    in_.getline(line_);
    std::vector<std::string> items_;
    std::stringstream check_(line_);

    while(std::getline(check_, val_, ' ')) {
      items_.push_back(val_);
    }
    for(int i=0; i<3; ++i) {
      ref_str(c_,i) = std::stold(items_[i]);
    }
    c_ += 1;
  }
  in_.close();

  //read precision
  in_.open(prec_f_name);

  std::string line, val;
  unsigned int c = 0;

  while(c < N) {
    in_.getline(line);

    // vector for storing the objects
    std::vector<std::string> items;

    // stringstream helps to treat a string like an ifstream!
    std::stringstream check(line);

    while (std::getline(check, val, ' ')) {
      items.push_back(val);
    }

    for(unsigned int i=0; i<N; ++i) {
      prec(c, i) = std::stold(items[i]);
    }

    c += 1;

  }
  in_.close();

  // read in the linear coeffs
  in_.open(coeffs_f_name);
  unsigned n_=0;
  std::string l_;
  while (n_ < N*3) {
    in_.getline(l_);
    linear_coeffs.push_back(std::stod(l_));
    n_ += 1;
  }
  linear_coeffs.resize(N*3);

  in_.close();

}



double position_linear_proj::determinant(int n, const std::vector<std::vector<double>>* B) {

  std::vector<std::vector<double>> A(n, std::vector<double>(n, 0));
  // make a copy first!
  for(int i=0; i<n; ++i) {
    for(int j=0; j<n; ++j) {
      A[i][j] = (*B)[i][j];
    }
  }


  //  It calculates determinant of a matrix using partial pivoting.

  double det = 1;

  // Row operations for i = 0, ,,,, n - 2 (n-1 not needed)
  for ( int i = 0; i < n - 1; i++ ) {
    // Partial pivot: find row r below with largest element in column i
    int r = i;
    double maxA = std::abs( A[i][i] );
    for ( int k = i + 1; k < n; k++ ) {
      double val = std::abs( A[k][i] );
      if ( val > maxA ) {
        r = k;
        maxA = val;
      }
    }
    if ( r != i ) {
      for ( int j = i; j < n; j++ ) {
        std::swap( A[i][j], A[r][j] );
      }
      det = -det;
    }

    // Row operations to make upper-triangular
    double pivot = A[i][i];
    if (std::abs( pivot ) < SMALL ) {
      return 0.0;  // Singular matrix
    }

    for ( int r = i + 1; r < n; r++ ) {                  // On lower rows
      double multiple = A[r][i] / pivot;                // Multiple of row i to clear element in ith column
      for ( int j = i; j < n; j++ ) {
        A[r][j] -= multiple * A[i][j];
      }
    }
    det *= pivot;                                        // Determinant is product of diagonal
  }

  det *= A[n-1][n-1];

  return det;
}

// kabsch rotation
void position_linear_proj::kabsch_rot_mat() {

  unsigned N=getNumberOfAtoms();

  Matrix<double> mobile_str_T(3,N);
  Matrix<double> prec_dot_ref_str(N,3);
  Matrix<double> correlation(3,3);


  transpose(mobile_str, mobile_str_T);
  mult(prec, ref_str, prec_dot_ref_str);
  mult(mobile_str_T, prec_dot_ref_str, correlation);


  int rw = correlation.nrows();
  int cl = correlation.ncols();
  int sz = rw*cl;

  // SVD part (taking from plu2/src/tools/Matrix.h: pseudoInvert function)

  std::vector<double> da(sz);
  unsigned k=0;

  // Transfer the matrix to the local array
  for (int i=0; i<cl; ++i)
    for (int j=0; j<rw; ++j) {
      da[k++]=static_cast<double>( correlation(j,i) );  // note! its [j][i] not [i][j]
    }

  int nsv, info, nrows=rw, ncols=cl;
  if(rw>cl) {
    nsv=cl;
  } else {
    nsv=rw;
  }

  // Create some containers for stuff from single value decomposition
  std::vector<double> S(nsv);
  std::vector<double> U(nrows*nrows);
  std::vector<double> VT(ncols*ncols);
  std::vector<int> iwork(8*nsv);

  // This optimizes the size of the work array used in lapack singular value decomposition
  int lwork=-1;
  std::vector<double> work(1);
  plumed_lapack_dgesdd( "A", &nrows, &ncols, da.data(), &nrows, S.data(), U.data(), &nrows, VT.data(), &ncols, work.data(), &lwork, iwork.data(), &info );
  //if(info!=0) return info;
  if(info!=0) {
    log.printf("info:", info);
  }

  // Retrieve correct sizes for work and rellocate
  lwork=(int) work[0];
  work.resize(lwork);

  // This does the singular value decomposition
  plumed_lapack_dgesdd( "A", &nrows, &ncols, da.data(), &nrows, S.data(), U.data(), &nrows, VT.data(), &ncols, work.data(), &lwork, iwork.data(), &info );
  //if(info!=0) return info;
  if(info!=0) {
    log.printf("info:", info);
  }


  // get U and VT in form of 2D vector (U_, VT_)
  std::vector<std::vector<double>> U_(nrows, std::vector<double>(nrows,0));
  std::vector<std::vector<double>> VT_(ncols, std::vector<double>(ncols,0));

  int  c=0;

  for(int i=0; i<nrows; ++i) {
    for(int j=0; j<nrows; ++j) {
      U_[j][i] = U[c];
      c += 1;
    }
  }
  c = 0; // note! its [j][i] not [i][j]
  for(int i=0; i<ncols; ++i) {
    for(int j=0; j<ncols; ++j) {
      VT_[j][i] = VT[c];
      c += 1;
    }
  }
  c=0; // note! its [j][i] not [i][j]


  // calculate determinants
  double det_u = determinant(nrows, &U_);
  double det_vt = determinant(ncols, &VT_);

  // check!
  if (det_u * det_vt < 0.0) {
    for(int i=0; i<nrows; ++i) {
      U_[i][nrows-1] *= -1;
    }
  }


  //Matrix<double> rotation(3,3);
  rotation.resize(3,3);
  Matrix<double> u(3,3), vt(3,3);
  for(int i=0; i<3; ++i) {
    for(int j=0; j<3; ++j) {
      u(i,j)=U_[i][j];
      vt(i,j)=VT_[i][j];
    }
  }

  // get rotation matrix
  mult(u, vt, rotation);

}


// calculates linear projection
double position_linear_proj::cal_position_linear_proj() {

  unsigned N=getNumberOfAtoms();

  Matrix<double> rotated_obj(N,3);
  // rotate the object
  mult(mobile_str, rotation, rotated_obj);

  // compute the displacement
  std::vector<double> disp(N*3);
  unsigned c=0;
  for(unsigned int i=0; i<N; ++i) {
    for(int j=0; j<3; ++j) {
      disp[c] = (rotated_obj(i,j)-ref_str(i,j));
      c+=1;
    }
  }

  //double proj_val = dotProduct(disp, linear_coeffs);
  double proj_val = 0.0;
  for(unsigned int i=0; i<N*3; ++i) {
    proj_val += (linear_coeffs[i]*disp[i]);
  }

  return proj_val;
}

// numeric gradient
void position_linear_proj::numeric_grad() {
  // This function performs numerical derivative.
  unsigned N=getNumberOfAtoms();

  unsigned stride;
  unsigned rank;
  if(serial) {
    // when using components the parallelisation do not work
    stride=1;
    rank=0;
  } else {
    stride=comm.Get_size();
    rank=comm.Get_rank();
  }

  for(unsigned i=rank; i<N; i+=stride) {
    for (unsigned j=0; j<3; ++j) {

      mobile_str(i,j) += delta;
      kabsch_rot_mat();
      derv_numeric(i,j) = ((cal_position_linear_proj() - proj)/delta);

      mobile_str(i,j) -= delta;
    }

  }

  if(!serial) {
    if(!derv_numeric.getVector().empty()) {
      comm.Sum(&derv_numeric(0,0), derv_numeric.getVector().size());
    }
  }


  for(unsigned i=0; i<N; ++i) {
    Vector vi(derv_numeric(i,0), derv_numeric(i,1), derv_numeric(i,2) );
    setAtomsDerivatives(i, vi);
  }

  // clear the matrix (very important step!!)
  derv_numeric *= 0;
}


// calculator
void position_linear_proj::calculate() {

  if(pbc) {
    makeWhole();
  }
  unsigned N=getNumberOfAtoms();

  mobile_str.resize(N,3);

  // load the mobile str
  for(unsigned int i=0; i<N; ++i) {
    Vector pos=getPosition(i);  // const PLMD::Vector
    for(unsigned j=0; j<3; ++j) {
      mobile_str(i,j) = pos[j];
    }
  }

  // translating the structure to center of geometry
  double center_of_geometry[3]= {0.0, 0.0, 0.0};

  for(unsigned int i=0; i<N; ++i) {
    center_of_geometry[0] += mobile_str(i,0);
    center_of_geometry[1] += mobile_str(i,1);
    center_of_geometry[2] += mobile_str(i,2);
  }

  for(unsigned int i=0; i<N; ++i) {
    for(int j=0; j<3; ++j) {
      mobile_str(i,j) -= (center_of_geometry[j]/N);
    }
  }

  kabsch_rot_mat();
  proj = cal_position_linear_proj();

  numeric_grad();
  setBoxDerivativesNoPbc();
  setValue(proj);


}

}
}



