/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2023 The plumed team
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
#ifndef __PLUMED_cuda_ndReduction_h
#define __PLUMED_cuda_ndReduction_h
#include "plumed/tools/Vector.h"
#include "plumed/tools/Tensor.h"
#include <vector>
namespace PLMD{

/// @brief reduce the component of a nat x 3 x nthreads vector
/// @return a std::vector that contains the the reduced nat Vectors
/// We are using nat as the number of vector because usually this is used to
/// reduce the per atom derivatives
/** @brief reduce the component of a 3 x nat x N vector
 *  @param cudaNVectorAddress the pointer to the memory in cuda
 *  @param N the number of Tensors to reduce
 *  @param nat the number 
 *  @param maxNumThreads limits the number of threads per block to be used 
 *  @return the reduced Vector
 * 
 * reduceTensor expects that cudaNVectorAddress is initializated with
 * cudaMalloc(&cudaVectorAddress,  N * nat * 3 * sizeof(double)).
 * The componensts of cudaVectorAddress in memory shoud be organized like
 *  [x0, y0, z0, x1, y1, z1 ... x(N-1), y(N-1), z(N-1)]
 * 
 * The algorithm will decide the best number of threads to be used in 
 * an extra nthreads * 3 * sizeof(double) will be allocated on the GPU,
 * where nthreads is the total number of threads that will be used
 */
std::vector<Vector> reduceNVectors(double* cudaNVectorAddress, unsigned N, unsigned nat, unsigned maxNumThreads=512);

/** @brief reduce the component of a 3 x N vector
 *  @param cudaVectorAddress the pointer to the memory in cuda
 *  @param N the number of Vectors to reduce
 *  @param maxNumThreads limits the number of threads per block to be used 
 *  @return the reduced Vector
 * 
 * reduceTensor expects that cudaVectorAddress is initializated with
 * cudaMalloc(&cudaVectorAddress,  N * 3 * sizeof(double)).
 * The componensts of cudaVectorAddress in memory shoud be organized like
 *  [x0, y0, z0, x1, y1, z1 ... x(N-1), y(N-1), z(N-1)]
 * 
 * The algorithm will decide the best number of threads to be used in 
 * an extra nthreads * 3 * sizeof(double) will be allocated on the GPU,
 * where nthreads is the total number of threads that will be used
 */
Vector reduceVector(double* cudaVectorAddress, unsigned N, unsigned maxNumThreads=512);

/** @brief reduces a Tensor from 3x3 x N to 3x3
 *  @param cudaTensorAddress the pointer to the memory in cuda
 *  @param N the number of Tensors to reduce
 *  @param maxNumThreads limits the number of threads per block to be used 
 *  @return the reduced Tensor
 * 
 * reduceTensor expects that cudaTensorAddress is initializated with
 * cudaMalloc(&cudaTensorAddress,  N * 9 * sizeof(double)).
 * The componensts of cudaVectorAddress in memory shoud be stred in organized in N sequential blocks 
 * whose compontents shall be [(0,0) (0,1) (0,2) (1,0) (1,1) (1,2) (2,0) (2,1) (2,2)]
 * 
 * The algorithm will decide the best number of threads to be used in 
 * an extra nthreads * 9 * sizeof(double) will be allocated on the GPU,
 * where nthreads is the total number of threads that will be used
 */ 
Tensor reduceTensor(double* cudaTensorAddress, unsigned N, unsigned maxNumThreads=512);

/** @brief reduce the component of a N vector to a scalar
 *  @param cudaScalarAddress the pointer to the memory in cuda
 *  @param N the number of scalar to reduce
 *  @param maxNumThreads limits the number of threads per block to be used 
 *  @return the reduced scalar
 * 
 * reduceTensor expects that cudaVectorAddress is initializated with
 * cudaMalloc(&cudaVectorAddress,  N * sizeof(double)).
 * The componensts of cudaScalarAddress in memory shoud be organized like
 *  [x0, x1, x2 ...,  x(N-1)]
 * 
 * The algorithm will decide the best number of threads to be used in 
 * an extra nthreads * sizeof(double) will be allocated on the GPU,
 * where nthreads is the total number of threads that will be used
 */ 
double reduceScalar(double* cudaScalarAddress, unsigned N, unsigned maxNumThreads=512);
} //namespace PLMD
#endif //__PLUMED_cuda_ndReduction_h