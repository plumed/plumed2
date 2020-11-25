/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2020 The plumed team
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
#include "FlexibleBin.h"
#include "ActionWithArguments.h"
#include <cmath>
#include <iostream>
#include <vector>
#include "tools/Matrix.h"

namespace PLMD {


FlexibleBin::FlexibleBin(int type, ActionWithArguments *paction, double const &d, std::vector<double> &smin, const std::vector<double> &smax):
  type(type),
  paction(paction),
  sigma(d),
  sigmamin(smin),
  sigmamax(smax)
{
  // initialize the averages and the variance matrices
  if(type==diffusion) {
    unsigned ncv=paction->getNumberOfArguments();
    std::vector<double> average(ncv*(ncv+1)/2);
    std::vector<double> variance(ncv*(ncv+1)/2);
  }
  paction->log<<"  Limits for sigmas using adaptive hills:  \n";
  for(unsigned i=0; i<paction->getNumberOfArguments(); ++i) {
    paction->log<<"   CV  "<<paction->getPntrToArgument(i)->getName()<<":\n";
    if(sigmamin[i]>0.) {
      limitmin.push_back(true);
      paction->log<<"       Min "<<sigmamin[i];
      sigmamin[i]*=sigmamin[i];	// this is because the matrix which is calculated is the sigmasquared
    } else {
      limitmin.push_back(false);
      paction->log<<"       Min No ";
    }
    if(sigmamax[i]>0.) {
      limitmax.push_back(true);
      paction->log<<"       Max "<<sigmamax[i];
      sigmamax[i]*=sigmamax[i];
    } else {
      limitmax.push_back(false);
      paction->log<<"       Max No ";
    }
    paction->log<<" \n";
  }

}

/// Constructure for 1D FB for PBMETAD
FlexibleBin::FlexibleBin(int type, ActionWithArguments *paction, unsigned iarg,
                         double const &d, std::vector<double> &smin, const std::vector<double> &smax):
  type(type),paction(paction),sigma(d),sigmamin(smin),sigmamax(smax)
{
  // initialize the averages and the variance matrices
  if(type==diffusion) {
    std::vector<double> average(1);
    std::vector<double> variance(1);
  }
  paction->log<<"  Limits for sigmas using adaptive hills:  \n";
  paction->log<<"   CV  "<<paction->getPntrToArgument(iarg)->getName()<<":\n";
  if(sigmamin[0]>0.) {
    limitmin.push_back(true);
    paction->log<<"       Min "<<sigmamin[0];
    sigmamin[0]*=sigmamin[0];
  } else {
    limitmin.push_back(false);
    paction->log<<"       Min No ";
  }
  if(sigmamax[0]>0.) {
    limitmax.push_back(true);
    paction->log<<"       Max "<<sigmamax[0];
    sigmamax[0]*=sigmamax[0];
  } else {
    limitmax.push_back(false);
    paction->log<<"       Max No ";
  }
  paction->log<<" \n";
}

/// Update the flexible bin
/// in case of diffusion based: update at every step
/// in case of gradient based: update only when you add the hill
void FlexibleBin::update(bool nowAddAHill) {
  unsigned ncv=paction->getNumberOfArguments();
  unsigned dimension=ncv*(ncv+1)/2;
  std::vector<double> delta;
  std::vector<double> cv;
  double decay=1./sigma;
  // this is done all the times from scratch. It is not an accumulator
  // here update the flexible bin according to the needs
  switch (type) {
  // This should be called every time
  case diffusion:
    // if you use this below then the decay is in time units
    //double decay=paction->getTimeStep()/sigma;
    // to be consistent with the rest of the program: everything is better to be in timesteps
    // THE AVERAGE VALUE
    // beware: the pbc
    delta.resize(ncv);
    for(unsigned i=0; i<ncv; i++) cv.push_back(paction->getArgument(i));
    if(average.size()==0) { // initial time: just set the initial vector
      average.resize(ncv);
      for(unsigned i=0; i<ncv; i++) average[i]=cv[i];
    } else { // accumulate
      for(unsigned i=0; i<ncv; i++) {
        delta[i]=paction->difference(i,average[i],cv[i]);
        average[i]+=decay*delta[i];
        average[i]=paction->bringBackInPbc(i,average[i]); // equation 8 of "Metadynamics with adaptive Gaussians"
      }
    }
    // THE VARIANCE
    if(variance.size()==0) {
      variance.resize(dimension,0.); // nonredundant members dimension=ncv*(ncv+1)/2;
    } else {
      unsigned k=0;
      for(unsigned i=0; i<ncv; i++) {
        for(unsigned j=i; j<ncv; j++) { // upper diagonal loop
          variance[k]+=decay*(delta[i]*delta[j]-variance[k]);
          k++;
        }
      }
    }
    break;
  case geometry:
    //this calculates in variance the \nabla CV_i \dot \nabla CV_j
    variance.resize(dimension);
    // now the signal for retrieving the gradients should be already given by checkNeedsGradients.
    // here just do the projections
    // note that the call  checkNeedsGradients() in BiasMetaD takes care of switching on the call to gradients
    if (nowAddAHill) { // geometry is sync with hill deposition
      unsigned k=0;
      for(unsigned i=0; i<ncv; i++) {
        for(unsigned j=i; j<ncv; j++) {
          // eq 12 of "Metadynamics with adaptive Gaussians"
          variance[k]=sigma*sigma*(paction->getProjection(i,j));
          k++;
        }
      }
    }
    break;
  default:
    plumed_merror("This flexible bin method is not recognized");
  }
}

std::vector<double> FlexibleBin::getMatrix() const {
  return variance;
}

/// Update the flexible bin for PBMetaD like FlexBin
/// in case of diffusion based: update at every step
/// in case of gradient based: update only when you add the hill
void FlexibleBin::update(bool nowAddAHill, unsigned iarg) {
  // this is done all the times from scratch. It is not an accumulator
  // here update the flexible bin according to the needs
  std::vector<double> cv;
  std::vector<double> delta;
  // if you use this below then the decay is in time units
  // to be consistent with the rest of the program: everything is better to be in timesteps
  double decay=1./sigma;
  switch (type) {
  // This should be called every time
  case diffusion:
    // THE AVERAGE VALUE
    delta.resize(1);
    cv.push_back(paction->getArgument(iarg));
    if(average.size()==0) { // initial time: just set the initial vector
      average.resize(1);
      average[0]=cv[0];
    } else { // accumulate
      delta[0]=paction->difference(iarg,average[0],cv[0]);
      average[0]+=decay*delta[0];
      average[0]=paction->bringBackInPbc(iarg,average[0]); // equation 8 of "Metadynamics with adaptive Gaussians"
    }
    // THE VARIANCE
    if(variance.size()==0) {
      variance.resize(1,0.); // nonredundant members dimension=ncv*(ncv+1)/2;
    } else {
      variance[0]+=decay*(delta[0]*delta[0]-variance[0]);
    }
    break;
  case geometry:
    //this calculates in variance the \nabla CV_i \dot \nabla CV_j
    variance.resize(1);
    // now the signal for retrieving the gradients should be already given by checkNeedsGradients.
    // here just do the projections
    // note that the call  checkNeedsGradients() in BiasMetaD takes care of switching on the call to gradients
    if (nowAddAHill) { // geometry is sync with hill deposition
      // eq 12 of "Metadynamics with adaptive Gaussians"
      variance[0]=sigma*sigma*(paction->getProjection(iarg,iarg));
    }
    break;
  default:
    plumed_merror("This flexible bin is not recognized");
  }
}

///
/// Calculate the matrix of  (dcv_i/dx)*(dcv_j/dx)^-1
/// that is needed for the metrics in metadynamics
///
///
std::vector<double> FlexibleBin::getInverseMatrix() const {
  unsigned ncv=paction->getNumberOfArguments();
  Matrix<double> matrix(ncv,ncv);
  unsigned i,j,k;
  k=0;
  // place the matrix in a complete matrix for compatibility
  for (i=0; i<ncv; i++) {
    for (j=i; j<ncv; j++) {
      matrix(j,i)=matrix(i,j)=variance[k];
      k++;
    }
  }
#define NEWFLEX
#ifdef NEWFLEX
  // diagonalize to impose boundaries (only if boundaries are set)
  Matrix<double>      eigenvecs(ncv,ncv);
  std::vector<double> eigenvals(ncv);

  //eigenvecs: first is eigenvec number, second is eigenvec component
  if(diagMat( matrix, eigenvals, eigenvecs )!=0) {plumed_merror("diagonalization in FlexibleBin failed! This matrix is weird\n");};

  for (i=0; i<ncv; i++) { //loop on the dimension
    if( limitmax[i] ) {
      //limit every  component that is larger
      for (j=0; j<ncv; j++) { //loop on components
        if(std::pow(eigenvals[j]*eigenvecs[j][i],2)>std::pow(sigmamax[i],2) ) {
          eigenvals[j]=std::sqrt(std::pow(sigmamax[i]/(eigenvecs[j][i]),2))*copysign(1.,eigenvals[j]);
        }
      }
    }
  }
  for (i=0; i<ncv; i++) { //loop on the dimension
    // find the largest one:  if it is smaller than min  then rescale
    if( limitmin[i] ) {
      unsigned imax=0;
      double fmax=-1.e10;
      for (j=0; j<ncv; j++) { //loop on components
        double fact=std::pow(eigenvals[j]*eigenvecs[j][i],2);
        if(fact>fmax) {
          fmax=fact; imax=j;
        }
      }
      if(fmax<std::pow(sigmamin[i],2) ) {
        eigenvals[imax]=std::sqrt(std::pow(sigmamin[i]/(eigenvecs[imax][i]),2))*copysign(1.,eigenvals[imax]);
      }
    }
  }

  // normalize eigenvecs
  Matrix<double> newinvmatrix(ncv,ncv);
  for (i=0; i<ncv; i++) {
    for (j=0; j<ncv; j++) {
      newinvmatrix[j][i]=eigenvecs[j][i]/eigenvals[j];
    }
  }

  std::vector<double> uppervec(ncv*(ncv+1)/2);
  k=0;
  for (i=0; i<ncv; i++) {
    for (j=i; j<ncv; j++) {
      double scal=0;
      for(unsigned l=0; l<ncv; ++l) {
        scal+=eigenvecs[l][i]*newinvmatrix[l][j];
      }
      uppervec[k]=scal; k++;
    }
  }
#else
  // get the inverted matrix
  Matrix<double> invmatrix(ncv,ncv);
  Invert(matrix,invmatrix);
  std::vector<double> uppervec(ncv*(ncv+1)/2);
  // upper diagonal of the inverted matrix (that is symmetric)
  k=0;
  for (i=0; i<ncv; i++) {
    for (j=i; j<ncv; j++) {
      uppervec[k]=invmatrix(i,j);
      k++;
    }
  }
#endif
  return uppervec;
}

///
/// Calculate the matrix of  (dcv_i/dx)*(dcv_j/dx)^-1
/// that is needed for the metrics in metadynamics
/// for PBMetaD like FlexBin
///
std::vector<double> FlexibleBin::getInverseMatrix(unsigned iarg) const {
  // diagonalize to impose boundaries (only if boundaries are set)
  std::vector<double> eigenvals(1, variance[0]);
  if( limitmax[0] ) {
    if(eigenvals[0]>sigmamax[0]) {
      eigenvals[0]=sigmamax[0];
    }
  }
  // find the largest one:  if it is smaller than min  then rescale
  if( limitmin[0] ) {
    double fmax=-1.e10;
    double fact=eigenvals[0];
    if(fact>fmax) {
      fmax=fact;
    }
    if(fmax<sigmamin[0]) {
      eigenvals[0]=sigmamin[0];
    }
  }
  std::vector<double> uppervec(1,1./eigenvals[0]);

  return uppervec;
}

}
