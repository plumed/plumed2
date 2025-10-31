/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2023 The plumed team
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
#include <iostream>
#include <complex>
#include "gridtools/ActionWithGrid.h"
#include "core/ActionRegister.h"
#ifdef __PLUMED_HAS_FFTW
#include <fftw3.h> // FFTW interface
#endif

namespace PLMD {
namespace fourier {

//+PLUMEDOC GRIDANALYSIS FOURIER_TRANSFORM
/*
Compute the Discrete Fourier Transform (DFT) by means of FFTW of data stored on a 2D grid.

This action takes a function on a two-dimensional grid as input and computes a fourier transform upon the input function using [FFTW](https://www.fftw.org).
Currently, this actions performs a complex fourier transition even if the input data is real.  The functionality here was developed and used in the paper cited
below. The following input was used in that paper:

```plumed
UNITS NATURAL

# These two commands calculate one symmetry function for each atom.  These
# symmetry functions tell us whether the environment around each atom resembles
# the environment in the solid or the environment in the liquid.
fcc: FCCUBIC SPECIES=1-20736 SWITCH={CUBIC D_0=1.2 D_MAX=1.5} ALPHA=27
smapfcc: MORE_THAN ARG=fcc SWITCH={SMAP R_0=0.5 A=8 B=8}
smapfcc_grp: GROUP ATOMS=1-20736
# This calculates the position of the center of the solid region of the simulation box.  What we are computing here a weighted average position
# the weights are the order parameters computed using the two commands above.
center: CENTER PHASES ATOMS=fcc WEIGHTS=smapfcc
# This calculates the phase field that tells us whether the structure in each part of the simulation box is solid-like or liquid like.
dens: MULTICOLVARDENS DATA=smapfcc ORIGIN=center DIR=xyz NBINS=50,80,80 BANDWIDTH=1.0,1.0,1.0 GRID_MIN=0.0,auto,auto GRID_MAX=20.0,auto,auto STRIDE=1 CLEAR=1
# This finds the instantaneous location of the interface between the solid and liquid phases
contour: FIND_CONTOUR_SURFACE ARG=dens CONTOUR=0.5 SEARCHDIR=dens_dist.x
DUMPGRID ARG=contour FILE=contour.dat
# This does the fourier transform of the location of the interface.  We can extract the interfacial stiffness from the average of this fourier transform
ft: FOURIER_TRANSFORM ARG=contour FT_TYPE=norm FOURIER_PARAMETERS=-1,1
DUMPGRID ARG=ft FILE=fourier.dat STRIDE=10
```

We do not think this action has been used in any other paper. If you are interested in using the functionality here we would recommend you carefully
read the documentation for the FFTW library.  You may find it necessary to modify the code in this action for your particular purpose.

To what FFTW computes we would recommend reading [this page](http://www.fftw.org/doc/What-FFTW-Really-Computes.html#What-FFTW-Really-Computes).

Notice that the keyword "FOURIER_PARAMETERS" specifies how the Fourier transform will be normalized. The keyword takes two numerical parameters $a$ and $b$ that are both set equal to one by default.
The normalization is then defined as:

$$
\frac{1}{n^{(1-a)/2}} \sum_{j=0}^{n-1} X_j e^{2\pi b\, j k \sqrt{-1}/n}
$$

*/
//+ENDPLUMEDOC


class FourierTransform : public gridtools::ActionWithGrid {
private:
  bool firsttime;
  std::string output_type;
  bool real_output, store_norm;
  std::vector<int> fourier_params;
  gridtools::GridCoordinatesObject gridcoords;
public:
  static void registerKeywords( Keywords& keys );
  explicit FourierTransform(const ActionOptions&ao);
  unsigned getNumberOfDerivatives() override ;
  const gridtools::GridCoordinatesObject& getGridCoordinatesObject() const override ;
  std::vector<std::string> getGridCoordinateNames() const override ;
  void calculate() override ;
};

PLUMED_REGISTER_ACTION(FourierTransform,"FOURIER_TRANSFORM")

void FourierTransform::registerKeywords( Keywords& keys ) {
  ActionWithGrid::registerKeywords( keys );
  keys.addInputKeyword("compulsory","ARG","grid","the label of the grid that you want to fourer transform");
  keys.add("optional","FT_TYPE","choose what kind of data you want as output on the grid. Possible values are: ABS = compute the complex modulus of Fourier coefficients (DEFAULT); NORM = compute the norm (i.e. ABS^2) of Fourier coefficients; COMPLEX = store the FFTW complex output on the grid (as a vector).");
  keys.add("compulsory","FOURIER_PARAMETERS","default","what kind of normalization is applied to the output and if the Fourier transform in FORWARD or BACKWARD. This keyword takes the form FOURIER_PARAMETERS=A,B, where A and B can be 0, 1 or -1. The default values are A=1 (no normalization at all) and B=1 (forward FFT). Other possible choices for A are: "
           "A=-1: normalize by the number of data, "
           "A=0: normalize by the square root of the number of data (one forward and followed by backward FFT recover the original data). ");
  keys.addOutputComponent("real","FT_TYPE","grid","the real part of the function");
  keys.addOutputComponent("imag","FT_TYPE","grid","the imaginary part of the function");
  keys.setValueDescription("grid","the fourier transform of the input grid");
  keys.addDOI("10.1088/1361-648X/aa893d");
}

FourierTransform::FourierTransform(const ActionOptions&ao):
  Action(ao),
  ActionWithGrid(ao),
  firsttime(true),
  real_output(true),
  store_norm(false),
  fourier_params(2) {
  if( getPntrToArgument(0)->getRank()!=2 ) {
    error("fourier transform currently only works with two dimensional grids");
  }

  // Get the type of FT
  parse("FT_TYPE",output_type);
  if (output_type.length()==0) {
    log<<"  keyword FT_TYPE unset. By default output grid will contain REAL Fourier coefficients\n";
  } else if ( output_type=="ABS" || output_type=="abs") {
    log << "  keyword FT_TYPE is '"<< output_type << "' : will compute the MODULUS of Fourier coefficients\n";
  } else if ( output_type=="NORM" || output_type=="norm") {
    log << "  keyword FT_TYPE is '"<< output_type << "' : will compute the NORM of Fourier coefficients\n";
    store_norm=true;
  } else if ( output_type=="COMPLEX" || output_type=="complex" ) {
    log<<"  keyword FT_TYPE is '"<< output_type <<"' : output grid will contain the COMPLEX Fourier coefficients\n";
    real_output=false;
  } else {
    error("keyword FT_TYPE unrecognized!");
  }

  // Normalize output?
  std::string params_str;
  parse("FOURIER_PARAMETERS",params_str);
  if (params_str=="default") {
    fourier_params.assign( fourier_params.size(), 1 );
    log.printf("  default values of Fourier parameters A=%i, B=%i : the output will NOT be normalized and BACKWARD Fourier transform is computed \n", fourier_params[0],fourier_params[1]);
  } else {
    std::vector<std::string> fourier_str = Tools::getWords(params_str, "\t\n ,");
    if (fourier_str.size()>2) {
      error("FOURIER_PARAMETERS can take just two values");
    }
    for (unsigned i=0; i<fourier_str.size(); ++i) {
      Tools::convert(fourier_str[i],fourier_params[i]);
      if (fourier_params[i]>1 || fourier_params[i]<-1) {
        error("values accepted for FOURIER_PARAMETERS are only -1, 1 or 0");
      }
    }
    log.printf("  Fourier parameters are A=%i, B=%i \n", fourier_params[0],fourier_params[1]);
  }

  std::vector<std::size_t> shape( getPntrToArgument(0)->getRank() );
  if (real_output) {
    addValueWithDerivatives( shape );
  } else {
    addComponentWithDerivatives( "real", shape );
    addComponentWithDerivatives( "imag", shape );
  }

  unsigned dimension = getPntrToArgument(0)->getRank();
  gridtools::ActionWithGrid* ag=dynamic_cast<gridtools::ActionWithGrid*>( getPntrToArgument(0)->getPntrToAction() );
  if( !ag ) {
    error("input action should be a grid");
  }
  const gridtools::GridCoordinatesObject & gcoords( ag->getGridCoordinatesObject() );
  if( gcoords.getGridType()=="fibonacci" ) {
    error("cannot fourier transform fibonacci grids");
  }
  std::vector<bool> ipbc( dimension );
  for(unsigned i=0; i<dimension; ++i) {
    ipbc[i] = gcoords.isPeriodic(i);
  }
  gridcoords.setup( "flat", ipbc, 0, 0.0 );
  checkRead();
#ifndef __PLUMED_HAS_FFTW
  error("this feature is only available if you compile PLUMED with FFTW");
#endif
}

unsigned FourierTransform::getNumberOfDerivatives() {
  return 2;
}

const gridtools::GridCoordinatesObject& FourierTransform::getGridCoordinatesObject() const {
  return gridcoords;
}

std::vector<std::string> FourierTransform::getGridCoordinateNames() const {
  gridtools::ActionWithGrid* ag=dynamic_cast<gridtools::ActionWithGrid*>( getPntrToArgument(0)->getPntrToAction() );
  return ag->getGridCoordinateNames();
}

void FourierTransform::calculate() {
  if( firsttime ) {
    gridtools::ActionWithGrid* ag=dynamic_cast<gridtools::ActionWithGrid*>( getPntrToArgument(0)->getPntrToAction() );
    const gridtools::GridCoordinatesObject & gcoords( ag->getGridCoordinatesObject() );
    std::vector<double> fspacing;
    std::vector<std::size_t> snbins( getGridCoordinatesObject().getDimension() );
    std::vector<std::string> smin( gcoords.getDimension() ), smax( gcoords.getDimension() );
    for(unsigned i=0; i<getGridCoordinatesObject().getDimension(); ++i) {
      smin[i]=gcoords.getMin()[i];
      smax[i]=gcoords.getMax()[i];
      // Compute k-grid extents
      double dmin, dmax;
      snbins[i]=gcoords.getNbin(false)[i];
      Tools::convert(smin[i],dmin);
      Tools::convert(smax[i],dmax);
      dmax=2.0*pi*snbins[i]/( dmax - dmin );
      dmin=0.0;
      Tools::convert(dmin,smin[i]);
      Tools::convert(dmax,smax[i]);
    }
    gridcoords.setBounds( smin, smax, snbins, fspacing );
    firsttime=false;
    for(unsigned i=0; i<getNumberOfComponents(); ++i) {
      getPntrToComponent(i)->setShape( gcoords.getNbin(true) );
    }
  }

#ifdef __PLUMED_HAS_FFTW
  // *** CHECK CORRECT k-GRID BOUNDARIES ***
  //log<<"Real grid boundaries: \n"
  //    <<"  min_x: "<<mygrid->getMin()[0]<<"  min_y: "<<mygrid->getMin()[1]<<"\n"
  //    <<"  max_x: "<<mygrid->getMax()[0]<<"  max_y: "<<mygrid->getMax()[1]<<"\n"
  //    <<"K-grid boundaries:"<<"\n"
  //    <<"  min_x: "<<ft_min[0]<<"  min_y: "<<ft_min[1]<<"\n"
  //    <<"  max_x: "<<ft_max[0]<<"  max_y: "<<ft_max[1]<<"\n";

  // Get the size of the input data arrays (to allocate FFT data)
  std::vector<std::size_t> N_input_data( gridcoords.getNbin(true) );
  size_t fft_dimension=1;
  for(unsigned i=0; i<N_input_data.size(); ++i) {
    fft_dimension*=static_cast<size_t>( N_input_data[i] );
  }
  // FFT arrays
  std::vector<std::complex<double> > input_data(fft_dimension), fft_data(fft_dimension);

  // Fill real input with the data on the grid
  Value* arg=getPntrToArgument(0);
  std::vector<unsigned> ind( arg->getRank() );
  for (unsigned i=0; i<arg->getNumberOfValues(); ++i) {
    // Get point indices
    gridcoords.getIndices(i, ind);
    // Fill input data in row-major order
    input_data[ind[0]*N_input_data[0]+ind[1]].real( arg->get( i ) );
    input_data[ind[0]*N_input_data[0]+ind[1]].imag( 0.0 );
  }

  // *** HERE is the only clear limitation: I'm computing explicitly a 2D FT. It should not happen to deal with other than two-dimensional grid ...
  fftw_plan plan_complex = fftw_plan_dft_2d(N_input_data[0], N_input_data[1], reinterpret_cast<fftw_complex*>(&input_data[0]), reinterpret_cast<fftw_complex*>(&fft_data[0]), fourier_params[1], FFTW_ESTIMATE);

  // Compute FT
  fftw_execute( plan_complex );

  // Compute the normalization constant
  double norm=1.0;
  for (unsigned i=0; i<N_input_data.size(); ++i) {
    norm *= pow( N_input_data[i], (1-fourier_params[0])/2 );
  }

  // Save FT data to output grid
  std::vector<std::size_t> N_out_data ( getGridCoordinatesObject().getNbin(true) );
  std::vector<unsigned> out_ind ( getPntrToArgument(0)->getRank() );
  for(unsigned i=0; i<getPntrToArgument(0)->getNumberOfValues(); ++i) {
    gridcoords.getIndices( i, out_ind );
    if (real_output) {
      double ft_value;
      // Compute abs/norm and fix normalization
      if (!store_norm) {
        ft_value=std::abs( fft_data[out_ind[0]*N_out_data[0]+out_ind[1]] / norm );
      } else {
        ft_value=std::norm( fft_data[out_ind[0]*N_out_data[0]+out_ind[1]] / norm );
      }
      // Set the value
      getPntrToComponent(0)->set( i, ft_value);
    } else {
      double ft_value_real, ft_value_imag;
      ft_value_real=fft_data[out_ind[0]*N_out_data[0]+out_ind[1]].real() / norm;
      ft_value_imag=fft_data[out_ind[0]*N_out_data[0]+out_ind[1]].imag() / norm;
      // Set values
      getPntrToComponent(0)->set( i, ft_value_real );
      getPntrToComponent(1)->set( i, ft_value_imag );
    }
  }

  // Free FFTW stuff
  fftw_destroy_plan(plan_complex);
#endif
}

} // end namespace 'gridtools'
} // end namespace 'PLMD'
