/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2023 The plumed team
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
#ifndef __PLUMED_gridtools_KDEGridTools_h
#define __PLUMED_gridtools_KDEGridTools_h

#include "KDE.h"

namespace PLMD {
namespace gridtools {

template <class K, class P>
class KDEGridTools {
public:
  double dp2cutoff;
  std::vector<double> gspacing;
  std::vector<std::size_t> nbin;
  std::vector<std::string> gmin, gmax;
  static void registerKeywords( Keywords& keys );
  static void readHeightKeyword( bool canusevol, std::size_t nargs, const std::vector<std::string>& bw, ActionWithArguments* action );
  static void readBandwidth( std::size_t nargs, ActionWithArguments* action, std::vector<std::string>& bw );
  static void readBandwidthKeyword( std::size_t nargs, ActionWithArguments* action, std::vector<std::string>& bw, std::vector<Value*>& bwargs );
  static void readBandwidthAndHeight( const P& params, ActionWithArguments* action );
  static void convertHeightsToVolumes( const std::size_t& nargs, const std::vector<std::string>& bw, const std::string& volstr, ActionWithArguments* action );
  static void readGridParameters( KDEGridTools<K,P>& g, ActionWithArguments* action, GridCoordinatesObject& gridobject, std::vector<std::size_t>& shape );
  static void setupGridBounds( KDEGridTools<K,P>& g, const Tensor& box, GridCoordinatesObject& gridobject, const std::vector<Value*>& args, Value* myval );
  static void getDiscreteSupport( const KDEGridTools<K,P>& g, P& p, const View<const double>& shape, std::vector<unsigned>& nneigh, GridCoordinatesObject& gridobject );
  static void getNeighbors( const P& p, View<double> at, const GridCoordinatesObject& gridobject, const std::vector<unsigned>& nneigh, unsigned& num_neighbors, std::vector<unsigned>& neighbors );
};

template <class K, class P>
void KDEGridTools<K,P>::registerKeywords( Keywords& keys ) {
  keys.add("optional","BANDWIDTH","the bandwidths for kernel density esimtation");
  keys.add("optional","VOLUMES","this keyword take the label of an action that calculates a vector of values.  The elements of this vector "
           "divided by the volume of the Gaussian are used as weights for the Gaussians");
  keys.add("optional","HEIGHTS","this keyword takes the label of an action that calculates a vector of values. The elements of this vector "
           "are used as weights for the Gaussians.");
  keys.add("compulsory","GRID_MIN","auto","the lower bounds for the grid");
  keys.add("compulsory","GRID_MAX","auto","the upper bounds for the grid");
  keys.add("compulsory","CUTOFF","6.25","the cutoff at which to stop evaluating the kernel functions is set equal to sqrt(2*x)*bandwidth in each direction where x is this number");
  keys.add("optional","GRID_SPACING","the approximate grid spacing (to be used as an alternative or together with GRID_BIN)");
  keys.add("optional","GRID_BIN","the number of bins for the grid");
}

template <class K, class P>
void KDEGridTools<K,P>::readHeightKeyword( bool canusevol, std::size_t nargs, const std::vector<std::string>& bw, ActionWithArguments* action ) {
  std::string weight_str;
  action->parse("HEIGHTS",weight_str);
  std::string str_nvals;
  if( (action->getPntrToArgument(0))->getRank()==2 ) {
    std::string nr, nc;
    Tools::convert( (action->getPntrToArgument(0))->getShape()[0], nr );
    Tools::convert( (action->getPntrToArgument(0))->getShape()[1], nc );
    str_nvals = nr + "," + nc;
  } else {
    Tools::convert( (action->getPntrToArgument(0))->getNumberOfValues(), str_nvals );
  }
  if( weight_str.length()>0 ) {
    KDEHelper<K,P,KDEGridTools<K,P>>::readKernelParameters( weight_str, action, "_heights", true );
  } else if( canusevol ) {
    action->plumed.readInputWords( Tools::getWords(action->getLabel() + "_volumes: ONES SIZE=" + str_nvals ), false );
    KDEGridTools<K,P>::convertHeightsToVolumes(nargs,bw,action->getLabel() + "_volumes",action);
  } else {
    action->plumed.readInputWords( Tools::getWords(action->getLabel() + "_heights: ONES SIZE=" + str_nvals ), false );
    KDEHelper<K,P,KDEGridTools<K,P>>::addArgument( action->getLabel() + "_heights", action );
  }
}

template<class K, class P>
void KDEGridTools<K, P>::readBandwidthKeyword( std::size_t nargs, ActionWithArguments* action, std::vector<std::string>& bw, std::vector<Value*>& bwargs ) {
  action->parseVector("BANDWIDTH",bw);
  if( nargs>1 && bw.size()==1 ) {
    ActionWithArguments::interpretArgumentList( bw, action->plumed.getActionSet(), action, bwargs );
    if( bwargs.size()!=1 ) {
      action->error("invalid bandwidth found");
    }
    // Create a vector of ones with the right size
    std::string nvals;
    if( (action->getPntrToArgument(0))->getRank()==2 ) {
      std::string nr, nc;
      Tools::convert( (action->getPntrToArgument(0))->getShape()[0], nr );
      Tools::convert( (action->getPntrToArgument(0))->getShape()[1], nc );
      nvals = nr + "," + nc;
    } else {
      Tools::convert( (action->getPntrToArgument(0))->getNumberOfValues(), nvals );
    }
    action->plumed.readInputWords( Tools::getWords(action->getLabel() + "_bwones: ONES SIZE=" + nvals ), false );
  }
}

template <class K, class P>
void KDEGridTools<K,P>::readBandwidth( std::size_t nargs, ActionWithArguments* action, std::vector<std::string>& bw ) {
  plumed_assert( typeid(K)==typeid(DiagonalKernelParams) );
  std::vector<Value*> bwargs;
  readBandwidthKeyword( nargs, action, bw, bwargs );
  if( nargs>1 && bw.size()==1 ) {
    if( bwargs[0]->getRank()!=1 || bwargs[0]->getNumberOfValues()!=nargs ) {
      action->error("invalid input for bandwidth parameter");
    }
    std::string str_i;
    bw.resize( nargs );
    for(unsigned i=0; i<nargs; ++i) {
      Tools::convert( i+1, str_i );
      action->plumed.readInputWords( Tools::getWords(action->getLabel() + "_scalar_bw" + str_i + ": SELECT_COMPONENTS ARG=" + bwargs[0]->getName() + " COMPONENTS=" + str_i ), false );
      bw[i] = action->getLabel() + "_bw" + str_i;
      action->plumed.readInputWords( Tools::getWords( bw[i] + ": CUSTOM ARG=" + action->getLabel() + "_bwones," + action->getLabel() + "_scalar_bw" + str_i + " FUNC=x*y PERIODIC=NO"), false );
    }
  } else if( bw.size()==nargs ) {
    double bwval;
    std::string str_i;
    for(unsigned i=0; i<nargs; ++i) {
      Tools::convert( i+1, str_i );
      if( Tools::convertNoexcept( bw[i], bwval ) && fabs(bwval)<epsilon ) {
        KDEHelper<K,P,KDEGridTools<K,P>>::readKernelParameters( bw[i], action, "_bwz" + str_i, true );
      } else {
        KDEHelper<K,P,KDEGridTools<K,P>>::readKernelParameters( bw[i], action, "_bw" + str_i, true );
      }
    }
  } else {
    action->error("wrong number of arguments specified in input to bandwidth parameter");
  }
}

template <class K, class P>
void KDEGridTools<K, P>::convertHeightsToVolumes( const std::size_t& nargs, const std::vector<std::string>& bw, const std::string& volstr, ActionWithArguments* action ) {
  if( bw.size()==nargs ) {
    unsigned nonzeroargs=0;
    for(unsigned i=0; i<nargs; ++i) {
      if( bw[i].find("_bwz")==std::string::npos ) {
        nonzeroargs++;
      }
    }
    std::string str_i, nargs_str;
    Tools::convert( nonzeroargs, nargs_str );
    std::string varstr = "VAR=h", funcstr = "FUNC=h/(sqrt((2*pi)^" + nargs_str + ")", argstr = "ARG=" + volstr;
    for(unsigned i=0; i<nargs; ++i) {
      if( bw[i].find("_bwz")!=std::string::npos ) {
        continue;
      }
      Tools::convert( i+1, str_i );
      varstr += ",b" + str_i;
      funcstr += "*b" + str_i;
      argstr += "," + bw[i];
    }
    funcstr += ")";
    if( (action->getPntrToArgument(0))->getNumberOfValues()==1 ) {
      action->plumed.readInputWords( Tools::getWords(action->getLabel() + "_heights: CUSTOM PERIODIC=NO " + argstr + " " + varstr + " " + funcstr), false );
    } else {
      action->plumed.readInputWords( Tools::getWords(action->getLabel() + "_heights: CUSTOM PERIODIC=NO " + argstr + " " + varstr + " " + funcstr + " MASK=" + volstr), false );
    }
    KDEHelper<K,P,KDEGridTools<K,P>>::addArgument( action->getLabel() + "_heights", action );
  } else if( bw.size()==nargs*nargs ) {
    action->error("have not implemented normalization parameter for non-diagonal kernels");
  }

}

template <class K, class P>
void KDEGridTools<K,P>::readGridParameters( KDEGridTools<K,P>& g, ActionWithArguments* action, GridCoordinatesObject& gridobject, std::vector<std::size_t>& shape ) {
  g.gmin.resize( shape.size() );
  g.gmax.resize( shape.size() );
  action->parseVector("GRID_MIN",g.gmin);
  action->parseVector("GRID_MAX",g.gmax);
  for(unsigned i=0; i<g.gmin.size(); ++i) {
    if( g.gmin[i]=="auto" ) {
      action->log.printf("  for %dth coordinate min and max are set automatically \n", (i+1) );
      if( g.gmax[i]!="auto" ) {
        action->error("if gmin is set automatically gmax must also be set automatically");
      }
      Value* myarg = action->getPntrToArgument(i);
      if( myarg->isPeriodic() ) {
        if( g.gmin[i]=="auto" ) {
          myarg->getDomain( g.gmin[i], g.gmax[i] );
        } else {
          std::string str_min, str_max;
          myarg->getDomain( str_min, str_max );
          if( str_min!=g.gmin[i] || str_max!=g.gmax[i] ) {
            action->error("all periodic arguments should have the same domain");
          }
        }
      } else if( myarg->getName().find(".")!=std::string::npos ) {
        std::size_t dot = myarg->getName().find_first_of(".");
        std::string name = myarg->getName().substr(dot+1);
        if( name!="x" && name!="y" && name!="z" ) {
          action->error("cannot set GRID_MIN and GRID_MAX automatically if input argument is not component of distance");
        }
      } else {
        action->error("cannot set GRID_MIN and GRID_MAX automatically if input argument is not component of distance");
      }
    } else {
      action->log.printf("  for %dth coordinate min is set to %s and max is set to %s \n", (i+1), g.gmin[i].c_str(), g.gmax[i].c_str() );
    }
  }

  action->parseVector("GRID_BIN",g.nbin);
  action->parseVector("GRID_SPACING",g.gspacing);
  action->parse("CUTOFF",g.dp2cutoff);

  if( g.nbin.size()!=shape.size() && g.gspacing.size()!=shape.size() ) {
    action->error("GRID_BIN or GRID_SPACING must be set");
  }
  // Create a value
  std::vector<bool> ipbc( shape.size() );
  for(unsigned i=0; i<shape.size(); ++i) {
    if( (action->getPntrToArgument( i ))->isPeriodic() || g.gmin[i]=="auto" ) {
      ipbc[i]=true;
      if( g.nbin.size()==shape.size() ) {
        shape[i] = g.nbin[i];
      }
    } else {
      ipbc[i]=false;
      if( g.nbin.size()==shape.size() ) {
        shape[i] = g.nbin[i]+1;
      }
    }
  }
  gridobject.setup( "flat", ipbc, 0, 0.0 );
}

template <class K, class P>
void KDEGridTools<K,P>::setupGridBounds( KDEGridTools<K,P>& g, const Tensor& box, GridCoordinatesObject& gridobject, const std::vector<Value*>& args, Value* myval ) {
  for(unsigned i=0; i<gridobject.getDimension(); ++i) {
    if( g.gmin[i]=="auto" ) {
      double lcoord, ucoord;
      std::size_t dot = args[i]->getName().find_first_of(".");
      std::string name = args[i]->getName().substr(dot+1);
      if( name=="x" ) {
        lcoord=-0.5*box(0,0);
        ucoord=0.5*box(0,0);
      } else if( name=="y" ) {
        lcoord=-0.5*box(1,1);
        ucoord=0.5*box(1,1);
      } else if( name=="z" ) {
        lcoord=-0.5*box(2,2);
        ucoord=0.5*box(2,2);
      } else {
        plumed_error();
      }
      // And convert to strings for bin and bmax
      Tools::convert( lcoord, g.gmin[i] );
      Tools::convert( ucoord, g.gmax[i] );
    }
  }
  // And setup the grid object
  gridobject.setBounds( g.gmin, g.gmax, g.nbin, g.gspacing );
  myval->setShape( gridobject.getNbin(true) );
}

template <> inline
void KDEGridTools<DiagonalKernelParams,DiscreteKernel>::getDiscreteSupport( const KDEGridTools<DiagonalKernelParams,DiscreteKernel>& g, DiscreteKernel& p, const View<const double>& shape, std::vector<unsigned>& nneigh, GridCoordinatesObject& gridobject ) {
  return;
}

template <class K, class P>
void KDEGridTools<K, P>::getDiscreteSupport( const KDEGridTools<K,P>& g, P& p, const View<const double>& shape, std::vector<unsigned>& nneigh, GridCoordinatesObject& gridobject ) {
  std::size_t ng = gridobject.getDimension();
  plumed_assert( nneigh.size()==ng );
  std::vector<double> support( ng );
  P::getSupport( p, shape, g.dp2cutoff, support );
  for(unsigned i=0; i<ng; ++i) {
    nneigh[i] = static_cast<unsigned>( ceil( support[i]/gridobject.getGridSpacing()[i] ));
  }
}

template <> inline
void KDEGridTools<DiagonalKernelParams,DiscreteKernel>::getNeighbors( const DiscreteKernel& p, View<double> at, const GridCoordinatesObject& gridobject, const std::vector<unsigned>& nneigh, unsigned& num_neighbors, std::vector<unsigned>& neighbors ) {
  num_neighbors=1;
  neighbors.resize(1);
  for(unsigned i=0; i<at.size(); ++i) {
    at[i] += 0.5*gridobject.getGridSpacing()[i];
  }
  neighbors[0]=gridobject.getIndex( View<const double>(at.data(),at.size()) );
}

template <class K, class P>
void KDEGridTools<K,P>::getNeighbors( const P& p, View<double> at, const GridCoordinatesObject& gridobject, const std::vector<unsigned>& nneigh, unsigned& num_neighbors, std::vector<unsigned>& neighbors ) {
  gridobject.getNeighbors( View<const double>(at.data(),at.size()), nneigh, num_neighbors, neighbors );
}

}
}
#endif
