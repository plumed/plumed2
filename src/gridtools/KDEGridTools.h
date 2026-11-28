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
#include "tools/GridExtension.h"

namespace PLMD {
namespace gridtools {

/// The component (0, 1 or 2) of the box for an argument x, y or z of a vector (e.g. DISTANCE ... COMPONENTS), -1 otherwise:
/// with no GRID_MIN and GRID_MAX the grid along these arguments covers the box
inline int getBoxComponent( const Value* arg ) {
  const std::string& name=arg->getName();
  const std::size_t dot=name.find_first_of(".");
  if( dot==std::string::npos ) {
    return -1;
  }
  const std::string c=name.substr(dot+1);
  if( c=="x" ) {
    return 0;
  } else if( c=="y" ) {
    return 1;
  } else if( c=="z" ) {
    return 2;
  }
  return -1;
}

/// The grid of KDE. As in METAD, the grid is extended along the arguments that are not periodic whenever the kernels
/// do not fit in it, so that it contains all the points where the density is not zero (plus one point on each side).
/// The spacing does not change and the new points are on the same lattice (GridExtension plans the extension as for METAD,
/// GridCoordinatesObject::extend applies it). With GRID_SPACING and no GRID_MIN and GRID_MAX the grid is placed around the
/// kernels of the first step that has some, on the lattice with a point at zero.
template <class K, class P>
class KDEGridTools {
public:
/// GRID_SPACING and GRID_BIN as given in input (empty if not given)
  std::vector<double> gspacing;
  std::vector<std::size_t> nbin;
/// GRID_MIN and GRID_MAX (auto: from the box or from the kernels)
  std::vector<std::string> gmin, gmax;
/// the range along this argument is set from the kernels (no GRID_MIN and GRID_MAX, not periodic, not a box component)
  std::vector<bool> automatic;
/// the grid is extended along this argument (not periodic, not a box component, no NO_EXTEND)
  std::vector<bool> extendable;
/// maximum number of points of the grid (GRID_SIZE_LIMIT, as given in input)
  double sizelimit{1e8};
/// the extensions of the grid (the limit can be raised by a large initial grid)
  GridExtension extension;
/// the range from the kernels has not been set yet (no kernels with a non-zero height so far)
  bool pristine{false};
/// the grid was restored from a checkpoint
  bool restored{false};
  static void registerKeywords( Keywords& keys );
  static void readHeightKeyword( bool canusevol, std::size_t nargs, const std::vector<std::string>& bw, ActionWithArguments* action );
  static void readBandwidth( std::size_t nargs, ActionWithArguments* action, std::vector<std::string>& bw );
  static void readBandwidthKeyword( std::size_t nargs, ActionWithArguments* action, std::vector<std::string>& bw, std::vector<Value*>& bwargs );
  static void readBandwidthAndHeight( const P& params, ActionWithArguments* action );
  static void convertHeightsToVolumes( const std::size_t& nargs, const std::vector<std::string>& bw, const std::string& volstr, ActionWithArguments* action );
  static void readGridParameters( KDEGridTools<K,P>& g, ActionWithArguments* action, GridCoordinatesObject& gridobject, std::vector<std::size_t>& shape );
/// the grid at the first step: periodic domains, the box, GRID_MIN and GRID_MAX, or a first range of one bin
  static void setupGridBounds( KDEGridTools<K,P>& g, P& p, const Tensor& box, GridCoordinatesObject& gridobject, const std::vector<Value*>& args, Value* myval, Action& action );
/// place the grid around the kernels (pristine grid) or extend it if the kernels do not fit in it; true if the grid changed
  static bool extendGridBounds( KDEGridTools<K,P>& g, P& p, GridCoordinatesObject& gridobject, const std::vector<Value*>& args, Value* myval, Action& action );
/// the range of the kernels with a non-zero height along each argument (centers plus and minus the support); false if there are none
  static bool getKernelsRange( P& p, const std::vector<Value*>& args, std::size_t ndim, std::vector<double>& lower, std::vector<double>& upper );
  static bool hasCheckpointState( const KDEGridTools<K,P>& g );
  static void writeCheckpoint( const KDEGridTools<K,P>& g, const GridCoordinatesObject& gridobject, CheckpointWriter& w );
  static void readCheckpoint( KDEGridTools<K,P>& g, GridCoordinatesObject& gridobject, Value* myval, CheckpointReader& r );
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
  keys.add("compulsory","GRID_MIN","auto","the lower bounds for the grid (auto: from the kernels, from the box for the components x, y and z of a vector, the domain for periodic arguments); the grid is extended when the kernels do not fit in it");
  keys.add("compulsory","GRID_MAX","auto","the upper bounds for the grid (auto: from the kernels, from the box for the components x, y and z of a vector, the domain for periodic arguments); the grid is extended when the kernels do not fit in it");
  keys.add("optional","GRID_SPACING","the grid spacing (with GRID_MIN and GRID_MAX, the grid spacing is not larger than this)");
  keys.add("optional","GRID_BIN","the number of bins for the grid (with no GRID_MIN and GRID_MAX, the range of the kernels of the first step is divided in this number of bins)");
  keys.add("compulsory","GRID_SIZE_LIMIT","100000000","the maximum number of points of the grid, which is extended when the kernels do not fit in it");
  keys.addFlag("NO_EXTEND",false,"do not extend the grid: the density of the kernels outside GRID_MIN and GRID_MAX is not computed (as in PLUMED 2.10)");
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
      if( Tools::convertNoexcept( bw[i], bwval ) && bwval<0 ) {
        action->error("BANDWIDTH should not be negative");
      }
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
  const std::size_t ndim=shape.size();
  g.gmin.resize( ndim );
  g.gmax.resize( ndim );
  action->parseVector("GRID_MIN",g.gmin);
  action->parseVector("GRID_MAX",g.gmax);
  action->parseVector("GRID_BIN",g.nbin);
  action->parseVector("GRID_SPACING",g.gspacing);
  action->parse("GRID_SIZE_LIMIT",g.sizelimit);
  bool noextend=false;
  action->parseFlag("NO_EXTEND",noextend);
  if( g.nbin.size()!=ndim && g.gspacing.size()!=ndim ) {
    action->error("GRID_BIN or GRID_SPACING must be set");
  }
  for(const auto s : g.gspacing) {
    if( !(s>0.0) ) {
      action->error("GRID_SPACING should be positive");
    }
  }
  const std::string limiterror=GridExtension::checkSizeLimit( g.sizelimit );
  if( !limiterror.empty() ) {
    action->error( limiterror );
  }
  g.extension.setSizeLimit( g.sizelimit, {} );
  g.automatic.assign( ndim, false );
  g.extendable.assign( ndim, false );
  std::vector<bool> ipbc( ndim );
  for(unsigned i=0; i<ndim; ++i) {
    Value* myarg = action->getPntrToArgument(i);
    const bool isauto=(g.gmin[i]=="auto");
    if( isauto!=(g.gmax[i]=="auto") ) {
      action->error("GRID_MIN and GRID_MAX should be either both given or both auto");
    }
    if( myarg->isPeriodic() ) {
      ipbc[i]=true;
      if( isauto ) {
        myarg->getDomain( g.gmin[i], g.gmax[i] );
      }
      action->log.printf("  for %dth coordinate (periodic) min is set to %s and max is set to %s \n", (i+1), g.gmin[i].c_str(), g.gmax[i].c_str() );
    } else if( isauto && getBoxComponent(myarg)>=0 ) {
      ipbc[i]=true;
      action->log.printf("  for %dth coordinate min and max are set from the box at the first step \n", (i+1) );
    } else {
      ipbc[i]=false;
      g.automatic[i]=isauto;
      g.extendable[i]=!noextend;
      if( isauto ) {
        g.pristine=true;
        action->log.printf("  for %dth coordinate min and max are set from the kernels \n", (i+1) );
      } else {
        action->log.printf("  for %dth coordinate min is set to %s and max is set to %s \n", (i+1), g.gmin[i].c_str(), g.gmax[i].c_str() );
      }
    }
    if( g.nbin.size()==ndim ) {
      shape[i] = ipbc[i] ? g.nbin[i] : g.nbin[i]+1;
    }
  }
  if( g.pristine && g.gspacing.size()!=ndim ) {
    action->log.printf("  with no GRID_SPACING, the range of the kernels of the first step is divided in GRID_BIN bins\n");
  }
  if( std::find(g.extendable.begin(),g.extendable.end(),true)!=g.extendable.end() ) {
    action->log.printf("  the grid is extended when the kernels do not fit in it, up to %s points (GRID_SIZE_LIMIT)\n", Tools::boundString(g.sizelimit).c_str() );
  } else if( noextend ) {
    action->log.printf("  the grid is not extended (NO_EXTEND): the density outside it is not computed\n");
  }
  gridobject.setup( "flat", ipbc, 0, 0.0 );
}

template <class K, class P>
void KDEGridTools<K,P>::setupGridBounds( KDEGridTools<K,P>& g, P& p, const Tensor& box, GridCoordinatesObject& gridobject, const std::vector<Value*>& args, Value* myval, Action& action ) {
  if( g.restored ) {
    return;
  }
  const std::size_t ndim = gridobject.getDimension();
  std::vector<std::string> smin( g.gmin ), smax( g.gmax );
  for(unsigned i=0; i<ndim; ++i) {
    if( g.automatic[i] ) {
      // a first range of one bin (or GRID_BIN bins): the grid is placed around the kernels by extendGridBounds
      smin[i]="0";
      if( g.gspacing.size()==ndim ) {
        smax[i]=Tools::boundString( (g.nbin.size()==ndim ? g.nbin[i] : 1)*g.gspacing[i] );
      } else {
        smax[i]="1";
      }
    } else if( smin[i]=="auto" ) {
      const int c=getBoxComponent( args[i] );
      plumed_assert( c>=0 );
      smin[i]=Tools::boundString( -0.5*box(c,c) );
      smax[i]=Tools::boundString( 0.5*box(c,c) );
    }
  }
  std::vector<double> spacing( g.gspacing );
  gridobject.setBounds( smin, smax, g.nbin, spacing );
  myval->setShape( gridobject.getNbin(true) );
  // as in METAD, an initial grid that is larger than GRID_SIZE_LIMIT raises the limit to twice its size
  if( std::find(g.extendable.begin(),g.extendable.end(),true)!=g.extendable.end() && g.extension.setSizeLimit( g.sizelimit, gridobject.getNbin(true) ) ) {
    action.log.printf("  The initial grid has more points than GRID_SIZE_LIMIT: the limit is set to %s points\n", Tools::boundString(g.extension.getSizeLimit()).c_str());
  }
}

template <class K, class P>
bool KDEGridTools<K,P>::getKernelsRange( P& p, const std::vector<Value*>& args, std::size_t ndim, std::vector<double>& lower, std::vector<double>& upper ) {
  lower.assign( ndim, std::numeric_limits<double>::max() );
  upper.assign( ndim, -std::numeric_limits<double>::max() );
  std::vector<double> argval( args.size() ), support( ndim, 0.0 ), kparams( K::getNumberOfParameters(ndim), 0.0 );
  View<double> kview( kparams.data(), kparams.size() );
  // with the same bandwidth for all the kernels the support is computed once
  const bool samewidth = K::bandwidthsAllSame( ndim, args );
  bool havesupport=false, found=false;
  auto addKernel = [&]() {
    if( !K::setKernelAndCheckHeight( kview, ndim, argval ) ) {
      return;
    }
    if( !havesupport ) {
      std::fill( support.begin(), support.end(), 0.0 );
      P::getSupport( p, View<const double>( kparams.data()+ndim, K::getNumberOfShapeParameters(ndim) ), support );
      havesupport=samewidth;
    }
    for(unsigned i=0; i<ndim; ++i) {
      lower[i]=std::min( lower[i], kparams[i]-support[i] );
      upper[i]=std::max( upper[i], kparams[i]+support[i] );
    }
    found=true;
  };
  // the kernels are the stored values of the heights (the last argument), as in KDEHelper::transferKernels
  const Value* heights=args[args.size()-1];
  if( heights->getRank()==2 ) {
    const unsigned nc    = heights->getShape()[1];
    const unsigned nrows = heights->getShape()[0];
    const unsigned ncs   = heights->getNumberOfColumns();
    for(unsigned i=0; i<nrows; ++i) {
      const unsigned ncols = heights->getRowLength(i);
      for(unsigned j=0; j<ncols; ++j) {
        const unsigned jind = heights->getRowIndex( i, j );
        for(unsigned k=0; k<args.size(); ++k) {
          if( jind==args[k]->getRowIndex( i, j ) ) {
            argval[k] = args[k]->get( i*ncs+ j, false );
          } else {
            argval[k] = args[k]->get( i*nc + jind );
          }
        }
        addKernel();
      }
    }
  } else {
    const std::size_t nkernels = heights->getNumberOfStoredValues();
    for(unsigned i=0; i<nkernels; ++i) {
      for(unsigned j=0; j<args.size(); ++j) {
        argval[j] = args[j]->get(i,false);
      }
      addKernel();
    }
  }
  return found;
}

template <class K, class P>
bool KDEGridTools<K,P>::extendGridBounds( KDEGridTools<K,P>& g, P& p, GridCoordinatesObject& gridobject, const std::vector<Value*>& args, Value* myval, Action& action ) {
  if( std::find(g.extendable.begin(),g.extendable.end(),true)==g.extendable.end() && !g.pristine ) {
    return false;
  }
  const std::size_t ndim = gridobject.getDimension();
  std::vector<double> lower, upper;
  if( !getKernelsRange( p, args, ndim, lower, upper ) ) {
    // no kernels with a non-zero height: nothing to fit in the grid
    return false;
  }
  for(unsigned i=0; i<ndim; ++i) {
    if( ((g.pristine && g.automatic[i]) || g.extendable[i]) && (!std::isfinite(lower[i]) || !std::isfinite(upper[i])) ) {
      action.error("cannot extend the grid along argument "+args[i]->getName()+" to contain the kernels: non finite value");
    }
  }
  const bool placed=g.pristine;
  g.pristine=false;
  if( placed ) {
    // the grid along the arguments with no GRID_MIN and GRID_MAX is placed around the kernels (a new lattice)
    std::vector<double> gmin( gridobject.getMinValues() ), dx( gridobject.getGridSpacing() );
    std::vector<std::size_t> nb( gridobject.getNbin(false) );
    std::vector<std::string> smin( gridobject.getMin() ), smax( gridobject.getMax() );
    double size=1.0;
    for(unsigned i=0; i<ndim; ++i) {
      if( g.automatic[i] ) {
        if( g.gspacing.size()==ndim ) {
          // on the lattice with a point at zero, plus one point on each side (as in METAD)
          dx[i]=g.gspacing[i];
          const double first=std::floor(lower[i]/dx[i])-1.0, last=std::ceil(upper[i]/dx[i])+1.0;
          if( std::fabs(first)>1e15 || std::fabs(last)>1e15 ) {
            action.error("cannot place the grid along argument "+args[i]->getName()+" around the kernels");
          }
          gmin[i]=first*dx[i];
          nb[i]=static_cast<std::size_t>( last-first );
        } else {
          // GRID_BIN only: the range of the kernels of the first step is divided in GRID_BIN bins, plus one on each side
          if( !(upper[i]>lower[i]) ) {
            action.error("the kernels of the first step have no extent along argument "+args[i]->getName()+": use GRID_SPACING to set the grid");
          }
          dx[i]=(upper[i]-lower[i])/g.nbin[i];
          gmin[i]=lower[i]-dx[i];
          nb[i]=g.nbin[i]+2;
        }
        smin[i]=Tools::boundString( gmin[i] );
        smax[i]=Tools::boundString( gmin[i]+nb[i]*dx[i] );
      }
      size*=static_cast<double>(nb[i]) + (gridobject.isPeriodic(i) ? 0.0 : 1.0);
    }
    if( size>g.extension.getSizeLimit() ) {
      action.error("the kernels require a grid with "+Tools::boundString(size)+" points, more than GRID_SIZE_LIMIT="
                   +Tools::boundString(g.extension.getSizeLimit())+": increase GRID_SIZE_LIMIT, or check the arguments and the bandwidth");
    }
    gridobject.setLattice( gmin, dx, nb, smin, smax );
  }
  // along the other arguments that can be extended: all the points where the density is not zero, plus one on each
  // side, with the growth of METAD (see GridExtension)
  std::vector<GridExtension::Mode> mode( ndim, GridExtension::Mode::fixed );
  std::vector<double> firstpoint( ndim, 0.0 ), lastpoint( ndim, 0.0 );
  for(unsigned i=0; i<ndim; ++i) {
    if( !g.extendable[i] || (placed && g.automatic[i]) ) {
      continue;
    }
    mode[i]=GridExtension::Mode::extend;
    firstpoint[i]=std::floor(gridobject.getLatticeCoordinate(i,lower[i]))-1.0;
    lastpoint[i]=std::ceil(gridobject.getLatticeCoordinate(i,upper[i]))+1.0;
  }
  const GridExtension::Plan plan=g.extension.plan( gridobject.getNbin(true), mode, firstpoint, lastpoint );
  if( plan.nonfinite!=std::string::npos ) {
    action.error("cannot extend the grid along argument "+args[plan.nonfinite]->getName()+" to contain the kernels");
  }
  if( plan.toolarge ) {
    action.error("the kernels require a grid with "+Tools::boundString(plan.npoints)+" points, more than GRID_SIZE_LIMIT="
                 +Tools::boundString(g.extension.getSizeLimit())+": increase GRID_SIZE_LIMIT, or check the arguments (e.g. a distance that grows without limit) "
                 "and the bandwidth, or use NO_EXTEND");
  }
  if( !placed && !plan.changed ) {
    return false;
  }
  if( plan.changed ) {
    gridobject.extend( plan.lower, plan.upper );
  }
  myval->setShape( gridobject.getNbin(true) );
  action.log.printf("  %s grid: min", placed ? "Placing the" : "Extending the");
  for(unsigned i=0; i<ndim; ++i) {
    action.log.printf(" %s",gridobject.getMin()[i].c_str());
  }
  action.log.printf(" max");
  for(unsigned i=0; i<ndim; ++i) {
    action.log.printf(" %s",gridobject.getMax()[i].c_str());
  }
  action.log.printf(" points");
  for(const auto n : gridobject.getNbin(true)) {
    action.log.printf(" %zu",n);
  }
  action.log.printf(" at time %f\n",action.getTime());
  return true;
}

template <class K, class P>
bool KDEGridTools<K,P>::hasCheckpointState( const KDEGridTools<K,P>& g ) {
  // the grid depends on the kernels or on the box of the first step
  for(unsigned i=0; i<g.gmin.size(); ++i) {
    if( g.extendable[i] || g.automatic[i] || g.gmin[i]=="auto" ) {
      return true;
    }
  }
  return false;
}

template <class K, class P>
void KDEGridTools<K,P>::writeCheckpoint( const KDEGridTools<K,P>& g, const GridCoordinatesObject& gridobject, CheckpointWriter& w ) {
  // nothing to store before the first step (the grid has not been set up)
  const bool set=gridobject.getNumberOfPoints()>0 && gridobject.getMin().size()==gridobject.getDimension();
  w.write("grid_set",(long long) set);
  if( !set ) {
    return;
  }
  w.write("pristine",(long long) g.pristine);
  gridobject.writeCheckpoint( w, "grid_" );
}

template <class K, class P>
void KDEGridTools<K,P>::readCheckpoint( KDEGridTools<K,P>& g, GridCoordinatesObject& gridobject, Value* myval, CheckpointReader& r ) {
  long long set;
  r.read("grid_set",set);
  if( !set ) {
    return;
  }
  long long pristine;
  r.read("pristine",pristine);
  g.pristine=(pristine!=0);
  gridobject.readCheckpoint( r, "grid_" );
  myval->setShape( gridobject.getNbin(true) );
  g.restored=true;
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
  P::getSupport( p, shape, support );
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
