/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
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
#include "LandmarkSelectionBase.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionSetup.h"
#include "core/CollectFrames.h"

namespace PLMD {
namespace analysis {

void LandmarkSelectionBase::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys ); ActionWithValue::registerKeywords( keys ); ActionWithArguments::registerKeywords( keys );
  keys.add("optional","DATA","the label of a COLLECT_FRAMES action from which landmarks should be selected");
  keys.use("ARG"); keys.add("compulsory","NLANDMARKS","the number of landmarks that you would like to select");
  keys.addOutputComponent("_rect","default","the rectangular matrix containing the distances between the landmarks and all the points");
  keys.addOutputComponent("_sqr","default","the square matrix containing the distances between the various landmark points");
  ActionWithValue::useCustomisableComponents( keys );
}

LandmarkSelectionBase::LandmarkSelectionBase( const ActionOptions& ao ):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao),
  jframe(0),
  nlandmarks(0),
  nvals(0)
{
  std::string data; parse("DATA",data);
  if( data.length()>0 ) {
      CollectFrames* myfram = plumed.getActionSet().selectWithLabel<CollectFrames*>(data);
      if( !myfram ) error( data + " is not a valid COLLECT_FRAMES action so cannot collect data");
      std::vector<Value*> vals( getArguments() );
      for(unsigned j=0;j<myfram->getNumberOfComponents();++j) { 
          vals.push_back( myfram->copyOutput(j) ); arg_ends.push_back(vals.size()); 
      }
      requestArguments( vals, false );
  }
 
  if( keywords.exists("NLANDMARKS") ) parse("NLANDMARKS",nlandmarks);
  log.printf("  selecting %u landmark points \n",nlandmarks);

  nvals = 0; landmarks.resize(nlandmarks);
  for(unsigned j=arg_ends[0];j<arg_ends[1];++j) {
      if( getPntrToArgument(j)->getRank()==1 ) nvals += getPntrToArgument(j)->getNumberOfValues( getLabel() );
      else if( getPntrToArgument(j)->getRank()==2 ) nvals += getPntrToArgument(j)->getShape()[1];
  }

  std::vector<unsigned> shape(1); std::vector<std::string> foundmat;
  for(unsigned i=0;i<arg_ends.size()-1;++i) { 
      unsigned tvals=0; 
      for(unsigned j=arg_ends[i];j<arg_ends[i+1];++j) {
          if( getPntrToArgument(j)->getRank()==1 ) tvals += getPntrToArgument(j)->getNumberOfValues( getLabel() );
          else if( getPntrToArgument(j)->getRank()==2 ) tvals += getPntrToArgument(j)->getShape()[1];
      }
      if( tvals!=nvals ) {
          if( tvals==0 ) {
              // Check here if input is a matrix that was read in during startup
              ActionSetup* as = dynamic_cast<ActionSetup*>( getPntrToArgument(arg_ends[0])->getPntrToAction() );
              if( !as ) error("mismatch between sizes of input positions");
          } else error("mismatch between sizes of input positions");
      }
      // Add suitable argument
      Value* argi = getPntrToArgument(arg_ends[i]);
      if( argi->hasDerivatives() ) {
          error("cannot select landmarks for value " + argi->getName() );
      } else if( argi->getRank()==2 ) {
          shape.resize(2); shape[0] = nlandmarks; shape[1] = nvals; argi->buildDataStore( getLabel() );
      } else if( argi->getRank()==1 ) {
          shape.resize(1); shape[0] = nlandmarks;
      } else error("cannot select landmarks for value " + argi->getName() );
      std::string vname = argi->getName(); 
      CollectFrames* cf = dynamic_cast<CollectFrames*>( argi->getPntrToAction() );
      if( cf ) { 
          if( shape.size()==2 ) error("input from COLLECT_FRAMES should be a vector -- how have you done this?");
          std::size_t dot=vname.find_first_of('.'); vname = vname.substr(dot+1); 
      } else if( shape.size()==2 ) { foundmat.push_back(vname); vname = vname + "_rect"; }
      addComponent( vname, shape );
      if( argi->isPeriodic() ) {
          std::string min, max; getPntrToArgument(arg_ends[i])->getDomain( min, max );
          componentIsPeriodic( vname, min, max );
      } else componentIsNotPeriodic( vname ); 
      if( argi->isTimeSeries() ) getPntrToOutput( getNumberOfComponents()-1 )->makeTimeSeries();
      getPntrToOutput( getNumberOfComponents()-1 )->alwaysStoreValues();
  }
  if( foundmat.size()>0 ) {
      shape.resize(2); shape[0]=shape[1]=nlandmarks; 
      for(unsigned i=0;i<foundmat.size();++i) {
          addComponent( foundmat[i] + "_sqr", shape ); componentIsNotPeriodic( foundmat[i] + "_sqr" ); 
          getPntrToOutput( getNumberOfComponents()-1 )->alwaysStoreValues(); 
          if( getPntrToOutput(0)->isTimeSeries() ) getPntrToOutput( getNumberOfComponents()-1 )->makeTimeSeries();
      }
  }
}

void LandmarkSelectionBase::getGridPointIndicesAndCoordinates( const unsigned& ind, std::vector<unsigned>& indices, std::vector<double>& coords ) const {
  if( indices[0]>landmarks.size() ) { 
      getPntrToArgument(0)->getPntrToAction()->getGridPointIndicesAndCoordinates( ind, indices, coords );
  } else if( getPntrToOutput(0)->isTimeSeries() ) { 
      ActionWithValue* cf = getPntrToArgument(0)->getPntrToAction(); 
      std::vector<double> first(1), second(1); cf->getGridPointIndicesAndCoordinates( 0, indices, first ); cf->getGridPointIndicesAndCoordinates( 1, indices, second ); 
      coords[0] = first[0] + landmarks[ind]*getTimeStep()*( second[0] - first[0] );
  } else getPntrToArgument(0)->getPntrToAction()->getGridPointIndicesAndCoordinates( landmarks[ind], indices, coords );
}

void LandmarkSelectionBase::selectFrame( const unsigned& iframe ) {
  for(unsigned i=0;i<arg_ends.size()-1;++i) {
      Value* val0=getPntrToOutput(i); 
      if( val0->getRank()==2 && val0->storingData() ) {
          Value* arg0 = getPntrToArgument(arg_ends[i]); unsigned ncols = arg0->getNumberOfColumns();
          for(unsigned j=0;j<arg0->getRowLength(iframe);++j) {
              val0->set( nvals*jframe + arg0->getRowIndex(iframe,j), arg0->get(ncols*iframe+j) );
          }
      } else if( val0->getRank()==1 && val0->storingData() ) val0->set( jframe, retrieveRequiredArgument( i, iframe ) );
  }
  landmarks[jframe]=iframe; jframe++; if( jframe==nlandmarks ) jframe=0;
}

void LandmarkSelectionBase::setLandmarkSeparations() {
  if( getNumberOfComponents()==(arg_ends.size()-1) ) return;
  for(unsigned k=arg_ends.size()-1;k<getNumberOfComponents();++k) {
      Value* val0 = getPntrToOutput( k );
      if( !val0->storingData() ) continue;

      Value* arg0; std::string fname = getPntrToOutput(k)->getName(); fname=fname.substr(0,fname.size()-4);
      for(unsigned i=0;i<arg_ends.size()-1;++i) {
          std::string tname = getPntrToOutput(i)->getName(); tname=tname.substr(0,tname.size()-5);
          if( tname==fname ) { arg0 = getPntrToArgument(arg_ends[i]); break; }
      }
      plumed_assert( arg0 ); 

      for(unsigned i=1;i<nlandmarks;++i) {
          for(unsigned j=0;j<i;++j) {
              double myval = arg0->get(nvals*landmarks[i] + landmarks[j] );
              val0->set( i*nlandmarks+j, myval ); val0->set( j*nlandmarks+i, myval );
          }
      }
  }
}

void LandmarkSelectionBase::calculate() {
  if( skipCalculate() ) return;
  plumed_dbg_assert( !actionInChain() );
  selectLandmarks(); setLandmarkSeparations();
}

void LandmarkSelectionBase::update() {
  if( skipUpdate() ) return;
  plumed_dbg_assert( !actionInChain() );
  selectLandmarks(); setLandmarkSeparations();
}

void LandmarkSelectionBase::runFinalJobs() {
  if( skipUpdate() ) {
      bool allsetup=true; 
      for(unsigned i=0;i<getNumberOfArguments();++i) {
          ActionSetup* as=dynamic_cast<ActionSetup*>( getPntrToArgument(i)->getPntrToAction() );
          if(!as) allsetup=false;
      }
      if( !allsetup ) return;
  }
  plumed_dbg_assert( !actionInChain() ); nvals=0;
  for(unsigned j=arg_ends[0];j<arg_ends[1];++j) {
      if( getPntrToArgument(j)->getRank()==1 ) nvals += getPntrToArgument(j)->getNumberOfValues( getLabel() );
      else if( getPntrToArgument(j)->getRank()==2 ) nvals += getPntrToArgument(j)->getShape()[1];
  }
  std::vector<unsigned> shape(2); shape[0]=nlandmarks; shape[1]=nvals;
  for(unsigned i=0;i<arg_ends.size()-1;++i) {
      unsigned tvals=0;
      for(unsigned j=arg_ends[i];j<arg_ends[i+1];++j) {
          if( getPntrToArgument(j)->getRank()==1 ) tvals += getPntrToArgument(j)->getNumberOfValues( getLabel() );
          else if( getPntrToArgument(j)->getRank()==2 ) tvals += getPntrToArgument(j)->getShape()[1];
      }
      if( tvals!=nvals ) error("mismatch between sizes of input positions");
      Value* val0 = getPntrToComponent(i); if( val0->getRank()==2 ) { val0->setShape( shape ); } 
  }
  selectLandmarks(); setLandmarkSeparations();
}

}
}
