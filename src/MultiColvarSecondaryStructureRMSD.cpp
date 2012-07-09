#include "MultiColvarSecondaryStructureRMSD.h"
#include "PlumedMain.h"
#include "Atoms.h"

namespace PLMD {

void MultiColvarSecondaryStructureRMSD::registerKeywords( Keywords& keys ){
  MultiColvar::registerKeywords( keys );
  ActionWithDistribution::autoParallelize( keys );
  keys.add("atoms","BACKBONE","");
  keys.add("compulsory","TYPE","DRMSD","the manner in which RMSD alignment is performed. Should be OPTIMAL, SIMPLE or DRMSD.");
  keys.remove("MORE_THAN");
}

MultiColvarSecondaryStructureRMSD::MultiColvarSecondaryStructureRMSD(const ActionOptions&ao):
Action(ao),
MultiColvar(ao)
{
  parse("TYPE",alignType);
  log.printf("  distances from secondary structure elements are calculated using %s algorithm\n",alignType.c_str() );
}

MultiColvarSecondaryStructureRMSD::~MultiColvarSecondaryStructureRMSD(){
  for(unsigned i=0;i<secondary_rmsd.size();++i) delete secondary_rmsd[i]; 
  for(unsigned i=0;i<secondary_drmsd.size();++i) delete secondary_drmsd[i];
}

void MultiColvarSecondaryStructureRMSD::setSecondaryStructure( std::vector<Vector>& structure, double bondlength, double units ){

  // Convert into correct units
  for(unsigned i=0;i<structure.size();++i){
     structure[i][0]*=units; structure[i][1]*=units; structure[i][2]*=units;
  }

  if( secondary_drmsd.size()==0 && secondary_rmsd.size()==0 ){ 
     int natoms; readAtoms(natoms); requestDistribution();
  }

  // If we are in natural units get conversion factor from nm into natural length units
  if( plumed.getAtoms().usingNaturalUnits() ){
      error("cannot use this collective variable when using natural units");
  }

  // Set the reference structure
  new_deriv.resize( structure.size() );
  if( alignType=="DRMSD" ){
    secondary_drmsd.push_back( new DRMSD() ); 
    secondary_drmsd[secondary_drmsd.size()-1]->setReference( bondlength, structure );
  } else {
    std::vector<double> align( structure.size(), 1.0 );
    secondary_rmsd.push_back( new RMSD(log) );
    secondary_rmsd[secondary_rmsd.size()-1]->setType( alignType );
    secondary_rmsd[secondary_rmsd.size()-1]->setReference( structure );
    secondary_rmsd[secondary_rmsd.size()-1]->setAlign( align ); 
    secondary_rmsd[secondary_rmsd.size()-1]->setDisplace( align );
  }
}

double MultiColvarSecondaryStructureRMSD::compute( const unsigned& j, const std::vector<Vector>& pos, std::vector<Vector>& deriv, Tensor& virial ){
  double r,nr; Tensor new_virial;

  if( secondary_drmsd.size()>0 ){
    r=secondary_drmsd[0]->calculate( pos, deriv, virial ); 
    for(unsigned i=1;i<secondary_drmsd.size();++i){
        nr=secondary_drmsd[i]->calculate( pos, new_deriv, new_virial );
        if(nr<r){
           r=nr;
           for(unsigned i=0;i<new_deriv.size();++i) deriv[i]=new_deriv[i];
           virial=new_virial; 
        }
    }
  } else {
    r=secondary_rmsd[0]->calculate( pos, deriv );
    for(unsigned i=1;i<secondary_rmsd.size();++i){
        nr=secondary_rmsd[i]->calculate( pos, new_deriv );
        if(nr<r){
           r=nr;
           for(unsigned i=0;i<new_deriv.size();++i) deriv[i]=new_deriv[i];
        }
    } 
    for(unsigned i=0;i<deriv.size();i++) virial=virial+(-1.0*Tensor(pos[i],deriv[i]));
  }
  return r;
}

unsigned MultiColvarSecondaryStructureRMSD::getNumberOfFieldDerivatives(){
  plumed_massert(0,"Fields are not allowed for secondary structure variables");
}

}
