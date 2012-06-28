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
MultiColvar(ao),
secondary_rmsd(NULL),
secondary_drmsd(NULL)
{
  std::string tstring;
  parse("TYPE",tstring);
  log.printf("  distances from secondary structure elements are calculated using %s algorithm\n",tstring.c_str() );
  if( tstring=="DRMSD" ){
     secondary_drmsd=new DRMSD();
  } else {
     secondary_rmsd=new RMSD(log);
     secondary_rmsd->setType( tstring );
  }
}

MultiColvarSecondaryStructureRMSD::~MultiColvarSecondaryStructureRMSD(){
  delete secondary_rmsd; delete secondary_drmsd;
}

void MultiColvarSecondaryStructureRMSD::setSecondaryStructure( std::vector<Vector>& structure, double bondlength, double units ){

  // Convert into correct units
  for(unsigned i=0;i<structure.size();++i){
     structure[i][0]*=units; structure[i][1]*=units; structure[i][2]*=units;
  }

  int natoms; readAtoms(natoms);
  requestDistribution();

  // If we are in natural units get conversion factor from nm into natural length units
  if( plumed.getAtoms().usingNaturalUnits() ){
      error("cannot use this collective variable when using natural units");
  }

  // Set the reference structure
  if( secondary_drmsd ){
    secondary_drmsd->setReference( bondlength, structure );
  } else {
    std::vector<double> align( structure.size(), 1.0 );
    secondary_rmsd->setReference( structure );
    secondary_rmsd->setAlign( align ); 
    secondary_rmsd->setDisplace( align );
  }
}

double MultiColvarSecondaryStructureRMSD::compute( const unsigned& j, const std::vector<Vector>& pos, std::vector<Vector>& deriv, Tensor& virial ){
  double r;
  if( secondary_drmsd ){
    r=secondary_drmsd->calculate( pos, deriv, virial ); 
  } else {
    r=secondary_rmsd->calculate( pos, deriv );
    for(unsigned i=0;i<deriv.size();i++) virial=virial+(-1.0*Tensor(pos[i],deriv[i]));
  }
  return r;
}

unsigned MultiColvarSecondaryStructureRMSD::getNumberOfFieldDerivatives(){
  plumed_massert(0,"Fields are not allowed for secondary structure variables");
}

}
