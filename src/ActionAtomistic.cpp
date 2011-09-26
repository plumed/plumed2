#include "ActionAtomistic.h"
#include "PlumedMain.h"
#include <vector>
#include <string>
#include <cassert>
#include "ActionWithVirtualAtom.h"
#include "GenericGroup.h"

using namespace std;
using namespace PLMD; 

ActionAtomistic::~ActionAtomistic(){
// forget the pending request
  plumed.getAtoms().remove(this);
// Get rid of this actions group
  if( atomGroupName==getLabel() ) plumed.getAtoms().removeGroup(atomGroupName);
}

ActionAtomistic::ActionAtomistic(const ActionOptions&ao):
ActionWithExternalArguments(ao),
atomGroupName("none"),
pbcOn(true),
lockRequestAtoms(false)
{
  plumed.getAtoms().add(this);

  registerKeyword(0, "PBC", "(default) use periodic boundary conditions when calculating the vectors that connect atoms");
  registerKeyword(0, "NOPBC", "do not use periodic boundary conditions when calculating the vectors that connect atoms");
  registerKeyword(2, "ATOMS", "specify the atoms involved in the colvar as a list of atoms. To specify multiple colvars of this type use a list of ATOMS keywords: ATOMS1, ATOMS2, ATOMS3..." );
  registerKeyword(2, "GROUP", "specify the atoms involved in the colvar using one or two groups (GROUP1, GROUP2) of atoms. These groups can either be defined in another action or the atoms in the group can be enumerated on after the GROUP keywords. If one group is present the value of the colvar is calculated for every conceivable combination of atoms in the group. If multiple groups are specified the colvar is calculated using every combination that contains at least one atom from every group"); 
}

void ActionAtomistic::readActionAtomistic( int& maxatoms, unsigned& maxgroups ){
   readAction();
   std::vector<std::string> strings; Atoms& atoms(plumed.getAtoms());

   if( testForKey("ATOMS") && testForKey("GROUP") ) error("you cannot mix ATOMS and GROUP keywords");

   if( testForKey("GROUP") ){
       assert( maxgroups>0 );
       std::vector< std::vector<unsigned> > groups; std::vector<unsigned> t;
       if( testForNumberedKeys("GROUP") ){
          atomGroupName=getLabel();
          for(int i=1;i<=maxgroups;++i ){
             if( !parseNumberedVector( "GROUP", i, strings ) ) break;
             atoms.readAtomsIntoGroup( atomGroupName, strings, t );
             groups.push_back( t );
          }
       } else {
          parseVector("GROUP",strings);
          GenericGroup* grp=plumed.getActionSet().selectWithLabel<GenericGroup*>( strings[0] );
          if( grp ){
              if(strings.size()!=1) error("you can only use one group at a time when you specify an GROUP");
              atomGroupName=strings[0];
              atoms.getGroupIndices( atomGroupName, t ); 
              groups.push_back( t );
          } else {
              atomGroupName=getLabel();
              atoms.readAtomsIntoGroup( atomGroupName, strings, t );
              groups.push_back( t ); 
          }
       }
       interpretGroupsKeyword( maxatoms, atomGroupName, groups );
   } 
   if( testForKey("ATOMS") ){
       assert( maxatoms!=0 );
       std::vector< std::vector<unsigned> > flist; std::vector<unsigned> t;
       atomGroupName=getLabel();
       if( testForNumberedKeys("ATOMS") ){
          for(int i=1;;++i ){
             if( !parseNumberedVector( "ATOMS", i, strings ) ) break;
             atoms.readAtomsIntoGroup( atomGroupName, strings, t );
             if( i==1 && maxatoms<0 ){ 
                 maxatoms=t.size(); 
             } else if ( t.size()!=maxatoms ){
                 std::string nn, ss, tt; Tools::convert( i, nn ); Tools::convert( t.size(), ss ); Tools::convert( maxatoms, tt );
                 error("in ATOMS" + nn + " keyword found " + ss + " atoms when there should only be " + tt + "atoms"); 
             }
             flist.push_back( t );
          }
       } else {
          parseVector("ATOMS",strings);
          atoms.readAtomsIntoGroup( atomGroupName, strings, t );
          if( maxatoms<0 ){ 
             maxatoms=t.size(); 
          } else if ( t.size()!=maxatoms ) {
             std::string ss, tt; Tools::convert( t.size(), ss ); Tools::convert( maxatoms, tt );
             error("in ATOMS keyword found " + ss + " atoms when there should only be " + tt + " atoms");
          }
          flist.push_back( t );
       }
       interpretAtomsKeyword( flist );
   }

   // Get the indices for this group
   std::vector<unsigned> indexes;
   atoms.getGroupIndices( atomGroupName, indexes );

   // Resize everything
   unsigned nn=indexes.size();
   positions.resize(nn); masses.resize(nn); 
   charges.resize(nn); skips.resize(nn); forces.resize(nn);

   // Now sort out the dependencies
   clearDependencies();  
   for(unsigned j=0;j<indexes.size();++j){
     skips[j]=false;
     if(indexes[j]>=atoms.getNatoms()) addDependency(atoms.virtualAtomsActions[indexes[j]-atoms.getNatoms()]);
   }

   // Read periodic boundary condition stuff
   bool nopbc=!pbcOn; parseFlag("NOPBC",nopbc);
   pbcOn=!nopbc; parseFlag("PBC",pbcOn);
   if(pbcOn) log.printf("  using periodic boundary conditions\n");
   else log.printf("  without periodic boundary conditions\n");
}

void ActionAtomistic::calculateNumericalDerivatives(){
  const int nval=getNumberOfValues();
  const int natoms=positions.size(); 
  std::vector<Vector> value(nval*natoms);
  std::vector<Tensor> valuebox(nval);
  std::vector<Vector> savedPositions(natoms);
  const double delta=sqrt(epsilon);

  for(int i=0;i<natoms;i++){
    for(int k=0;k<3;k++){
       savedPositions[i][k]=positions[i][k];
       positions[i][k]=positions[i][k]+delta;
       calculate();
       positions[i][k]=savedPositions[i][k];
       for(int j=0;j<nval;j++){
         value[j*natoms+i][k]=getValue(j); 
       }
    }
 }
 for(int i=0;i<3;i++){ 
    for(int k=0;k<3;k++){
       double arg0=box(i,k);
       for(int j=0;j<natoms;j++) positions[j]=pbc.realToScaled(positions[j]);
       box(i,k)=box(i,k)+delta;
       pbc.setBox(box);
       for(int j=0;j<natoms;j++) positions[j]=pbc.scaledToReal(positions[j]);
       calculate();
       box(i,k)=arg0;
       pbc.setBox(box);
       for(int j=0;j<natoms;j++) positions[j]=savedPositions[j];
       for(int j=0;j<nval;j++) valuebox[j](i,k)=getValue(j); 
    }
 }

 calculate();
 clearDerivatives();
 for(int j=0;j<nval;j++){
    double ref=getValue(j);
    for(int i=0;i<natoms;i++){
      for(int k=0;k<3;k++){
         addDerivative( j, 3*i+k, (value[j*natoms+i][k]-ref)/delta );
      } 
    }
    Tensor virial;
    for(int i=0;i<3;i++) for(int k=0;k<3;k++) virial(i,k) = (valuebox[j](i,k)-ref)/delta;
// BE CAREFUL WITH NON ORTHOROMBIC CELL
    virial=-1.0*matmul(box.transpose(),virial.transpose());
    for(int i=0;i<3;i++) for(int k=0;k<3;k++) addDerivative( j, 3*natoms+3*k+i, virial(i,k) ); 
  }
}

void ActionAtomistic::retrieveData(){
  box=plumed.getAtoms().box; 
  pbc.setBox(box);
  plumed.getAtoms().getAtomsInGroup( atomGroupName, positions, charges, masses );
}

void ActionAtomistic::applyForces( const std::vector<Vector>& forces, const Tensor& virial ){
  plumed.getAtoms().applyForceToAtomsInGroup( atomGroupName, forces, virial );
}
