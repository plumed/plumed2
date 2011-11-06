#include "ActionAtomistic.h"
#include "PlumedMain.h"
#include <vector>
#include <string>
#include <cassert>
#include "ActionWithVirtualAtom.h"
#include "Group.h"

using namespace std;
using namespace PLMD; 

ActionAtomistic::~ActionAtomistic(){
// forget the pending request
  plumed.getAtoms().remove(this);
// Get rid of this actions group
  if( atomGroupName==getLabel() ) plumed.getAtoms().removeGroup(atomGroupName);
  if( forcefile ) fclose(forcefile);
}

ActionAtomistic::ActionAtomistic(const ActionOptions&ao):
ActionWithExternalArguments(ao),
atomGroupName("none"),
doneRead(false),
pbcOn(true),
forcefile(NULL),
lockRequestAtoms(false),
updateFreq(0),
lastUpdate(0),
nl_cut(0)
{
  plumed.getAtoms().add(this);

  registerKeyword(0, "PBC", "(default) use periodic boundary conditions when calculating the vectors that connect atoms");
  registerKeyword(0, "NOPBC", "do not use periodic boundary conditions when calculating the vectors that connect atoms");
  registerKeyword(2, "ATOMS", "specify the atoms involved in the colvar as a list of atoms. To specify multiple colvars of this type use a list of ATOMS keywords: ATOMS1, ATOMS2, ATOMS3..." );
  registerKeyword(2, "GROUP", "specify the atoms involved in the colvar using one or two groups (GROUP1, GROUP2) of atoms. These groups can either be defined in another action or the atoms in the group can be enumerated on after the GROUP keywords. If one group is present the value of the colvar is calculated for every conceivable combination of atoms in the group. If multiple groups are specified the colvar is calculated using every combination that contains at least one atom from every group"); 
  registerKeyword(0, "UPDATE", "frequency for updates of neighbour lists and dynamic groups");
  registerKeyword(0, "NL_CUTOFF", "the cutoff for distances inside the neighbour list");
  registerKeyword(0, "DG_CUTOFF", "dynamic groups have a quanity between 1 and 0 which measures the extent to which a quantity can be said to be a member.  This states that atoms whole assignment quantity is less than this value can safely be ignored");
  registerKeyword(0, "DUMPFORCES", "(debug only) this keyword can be used to force action atomistic to print out the forces on it every step so that they can be debugged");
}

void ActionAtomistic::readActionAtomistic( int& maxatoms, unsigned& maxgroups ){
   std::vector<std::string> strings; Atoms& atoms(plumed.getAtoms());

   if( testForKey("ATOMS") && testForKey("GROUP") ) error("you cannot mix ATOMS and GROUP keywords");
   parse("UPDATE",updateFreq);

   parse("NL_CUTOFF",nl_cut);
   if(nl_cut==0) parse("DG_CUTOFF",nl_cut);
   else error("no action should define both a NL_CUTOFF and a DG_CUTOFF");

   if(!doneRead){
      readAction();
      if( testForKey("GROUP") ){
          assert( maxgroups>0 ); assert( atomGroupName=="none" );
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
             Group* grp=plumed.getActionSet().selectWithLabel<Group*>( strings[0] );
             if( grp ){
                 grp->setUpdateFreq( updateFreq ); // Do we really need this - i.e. to force neighbour lists to update at the same time as dynamic groups
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
          interpretGroupsKeyword( maxatoms, atomGroupName, groups ); doneRead=true;
      } 
      if( testForKey("ATOMS") ){
          assert( maxatoms!=0 ); assert( atomGroupName=="none" );
          std::vector< std::vector<unsigned> > flist; std::vector<unsigned> t;
          atomGroupName=getLabel();
          if( testForNumberedKeys("ATOMS") ){
             for(int i=1;;++i ){
                if( !parseNumberedVector( "ATOMS", i, strings ) ) break;
                atoms.readAtomsIntoGroup( atomGroupName, strings, t );
                if( i==1 && maxatoms<0 ){ 
                    maxatoms=t.size(); 
                } else if ( t.size()!=maxatoms ){
                    std::string nn, ss, tt; Tools::convert( i, nn ); Tools::convert( static_cast<int>( t.size() ), ss ); Tools::convert( maxatoms, tt );
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
                std::string ss, tt; Tools::convert( static_cast<int>( t.size() ), ss ); Tools::convert( maxatoms, tt );
                error("in ATOMS keyword found " + ss + " atoms when there should only be " + tt + " atoms");
             }
             flist.push_back( t );
          }
          interpretAtomsKeyword( flist ); doneRead=true;
      }
   }
   if( !doneRead ) error("I have not read in the list of atoms involved in this action");

   // Get the indices for this group
   std::vector<unsigned> indexes;
   atoms.getGroupIndices( atomGroupName, indexes );

   // Resize everything
   unsigned nn=indexes.size();
   positions.resize(nn); masses.resize(nn); gderivs.resize(nn); 
   charges.resize(nn); skips.resize(nn); forces.resize(nn);
   // This is stuff for super groups
   group_f.resize(nn); group_df.resize(nn);

   // Now sort out the dependencies
   clearDependencies();  
   // Virtual atoms
   for(unsigned j=0;j<indexes.size();++j){
     skips[j]=false;
     if(indexes[j]>=atoms.getNatoms()) addDependency(atoms.virtualAtomsActions[indexes[j]-atoms.getNatoms()]);
   }
   // Groups
   if( atomGroupName!=getLabel() ){
       Group* grp=plumed.getActionSet().selectWithLabel<Group*>( atomGroupName );
       addDependency( grp );
   }

   // Read periodic boundary condition stuff
   bool nopbc=!pbcOn; parseFlag("NOPBC",nopbc);
   pbcOn=!nopbc; parseFlag("PBC",pbcOn);
   if(pbcOn) log.printf("  using periodic boundary conditions\n");
   else log.printf("  without periodic boundary conditions\n");

  // Read force dumping stuff
  std::string forcefilename; parse("DUMPFORCES",forcefilename); 
  if( forcefilename.length()>0 ){
    if(comm.Get_rank()==0){
      forcefile=fopen(forcefilename.c_str(),"wa");
      log.printf("  writing forces on file %s to debug\n",forcefilename.c_str());
      fprintf(forcefile,"#! FIELDS time parameter force");
      fprintf(forcefile,"\n");
    }
  }
}

bool ActionAtomistic::readBackboneAtoms( const std::string& type, std::vector< std::vector<unsigned> >& backbone ){
   readAction(); 
   if ( testForKey("BACKBONE") ){
        assert( atomGroupName=="none" ); atomGroupName=getLabel(); doneRead=true;
        std::vector<std::string> residues; parseVector("BACKBONE",residues);
        plumed.getAtoms().putBackboneInGroup( atomGroupName, type, residues, backbone ); 
        return true;
   }
   return false;
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
       if( atomGroupName!=getLabel() ){ 
           Group* grp=plumed.getActionSet().selectWithLabel<Group*>( atomGroupName ); grp->calculate();
       }
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
       if( atomGroupName!=getLabel() ){
           Group* grp=plumed.getActionSet().selectWithLabel<Group*>( atomGroupName ); grp->calculate();
       }
       calculate();
       box(i,k)=arg0;
       pbc.setBox(box);
       for(int j=0;j<natoms;j++) positions[j]=savedPositions[j];
       for(int j=0;j<nval;j++) valuebox[j](i,k)=getValue(j); 
    }
 }

 if( atomGroupName!=getLabel() ){
     Group* grp=plumed.getActionSet().selectWithLabel<Group*>( atomGroupName ); grp->calculate();
 }
 calculate();
 clearDerivatives(); bool isvatom=false;
 ActionWithVirtualAtom* vatom=dynamic_cast<ActionWithVirtualAtom*>(this);
 if ( vatom ){ isvatom=true; } 
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
    if (!isvatom) {
       for(int i=0;i<3;i++) for(int k=0;k<3;k++) addDerivative( j, 3*natoms+3*k+i, virial(i,k) ); 
    }
  }
}

void ActionAtomistic::prepare(){
  if( updateFreq>0 && (getStep()-lastUpdate)>=updateFreq ){
     for(unsigned i=0;i<skips.size();++i) skips[i]=false;
     plumed.getAtoms().updateSkipsForGroup( atomGroupName, skips );
  }
}

void ActionAtomistic::retrieveData(){
  box=plumed.getAtoms().box; 
  pbc.setBox(box);
  plumed.getAtoms().getAtomsInGroup( atomGroupName, positions, charges, masses );
}

void ActionAtomistic::calculateAtomisticActions(){
  if( updateFreq>0 && (getStep()-lastUpdate)>=updateFreq ){
      skips.assign(skips.size(),false); 
      if( atomGroupName!=getLabel() ){
         Group* grp=plumed.getActionSet().selectWithLabel<Group*>( atomGroupName );
         grp->retrieveSkips( skips );  
      }
      updateDynamicContent( nl_cut, skips );   
      plumed.getAtoms().updateSkipsForGroup( atomGroupName, skips );
      lastUpdate=getStep();
  }

  if( atomGroupName!=getLabel() ){
     Group* grp=plumed.getActionSet().selectWithLabel<Group*>( atomGroupName );
     group_val=grp->getGroupData( group_f, group_df, group_vir );   
  }
}

void ActionAtomistic::retrieveSkips( std::vector<bool>& s ) const {
  assert( s.size()==skips.size() );
  for(unsigned i=0;i<skips.size();++i){
     if( !s[i] && skips[i] ) s[i]=true;
  }
}

void ActionAtomistic::applyForces( const std::vector<Vector>& forces, const Tensor& virial ){
  plumed.getAtoms().applyForceToAtomsInGroup( atomGroupName, forces, virial );
  if( forcefile ){
      for(unsigned i=0;i<getNumberOfAtoms();++i){
          fprintf( forcefile, "%f %d %f \n", getTime(), 3*i+0, forces[i][0] );
          fprintf( forcefile, "%f %d %f \n", getTime(), 3*i+1, forces[i][1] );
          fprintf( forcefile, "%f %d %f \n", getTime(), 3*i+2, forces[i][2] );
      }
      unsigned natoms=getNumberOfAtoms();
      for(unsigned j=0;j<3;++j) for(unsigned k=0;k<3;++k){
         fprintf( forcefile, "%f %d %f \n", getTime(), 3*natoms+3*j+k, virial(j,k) );
      }
  } 
}
