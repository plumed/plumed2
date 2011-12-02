#include "ActionAtomistic.h"
#include "PlumedMain.h"
#include <vector>
#include <string>
#include <cassert>
#include "ActionWithVirtualAtom.h"

using namespace std;
using namespace PLMD; 

ActionAtomistic::~ActionAtomistic(){
// forget the pending request
  plumed.getAtoms().remove(this);
// Get rid of this actions group
  assert( atomGroupName==getLabel() );
  plumed.getAtoms().removeGroup(atomGroupName);
  if( forcefile ) fclose(forcefile);
}

ActionAtomistic::ActionAtomistic(const ActionOptions&ao):
ActionWithExternalArguments(ao),
atomGroupName("none"),
pbcOn(true),
forcefile(NULL),
lockRequestAtoms(false),
nliststyle(3),
updateFreq(0),
nl_cut(0),
lastUpdate(0)
{
  plumed.getAtoms().add(this);

  registerKeyword(0, "PBC", "(default) use periodic boundary conditions when calculating the vectors that connect atoms");
  registerKeyword(0, "NOPBC", "do not use periodic boundary conditions when calculating the vectors that connect atoms");
  registerKeyword(0, "UPDATE", "frequency with which to update all neighbour lists");
  registerKeyword(0, "NL_CUT", "cutoff for neighbour lists");
  registerKeyword(0, "DUMPFORCES", "(debug only) this keyword can be used to force action atomistic to print out the forces on it every step so that they can be debugged");
}

void ActionAtomistic::setNeighbourListStyle( const std::string style ){
  if ( style=="none" ){ 
     nliststyle=0; forbidKeyword("UPDATE"); forbidKeyword("NL_CUT"); 
  } else  if ( style=="skipAll" ){
     nliststyle=1;
  } else if ( style=="skipBigOnly" ){
     nliststyle=2;
  }
}

void ActionAtomistic::addColvar( AtomicNeighbourList& newlist ){
   assert( newlist.neighbours.size()==0 );   // Neighbour lists should be set up at the end using setupNeighbourList
   newlist.style=nliststyle; newlist.rcut=nl_cut;
   nlists.push_back( newlist );

   // Report the atoms in the log file
   Atoms& atoms(plumed.getAtoms());
   log.printf("  Creating %d th cv from atoms : %s", nlists.size(), atoms.interpretIndex( atomGroupName, newlist.all_atoms[0] ).c_str() );
   for(unsigned j=1;j<newlist.all_atoms.size();++j){ log.printf(", %s", atoms.interpretIndex( atomGroupName, newlist.all_atoms[j] ).c_str() ); }
   log.printf("\n");
}

void ActionAtomistic::setupNeighbourList( const std::vector< std::pair<unsigned, unsigned> >& nlist_template ){
   assert( nliststyle>0 );
   for(unsigned i=0;i<nlists.size();++i){
      for(unsigned j=0;j<nlist_template.size();++j){ 
          nlists[i].addPair( nlists[i].all_atoms[ nlist_template[j].first ], nlists[i].all_atoms[ nlist_template[j].second ] );
      }
   }
}

void ActionAtomistic::checkNeighbourLists() const {
   if( nliststyle>0 ){
       for(unsigned i=0;i<nlists.size();++i) {
          if ( nlists[i].neighbours.size()==0 ) error("You have specified that you want to use a neighbour list but you haven't set it up");
       } 
   }
}

void ActionAtomistic::checkForBadNeighbourLists() const {
   if( nliststyle==1 && updateFreq>0 ) error("neighbour lists can only be used with this colvar when you are using modifiers");
}

void ActionAtomistic::readActionAtomistic(){
   assert( nliststyle==0 || nliststyle==1 || nliststyle==2 ); 

   // Read in action stuff
   readAction();

   // Read in neighbour list stuff
   if ( nliststyle!=0 ){
        parse("UPDATE",updateFreq); 
        if(updateFreq>0){
           parse("NL_CUT",nl_cut);
           if( nl_cut==0 ) error("you have not specified a cutoff for neighbour lists use NL_CUT");

           log.printf("  updating neighbour lists every %d steps. ", updateFreq);
           if(nliststyle==1){
               log.printf("Skipping calculation if any pair is greater than %f\n",nl_cut);
           } else if(nliststyle==2){
               log.printf("Excluding pairs of atoms separated by more than %f\n",nl_cut);
           } else {
               assert(false);   // You have not specified a valid neighbour list style
           }
        }
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

bool ActionAtomistic::parseAtomList( const std::string& key, const unsigned num ){
  std::vector<std::string> strings; 
 
  bool readStuff=true;
  if(num==0){ 
     parseVector( key, strings ); 
  } else {
     readStuff=parseNumberedVector( key, num, strings );
  }
  if(!readStuff) return false;

  if( atomGroupName=="none" ) atomGroupName=getLabel();
  assert( atomGroupName==getLabel() );

  // Read the atoms into the colvar
  Atoms& atoms(plumed.getAtoms()); std::vector<unsigned> t;
  atoms.readAtomsIntoGroup( atomGroupName, strings, t );

  // Get the indices for this group
  std::vector<unsigned> indexes;
  atoms.getGroupIndices( atomGroupName, indexes );

   // Resize everything
   unsigned nn=indexes.size();
   positions.resize(nn); masses.resize(nn); 
   charges.resize(nn); skips.resize(nn);

   // Now sort out the dependencies
   clearDependencies();
   // Virtual atoms
   for(unsigned j=0;j<indexes.size();++j){
     skips[j]=false;
     if(indexes[j]>=atoms.getNatoms()) addDependency(atoms.virtualAtomsActions[indexes[j]-atoms.getNatoms()]);
   }
   return true;
}

void ActionAtomistic::printAllAtoms( const std::string report ){
   // Report the atoms in the log file
   Atoms& atoms(plumed.getAtoms());
   log.printf("  %s : %s",report.c_str(), atoms.interpretIndex( atomGroupName, 0 ).c_str() );
   for(unsigned j=1;j<getNumberOfAtoms();++j){ log.printf(", %s", atoms.interpretIndex( atomGroupName, j ).c_str() ); }
   log.printf("\n");
}

bool ActionAtomistic::readBackboneAtoms( const std::string& type, std::vector< std::vector<unsigned> >& backbone ){
   if ( testForKey("BACKBONE") ){
        assert( atomGroupName=="none" ); atomGroupName=getLabel(); 
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
  assert( atomGroupName==getLabel() );

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
  assert( atomGroupName==getLabel() );
  if( updateFreq>0 && (getStep()-lastUpdate)>=updateFreq ){
     for(unsigned i=0;i<skips.size();++i) skips[i]=false;
     plumed.getAtoms().updateSkipsForGroup( atomGroupName, skips );
  }
}

void ActionAtomistic::updateNeighbourLists(){
  assert(nliststyle>0);

  if( updateFreq>0 && (getStep()-lastUpdate)>=updateFreq ){
      for(unsigned i=0;i<skips.size();++i) skips[i]=true;
      for(unsigned j=0;j<nlists.size();++j){ nlists[j].update( skips ); }
      plumed.getAtoms().updateSkipsForGroup( atomGroupName, skips );
      lastUpdate=getStep();
  } 
}

void ActionAtomistic::retrieveData(){
  assert( atomGroupName==getLabel() );
  box=plumed.getAtoms().box; 
  pbc.setBox(box);
  plumed.getAtoms().getAtomsInGroup( atomGroupName, positions, charges, masses );
}

void ActionAtomistic::applyForces( const std::vector<Vector>& forces, const Tensor& virial ){
  assert( atomGroupName==getLabel() );
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
