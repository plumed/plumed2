#include "plumed/reference/ReferenceConfiguration.h"
#include "plumed/reference/MultiReferenceBase.h"
#include "plumed/reference/PathReparameterization.h"
#include "plumed/tools/AtomNumber.h"
#include "plumed/core/Value.h"
#include "plumed/tools/PDB.h"
#include "plumed/tools/Pbc.h"
#include "plumed/tools/OFile.h"

bool end( const int& ii, const int& last, const int& inc ){
  if( inc>0 && ii<last ) return false;
  if( inc<0 && ii>last ) return false;
  return true;
}

int main(){
  // Test path reparamterization with euclidean distances between frames
  std::string mtype="EUCLIDEAN";
  std::string inputfilename="epath.pdb"; 
  PLMD::MultiReferenceBase mymap( mtype, false ); 

  FILE* fp=fopen(inputfilename.c_str(),"r");
  bool do_read=true; unsigned nfram=0;
  while (do_read){
     // Read the pdb file
     PLMD::PDB mypdb; do_read=mypdb.readFromFilepointer(fp,false,1.0);
     if( do_read ){
         mymap.readFrame( mypdb ); nfram++;
     } else {
         break;
     }
  }
  fclose(fp);
  std::vector<PLMD::AtomNumber> atoms; std::vector<std::string> args;
  mymap.getAtomAndArgumentRequirements( atoms, args );

  // Make a fake vector of Values 
  std::vector<PLMD::Value*> fvals; PLMD::Pbc fpbc;
  for(unsigned i=0;i<(mymap.getFrame( 0 ))->getNumberOfReferenceArguments();++i){
      fvals.push_back( new PLMD::Value() ); fvals[i]->setNotPeriodic();
  } 
  // Reparameterize the path
  PLMD::PathReparameterization myreparam = PLMD::PathReparameterization( fpbc, fvals, mymap.getReferenceConfigurations() );
  myreparam.reparameterize( 0, nfram-1, 1.E-6 );
  // Needs something to print reparameterized path
  PLMD::OFile ofile; ofile.open("epath-out.pdb");
  std::vector<PLMD::ReferenceConfiguration*>& oframes=mymap.getReferenceConfigurations();
  for(unsigned i=0;i<oframes.size();++i){ oframes[i]->print( ofile, "%8.4f", 1.0 ); }
  ofile.close();


  // Delete pointers
  for(unsigned i=0;i<fvals.size();++i) delete fvals[i];
  fvals.resize(0);

  // Test path reparameterization with RMSD distances between frames
  mtype="OPTIMAL"; PLMD::MultiReferenceBase mymap2( mtype, false );
  FILE* fp2=fopen("all.pdb","r");
  do_read=true; nfram=0;
  while (do_read){
     // Read the pdb file
     PLMD::PDB mypdb; do_read=mypdb.readFromFilepointer(fp2,false,1.0);
     if( do_read ){
         mymap2.readFrame( mypdb ); nfram++;
     } else {
         break;
     }
  }
  fclose(fp2); atoms.resize(0); args.resize(0);
  mymap2.getAtomAndArgumentRequirements( atoms, args );

  // Reparaemterize the path
  PLMD::PathReparameterization myreparam2 = PLMD::PathReparameterization( fpbc, fvals, mymap2.getReferenceConfigurations() );
  myreparam2.reparameterize( 3, 6, 1.E-6 );

  // Needs something to print reparameterized path
  PLMD::OFile ofile2; ofile2.open("path-out.pdb");
  std::vector<PLMD::ReferenceConfiguration*>& oframes2=mymap2.getReferenceConfigurations();
  for(unsigned i=0;i<oframes2.size();++i){ oframes2[i]->print( ofile2, "%8.4f", 1.0 ); }
  ofile2.close();

  return 0;
}
