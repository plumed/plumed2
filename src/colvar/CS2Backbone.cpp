#ifdef __PLUMED_HAS_ALMOST

#include "Colvar.h"
#include "ActionRegister.h"
#include "tools/Matrix.h"
#include "tools/Communicator.h"

#include <string>
#include <cmath>
#include <cassert>

#include <almost/mdb.h>
#include <almost/pdb.h>
#include <almost/forcefield/const/camshift2.h>
#include <almost/io/formostream.h>

using namespace std;
using namespace Almost;

namespace PLMD{

//+PLUMEDOC COLVAR CS2-BACKBONE 
/*
*/
//+ENDPLUMEDOC

class ColvarCS2Backbone : public Colvar {
  vector<CamShift2> cam_list;
  Molecules molecules;
  int  numResidues;
  int  pperiod;
  int  ens_dim;
  bool ensemble;
  bool serial;
  double **sh;
public:
  ColvarCS2Backbone(const ActionOptions&);
  ~ColvarCS2Backbone();
  static void registerKeywords( Keywords& keys );
  virtual void calculate();
};

PLUMED_REGISTER_ACTION(ColvarCS2Backbone,"CS2BACKBONE")

void ColvarCS2Backbone::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords( keys );
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.add("atoms","ATOMS","The atoms to be included in the calculatios, e.g. the whole protein.");
  keys.add("compulsory","DATA","data/","The folder with the experimental chemical shifts.");
  keys.add("compulsory","FF","a03_gromacs.mdb","The ALMOST force-field to map the atoms' names.");
  keys.add("compulsory","FLAT","1.0","Flat region in the scoring function.");
  keys.add("compulsory","NEIGH_FREQ","10","Period in step for neighbour list update.");
  keys.add("compulsory","WRITE_CS","0","Write chemical shifts period.");
  keys.add("compulsory","NRES","Number of residues, corresponding to the number of chemical shifts.");
  keys.addFlag("CYS-DISU",false,"Set to TRUE if your system has disulphide bridges");  
  keys.addFlag("ENSEMBLE",false,"Set to TRUE if you want to average over multiple replicas");  
}

ColvarCS2Backbone::ColvarCS2Backbone(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao)
{
  string stringadb;
  string stringamdb;
  string stringapdb;

  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  checkRead();

  serial=false;
  parseFlag("SERIAL",serial);

  string stringa_data;
  parse("DATA",stringa_data);

  string stringa_forcefield;
  parse("FF",stringa_forcefield);

  bool disu=false;
  parseFlag("CYS-DISU",disu);

  double grains=1.0;
  parse("FLAT", grains);

  int neigh_f=10;
  parse("NEIGH_FREQ", neigh_f);

  int w_period=0;
  parse("WRITE_CS", w_period);
  pperiod=w_period;

  parse("NRES", numResidues);

  ensemble=false;
  parseFlag("ENSEMBLE",ensemble);
  if(multi_sim_com.Get_size()<2) plumed_merror("You CANNOT run Replica-Averaged simulations without running multiple replicas!\n");

  stringadb  = stringa_data + string("/camshift.db");
  stringamdb = stringa_data + string("/") + stringa_forcefield;
  stringapdb = stringa_data + string("/template.pdb");

  log.printf("  loading force-field %s\n", stringamdb.c_str()); log.flush();
  Almost::MDB mdb((char*)stringamdb.c_str());
  log.printf("  loading template %s\n", stringapdb.c_str()); log.flush();
  Almost::PDB pdb((char*)stringapdb.c_str());

  log.printf("  building molecule ..."); log.flush();
  for(unsigned i=0;i<pdb[0].size();i++){
    string str;
    str +='A'+i;
    Protein p(str);
    p.build_missing(pdb[0][i],mdb,"DEFAULT","NONE");
    if(disu) p.auto_disu_bonds(2.9,mdb);
    molecules.add_protein(p);
  }
  log.printf(" done!\n"); log.flush(); 

  log.printf("  Initialization of the predictor ...\n"); log.flush();
  CamShift2 a = CamShift2(molecules, stringadb);

  log.printf("  Reading experimental data ...\n"); log.flush();
  stringadb = stringa_data + string("/CAshifts.dat");
  log.printf("  Initializing CA shifts %s\n", stringadb.c_str()); log.flush();
  a.read_cs(stringadb, "CA");
  stringadb = stringa_data + string("/CBshifts.dat");
  log.printf("  Initializing CB shifts %s\n", stringadb.c_str()); log.flush();
  a.read_cs(stringadb, "CB");
  stringadb = stringa_data + string("/Cshifts.dat");
  log.printf("  Initializing C' shifts %s\n", stringadb.c_str()); log.flush();
  a.read_cs(stringadb, "C");
  stringadb = stringa_data + string("/HAshifts.dat");
  log.printf("  Initializing HA shifts %s\n", stringadb.c_str()); log.flush();
  a.read_cs(stringadb, "HA");
  stringadb = stringa_data + string("/Hshifts.dat");
  log.printf("  Initializing H shifts %s\n", stringadb.c_str()); log.flush();
  a.read_cs(stringadb, "H");
  stringadb = stringa_data + string("/Nshifts.dat");
  log.printf("  Initializing N shifts %s\n", stringadb.c_str()); log.flush();
  a.read_cs(stringadb, "N");
  /* this is a workaround for those chemical shifts that can result in too large forces */
  a.remove_problematic("GLN","CB");
  a.remove_problematic("ILE","CB");
  /* this is a workaround for those chemical shifts that are not parameterized */
  a.remove_problematic("HIE", "HA"); a.remove_problematic("HIP", "HA"); a.remove_problematic("HSP", "HA");
  a.remove_problematic("HIE", "H");  a.remove_problematic("HIP", "H");  a.remove_problematic("HSP", "H"); 
  a.remove_problematic("HIE", "N");  a.remove_problematic("HIP", "N");  a.remove_problematic("HSP", "N"); 
  a.remove_problematic("HIE", "CA"); a.remove_problematic("HIP", "CA"); a.remove_problematic("HSP", "CA");
  a.remove_problematic("HIE", "CB"); a.remove_problematic("HIP", "CB"); a.remove_problematic("HSP", "CB");
  a.remove_problematic("HIE", "C");  a.remove_problematic("HIP", "C");  a.remove_problematic("HSP", "C"); 
  a.remove_problematic("GLH", "HA"); a.remove_problematic("ASH", "HA"); a.remove_problematic("HSE", "HA");
  a.remove_problematic("GLH", "H");  a.remove_problematic("ASH", "H");  a.remove_problematic("HSE", "H");
  a.remove_problematic("GLH", "N");  a.remove_problematic("ASH", "N");  a.remove_problematic("HSE", "N");
  a.remove_problematic("GLH", "CA"); a.remove_problematic("ASH", "CA"); a.remove_problematic("HSE", "CA");
  a.remove_problematic("GLH", "CB"); a.remove_problematic("ASH", "CB"); a.remove_problematic("HSE", "CB");
  a.remove_problematic("GLH", "C");  a.remove_problematic("ASH", "C");  a.remove_problematic("HSE", "C");
  if(disu) { 
    a.remove_problematic("CYS", "HA");
    a.remove_problematic("CYS", "H");
    a.remove_problematic("CYS", "N");
    a.remove_problematic("CYS", "CA");
    a.remove_problematic("CYS", "CB");
    a.remove_problematic("CYS", "C");
  }
  /* done */

  log.printf("  Setting parameters ...\n"); log.flush();
  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  if(serial) {stride=1; rank=0;}
  if(stride>1) log.printf("  Parallelized over %i processors\n", stride);
  a.set_mpi(stride, rank);
  
  if(ensemble) { ens_dim=multi_sim_comm.Get_size(); log.printf("  ENSEMBLE averaging over %i replicas\n", ens_dim); }

  a.set_flat_bottom_const(grains);
  a.set_box_nupdate(neigh_f);
  a.set_lambda(1);
  cam_list.push_back(a);
  log.printf("  Writing converted template.pdb ...\n"); log.flush();
  mol2pdb(molecules,"converted-template.pdb");

  sh = new double*[numResidues];
  sh[0] = new double[numResidues*6];
  for ( int i = 1 ; i < numResidues ; i++)  sh[i] = sh[i-1] + 6; 

  addValueWithDerivatives();
  setNotPeriodic();
  requestAtoms(atoms);
  log.printf("  DONE!\n"); log.flush();
}

ColvarCS2Backbone::~ColvarCS2Backbone()
{
  delete[] sh[0];
  delete[] sh;
}


void ColvarCS2Backbone::calculate()
{
  double energy=0.;
  Tensor virial;
  virial.zero();
  vector<Vector> deriv(getNumberOfAtoms());
  int N = getNumberOfAtoms();
  Coor<double> coor(N); 
  Coor<double> forces(N);

  forces.clear();
  for(unsigned i=0; i<numResidues; i++) for(unsigned j=0; j<6; j++) sh[i][j]=0.;

  for (int i = 0; i < N; i++) {
     int ipos = 4 * i;
     Vector Pos = getPosition(i);
     coor.coor[ipos]   = 10.*Pos[0];
     coor.coor[ipos+1] = 10.*Pos[1];
     coor.coor[ipos+2] = 10.*Pos[2];
  }
  cam_list[0].ens_return_shifts(coor, sh);
  if(!serial) comm.Sum(&sh[0][0], numResidues*6);

  int printout=0;
  if(pperiod>0&&comm.Get_rank()==0) printout = (!(getStep()%pperiod));
  if(printout) {
    string csfile;
    char tmps1[21], tmps2[21];
    sprintf(tmps1, "%li", getStep());
    if(ensemble) {
      sprintf(tmps2, "%i", multi_sim_comm.Get_rank());
      csfile = string("cs")+tmps2+"-"+tmps1+string(".dat");
    } else csfile = string("cs")+tmps1+string(".dat");
    cam_list[0].printout_chemical_shifts(csfile.c_str());
  }

  if(ensemble) {
    if(comm.Get_rank()==0) { // I am the master of my replica
      // among replicas
      double fact = 1./((double) ens_dim);
      multi_sim_comm.Sum(&sh[0][0], numResidues*6);
      multi_sim_comm.Barrier(); 
      for(unsigned i=0;i<6;i++) for(unsigned j=0;j<numResidues;j++) sh[j][i] *= fact; 
    } else for(unsigned i=0;i<6;i++) for(unsigned j=0;j<numResidues;j++) sh[j][i] = 0.;
    // inside each replica
    comm.Sum(&sh[0][0], numResidues*6);
  }

  energy = cam_list[0].ens_energy_force(coor, forces, sh);
  if(!serial) comm.Sum(&forces[0][0], N*4);

  for (int i = 0; i < N; i++)
  {
    Vector For;
    int ipos = 4 * i;
    For[0] = forces.coor[ipos];
    For[1] = forces.coor[ipos+1];
    For[2] = forces.coor[ipos+2];
    deriv[i] = 41.86*For;
    if(ensemble) {double fact = 1./((double) ens_dim); deriv[i]*=fact;}
    virial=virial+(-1.*Tensor(getPosition(i),deriv[i]));
  }

  for(unsigned i=0;i<getNumberOfAtoms();++i) setAtomsDerivatives(i,deriv[i]);
  setValue           (4.186*energy);
  setBoxDerivatives  (virial);
}

}
#endif
