#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include "plumed/tools/Exception.h"
#include "plumed/tools/RMSD.h"
#include "plumed/tools/Stopwatch.h"
#include "plumed/tools/PDB.h"
#include "plumed/tools/Log.h"
#include "plumed/tools/Vector.h"
#include "plumed/tools/Matrix.h"
#include <fstream>
#include "plumed/tools/AtomNumber.h"
#include "plumed/tools/RMSD.h"


using namespace std ;
using namespace PLMD;

int main(int argc, char* argv[]) {

  // findiff step: eps
  double eps=1.e-6;

  // various setups: com needs to be reset
  bool remove_com=true;
  // normalize weights or not (to check if the rmsd is proportional to weights) 
  bool normalize_weights=true;
  // use msd instead of rmsd to check consistency
  bool squared=false;
  // enhance com to enphasize the issues with COM treatment
  bool enhance_com=false;

  // parse input instructions
  // default task, calculate the RMSD of one frame respect to a set of others 
  vector<int> task(1);task[0]=0;

  // this test wants to be only for OPTIMAL case
  string type; type.assign("OPTIMAL"); 

  // first parse the task: in this applications the tasks are simple integers that pass the action that is required
  for(int i = 1; i < argc; i++){ 
      task.push_back(atoi(argv[i]));
  }  
  if(std::find(task.begin(), task.end(), -1)!=task.end()){cout<<"squared=true (default false)"<<endl;squared=true;}
  if(std::find(task.begin(), task.end(), -2)!=task.end()){cout<<"normalize_weights=false (default true)"<<endl;normalize_weights=false;}
  if(std::find(task.begin(), task.end(), -3)!=task.end()){cout<<"remove_com=false (default true)"<<endl; remove_com=false;}
  if(std::find(task.begin(), task.end(), -4)!=task.end()){cout<<"OPTIMAL-FAST (default OPTIMAL)"<<endl; type.assign("OPTIMAL-FAST");}
  if(std::find(task.begin(), task.end(), -5)!=task.end()){cout<<"enhance_com=true (default false) include option -3 (no com removal) "<<endl;enhance_com=true ; }
  if(enhance_com)remove_com=false;


  cout<<"ARGUMENTS: \n";
  cout<<"OPTIONS that go on top of tasks:\n";
  cout<<" -1 : squared=true (default=false)\n";  
  cout<<" -2 : normalize_weights=false (default=true)\n";
  cout<<" -3 : remove_com=false (default=true) \n";
  cout<<" -4 : OPTIMAL-FAST (default=OPTIMAL) \n";
  cout<<" -5 : enhance_com=true (default false) automatically set option -3 (i.e. check the trouble with com when you do not remove it)\n";
  cout<<"TASKS (can choose more than one):\n";
  cout<<"  0 : normal rmsd/msd calculation  and derivative dumps (default: always done)\n";
  cout<<"  1 : findiff test for  d msd / d position  (inhomogenehous weights)\n";
  cout<<"  2 : findiff test for  d msd / d reference (inhomogenehous weights)\n";
  cout<<"  3 : findiff test for  d msd / d position  (homogenehous weights)\n";
  cout<<"  4 : findiff test for  d msd / d reference (homogenehous weights)\n";
  cout<<"  5 : findiff test for  d Rot / d position  (inhomogenehous weights)  (reference->position)\n";
  cout<<"  6 : findiff test for  d Rot / d reference  (inhomogenehous weights) (reference->position)\n";
  cout<<"  7 : consistency check for MSD proportionality (works with squared=true through option -1 )\n";
  cout<<"  8 : do some timings for all the above routines and for a growing number of atoms\n";
  cout<<"  9 : test the rotation order: print position.pdb reference.pdb aligned.pdb and check that it makes sense(should be reference aligned onto positions)\n";
  cout<<" 10 : findiff test for  d Rot / d position  ( position -> reference ) \n";
  cout<<" 11 : findiff test for  d Rot / d position  (homogenehous weights)  (reference->position)\n";
  cout<<" 12 : findiff test for  d Rot / d reference  (homogenehous weights) (reference->position)\n";
  cout<<" 13 : do timings only for the most common use (aligment +derivatives) and repeat for homogeneous weights\n";
  cout<<" 20 : a simple test calculating only derivative of positions which is needed to compare new version and old version (need to compile with -DOLDRMSD both plumed and the test)\n";

  PDB pdbref;

  // now create the object: does not do anything but set the typer to SIMPLE
  PLMD::RMSD* rmsd=new RMSD();

  // set the reference pdb 
  string reference; reference.assign("1GB1_mdl1_rototranslated.pdb");
  PDB pdb;
  if( !pdb.read(reference,false,0.1) ) 
        cout<<"missing input file 1GB1_mdl1_rototranslated.pdb "<<"\n";
  if (enhance_com){
	vector<Vector> v=pdb.getPositions();
	for(unsigned i=0;i<v.size();i++){v[i][0]+=10.;v[i][1]+=20.;v[i][2]+=30.;}
	pdb.setPositions(v);
  }

  
  cout<<"NOW CREATING THE RMSD OBJECT... with set() method";  
  // remember that "set" method parses the reference, align and displace from the PDB by calling the following methods 
  //      setReference(pdb.getPositions());  -> set the reference and put the weights as 1/n, remove the com according to such weight 
  //      setAlign(pdb.getOccupancy());    -> normalizes, set the alignment vector and remove the Com using such weights 
  //      setDisplace(pdb.getBeta());    -> normalizes and set the displacement vector 
  //      setType(mytype);
  rmsd->set(pdb,type,remove_com,normalize_weights);
  // store other vectors from further manipulation
  std::vector<Vector> ref ;  ref=pdb.getPositions() ;
  std::vector<double> align; align=pdb.getOccupancy(); // non-normalized !
  std::vector<double> displace;  displace=pdb.getBeta(); // non-normalized ! 

  cout<<"DONE!"<<endl;  

  // now take another conformation to compare with: the running frame
  PDB pdbrun; 
  // mimic gromacs: do it in nm
  if( !pdbrun.read("1GB1_mdl2.pdb",false,0.1) )
        cout<<"missing input file 1GB1_mdl2.pdb\n" ;
  std::vector<Vector> run ;  run=pdbrun.getPositions() ;
  std::vector<Vector> derivatives ; 
  if (enhance_com){
        for(unsigned i=0;i<run.size();i++){run[i][0]-=10.;run[i][1]-=20.;run[i][2]-=30.;}
  }




// Task 0: calculate the alignment and dump some data
  if(std::find(task.begin(), task.end(), 0)!=task.end()){

      cout<<"Task 0: calculates the alignment and retrieve some data"<<endl;
      double r=rmsd->calculate( run, derivatives, squared );
      cout<<"RMSD IS "<<r<<"\n";
      // now dump some more information
      ofstream myfile;
      myfile.open ("output_rmsd");
      myfile<<"RMSD ";
      myfile.setf( std::ios::fixed, std:: ios::floatfield );
      myfile.width(12);myfile.precision(6); myfile<<std::right<<r<<"\n";
      // the derivatives
      for(unsigned int i=0;i<run.size();++i){
       myfile<<"DDIST_DPOS "<<std::setw(5)<<i<<" X= "<<std::setw(18)<<std::right<<derivatives[i][0]<<" Y= "<<std::setw(18)<<std::right<<derivatives[i][1]<<" Z= "<<std::setw(18)<<std::right<<derivatives[i][2]<<"\n";
      }

      std::vector<Vector> DDistDRef; DDistDRef.resize(align.size());
      r=rmsd->calc_DDistDRef( run,derivatives,DDistDRef,squared);
      for(unsigned int i=0;i<run.size();++i){
       myfile<<"DDIST_DREF "<<std::setw(5)<<i<<" X= "<<std::setw(18)<<std::right<<DDistDRef[i][0]<<" Y= "<<std::setw(18)<<std::right<<DDistDRef[i][1]<<" Z= "<<std::setw(18)<<std::right<<DDistDRef[i][2]<<"\n";
      }
      myfile.close();
  }

  // Task 1: calculate findiff of running frame
  if(std::find(task.begin(), task.end(), 1)!=task.end()){
	cout<<"Task 1: calculates the finite difference for the running frame"<<endl;
  	double r_old=rmsd->calculate( run, derivatives, squared ); 
	std::vector<Vector> run_save=run;	
	for(unsigned int comp=0;comp<3;comp++){	
		for(unsigned int i=0;i<run.size();++i){
			//change the position		
			run[i][comp]+=eps;
  			double r=rmsd->calculate( run, derivatives, squared ); 
			cout<<"DDIST_DPOS COMPONENT "<<comp<<" "<<(r-r_old)/(run[i][comp]-run_save[i][comp])<<" "<<derivatives[i][comp]<<"\n";
			// restore the old position
			run=run_save;
		}
	}
  }

#ifndef OLDRMSD
  // Task 2: calculate findiff of reference frame
  if(std::find(task.begin(), task.end(), 2)!=task.end()){
	cout<<"Task 2: calculates the finite difference for the reference frame"<<endl;
  	double r_old=rmsd->calculate( run, derivatives, squared ); 
	std::vector<Vector> ref_save=ref;	
	std::vector<Vector> DDistDRef;	
        rmsd->clear();
        rmsd->set(align, displace, ref,type,remove_com,normalize_weights  );	
	for(unsigned int comp=0;comp<3;comp++){	
		for(unsigned int i=0;i<run.size();++i){
			//change the position		
			ref[i][comp]+=eps;
			// the asymmetry between reference and positions requires that 
			// in order to modify the reference one has to reset the rmsd object  
		        rmsd->clear();
			rmsd->set(align, displace, ref,type,remove_com,normalize_weights  );
  			double r=rmsd->calc_DDistDRef( run, derivatives, DDistDRef, squared); 
			cout<<"DDIST_DREF COMPONENT "<<comp<<" "<<(r-r_old)/(ref[i][comp]-ref_save[i][comp])<<" "<<DDistDRef[i][comp]<<"\n";
			// restore the old position
			ref=ref_save;
		}
	}
  }
 // Task 3 calculate findiff of running frame for alEqDis  version (align=displace)
  if(std::find(task.begin(), task.end(), 3)!=task.end()){
	std::vector<double> newalign(run.size(),1.); 
	std::vector<double> newdisplace(run.size(),1.); 
	rmsd->clear();
	rmsd->set(newalign, newdisplace, ref,type,remove_com,normalize_weights  );
	cout<<"Task 3: calculates the finite difference for the running frame (fast version)"<<endl;
  	double r_old=rmsd->calculate( run, derivatives, squared ); 
	std::vector<Vector> run_save=run;	
	for(unsigned int comp=0;comp<3;comp++){	
		for(unsigned int i=0;i<run.size();++i){
			//change the position		
			run[i][comp]+=eps;
  			double r=rmsd->calculate( run, derivatives, squared ); 
			cout<<"DDIST_DPOS_FAST COMPONENT "<<comp<<" "<<(r-r_old)/(run[i][comp]-run_save[i][comp])<<" "<<derivatives[i][comp]<<"\n";
			// restore the old position
			run=run_save;
		}
	}
  }

  // Task 4: calculate findiff of reference frame for alEqDis version (align=displace)
  if(std::find(task.begin(), task.end(), 4)!=task.end()){
	cout<<"Task 4: calculates the finite difference for the reference frame"<<endl;
	std::vector<double> newalign(run.size(),1.); 
	std::vector<double> newdisplace(run.size(),1.); 
	rmsd->clear();
	rmsd->set(newalign, newdisplace, ref,type,remove_com,normalize_weights);
  	double r_old=rmsd->calculate( run, derivatives, squared ); 
	std::vector<Vector> DDistDRef;	
	std::vector<Vector> ref_save=ref;	
	for(unsigned int comp=0;comp<3;comp++){	
		for(unsigned int i=0;i<run.size();++i){
			//change the position		
			ref[i][comp]+=eps;
			// this function below also reset the com of the reference (but not of the running frame)
			rmsd->clear();
			rmsd->set(newalign, newdisplace, ref,type,remove_com,normalize_weights  );
  			double r=rmsd->calc_DDistDRef( run, derivatives,DDistDRef, squared); 
			cout<<"DDIST_DREF_FAST COMPONENT "<<comp<<" "<<(r-r_old)/(ref[i][comp]-ref_save[i][comp])<<" "<<DDistDRef[i][comp]<<" "<<r<<" "<<r_old<<"\n";
			// restore the old position
			ref=ref_save;
		}
	}
  }
  
  // Task 5: calculate findiff of derivative of the rotation matrix respect to running frame
  if(std::find(task.begin(), task.end(), 5)!=task.end()){
	cout<<"Task 5: calculates the finite difference for derivative of the rotation matrix respect to the the running frame"<<endl;
        rmsd->clear();
        rmsd->set(align, displace, ref,type,remove_com,normalize_weights  );	

	Tensor Rotation,OldRotation;
	Matrix<std::vector<Vector> > DRotDPos(3,3);
        std::vector<Vector> DDistDRef;
	rmsd->calc_DDistDRef_Rot_DRotDPos( run, derivatives, DDistDRef, OldRotation, DRotDPos, squared ); 
	std::vector<Vector> run_save=run;	
	for(unsigned int a=0;a<3;a++){	
		for(unsigned int b=0;b<3;b++){	
			for(unsigned int comp=0;comp<3;comp++){	
				for(unsigned int i=0;i<run.size();++i){
					//change the position		
					run[i][comp]+=eps;
					rmsd->calc_DDistDRef_Rot_DRotDPos( run, derivatives ,DDistDRef, Rotation, DRotDPos, squared ); 
					cout<<"DROT_DPOS COMPONENT "<<comp<<" "<<(Rotation[a][b]-OldRotation[a][b])/(run[i][comp]-run_save[i][comp])<<" "<<DRotDPos[a][b][i][comp]<<"\n";
					// restore the old position
					run=run_save;
				}
			}
		}
	}
  }
  // Task 6: calculate findiff of derivative of the rotation matrix respect to reference frame 
  if(std::find(task.begin(), task.end(), 6)!=task.end()){
	cout<<"Task 6: calculates the finite difference for derivative of the rotation matrix respect to the the reference frame"<<endl;
        rmsd->clear();
        rmsd->set(align, displace, ref,type,remove_com,normalize_weights  );	
	Tensor Rotation,OldRotation;
	Matrix<std::vector<Vector> > DRotDPos(3,3),DRotDRef(3,3);
        std::vector<Vector> DDistDRef;
	rmsd->calc_DDistDRef_Rot_DRotDPos_DRotDRef( run, derivatives,  DDistDRef, OldRotation , DRotDPos , DRotDRef , squared); 
	std::vector<Vector> ref_save=ref;	
	for(unsigned int a=0;a<3;a++){	
		for(unsigned int b=0;b<3;b++){	
			for(unsigned int comp=0;comp<3;comp++){	
				for(unsigned int i=0;i<run.size();++i){
					//change the position		
					ref[i][comp]+=eps;
					// this function below also reset the com of the reference (but not of the running frame)
		                        rmsd->clear();
       			                rmsd->set(align, displace, ref,type,remove_com,normalize_weights  );
					rmsd->calc_DDistDRef_Rot_DRotDPos_DRotDRef( run, derivatives, DDistDRef, Rotation , DRotDPos, DRotDRef ,squared ); 
					cout<<"DROT_DREF COMPONENT "<<comp<<" "<<(Rotation[a][b]-OldRotation[a][b])/(ref[i][comp]-ref_save[i][comp])<<" "<<DRotDRef[a][b][i][comp]<<"\n";
					// restore the old position
					ref=ref_save;
				}
			}
		}
	}
  }
  // Task 7:  check weight consistency 

  if(std::find(task.begin(), task.end(), 7)!=task.end()){
	cout<<"Task 7: calculates the weight (displacement) consistency: all these should same result when weights are normalized in input by setReferenceAtoms otherwise they should be proportional when squared=true \n When squared=false, each factor of 2 in weights should produce a factor of sqrt(2) in the total value  "<<endl;
	rmsd->clear();
        rmsd->set(align, displace, ref,type,remove_com, false   );
  	double r=rmsd->calculate( run, derivatives, squared ); 
	cout<<"STANDARD WEIGHT "<<r<<"\n"; 

        std::vector<double> newdisplace=displace;for(std::vector<double>::iterator p=newdisplace.begin();p!=newdisplace.end();++p ){(*p)*=2.;} 
	rmsd->clear();
        rmsd->set(align, newdisplace, ref,type,remove_com, false   );
  	r=rmsd->calculate( run, derivatives,  squared ); 
	cout<<"DOUBLE WEIGHT "<<r<<"\n"; 

        newdisplace=displace;for(std::vector<double>::iterator p=newdisplace.begin();p!=newdisplace.end();++p ){(*p)*=4.;} 
	rmsd->clear();
        rmsd->set(align,newdisplace, ref,type,remove_com, false );
  	r=rmsd->calculate( run, derivatives, squared ); 
	cout<<"FOUR WEIGHT "<<r<<"\n"; 
  }

  // Task 8: do some timings to get a flavor
  if(std::find(task.begin(), task.end(), 8)!=task.end()){
      cout<<"Task 8: makes some timings for increasing atoms and different routines "<<endl;
      rmsd->clear();
      rmsd->set(align, displace, ref,type,remove_com,normalize_weights  );	
      vector<Vector> r_run,r_ref;
      vector<double> r_al,r_disp;
      for (unsigned int i=0;i<10;i++){r_run.push_back(run[i]);r_ref.push_back(ref[i]);r_al.push_back(align[i]);r_disp.push_back(displace[i]);}

      for(unsigned int i=0;i<10;i++){
      	cout<<"NUMBER OF ATOMS : "<<r_run.size()<<endl;
      	unsigned ntest; ntest=100;
      	// test the fast routine
	rmsd->clear();
        rmsd->set(r_al,r_disp, r_ref, type,remove_com, normalize_weights );
 
      	Stopwatch sw;
      	sw.start();	
        for(unsigned int j=0;j<ntest;j++)rmsd->calculate( r_run, derivatives, squared );
      	sw.stop();	
      	cout<<"SIMPLE ROUTINE \n"<<sw<<endl;

        std::vector<Vector> DDistDRef;
      	Stopwatch sw2;
      	sw2.start();	
        for(unsigned int j=0;j<ntest;j++)rmsd->calc_DDistDRef( r_run, derivatives,DDistDRef, squared); 
      	sw2.stop();	
      	cout<<"WITH REFERENCE FRAME: \n"<<sw2<<endl;

      	Tensor Rotation;
      	Matrix<std::vector<Vector> > DRotDPos(3,3);	
      	Stopwatch sw3;
      	sw3.start();	
        for(unsigned int j=0;j<ntest;j++)rmsd->calc_DDistDRef_Rot_DRotDPos( r_run, derivatives,DDistDRef, Rotation , DRotDPos, squared); 
      	sw3.stop();	
      	cout<<"WITH ROTATION MATRIX DERIVATIVE: \n"<<sw3<<endl;

        Matrix<std::vector<Vector> > DRotDRef(3,3);
      	Stopwatch sw4;
      	sw4.start();	
        for(unsigned int j=0;j<ntest;j++)rmsd->calc_DDistDRef_Rot_DRotDPos_DRotDRef( r_run, derivatives ,DDistDRef, Rotation , DRotDPos, DRotDRef, squared); 
      	sw4.stop();	
      	cout<<"WITH ROTATION MATRIX DERIVATIVE OF REEFERENCE: \n"<<sw4<<endl;
      	// duplicate the atoms
      	unsigned s=r_run.size();
      	for (unsigned int i=0;i<s;i++){r_run.push_back(r_run[i]);r_ref.push_back(r_ref[i]);r_al.push_back(r_al[i]);r_disp.push_back(r_disp[i]);}
      
      }
  } 
  // Task 9: check the rotation
  if(std::find(task.begin(), task.end(), 9)!=task.end()){
      cout<<"Task 9: dump some pdbs so to check if all makes sense when using inverse transform. In particular positions_aligned.pdb should overlap with reference_centered.pdb "<<endl;
	rmsd->clear();
        rmsd->set(align, displace, ref,type,remove_com,normalize_weights  );	
	// dump the reference
	ofstream myfile;
	myfile.open ("reference.pdb");		
	std::vector<AtomNumber> at=pdb.getAtomNumbers();
	std::vector<Vector>   pos=pdb.getPositions();
	unsigned k=0;
	for(std::vector<AtomNumber>::iterator i=at.begin(); i!=at.end(); i++){
		myfile<<"ATOM";
                myfile.width(7);myfile<<std::right<<(*i).serial()<<" "; 
                myfile.width(4);myfile<<std::left<<pdb.getAtomName(*i); 
                myfile.width(4);myfile<<std::right<<pdb.getResidueName(*i)<<" A"; 
                myfile.width(4);myfile<<std::right<<pdb.getResidueNumber(*i)<<"    "; 
		myfile.setf( std::ios::fixed, std:: ios::floatfield );
                myfile.width(8);myfile.precision(3); myfile<<std::right<<pos[k][0]*10; 
                myfile.width(8);myfile.precision(3); myfile<<std::right<<pos[k][1]*10; 
                myfile.width(8);myfile.precision(3); myfile<<std::right<<pos[k][2]*10<<"  1.00  1.00\n"; 
		k++;	
	}	
	myfile.close();			
	// dump the position
	myfile.open ("positions.pdb");		
	at=pdbrun.getAtomNumbers();
	std::vector<Vector>   runpos=pdbrun.getPositions();
	k=0;
	for(std::vector<AtomNumber>::iterator i=at.begin(); i!=at.end(); i++){
		myfile<<"ATOM";
                myfile.width(7);myfile<<std::right<<(*i).serial()<<" "; 
                myfile.width(4);myfile<<std::left<<pdbrun.getAtomName(*i); 
                myfile.width(4);myfile<<std::right<<pdbrun.getResidueName(*i)<<" A"; 
                myfile.width(4);myfile<<std::right<<pdbrun.getResidueNumber(*i)<<"    "; 
		myfile.setf( std::ios::fixed, std:: ios::floatfield );
                myfile.width(8);myfile.precision(3); myfile<<std::right<<runpos[k][0]*10; 
                myfile.width(8);myfile.precision(3); myfile<<std::right<<runpos[k][1]*10; 
                myfile.width(8);myfile.precision(3); myfile<<std::right<<runpos[k][2]*10<<"  1.00  1.00\n"; 
		k++;	
	}	
	myfile.close();			
	// now do the alignment
	Tensor Rotation;
	Matrix<std::vector<Vector> > DRotDPos(3,3);
        std::vector<Vector> DDistDRef;
	std::vector<Vector> alignedpos;
	std::vector<Vector> centeredpos;
	std::vector<Vector> centeredref;
	std::vector<Vector> ddistdpos;
	rmsd->calc_PCAelements( run, derivatives, Rotation ,  DRotDPos , alignedpos ,centeredpos, centeredref ,squared); 
	myfile.open ("positions_aligned.pdb");		
	k=0;
	for(std::vector<AtomNumber>::iterator i=at.begin(); i!=at.end(); i++){
		myfile<<"ATOM";
                myfile.width(7);myfile<<std::right<<(*i).serial()<<" "; 
                myfile.width(4);myfile<<std::left<<pdbrun.getAtomName(*i); 
                myfile.width(4);myfile<<std::right<<pdbrun.getResidueName(*i)<<" A"; 
                myfile.width(4);myfile<<std::right<<pdbrun.getResidueNumber(*i)<<"    "; 
		myfile.setf( std::ios::fixed, std:: ios::floatfield );
                myfile.width(8);myfile.precision(3); myfile<<std::right<<alignedpos[k][0]*10; 
                myfile.width(8);myfile.precision(3); myfile<<std::right<<alignedpos[k][1]*10; 
                myfile.width(8);myfile.precision(3); myfile<<std::right<<alignedpos[k][2]*10<<"  1.00  1.00\n"; 
		k++;	
	}	
	myfile.close();			
	// dump the aligned	
	myfile.open ("reference_centered.pdb");		
	k=0;
	for(std::vector<AtomNumber>::iterator i=at.begin(); i!=at.end(); i++){
		myfile<<"ATOM";
                myfile.width(7);myfile<<std::right<<(*i).serial()<<" "; 
                myfile.width(4);myfile<<std::left<<pdbrun.getAtomName(*i); 
                myfile.width(4);myfile<<std::right<<pdbrun.getResidueName(*i)<<" A"; 
                myfile.width(4);myfile<<std::right<<pdbrun.getResidueNumber(*i)<<"    "; 
		myfile.setf( std::ios::fixed, std:: ios::floatfield );
                myfile.width(8);myfile.precision(3); myfile<<std::right<<centeredref[k][0]*10; 
                myfile.width(8);myfile.precision(3); myfile<<std::right<<centeredref[k][1]*10; 
                myfile.width(8);myfile.precision(3); myfile<<std::right<<centeredref[k][2]*10<<"  1.00  1.00\n"; 
		k++;	
	}	
	myfile.close();			
  }
  // Task 10: derivative of the rotation matrix (in case of reverse transition) 
  if(std::find(task.begin(), task.end(), 10)!=task.end()){
	cout<<"Task 10: calculates the finite difference for derivative of the rotation matrix respect to the the running frame"<<endl;
	rmsd->clear();
        rmsd->set(align, displace, ref,type,remove_com,normalize_weights  );	
	Tensor Rotation,OldRotation;
	Matrix<std::vector<Vector> > DRotDPos(3,3);
        std::vector<Vector> DDistDRef;
        std::vector<Vector> alignedpos;
        std::vector<Vector> centeredpos;
        std::vector<Vector> centeredref;
	std::vector<Vector> ddistdpos;
	rmsd->calc_PCAelements( run, derivatives, OldRotation , DRotDPos , alignedpos ,centeredpos, centeredref, squared ); 
	std::vector<Vector> run_save=run;	
	for(unsigned int a=0;a<3;a++){	
		for(unsigned int b=0;b<3;b++){	
			for(unsigned int comp=0;comp<3;comp++){	
				for(unsigned int i=0;i<run.size();++i){
					//change the position		
					run[i][comp]+=eps;
					rmsd->calc_PCAelements( run, derivatives, Rotation ,DRotDPos , alignedpos ,centeredpos, centeredref , squared); 
					cout<<"DROT_DPOS_INVERSE_TRANSFORM COMPONENT "<<comp<<" "<<(Rotation[a][b]-OldRotation[a][b])/(run[i][comp]-run_save[i][comp])<<" "<<DRotDPos[a][b][i][comp]<<"\n";
					// restore the old position
					run=run_save;
				}
			}
		}
	}
  }
  // Task 11: calculate findiff of derivative of the rotation matrix respect to running frame (homogenehous weights)
  if(std::find(task.begin(), task.end(), 11)!=task.end()){
	cout<<"Task 11: calculates the finite difference for derivative of the rotation matrix respect to the the running frame (homogeneous weights)"<<endl;
	std::vector<double> newalign(run.size(),1.); 
	std::vector<double> newdisplace(run.size(),1.); 
	rmsd->clear();
        rmsd->set(newalign, newdisplace, ref,type,remove_com,normalize_weights  );	
	Tensor Rotation,OldRotation;
	Matrix<std::vector<Vector> > DRotDPos(3,3);
        std::vector<Vector> DDistDRef;
	std::vector<Vector> alignedpos;
	std::vector<Vector> centeredpos;
	std::vector<Vector> centeredref;
	rmsd->calc_DDistDRef_Rot_DRotDPos( run,derivatives, DDistDRef, OldRotation , DRotDPos, squared  ); 
	//rmsd->calc_PCAelements( run, derivatives, OldRotation ,  DRotDPos , alignedpos ,centeredpos, centeredref ,squared); 
	std::vector<Vector> run_save=run;	
	for(unsigned int a=0;a<3;a++){	
		for(unsigned int b=0;b<3;b++){	
			for(unsigned int comp=0;comp<3;comp++){	
				for(unsigned int i=0;i<run.size();++i){
					//change the position		
					run[i][comp]+=eps;
					rmsd->calc_DDistDRef_Rot_DRotDPos( run, derivatives , DDistDRef, Rotation , DRotDPos, squared  ); 
					//rmsd->calc_PCAelements( run, derivatives, Rotation ,  DRotDPos , alignedpos ,centeredpos, centeredref ,squared); 
					cout<<"DROT_DPOS COMPONENT "<<comp<<" "<<(Rotation[a][b]-OldRotation[a][b])/(run[i][comp]-run_save[i][comp])<<" "<<DRotDPos[a][b][i][comp]<<"\n";
					// restore the old position
					run=run_save;
				}
			}
		}
	}
  }
  // Task 12: calculate findiff of derivative of the rotation matrix respect to reference frame 
  if(std::find(task.begin(), task.end(), 12)!=task.end()){
	cout<<"Task 12: calculates the finite difference for derivative of the rotation matrix respect to the the reference frame (homogeneous weights)"<<endl;
	std::vector<double> newalign(run.size(),1.); 
	std::vector<double> newdisplace(run.size(),1.); 
	rmsd->clear();
        rmsd->set(newalign, newdisplace, ref,type,remove_com,normalize_weights  );	
	Tensor Rotation,OldRotation;
	Matrix<std::vector<Vector> > DRotDPos(3,3),DRotDRef(3,3);
        std::vector<Vector> DDistDRef;
	rmsd->calc_DDistDRef_Rot_DRotDPos_DRotDRef( run, derivatives, DDistDRef, OldRotation , DRotDPos , DRotDRef ,squared); 
	std::vector<Vector> ref_save=ref;	
	for(unsigned int a=0;a<3;a++){	
		for(unsigned int b=0;b<3;b++){	
			for(unsigned int comp=0;comp<3;comp++){	
				for(unsigned int i=0;i<run.size();++i){
					//change the position		
					ref[i][comp]+=eps;
					// this function below also reset the com of the reference (but not of the running frame)
					rmsd->clear();
				        rmsd->set(newalign, newdisplace, ref,type,remove_com,normalize_weights  );	
					rmsd->calc_DDistDRef_Rot_DRotDPos_DRotDRef( run, derivatives, DDistDRef, Rotation , DRotDPos, DRotDRef ,squared ); 
					cout<<"DROT_DREF COMPONENT "<<comp<<" "<<(Rotation[a][b]-OldRotation[a][b])/(ref[i][comp]-ref_save[i][comp])<<" "<<DRotDRef[a][b][i][comp]<<"\n";
					// restore the old position
					ref=ref_save;
				}
			}
		}
	}
  }
   // Task 13: do some timings to get a flavor (only on alignment+derivatives)
  if(std::find(task.begin(), task.end(), 13)!=task.end()){
	cout<<"Task 13: timings for the most common use (only on alignment+derivatives)  "<<endl;
	vector<Vector> r_run,r_ref;
	vector<double> r_al,r_disp;
	for (unsigned int i=0;i<10;i++){r_run.push_back(run[i]);r_ref.push_back(ref[i]);r_al.push_back(align[i]);r_disp.push_back(displace[i]);}
	for(unsigned int i=0;i<10;i++){
		cout<<"NUMBER OF ATOMS : "<<r_run.size()<<endl;
		unsigned ntest; ntest=100;
		// test the fast routine
		rmsd->clear();
	        rmsd->set(r_al, r_disp, r_ref, type,remove_com,normalize_weights  );	
		Stopwatch sw;
		sw.start();	
	        for(unsigned int j=0;j<ntest;j++)rmsd->calculate( r_run, derivatives, squared );
		sw.stop();	
                cout<<"TIME \n"<<sw<<endl;
		// duplicate the atoms
		unsigned s=r_run.size();
		for (unsigned int i=0;i<s;i++){r_run.push_back(r_run[i]);r_ref.push_back(r_ref[i]);r_al.push_back(r_al[i]);r_disp.push_back(r_disp[i]);}
	
	}
	cout <<"NOW HOMOGENEOUS WEIGHTS\n";
	std::vector<double> newalign(10,1.); 
	std::vector<double> newdisplace(10,1.); 
	r_run.resize(0);
	r_ref.resize(0);
	r_al.resize(0);
	r_disp.resize(0);
	for (unsigned int i=0;i<10;i++){r_run.push_back(run[i]);r_ref.push_back(ref[i]);r_al.push_back(newalign[i]);r_disp.push_back(newdisplace[i]);}
        rmsd->clear();
        rmsd->set(r_al,r_disp , r_ref ,type,remove_com,normalize_weights  );
	for(unsigned int i=0;i<10;i++){
		cout<<"NUMBER OF ATOMS : "<<r_run.size()<<endl;
		unsigned ntest; ntest=100;
		// test the fast routine
                rmsd->clear();
                rmsd->set(r_al, r_disp, r_ref, type,remove_com,normalize_weights  );
		derivatives.resize(r_al.size());
		Stopwatch sw;
		sw.start();	
	        for(unsigned int j=0;j<ntest;j++)rmsd->calculate( r_run, derivatives, squared );
		sw.stop();	
                cout<<"TIME \n"<<sw<<endl;
		// duplicate the atoms
		unsigned s=r_run.size();
		for (unsigned int i=0;i<s;i++){r_run.push_back(r_run[i]);r_ref.push_back(r_ref[i]);r_al.push_back(r_al[i]);r_disp.push_back(r_disp[i]);}
	
	}
	
	

  } 


#endif
  // Task 20: do some timings to get a flavor
  if(std::find(task.begin(), task.end(), 20)!=task.end()){
      cout<<"Task 8: makes some timings for increasing atoms to compare new and old rmsd (need to be recompiled with -DOLDRMSD) "<<endl;
      vector<Vector> r_run,r_ref;
      vector<double> r_al,r_disp;
      for (unsigned int i=0;i<10;i++){r_run.push_back(run[i]);r_ref.push_back(ref[i]);r_al.push_back(align[i]);r_disp.push_back(displace[i]);}

      for(unsigned int i=0;i<10;i++){
      	cout<<"NUMBER OF ATOMS : "<<r_run.size()<<endl;
      	unsigned ntest; ntest=100;
      	// test the fast routine
	rmsd->clear();
        rmsd->set(r_al,r_disp, r_ref, type,remove_com, normalize_weights );
 
      	Stopwatch sw;
      	sw.start();	
        for(unsigned int j=0;j<ntest;j++)rmsd->calculate( r_run, derivatives, squared );
      	sw.stop();	
      	cout<<"SIMPLE ROUTINE \n"<<sw<<endl;

      	unsigned s=r_run.size();
      	for (unsigned int i=0;i<s;i++){r_run.push_back(r_run[i]);r_ref.push_back(r_ref[i]);r_al.push_back(r_al[i]);r_disp.push_back(r_disp[i]);}
      
      }
  } 


  return 0;
}
