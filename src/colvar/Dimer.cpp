// LEGAL CRAP HERE

#include "Colvar.h"
#include "ActionRegister.h"
#include <string>
#include <cmath>
#include <cassert>
#include <iostream>
#include <vector>

using namespace std;

namespace PLMD{
namespace colvar{

//+PLUMEDOC COLVAR DIMER
/*
WRITE DOCUMENTATION + INPUT EXAMPLE HERE
I DONT WANT TO DIE SO I'LL WRITE IT BEFORE PUSHING
*/
//+ENDPLUMEDOC

class Dimer : public Colvar {
	public:
		static void registerKeywords( Keywords& keys);
		Dimer(const ActionOptions&);
		virtual void calculate();
	protected:
		bool trimer,useall;
		int myrank, nranks, natoms;
	private:
		void consistencyCheck();
		double qexp,temperature,beta,dsigma;
		vector<double> dsigmas;
		vector<int> usedatoms;
};

PLUMED_REGISTER_ACTION(Dimer, "DIMER")



void Dimer::registerKeywords( Keywords& keys){
	Colvar::registerKeywords(keys);
	
	keys.add("compulsory","DSIGMA","The interaction strength of the dimer bond.");
	keys.add("compulsory", "Q", "The exponent of the dimer potential.");
	keys.add("compulsory", "TEMP", "The temperature (in Kelvin) of the simulation.");
	keys.add("compulsory", "ATOMS", "The list of atoms being considered by this CV. Used if ALLATOMS flag is missing");
	keys.add("compulsory","NATOMS","The number of dimerized atoms. Used with ATOMS list");
	keys.addFlag("ALLATOMS", false, "Use EVERY atom of the system. Overrides ATOMS keyword.");
	keys.addFlag("VSITES", false, "If present the configuration is with virtual sites at the centroids.");
	
}



Dimer::Dimer(const ActionOptions& ao):
	PLUMED_COLVAR_INIT(ao)
{
	
	log<<" Bibliography here...";
	parseVector("DSIGMA",dsigmas);
	parse("Q",qexp);
	
	parse("TEMP",temperature);
	
	
	vector<AtomNumber> atoms;   
	parseFlag("ALLATOMS",useall);
	trimer=false;
	parseFlag("VSITES",trimer);
	
	nranks=multi_sim_comm.Get_size();
	myrank=multi_sim_comm.Get_rank();
	if(dsigmas.size()==1)
		dsigma=dsigmas[0];
	else
		dsigma=dsigmas[myrank];
	
	
	
	
	if(useall)
	{
		// go with every atom in the system but not the virtuals...
		if(trimer)
			natoms= 2*getTotAtoms()/3;
		else
			natoms=getTotAtoms()/2;	

		for(unsigned int i=0;i<((unsigned int)natoms);i++)
		{
			AtomNumber ati;
			ati.setIndex(i);
			atoms.push_back(ati);
		}
	}
	else  // comment on this  
	{
		parseVector("ATOMS",usedatoms);
		double ntm;
		parse("NATOMS",ntm);
		natoms=ntm;
		unsigned int isz = usedatoms.size();
		for(unsigned int i=0;i<isz;i++)
		{
			AtomNumber ati;
			ati.setIndex(usedatoms[i]-1);
			atoms.push_back(ati);
			
		}
		
		for(unsigned int i=0;i<isz;i++)
		{
			AtomNumber atip2;
			atip2.setIndex(usedatoms[i]+natoms-1);
			if(usedatoms[i]>natoms)
				error("The Dimer CV requires that when choosing atoms you refere only to the first beads.");	
			atoms.push_back(atip2);
		}
		
	}
	consistencyCheck();
	checkRead(); 
	beta = 1./(kBoltzmann*temperature);
	
	addValueWithDerivatives();  // allocate
	
	requestAtoms(atoms);
	
	setNotPeriodic();
	
}

void Dimer::calculate()
{
	double cv_val=0;
	vector<Vector> derivatives;
	vector<Vector> my_pos=getPositions();
	int atms = my_pos.size();
	vector<Vector> der_b2;
	for(int i=0;i<atms/2;i++)
	{
		Vector dist;
		dist = pbcDistance(my_pos[i],my_pos[i+atms/2]);
		double distquad=0;
		for(int j=0;j<3;j++)
			distquad += dist[j]*dist[j];
		
		double dsigquad = dsigma*dsigma;
		double fac1 = 1.0 + distquad/(2*qexp*dsigquad);
		double fac1qm1 = pow(fac1,qexp-1);
		
		
		cv_val += (fac1*fac1qm1-1.0)/beta;
		Vector der_val;
		Vector mder_val;
		for(int j=0;j<3;j++)
		{
			der_val[j] = -fac1qm1*dist[j]/(dsigquad*beta);
			mder_val[j]=-der_val[j];
		}
		derivatives.push_back(der_val);
		der_b2.push_back(mder_val);
	}
	
	derivatives.insert(derivatives.end(), der_b2.begin(), der_b2.end());
	
	for(unsigned int i=0;i<derivatives.size();i++)
		setAtomsDerivatives(i,derivatives[i]);
	
	// Don't use variable box dimensions:
	// length of the dimer in pbc is ill posed if both replicas are doing 
	// different NPT simulations and I don't see any easy way of syncing them.
	// Well... if  we manage to make a dimer in a single replica probably we 
	// will also have this.
	Tensor virial; 
	plumed_dbg_assert( !mypack.virialWasSet() );
  	setBoxDerivativesNoPbc();
	setValue(cv_val);
	
}



/*****************
There are some conditions that a valid input should satisfy.
These are checked here and PLUMED error handlers are (eventually) called.
******************/
void Dimer::consistencyCheck()
{
	if(useall==false && natoms==0)
		error("With ATOMS also NATOMS is required to specify the number of dimerized atoms.");
		
	if(qexp<0.5 || qexp>1)
		warning("Dimer CV is meant to be used with q-exponents between 0.5 and 1. We are not responsible for any black hole. :-)");
	if(dsigma<0)
		error("Please use positive sigma values for the Dimer strength constant");	
	if(temperature<0)
		error("Please, use a positive value for the temperature...");
	
	// if dsigmas has only one element means that you are using different plumed.x.dat files
	if(dsigmas.size()!=nranks && dsigmas.size()!=1) 
		error("Mismatch between provided sigmas and number of replicas");

}


}}

