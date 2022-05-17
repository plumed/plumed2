/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2014 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#include "Colvar.h"
#include "ActionRegister.h"

#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

using namespace std;

namespace PLMD{
  namespace colvar{

    //+PLUMEDOC COLVAR MOVSHELL
    /*
      Calculate the solute and solvent concentration in a planar shell of the box
      !!! The derivative contains only the outer boundary terms, the resulting bias force is non-conservative !!! Only for CmuMD

      \par Examples


      \verbatim
      d: DIPOLE GROUP=1-10
      PRINT FILE=output STRIDE=5 ARG=5
      \endverbatim
      (see also \ref PRINT)

      \attention
      If the total charge Q of the group in non zero, then a charge Q/N will be subtracted to every atom,
      where N is the number of atoms. This implies that the dipole (which for a charged system depends
      on the position) is computed on the geometric center of the group.


    */
    //+ENDPLUMEDOC

    class Cmumd : public Colvar {
      //      SwitchingFunction switchingFunction; //instance sw.f
      bool print_intf;
      bool issolute,isdelta;
      bool isnotscaled;
      int  N_st, N_sv, Na_sv_permol, Na_st_permol, Na_sv, Na_st, N_mol, com_sv, com_st, nbin, asymm;
      double  iD_CR, iD_F, iCR_Size, iw_force, iw_in, iw_out, co_out, co_in, co_f, nint, fixi;
      //bool ismoving;
    public:
      Cmumd(const ActionOptions&);
      virtual void calculate();
      static void registerKeywords(Keywords& keys);
      double sigmon(double z, double Coff);
      double sigmoff(double z, double Coff);
      double dsig(double z, double Coff);
      ofstream fdbg,fdbg2;
    };

    PLUMED_REGISTER_ACTION(Cmumd,"CMUMD")

    void Cmumd::registerKeywords(Keywords& keys){
      Colvar::registerKeywords(keys);
      keys.add("atoms","GROUP","the group of atoms involved in the calculation");
      keys.add("compulsory","NSV","Solvent atoms");
      keys.add("optional","SOLUTE","Solute tot atoms");
      keys.add("optional","NST","Solute atoms");
      keys.add("compulsory","DCR","CR distance");
      keys.add("compulsory","CRSIZE","CR size");
      keys.add("optional","DF","Force distance");
      keys.add("compulsory","WF","force sigma length");
      keys.add("optional","COF","force sigma cutoff");
      keys.add("optional","WIN","in sigma length");
      keys.add("optional","COIN","in sigma cutoff");
      keys.add("optional","WOUT","out sigma length");
      keys.add("optional","COOUT","out sigma cutoff");
      keys.add("optional","FIXED","fixed interface");
      keys.add("optional","NZ","Interface localization: zbin");
      keys.add("optional","NINT","Interface localization: int density");
      keys.add("optional","COMST","solute COM");
      keys.add("optional","COMSV","solvent COM");
      keys.add("optional","ASYMM","only left(smaller z) or right (larger z) considered");
      keys.addFlag("INTFILE",false,"print interface file (not available!!!)");
      keys.addFlag("DELTA",false,"concentration gradient");
      keys.addFlag("NOSCALE",false,"use absolute length units");
      keys.remove("NOPBC");
    }

    Cmumd::Cmumd(const ActionOptions&ao):
      PLUMED_COLVAR_INIT(ao),
      //init bool parameters
      print_intf(false),
      isnotscaled(false)
      //isnotmoving(false),
      //serial(false)
    {

      //Read atom group
      vector<AtomNumber> at_list;

      parseAtomList("GROUP",at_list);

      Na_sv_permol=1; //default
      parse("NSV",Na_sv_permol); //get number of atoms per molecule

      N_st=0; //default
      Na_st=0;
      Na_st_permol=1;

      parse("SOLUTE",Na_st); //get number of solute atoms
      parse("NST",Na_st_permol);

      //Solution numbers
      N_st=(int)(Na_st/Na_st_permol); //Number of solute atoms
      Na_sv=at_list.size()-Na_st; //Number of solvent atoms
      N_sv=(int)(Na_sv/Na_sv_permol); //Number of solvent molecules
      N_mol=N_sv+N_st; //Number of total molecules

      log.printf("Number of atoms:\tw %d\tu %d\n",Na_sv,Na_st);
      log.printf("Number of molecules:\ttot %d\t w %d\tu %d\n",N_mol,N_sv,N_st);

      //Parameters (force position and switching function temperature)

      parse("DCR",iD_CR); //CR distance from interface
      parse("CRSIZE",iCR_Size); //CR Size
      iD_F=iD_CR+iCR_Size; //initialize D_F: force distance from interface
      parse("DF",iD_F);
      if(iD_F<iD_CR+iCR_Size){ //re-initialize D_F if inside CR
	iD_F=iD_CR+iCR_Size;
	log.printf("D_F inside CR region, reset at the boundary");
      }
      parse("WF",iw_force); //Fermi Fun T at DF
      co_f=20.0; //initialize cut-off in
      parse("COF",co_f); //cut-off for Fermi f
      iw_in=iw_force; //initialize w_in
      parse("WIN",iw_in); //Fermi Fun T at CRin
      co_in=co_f; //initialize cut-off in
      parse("COIN",co_in); //cut-off for Fermi f
      iw_out=iw_force; //initialize w_out
      parse("WOUT",iw_out); //Fermi Fun T at CRout
      co_out=co_f; //initialize cut-off in
      parse("COOUT",co_out); //cut-off for Fermi f

      log.printf("Geometry:\tD_CR %lf\tCR_size %lf\tD_F %lf\n",iD_CR,iCR_Size,iD_F);
      log.flush();

      fixi=-1.; //default fixed inactive
      parse("FIXED",fixi); //fixed interface coordinate (always scaled)
      if(fixi>=0){
	log.printf("Fixed interface at:\t %lf\n L_box",fixi);
      }

      //Asymmetry

      asymm=0; //default no asymmetry
      parse("ASYMM",asymm); //cut-off for Fermi f
      if(asymm<0){
	log.printf("Only left CR considered");
      }else if(asymm>0){
	log.printf("Only right CR considered");
      }

      parseFlag("DELTA",isdelta);
      if(isdelta){
	log.printf("Difference between right and left CR calculated");
      }

      //COM flags
      com_sv=-1;
      com_st=-1;
      parse("COMSV",com_sv);
      parse("COMST",com_st);

      nint=0.0; //default values
      nbin=100;
      parse("NINT",nint); //interface boundary concentration
      parse("NZ",nbin); //z histogram bins
      if(fixi<0){
	log.printf("Histogram:\tnint %lf\tnbin %d\n",nint,nbin);
      }
      log.flush(); //DBG
      //other bool parameters
      parseFlag("INTFILE",print_intf);
      parseFlag("NOSCALE",isnotscaled);
      //parseFlag("NOMOVE",isnotmoving);
      //parseFlag("SERIAL",serial);
      log.printf("after all parsing\n");       log.flush(); //DBG
      checkRead();
      addValueWithDerivatives();
      setNotPeriodic();

      //log atom lists
      log.printf("  of %d atoms\n",at_list.size());
      for(unsigned int i=0;i<at_list.size();++i){
	//log.printf("  %d", at_list[i].serial());
      }
      log.printf("  \n");
      if(N_st>0){
	log.printf("of which the first %d are solute atoms\n",N_st);
      }
      requestAtoms(at_list);
      log.printf("  \n");
      log.flush(); //DBG

      //open debug file DBG
      
      //fdbg.open("dbg.dat");
      //fdbg.flush(); //DBG
    }

    double Cmumd::sigmon(double z, double Coff){
      double sig;
      if( z < -Coff){
	sig=0.0;
      }else if(z > Coff){
	sig=1.0;
      }else{
	sig=1.0/(exp(-z)+1.0);
      }
      return(sig);
    }


    double Cmumd::sigmoff(double z, double Coff){
      double sig;
      if( z < -Coff){
	sig=1.0;
      }else if(z > Coff){
	sig=0.0;
      }else{
	sig=1.0/(exp(z)+1.0);
      }
      return(sig);
    }

    double Cmumd::dsig(double z, double Coff){
      double dsig;
      if(fabs(z) > Coff){
	dsig=0.0;
      }else{
	dsig=0.5/(1.0+cosh(z));
      }
      return(dsig);
    }



    // calculator
    void Cmumd::calculate()
    {

      double n_CR;
      Tensor virial;

      virial.zero();  //no virial contribution

      //Vector deriv;

      vector<Vector> deriv(getNumberOfAtoms());
      vector<Vector> com_solv(N_sv);
      Vector diff;

      Vector ze;
      ze.zero();
      //init derivatives
      fill(deriv.begin(), deriv.end(), ze);

      //Parallel parameters

      unsigned int stride;
      unsigned int rank;

      stride=comm.Get_size();  //Number of processes
      rank=comm.Get_rank(); //Rank of present process

      //fdbg<<fixed<< "stride\t" <<stride<< "\t rank\t" <<rank<<"\t"<<endl;

      //Solvent position matrix allocation
      vector<Vector> solve_x(Na_sv_permol);
      //Solvent mass array allocation
      vector<double> solve_m(Na_sv_permol);

      //Solvent masses and total mass
      double M_sv=0.0;
      for(int i=0;i<Na_sv_permol;++i){
	solve_m[i]=getMass(Na_st+i); //the first Na_st are skipped
	M_sv += solve_m[i];
      }

      //fdbg<<fixed<< "Na_st\t" <<Na_st<<"\t"<<endl;
      //fdbg<<fixed<< "Masses\t" <<M_sv<<"\t"<<solve_m[0]<<"\t"<<solve_m[1]<<endl;

      //Box dimensions

      double LBC[3];

      for(int i=0;i<3;++i) LBC[i]=getBox()[i][i];
      //LBC[2]=3.74769; //DBG
      //double Vbox=LBC[0]*LBC[1]*LBC[2]; //box volume

      //log.printf("Z size: %lf\n",LBC[2]);

      double D_CR,CR_Size,D_F,w_force,w_in,w_out;

      double fix_int=0;
      if(fixi>=0.0) fix_int=LBC[2]*fixi; //fixed interface

      if(!isnotscaled){ //rescale input distances

	//log.printf("scaled units\n");

	D_CR=LBC[2]*iD_CR;
	CR_Size=LBC[2]*iCR_Size;
	D_F=LBC[2]*iD_F;
	w_force=LBC[2]*iw_force;
	w_in=LBC[2]*iw_in;
	w_out=LBC[2]*iw_out;
	//co_f=LBC[2]*co_f;
	//co_in=LBC[2]*co_in;
	//co_out=LBC[2]*co_out;
      }

      //fdbg<<setprecision(5); //DBG
      //fdbg<<fixed<<"D_CR CR_Size D_F w_force w_in w_out"<<"\t"<<D_CR<<"\t"<<CR_Size<<"\t"<<D_F<<"\t"<<w_force<<"\t"<<w_in<<"\t"<<w_out<<"\t"<<w_force<<endl;

      //rescale the cut-offs
      //co_f=co_f/w_force;
      //co_in=co_in/w_in;
      //co_out=co_out/w_out;
      double Vbox=getBox().determinant();
      double VCR;
      if(asymm==0){
	VCR=2*iCR_Size*Vbox; //CR volume
      }else{
	VCR=iCR_Size*Vbox; //CR volume
      }
      //VCR=3.74769*3.74769*3.74769*iCR_Size; //
      //Histogram settings (for interface localization)

      //histz-array allocation
      vector<int> histz(nbin,0.0);
      int nz=0;
      //bins width (time dependent)
      double dz=LBC[2]/nbin;
      double Vbin=Vbox/nbin; //Bin volume [nm^3]

      //center of mass vector
      for(int i=rank; i<N_sv; i+=stride){
	com_solv[i].zero();
	//center of mass
	if (com_sv<0) {
	  solve_x[0] = getPosition(Na_st+i*Na_sv_permol);
	  ////fdbg<<i<<"\t"<<setprecision(5)<<fixed<< solve_x[0][0]<<"\t" << solve_x[0][1] << "\t" << solve_x[0][2] << endl;
	  for(int j=1;j<Na_sv_permol;++j){
	    solve_x[j] = getPosition(Na_st+i*Na_sv_permol+j);
	    diff = pbcDistance(solve_x[0],solve_x[j]);
	    //fdbg<<i<<"\t"<<setprecision(5)<<fixed<< solve_x[j][0]<<"\t" << solve_x[j][1] << "\t" << solve_x[j][2] << "\t" << diff.modulo() << endl;
	    com_solv[i] += solve_m[j]*diff;
	  }
	  com_solv[i] = com_solv[i] / M_sv + solve_x[0];
	  //fdbg2<<i<<"\t"<<setprecision(5)<<fixed<< com_solv[i][0]<<"\t" << com_solv[i][1] << "\t" << com_solv[i][2] << endl;
	  //fdbg<< com_solv[i][2]<<endl;
	  //impose PBC on com (useless on x and y for now!!!)
	  //for(int k=0; k<3; k++){
	  if(com_solv[i][2]<0) com_solv[i][2]=com_solv[i][2]+LBC[2];
	  if(com_solv[i][2]>=LBC[2]) com_solv[i][2]=com_solv[i][2]-LBC[2];
	  //}
	}else{
	  //no com
	  com_solv[i]=getPosition(Na_st+i*Na_sv_permol+com_sv);
	}
	//fdbg<<i<<"\t"<<setprecision(5)<<fixed<<com_solv[i][0]<<"\t"<<com_solv[i][1]<<"\t"<<com_solv[i][2]<<endl;
	if(fixi<0){
	  nz=(int)(com_solv[i][2]/dz); //fill histogram
	  histz[nz]+=1;
	}
      }
      /*for(int i=0; i<nbin; ++i){
            //fdbg<<histz[i]<<"\t";
	    }
	    //fdbg<<endl;*/
      //communicate
      comm.Sum(histz);
      comm.Sum(com_solv);

      //log.printf("solvent positions acquired \n");              log.flush(); //DBG
      //Get the liquid-crystal interfaces
      double halfbin, ileft, iright, zleft, zright;
      halfbin=(int)(LBC[2]/(2*dz));
      int p=0;
      int pmone=0;

      //interface finder
      if(fixi<0){
	//find the crystal if it's not at the half, it finds the crystal before halfbin exceeds the limits
	//3 adjacent bins with water concentration < than nint/3
	while((histz[halfbin]+histz[halfbin+1]+histz[halfbin-1]) > nint*Vbin){
	  p++;
	  pmone=2*(p%2)-1;
	  halfbin=halfbin+p*pmone; //Move through the bins
	}

	//put halfbin inside the crystal volume (3 bins, WARNING parameter dependent)

	/*if(j!=0){
	  halfbin=halfbin+3*pmone;
	  if(halfbin<0 || halfbin>=nbin) halfbin=halfbin-nbin*pmone; //set pbc on halfbin
	  }*/

	ileft=halfbin;
	while(histz[ileft] < nint*Vbin){
	  ileft=ileft-1;
	  if(ileft<0) ileft=ileft+nbin; //pbc on left
	}

	iright=ileft+10; //WARNING parameter dependent
	if(iright>=nbin) iright=iright-nbin; //pbc on right
	while(histz[iright]< nint*Vbin){
	  iright=iright+1;
	  if(iright>=nbin) iright=iright-nbin; //pbc on right
	}

	zleft=dz*(ileft+1); //left interface coordinate
	zright=dz*(iright); //right interface coordinate
      }else{
	zleft=fix_int;
	zright=fix_int;
      }
      //Fermi function parameters

      double ZCRrin, ZCRrout, ZCRlin, ZCRlout, ZFright, ZFleft;
      //log.printf("interface positions:\t%lf\t%lf\n",zleft,zright);            log.flush(); //DBG
      //log.printf("scaled parameters:\t%lf\t%lf\t%lf\n",D_CR,CR_Size,D_F);          log.flush(); //DBG
      ZCRlin=zleft-D_CR;
      ZCRlout=zleft-D_CR-CR_Size;
      ZFleft=zleft-D_F;

      ZCRrin=zright+D_CR;
      ZCRrout=zright+D_CR+CR_Size;
      ZFright=zright+D_F;

      //log.printf("CR and Force positions:\nLeft\t%lf\t%lf\t%lf\n",ZCRlout,ZCRlin,ZFleft);
      //log.printf("Right\t%lf\t%lf\t%lf\n",ZCRrin,ZCRrout,ZFright);
      //Evaluate concentration and derivatives
      //if isolute is true C counts the solute molecules, else the solvent ones
      //initializing


      //fdbg<<setprecision(5); //DBG
      //fdbg<<scientific;  //DBG
      //fdbg << fixed<< ZFleft << "\t" << ZCRlout << "\t" << ZCRlin << "\t"  << zleft << "\t"  << zright << "\t" << ZCRrin << "\t" << ZCRrout << "\t" << ZFright <<endl; //DBG
      //fdbg.flush(); //DBG
      ////fdbg << "gridsize: " << gridsize << endl;
      ////fdbg << "listsize: " << listsize << endl;


      n_CR=0.0;
      //deriv.zero();
      double zin,zout,n_lx,n_rx,n_x,zl,zr,dfunc,dl,dr;
      int k;
      if(N_st == 0){ //if solvent specie is restrained
	Vector dleft, dright;
	for(int i=rank; i<N_sv; i+=stride){
	  //for(int i=0; i<N_sv; ++i){
	  //Fermi-like weighting
	  dleft.zero();
	  dright.zero();
	  dfunc=0;
	  dl=0;
	  dr=0;
	  n_lx=0;
	  n_rx=0;
	  //left-side sigma
	  if(asymm<=0){
	    zin=(com_solv[i][2]-ZCRlin)/w_in;
	    zout=(com_solv[i][2]-ZCRlout)/w_out;
	    //with periodic image, sigma on at zout, off at zin
	    n_lx=sigmon(zout,co_out)*sigmoff(zin,co_in)+sigmon(zout-LBC[2]/w_out,co_out)*sigmoff(zin-LBC[2]/w_in,co_in);

	    //Derivatives (only outer boundary derivatives!!!)
	    dleft[2]=com_solv[i][2]-ZFleft;
	    zl=dleft[2]/w_force;
	    dl=(dsig(zl,co_f)+dsig(zl-LBC[2]/w_force,co_f))/w_force;

	  }
	  //right-side sigma
	  if(asymm>=0){
	    zin=(com_solv[i][2]-ZCRrin)/w_in;
	    zout=(com_solv[i][2]-ZCRrout)/w_out;
	    //fdbg<<setprecision(5); //DBG
	    //fdbg<<scientific;  //DBG
	    //with periodic image, sigma on at zin, off at zout
	    n_rx=sigmon(zin,co_in)*sigmoff(zout,co_out)+sigmon(zin+LBC[2]/w_in,co_in)*sigmoff(zout+LBC[2]/w_out,co_out);
	    //fdbg << com_solv[i][2] << "\t" << zin*w_in << "\t" << zout*w_out << "\t" << n_rx << "\t" << sigmon(zin,co_in)*sigmoff(zout,co_out) ; //DBG
	    dright[2]=com_solv[i][2]-ZFright;
	    zr=dright[2]/w_force;
	    dr=(-dsig(zr,co_f)-dsig(zr+LBC[2]/w_force,co_f))/w_force;
	  }

	  if(isdelta){
	    n_x=n_rx-n_lx;
	    //sum the two densities (for ASYMMETRIC, change here!!!)
	    dfunc=dr-dl;
	  }else{
	    //fdbg << "\t"  << n_rx << "\t" << n_lx << endl; //DBG
	    n_x=n_rx+n_lx;


	    dfunc=dr+dl;
	  }

	  //update CV (for now this is the number of molcules)
	  n_CR+=n_x;

	  //fdbg<<setprecision(10); //DBG
	  //fdbg<<i<<"\t"<<com_solv[i][2]<<"\t"<<zl<<"\t"<<zr<<"\t"<<dfunc<<"\t"; //DBG
	  //fdbg<<dfunc<<endl; //DBG
	  if(com_sv<0){ //com coordinates
	    for(int l=0; l<Na_sv_permol; ++l){
	      k=Na_st+i*Na_sv_permol+l; //atom counter
	      deriv[k][2]=getMass(k)/M_sv*(dfunc/VCR); //com affects the derivatives
	      //setAtomsDerivatives(k, deriv );
	      /*
	      //fdbg<<setprecision(8); //DBG
	      //fdbg << fixed << k << "\t" << deriv[0] << endl; //DBG
	      //fdbg << fixed << k << "\t" << deriv[1] << endl; //DBG
	      //fdbg << fixed << k << "\t" << deriv[2] << endl; //DBG
	      //fdbg.flush(); //DBG*/
	      virial-=transpose(Tensor(deriv[k],matmul(getBox(),getPbc().realToScaled(getPosition(k))))); //Virial component;
	    }
	  }else{//single atom coordinates
	    k=Na_st+i*Na_sv_permol+com_sv ; //atom counter (just the derivatives with respect to "com" atom coordinates)
	    deriv[k][2]=dfunc/VCR;

	    virial-=matmul(Tensor(deriv[k],getPbc().realToScaled(getPosition(k))),transpose(getBox())); //Virial component;
	    //setAtomsDerivatives(k, deriv );
	  }
	  //virial-=Tensor(deriv[k],getPosition(k)); //Virial component;
	}
	/*comm.Sum(deriv);
	for(int i=rank; i<Na_sv; i+=stride){
	  k=Na_st+i;
	  setAtomsDerivatives(k, deriv[k]);
	  }*/
	//	log.printf("Derivatives and CV evaluated \n");              log.flush(); //DBG
	vector<Vector>().swap(com_solv);
      }else{ //if solute specie is restrained
	Vector dleft, dright;
      	vector<Vector> com_solut(N_st);
	//Solute position matrix allocation
	vector<Vector> solut_x(Na_st_permol);
	//Solute mass array allocation
	vector<double> solut_m(Na_st_permol);
	//Solute masses and total mass
	double M_st=0.0;

	for(int i=0;i<Na_st_permol;++i){
	  solut_m[i]=getMass(i);
	  M_st += solut_m[i];
	}

	for(int i=rank; i<N_st; i+=stride){

	  dfunc=0;
	  dleft.zero();
	  dright.zero();
	  dl=0;
	  dr=0;
	  n_lx=0;
	  n_rx=0;
	  //for(int i=0; i<N_st; ++i){
	  com_solut[i].zero();
	  //center of mass
	  if (com_st<0) {
	    solut_x[0] = getPosition(i*Na_st_permol);
	    for(int j=1; j<Na_st_permol; ++j){
	      solut_x[j] = getPosition(i*Na_st_permol+j);
	      diff = pbcDistance(solut_x[0],solut_x[j]);
	      com_solut[i] += solut_m[j]*diff;
	    }
	    com_solut[i] = com_solut[i] / M_st + solut_x[0];
	    //PBC (only orthorhombic!!!)
	    //Only on z
	    if(com_solut[i][2]<0) com_solut[i][2]=com_solut[i][2]+LBC[2];
	    if(com_solut[i][2]>=LBC[2]) com_solut[i][2]=com_solut[i][2]-LBC[2];
	  }else{
	    com_solut[i]=getPosition(i*Na_st_permol+com_st);
	  }

	  //log.printf("\n");              log.flush(); //DBG
	  if(i<=2){
	    //log.printf("solute %d\t %lf\n",i,com_solut[i][2]);              log.flush(); //DBG
	  }
	  //Fermi-like weighting


	  if(asymm<=0){
	    //left-side sigma

	    zin=(com_solut[i][2]-ZCRlin)/w_in;
	    zout=(com_solut[i][2]-ZCRlout)/w_out;
	    //with periodic image
	    n_lx=sigmon(zout,co_out)*sigmoff(zin,co_in)+sigmon(zout-LBC[2]/w_out,co_out)*sigmoff(zin-LBC[2]/w_in,co_in);
	    dleft[2]=com_solut[i][2]-ZFleft;
	    zl=dleft[2]/w_force;
	    dl=(dsig(zl,co_f)+dsig(zl-LBC[2]/w_force,co_f))/w_force;

	    if(i<=2){
	       //log.printf("ZCRlin: %f\t ZCRlout: %f\t w_in: %f w_out: %f \n",ZCRlin,ZCRlout,w_in,w_out);  log.flush(); //DBG
	      //log.printf("zin: %f\t zout: %f\t n_lx: %f\n",zin,zout,n_lx);              log.flush(); //DBG
	      //	    log.printf("l: %f\t %f\t %f\n",zin,zout,n_lx);              log.flush(); //DBG
	    }
	  }

	  if(asymm>=0){
	    //right-side sigma
	    zin=(com_solut[i][2]-ZCRrin)/w_in;
	    zout=(com_solut[i][2]-ZCRrout)/w_out;

	    //with periodic image
	    n_rx=sigmon(zin,co_in)*sigmoff(zout,co_out)+sigmon(zin+LBC[2]/w_in,co_in)*sigmoff(zout+LBC[2]/w_out,co_out);
	    dright[2]=com_solut[i][2]-ZFright;
	    zr=dright[2]/w_force;
	    dr=(-dsig(zr,co_f)-dsig(zr+LBC[2]/w_force,co_f))/w_force;
	  }

	  if(isdelta){
	    n_x=n_rx-n_lx;
	    //sum the two densities (for ASYMMETRIC, change here!!!)
	    dfunc=dr-dl;
	  }else{
	    n_x=n_rx+n_lx;

	    dfunc=dr+dl;
	  }

	  //fdbg<<setprecision(10); //DBG
	  //fdbg<<i<<"\t"<<com_solut[i][2]<<"\t"<<zl<<"\t"<<zr<<"\t"<<dfunc<<"\t"; //DBG
	  //dfunc=(dsig(zl,co_f)-dsig(zr,co_f)+dsig(zl-LBC[2]/w_force,co_f)-dsig(zr+LBC[2]/w_force,co_f))/w_force;//)
	  //fdbg<<dfunc<<endl; //DBG
	  if(com_st<0){ //com coordinates
	    for(int l=0; l<Na_st_permol; ++l){
	      k=i*Na_st_permol+l; //atom counter
	      deriv[k][2] = getMass(k)/M_st*(dfunc/VCR);
	      //setAtomsDerivatives(k, deriv );
	      //fdbg<<setprecision(5); //DBG
	      //fdbg << scientific << k << "\t" << deriv[2] << endl; //DBG
	      //fdbg.flush(); //DBG
	      virial-=matmul(Tensor(deriv[k],getPbc().realToScaled(getPosition(k))),transpose(getBox())); //Virial component;
	    }
	  }else{//single atom coordinates
	    k=i*Na_st_permol+com_st ; //atom counter (just the derivatives with respect to "com" atom coordinates)
	    deriv[k][2] = dfunc/VCR;
	    virial-=matmul(Tensor(deriv[k],getPbc().realToScaled(getPosition(k))),transpose(getBox())); //Virial component;
	  }


	}
	vector<Vector>().swap(com_solut);
      }

      comm.Sum(deriv);
      comm.Sum(n_CR);
      comm.Sum(virial);
      int Natot=Na_st+Na_sv;
      for(int i=0; i< Natot; ++i){
	setAtomsDerivatives(i, deriv[i]);
      }

      vector<Vector>().swap(deriv);
      //setValue         (n_CR); //DBG
      setValue           (n_CR/VCR);
      setBoxDerivatives(virial+n_CR/VCR/getBox().determinant()*matmul(getBox(),transpose(getPbc().getInvBox())));
      //setBoxDerivativesNoPbc();
    }
  }
}
