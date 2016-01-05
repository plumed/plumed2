! <compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"
#include "../include/assert.fh"
#include "ncsu-config.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ driver routine for molecular dynamics
subroutine runmd(xx,ix,ih,ipairs,x,winv,amass,f, &
      v,vold,xr,xc,conp,skip,nsp,tma,erstop, qsetup)

   !  Runmd operates in kcal/mol units for energy, amu for masses,
   !     and angstroms for distances.  To convert the input time parameters
   !     from picoseconds to internal units, multiply by 20.455
   !     (which is 10.0*sqrt(4.184)).

   use state
   
#if !defined(DISABLE_NCSU) && defined(NCSU_ENABLE_BBMD)
   use ncsu_sander_hooks, only : ncsu_on_mdstep => on_mdstep
#endif

   use molecule, only: n_iwrap_mask_atoms, iwrap_mask_atoms
   use cmd_vars, only: activate, file_pos_cmd, file_vel_cmd, file_nrg_cmd,  &
                       nstep_cmd, t_cmd, eq_cmd, restart_cmd,  &
                       etot_cmd, eke_cmd, temp_cmd

   use pimd_vars, only: ipimd, nbead, natomCL, &
                        bnd_vir, Eimp_virial, equal_part, Epot_deriv,  &
                        tau_vol, Epot_spring, NMPIMD, CMD, cartpos, cartvel, &
                        itimass, real_mass

   use neb_vars, only: ineb, neb_nbead

   use lscivr_vars, only: ilscivr, ndof_lsc, natom_lsc, mass_lsc, v2_lsc, &
                          ilsc, x_lsc, f_lsc, dx_lsc

   use nose_hoover_module, only : thermo_lnv, x_lnv, x_lnv_old, v_lnv,  &
                                  f_lnv_p, f_lnv_v, c2_lnv, mass_lnv,  &
                                  Thermostat_init 
 
#ifdef RISMSANDER
   use sander_rism_interface, only: rismprm,rism_3d, RISM_NONE, RISM_FULL, RISM_INTERP, &
        rism_calc_type, rism_solvdist_thermo_calc, mylcm
#endif

   use full_pimd_vars, only: totener,totenert,totenert2,mybeadid

   use qmmm_module, only : qmmm_nml,qmmm_struct, qmmm_mpi, qm2_struct, &
                           qmmm_vsolv
   use file_io_dat
   use constants, only : third, ten_to_minus3, plumed, plumedfile
   use trace
   use stack
   use decomp, only : nat, nrs, decpr, jgroup, indx, irespw, &
#ifdef MPI
   ! -- ti decomp
                      collect_dec, &
#endif
                      checkdec, printdec
   use fastwt
   use bintraj, only: end_binary_frame
   use nblist,only: fill_tranvec,volume,oldrecip,ucell

   use nose_hoover_module, only: Thermostat_switch,  &
                                 Thermostat_integrate_1, Thermostat_integrate_2, & ! APJ
                                 Thermostat_hamiltonian, &                         ! APJ
                                 Adaptive_Thermostat_integrate, &                  ! APJ
                                 Adaptive_Thermostat_hamiltonian, &                ! APJ
                                 file_nhc, nchain, thermo, nthermo, Econserved     ! APJ

#ifdef MPI
   use evb_parm,  only: evb_dyn, nbias
   use evb_data,  only: evb_frc, evb_vel0, evb_bias, evb_nrg, evb_nrg_ave &
                      , evb_nrg_rms, evb_nrg_tmp, evb_nrg_old, evb_nrg_tmp2 &
                      , evb_nrg_old2
   use wigner,    only: rflux
   use remd, only : rem, mdloop, remd_ekmh, repnum, stagid, my_remd_data, &
                    hybrid_remd_ene, next_rem_method
#  ifdef LES
   use evb_pimd,  only: evb_pimd_dealloc
   use miller,    only: i_qi
#  endif
   use softcore, only: ifsc, sc_dvdl, sc_tot_dvdl, sc_tot_dvdl_partner, &
                       sc_dvdl_ee, sc_tot_dvdl_ee, sc_tot_dvdl_partner_ee, &
                       extra_atoms, mix_temp_scaling, sc_pscale, &
                       adj_dvdl_stat, sc_mix_velocities, &
                       sc_nomix_frc, sc_sync_x, sc_print_energies, &
                       calc_softcore_ekin, &
                       sc_ener, sc_ener_ave, sc_ener_rms, sc_lngdyn, &
                       sc_ener_tmp, sc_ener_tmp2, sc_ener_old, sc_ener_old2, &
                       sc_mix_position, sc_print_dvdl_values, &
                       sc_degrees_o_freedom, dynlmb, sc_change_clambda, ti_ene_cnt, &
                       sc_compare
   use mbar, only : ifmbar, bar_intervall, calc_mbar_energies, &
                       bar_collect_cont, do_mbar
#endif

   use amoeba_mdin, only: iamoeba
   use amoeba_runmd, only: AM_RUNMD_scale_cell
   use constantph, only: cnstphinit, cnstphwrite, cnstphupdatepairs, &
                         cnstphbeginstep, cnstphendstep, chrgdat,    &
                         cnstph_explicitmd, cnstphwriterestart, cphfirst_sol
   use emap, only:temap,emap_move
   use barostats, only : mcbar_trial, mcbar_summary
#ifdef EMIL
   use emil_mod,     only : emil_do_calc, emil_calc_AMBER, &
                            emil_save_pme, emil_save_gb, &
                            emil_init, emil_step
#endif
   use memory_module, only: coordinate, velocity, mass

! Self-Guided molecular/Langevin Dynamics (SGLD)
   use sgld, only : isgld, isgsta,isgend,trxsgld, &
                    sgenergy,sgldw,sgmdw,sgfshake, sg_fix_degree_count

   !AWG adaptive QM/MM
   use qmmm_adaptive_module, only: adaptive_qmmm

   use crg_reloc, only: ifcr, crprintcharges, cr_print_charge

   use abfqmmm_module, only: abfqmmm_param, abfqmmm_combine_forces  ! lam81
!AMD
   use amd_mod

!scaledMD
   use scaledMD_mod

!SEBOMD
   use sebomd_module, only : sebomd_obj, sebomd_gradient_write, sebomd_hessian_compute

!     Variable Descriptions
!
! Passed variables
!  xx          : global real array. See locmem.f for structure/pointers
!  ix          : global integer array. See locmem.f for structure/pointers
!  ih          : global hollerith array. See locmem.f for structure/pointers
!  ipairs      : ?? Global pairlist ?? --add description (JMS 11/2010)
!  x           : global position array *
!  winv        : array with inverse masses *
!  amass       : mass array *
!  f           : force array, used to hold old coordinates temporarily, too
!  v           : velocity array
!  vold        : old velocity array, from the previous step
!  xr          : coordinates with respect to COM of molecule
!  conp        : bond parameters for SHAKE
!  skip        : logical skip array for SHAKE (and QM/MM too, I think)
!  nsp         : submolecule index array (?)
!  tma         : submolecular weight array (?)
!  erstop      : should we stop in error (?)
!  qsetup      : Not quite sure what this does, if anything anymore.
!
! Local variables
!  factt       : degree-of-freedom correction factor for temperature scaling
!  nr          : local copy of nrp, number of atoms
!  nr3         : 3 * nr, used for runtime efficiency
!
! Common memory variables
!  nrp         : number of atoms, adjusted for LES copies

   implicit none
   character(kind=1,len=5) :: routine="runmd"
   integer   ipairs(*), ix(*)
   _REAL_ xx(*)
   character(len=4) ih(*)
   _REAL_ combination, rem_val

#ifdef MPI
#  include "parallel.h"
   include 'mpif.h'
#  ifdef LES
   _REAL_  :: fbead(3,natomCL), xbead(3,natomCL)
   integer :: mm, n
#  endif   
   _REAL_ mpitmp(8) !Use for temporary packing of mpi messages.
   integer ist(MPI_STATUS_SIZE), partner, ierr
#else
   ! mdloop and REM is always 0 in serial
   integer, parameter :: mdloop = 0, rem = 0
#endif

! The following variables are needed since nstep and nstlim
!  behave differently in a REMD run.
! In certain places where output is required, total_nstep and total_nstlim 
!  take the place of nstep and nstlim. This allows replica runs to output
!  less often than every exchange.
! They are the absolute step # of the REMD or MD simulation.
   integer total_nstep, total_nstlim

#include "../include/md.h"
#include "box.h"
#include "nmr.h"
#include "tgtmd.h"
#include "multitmd.h"
#include "../include/memory.h"
#include "extra.h"
#include "ew_frc.h"
#include "ew_cntrl.h"
#include "ew_mpole.h"
#include "def_time.h"
#include "extra_pts.h"
#if defined(LES)
#  include "les.h"
#endif
#include "../pbsa/pb_md.h"
#include "../lib/random.h"   
      
! additional variables for PIMD output
   _REAL_  :: xcmd(3*natomCL),vcmd(3*natomCL)
   integer :: ncmd
   ! for const press PIMD
   _REAL_ tmpvir(3,3),atomvir

   _REAL_ sgsta_rndfp, sgend_rndfp, ignore_solvent
   _REAL_ sysx,sysy,sysz,sysrange(3,2)
   logical mv_flag

#ifdef MMTSB
#  include "mmtsb.h"
   logical is_done_mmtsb      ! MMTSB replica exchange calculation completed
   _REAL_  lambda_mmtsb       ! MMTSB replica exchange new lambda
   _REAL_  pert_pe_mmtsb      ! MMTSB lambda replica exchange perturbed PE
   _REAL_  temp_mmtsb         ! MMTSB replica exchange new temperature
   _REAL_  unpert_pe_mmtsb    ! MMTSB lambda replica exchange unperturbed PE
#endif

   _REAL_ , dimension(1) :: shkh
   integer, dimension(1) :: ifstwr2
   integer :: nshkh

   integer idx, iatom, iatomCL,m
   _REAL_  Ekin2_tot,tmp  ! APJ
   integer :: idim, ithermo
   _REAL_ :: E_nhc, exp1, exp2, v_sum
   
   logical ivscm
   logical qspatial
   character(len=6)fnam

   logical resetvelo
   integer nshak
   _REAL_ ekgs,eold3,eold4,etot_save,ekpbs
   
   logical do_list_update
   logical skip(*),belly,lout,loutfm,erstop,vlim,onstep
   _REAL_ x(*),winv(*),amass(*),f(*),v(*),vold(*), &
         xr(*),xc(*),conp(*)
   type(state_rec) :: ener   ! energy values per time step
   type(state_rec) :: enert  ! energy values tallied over the time steps
   type(state_rec) :: enert2 ! energy values squared tallied over the time steps
   type(state_rec) :: enert_old, enert2_old
   type(state_rec) :: enert_tmp, enert2_tmp
   type(state_rec) :: ecopy, edvdl
   type(state_rec) :: edvdl_r
   _REAL_ rmu(3),fac(3),onefac(3),clfac, etot_start
   _REAL_ tma(*)

   _REAL_ tspan,atempdrop,fln,scaltp,scaltpo
   _REAL_ vel,vel2,vcmx,vcmy,vcmz,vmax,vx,vy,vz
   _REAL_ winf,aamass,rterm,ekmh,ekph,ekpht,wfac,rsd,ekav
   _REAL_ fit,fiti,fit2,vscalt
   logical is_langevin  ! Is this a Langevin dynamics simulation
   _REAL_ gammai,c_implic,c_explic,c_ave,sdfac,ekins0
   _REAL_ dtx,dtxinv,dt5,factt,ekin0,ekinp0,dtcp,dttp
   _REAL_ rndf,rndfs,rndfp,boltz2,pconv,tempsu
   _REAL_ xcm(3),acm(3),ocm(3),vcm(3),ekcm,ekrot
   _REAL_ emtmd

! Variables and parameters for constant surface tension:
   _REAL_, parameter :: ten_conv = 100.0d0 !ten_conv - converts
                                                    !dyne/cm to bar angstroms
   _REAL_  :: pres0x
   _REAL_  :: pres0y
   _REAL_  :: pres0z
   _REAL_  :: gamma_ten_int
   _REAL_  :: press_tan_ave

   integer nsp(*)
   integer idumar(4)
   integer l_temp
   integer i,j,im,i3,nitp,nits, iskip_start,iskip_end  ! APJ
   integer nstep,nrep,nrek,nren,iend,istart3,iend3
   integer nrx,nr,nr3,ntcmt,izero,istart
   logical ixdump,ivdump,itdump,ifdump
   logical qsetup
   _REAL_, allocatable, dimension(:) :: for ! lam81
#ifdef RISMSANDER
   logical irismdump
   _REAL_  cm(3),angvel(3),r(3),rxv(3),proj(3),moi,erot
#endif
   
   integer nvalid, nvalidi
   _REAL_ eke,eket
   _REAL_ extent

   _REAL_ xcen,ycen,zcen,extents(3,2)
   _REAL_, allocatable, dimension(:) :: frcti
   integer ier

   _REAL_ small
   data small/1.0d-7/
   data nren/51/

   !--- VARIABLES FOR DIPOLE PRINTING ---
   integer prndipngrp
   integer prndipfind
   character(len=4) prndiptest

   _REAL_,parameter :: pressure_constant = 6.85695d+4
   ! variables used in constant pressure PIMD
   _REAL_ :: Nkt,centvir,pressure, aa, arg2, poly, e2, e4, e6, e8 
   ! variable used in CMD
   real(8) :: tmp_eke_cmd !Use for temporary packing of mpi messages.

   _REAL_ :: box_center(3)

! for adaptive qm/mm runs

   _REAL_ :: adqmmm_first_energy, etotcorr, tadc
   integer :: nstepadc
   logical :: flag_first_energy = .true.

   _REAL_ :: xold(3*natom)
   _REAL_ :: corrected_energy
   _REAL_ :: kinetic_E_save(2)
   integer :: aqmmm_flag 

! variables for plumed
   _REAL_ :: plumed_box(3,3),plumed_virial(3,3), plumed_kbt
   integer :: plumed_version,plumed_stopflag,plumed_ms
   _REAL_ :: plumed_energyUnits,plumed_timeUnits,plumed_lengthUnits,plumed_chargeUnits

   !==========================================================================
  
   call trace_enter( 'runmd' )

   !     ----- INITIALIZE SOME VARIABLES -----
   
#ifdef MPI
   if( master ) then
      ! If remd, runmd will be called many times, so we dont want to open every
      !  time. For normal md, mdloop will just be 0.
      if (mdloop.eq.0) call amopen(7,mdinfo,'U','F',facc)
   endif

   if (rem < 3) then
      rem_val = temp0
   else if (rem == 4) then
      rem_val = solvph
   else
      rem_val = 0.d0
   end if
#else
   if( master ) call amopen(7,mdinfo,'U','F','W')
#endif
   vlim = vlimit > small
   ntcmt = 0
   izero = 0
   belly = ibelly > 0
   lout = .true.
   loutfm = ioutfm <= 0
   nr = nrp
   nr3 = 3*nr
   ekmh = 0.d0

   aqmmm_flag = 0

#ifdef LES
   ekmhles = 0.d0
#endif
      do_list_update=.false.
#ifdef MPI
   if ( mpi_orig ) then
      istart = 1
      iend = natom
   else
      istart = iparpt(mytaskid) + 1
      iend = iparpt(mytaskid+1)
   end if
#else
   istart = 1
   iend = nr
#endif
   istart3 = 3*istart -2
   iend3 = 3*iend

#ifdef MPI
   if( icfe /= 0 ) then
      allocate( frcti( nr3+3*extra_atoms ), stat = ier )
      REQUIRE( ier == 0 )
   end if
#endif

   ! If NTWPRT.NE.0, only print the atoms up to this value
   nrx  = nr3
   if (ntwprt > 0) nrx = ntwprt*3

   if (.not. allocated(for)) allocate(for(nr3))  ! lam81

   if(abfqmmm_param%abfqmmm == 1) then                  ! lam81
#ifdef MPI
     call xdist(v, xx(lfrctmp), natom)                  ! lam81 
#endif
     if(abfqmmm_param%system == 1) then                 ! lam81
       if(abfqmmm_param%qmstep == 1) abfqmmm_param%v(1:nr3+iscale) = v(1:nr3+iscale)  ! lam81
       v(1:nr3+iscale) = 0.d0                           ! lam81
       t = t+dt                                         ! lam81
       if(abfqmmm_param%maxqmstep == 0) t = 0           ! lam81
     else                                               ! lam81
       v(1:nr3+iscale) = abfqmmm_param%v(1:nr3+iscale)  ! lam81
     endif                                              ! lam81
   endif                                                ! lam81
   
   ! Cleanup the velocity if belly run
   if(belly) call bellyf(nr,ix(ibellygp),v)
   
   !=======================================================================
   ! Determine system degrees of freedom (for T scaling, reporting)
   
   ! Call DEGCNT to get the actual number of degrees of freedom for the
   ! solute and solvent. This call returns the correct numbers for belly
   ! simulations and simulations with separate solute/solvent scaling -- dap
   ! "IDUMAR" is dummy array. Used since this routine was also used w/ GIBBS.
   
#ifdef LES
   ! return LES and non-LES degrees,
   ! since separate solvent coupling no longer used
   ! large changes to degcnt were made
   ! cnum is now passed (LES copy number of each atom)
   call degcnt(ibelly,nr,ix(ibellygp),nsolut,nbonh,nbona,0, &
         ix(iibh),ix(ijbh),ix(iiba),ix(ijba),idumar, &
         idumar,ntc,idumar,0,0,0, &
         idumar,rndfp,rndfles,cnum,temp0les)
   
   ! RNDFP = # degrees of freedom for solute
   ! RNDFS = # degrees of freedom for solvent
   ! RNDF = total number of degrees of freedom.
   ! RNDFLES = # degrees of freedom for LES groups
   
   ! temp0les was init to negative number to signify not to use a LES bath
   ! just do standard code (meaning use solute/solvent baths)
   ! any positive (or zero) means to use LES bath with that target
   
   ! degcnt returns rndfs or rndfles in the rndfles variable
   ! depending on whether a LES bath was specified
   ! do this instead of duplicating call with rndfs or rndfles
   
   if (temp0les < 0.d0) then
      rndfs=rndfles
      rndfles=0.d0
   else
      rndfs=0.d0
   end if

   if (master) then
      write (6,'(a,f8.0)') &
            "# degrees of freedom in non-LES region: ",rndfp
      write (6,'(a,f8.0)') &
            "# degrees of freedom in     LES region: ",rndfles
   end if
   
   !    modify RNDFP to reflect NDFMIN (set in mdread)
   
   rndfp = rndfp - ndfmin

   if (temp0les < 0.d0) then
      rndf = rndfp+rndfs
   else
      rndf = rndfp+rndfles
   end if
   
#else
   
   call degcnt(ibelly,nr,ix(ibellygp),nsolut,nbonh,nbona,0, &
         ix(iibh),ix(ijbh),ix(iiba),ix(ijba),idumar, &
         idumar,ntc,idumar,0,0,0, &
         idumar,rndfp,rndfs)
   
   ! RNDFP = # degrees of freedom for solute
   ! RNDFS = # degrees of freedom for solvent
   ! RNDF = total number of degrees of freedom.

#ifdef MPI
   if (mdloop .eq. 0 .and. master) then
#else
   if (master) then
#endif
   if (abfqmmm_param%abfqmmm /= 1 .or. (abfqmmm_param%qmstep == 1 .and. abfqmmm_param%system == 1)) then  ! lam81
      write (6,'(a,f8.0)') &
            "|  # of SOLUTE  degrees of freedom (RNDFP): ",rndfp
      write (6,'(a,f8.0)') &
            "|  # of SOLVENT degrees of freedom (RNDFS): ",rndfs
   end if ! lam81
   end if
   ! qtw - substract the number of overlapping noshake QM atoms in noshakemask
   rndfp = rndfp - qmmm_struct%noshake_overlap
   !    modify RNDFP to reflect NDFMIN (set in mdread) and num_noshake
   rndfp = rndfp - ndfmin + num_noshake
   rndf = rndfp + rndfs
#ifdef MPI
   if (mdloop .eq. 0 .and. master) then
#else
   if (master) then
#endif
   if (abfqmmm_param%abfqmmm /= 1 .or. (abfqmmm_param%qmstep == 1 .and. abfqmmm_param%system == 1)) then  ! lam81
      if (qmmm_nml%ifqnt) then
         write (6,'(a,i6)') &
               "|  QMSHAKE_NOSHAKEMASK_OVERLAP = ", qmmm_struct%noshake_overlap
      endif
      write (6,'(a,f8.0,a,i6,a,f8.0)') &
            "|  NDFMIN = ",rndfp, "     NUM_NOSHAKE = ",num_noshake, "     CORRECTED RNDFP = ", rndfp
      write (6,'(a,f8.0)') &
            "|  TOTAL # of degrees of freedom (RNDF) = ", rndf
   end if  ! lam81
   end if

#endif

   call fix_degree_count(rndf) ! correct for extra points
! Warning - NOTE that rndfp, rndfs are uncorrected in an extra points context!
   
#ifndef LES
   if (isgld > 0) then
      ! number of degrees of freedom in the SGLD part
      if (isgsta == 1) then
         sgsta_rndfp = 0
      else
         call degcnt(ibelly,nr,ix(ibellygp),isgsta-1,nbonh,nbona,0, &
               ix(iibh),ix(ijbh),ix(iiba),ix(ijba),idumar, &
               idumar,ntc,idumar,0,0,0,idumar,sgsta_rndfp,ignore_solvent)
      end if
      if (isgend == nr) then
         sgend_rndfp = rndf
      else
         call degcnt(ibelly,nr,ix(ibellygp),isgend,nbonh,nbona,0, &
               ix(iibh),ix(ijbh),ix(iiba),ix(ijba),idumar, &
               idumar,ntc,idumar,0,0,0,idumar,sgend_rndfp,ignore_solvent)
      end if
! Warning - NOTE that the solute ndf outputs above from degcnt are uncorrected
! for qmmm_struct%noshake_overlap, num_noshake, and extra points;
! also ndfmin is not always being handled.
      call sg_fix_degree_count(sgsta_rndfp, sgend_rndfp, ndfmin, rndf)
   end if
#endif

#ifdef MPI /* SOFT CORE */
   if (ifsc /=0 ) call sc_degrees_o_freedom(ndfmin)
#endif

   ! End of degrees of freedom setup
   !=======================================================================
   
   boltz2 = 8.31441d-3 * 0.5d0
   pconv = 1.6604345d+04  ! factor to convert the pressure kcal/mole to bar
   
   !     ---convert to kcal/mol units
   
   boltz2 = boltz2/4.184d0   ! k-sub-B/2
   dtx = dt*20.455d+00
   dtxinv = 1.0d0 / dtx
   dt5 = dtx * 0.5d0
   pconv = pconv*4.184d0
   
   ! FAC() are #deg freedom * kboltz / 2
   ! multiply by T to get expected kinetic energy
   ! FAC(1) is for total system
   
   fac(1) = boltz2*rndf
   fac(2) = boltz2*rndfp

   if(rndfp < 0.1d0) fac(2) = 1.d-6

#ifdef LES
   ! replaced solvent variables with LES ones
   ! since separate solvent coupling no longer used
   ! ASSUME SAME COUPLING CONSTANT FOR BOTH BATHS, just different target T
   
   ! will also have to accumulate LES and non-LES kinetic energies separately
   
   if (temp0les < 0.d0) then
      fac(3) = boltz2*rndfs
      if(rndfs < 0.1d0) fac(3) = 1.d-6
   else
      fac(3) = boltz2*rndfles
      if(rndfles < 0.1d0) fac(3) = 1.d-6
   end if
#else
   fac(3) = boltz2*rndfs
   if(rndfs < 0.1d0) fac(3) = 1.d-6
#endif
   if ( ipimd==CMD ) then
      if ( eq_cmd ) then
         fac(1) = boltz2 * dble( 3*natomCL )
      else
         fac(1) = boltz2 * dble( 3*(natomCL-1) )
      endif
   endif
   onefac(1) = 1.0d0/fac(1)
   onefac(2) = 1.0d0/fac(2)
   onefac(3) = 1.0d0/fac(3)
   factt = rndf/(rndf+ndfmin)
   
   ! these are "desired" kinetic energies based on
   ! # degrees freedom and target temperature
   ! they will be used for calculating the velocity scaling factor
   
   ekinp0 = fac(2)*temp0
#ifdef LES
   
   ! modified for LES temperature
   
   ekins0=0.d0
   ekinles0=0.d0
   if (temp0les < 0.d0) then
      ekins0 = fac(3) * temp0
      ekin0  = fac(1) * temp0
      if (master) &
            write (6,*) "Single temperature bath for LES and non-LES"
   else
      ekinles0 = fac(3)*temp0les
      ekin0  = ekinp0 + ekinles0
      if (master) then
         write (6,*) "LES particles coupled to separate bath"
         write (6,'(a,f8.2)')"    LES target temperature:    ",temp0les
         write (6,'(a,f8.2)')"    LES target kinetic energy: ",ekinles0
         write (6,'(a,f8.2)')"non-LES target temperature:    ",temp0
         write (6,'(a,f8.2)')"non-LES target kinetic energy: ",ekinp0
      end if
   end if
#else
   ekins0 = fac(3)*temp0
   ekin0  = fac(1)*temp0
#endif

#ifdef LES
   if(abfqmmm_param%abfqmmm /= 1) then                  ! lam81
    if ( ntt==4 ) call nose_hoover_init_LES(amass,v,f)  ! APJ
   else                                                 ! lam81
    if ( ntt==4 .and. abfqmmm_param%qmstep == 1 .and. abfqmmm_param%system == 1 ) &  ! lam81
      call nose_hoover_init_LES(amass,abfqmmm_param%v,abfqmmm_param%f)               ! lam81
   endif                                                ! lam81
#else
   if(abfqmmm_param%abfqmmm /= 1) then                  ! lam81
    if ( ntt>=4 .and. ntt<=8 ) call nose_hoover_init(amass,v,f)  ! APJ
   else                                                 ! lam81
    if ( ntt>=4 .and. ntt<=8 .and. abfqmmm_param%qmstep == 1 .and. abfqmmm_param%system == 1 ) &  ! lam81
      call nose_hoover_init(amass,abfqmmm_param%v,abfqmmm_param%f)                                ! lam81
   endif
#endif

   ! Langevin dynamics setup:

   is_langevin = gamma_ln > 0.0d0
   gammai = gamma_ln/20.455d0
   c_implic = 1.d0/(1.d0+gammai*dt5)
   c_explic = 1.d0 - gammai*dt5
   c_ave    = 1.d0+gammai*dt5
   sdfac = sqrt( 4.d0*gammai*boltz2*temp0/dtx )
#ifdef LES
   if( temp0les < 0.d0 ) then
      sdfacles = sqrt( 4.d0*gammai*boltz2*temp0/dtx )
   else
      sdfacles = sqrt( 4.d0*gammai*boltz2*temp0les/dtx )
   endif
#endif
   if (is_langevin .and. ifbox==0) then
      call get_position(nr,x,sysx,sysy,sysz,sysrange,0)
#ifdef MPI /* SOFT CORE */
      if (ifsc == 1) call sc_mix_position(sysx,sysy,sysz,clambda)
#endif
   end if

   !     Constant pH setup 
   !
   if (icnstph /= 0 .and. mdloop .eq. 0) &
      call cnstphinit(x, ig)
   
   if (ntt == 1) dttp = dt/tautp
   if (ntp > 0) dtcp = comp * 1.0d-06 * dt / taup

   ! Constant surface tension setup:

   if (csurften > 0) then

      ! Set pres0 in direction of surface tension.
      ! The reference pressure is held constant in on direction dependent
      ! on what the surface tension direction is set to.
      if (csurften .eq. 1) then           ! pres0 in the x direction
        pres0x = pres0

      else if (csurften .eq. 2) then      ! pres0 in the y direction
        pres0y = pres0

      !else if (csurften .eq. 3) then      ! pres0 in the z direction
      else
        pres0z = pres0

      end if

      ! Multiply surface tension by the number of interfaces
      gamma_ten_int = dble(ninterface) * gamma_ten
 
   end if

   nrek = 4
   nrep = 15
   
   nvalid = 0
   nvalidi = 0
   nstep = 0
   total_nstep = 0
#ifdef MPI
   ! For REMD, total_nstep is the number of steps * the number of exchanges
   ! we've already attempted
   if (rem /= 0) &
      total_nstep = (mdloop - 1) * nstlim
#endif
   fit = 0.d0
   fiti = 0.d0
   fit2 = 0.d0

   ! Zero all elements of these sequence types
   ener       = null_state_rec
   enert      = null_state_rec
   enert2     = null_state_rec
   enert_old  = null_state_rec
   enert2_old = null_state_rec
   edvdl      = null_state_rec
   edvdl_r    = null_state_rec
   ! for PIMD/NMPIMD/CMD/RPMD:
   totenert   = null_state_rec
   totenert2  = null_state_rec

   ener%kin%pres_scale_solt = 1.d0
   ener%kin%pres_scale_solv = 1.d0
   ener%box(1:3) = box(1:3)

   
   ener%cmt(1:4) = 0.d0
   nitp = 0
   nits = 0

! init PLUMED
   if(plumed.eq.1) then
     call plumed_f_gcreate()
#ifdef DPREC
     call plumed_f_gcmd("setRealPrecision"//char(0),8)
#else
     call plumed_f_gcmd("setRealPrecision"//char(0),4)
#endif
     call plumed_f_gcmd("getApiVersion"//char(0),plumed_version)
     if(plumed_version>1)then
       plumed_kbt=2.0*temp0*boltz2
       call plumed_f_gcmd("setKbT"//char(0),plumed_kbt)
     endif
     plumed_energyUnits=4.184
     plumed_lengthUnits=0.1
     plumed_timeUnits=1.0
     plumed_chargeUnits=1.0/18.2223
     call plumed_f_gcmd("setMDEnergyUnits"//char(0),plumed_energyUnits)
     call plumed_f_gcmd("setMDLengthUnits"//char(0),plumed_lengthUnits)
     call plumed_f_gcmd("setMDTimeUnits"//char(0),plumed_timeUnits)
     if(plumed_version>3)then
       call plumed_f_gcmd("setMDChargeUnits"//char(0),plumed_chargeUnits)
     endif
     call plumed_f_gcmd("setPlumedDat"//char(0),trim(adjustl(plumedfile))//char(0))
     call plumed_f_gcmd("setNatoms"//char(0),nr)
     call plumed_f_gcmd("setMDEngine"//char(0),"amber")
     call plumed_f_gcmd("setTimestep"//char(0),dt)
#  ifdef MPI
     call plumed_f_gcmd("setMPIFComm"//char(0),commsander)
     if(numgroup>1)then
       call plumed_f_gcmd("GREX setMPIFIntracomm"//char(0),commsander)
       if(master) then
         call plumed_f_gcmd("GREX setMPIFIntercomm"//char(0),commmaster)
       endif
       call plumed_f_gcmd("GREX init"//char(0),0)
     endif
#  endif
     call plumed_f_gcmd("init"//char(0),0);


!     if(ifbox/=0 .and. ifbox/=1 .and. ifbox/=2) then
!      write (6,*) "!!!!! PLUMED ERROR: Only orthorhombic and truncted octahedron cells are supported in this release."
!      write (6,*) "!!!!! ABORTING RUN"
!      stop
!     endif
!     call init_metadyn(nr,dt,amass,xx(l15),ifbox,0,trim(adjustl(plumedfile))//char(0))
     continue
   endif
 ! end init PLUMED



   !=======================================================================
   !     ----- MAKE A FIRST DYNAMICS STEP -----
   !=======================================================================
   !  init = 3:  general startup if not continuing a previous run

   if( ipimd.eq.NMPIMD .or. ipimd.eq.CMD) then
      call trans_pos_cart_to_nmode( x )
   end if
  
   if( init == 3 .or. nstlim == 0 .or. (abfqmmm_param%abfqmmm == 1 .and. abfqmmm_param%system == 1) ) then  ! lam81
      if (ntp > 0 .and. iamoeba==0 .and. ipimd==0) then
         xr(1:nr3) = x(1:nr3)
         
         ! ----- CALCULATE THE CENTER OF MASS ENERGY AND THE COORDINATES
         !       OF THE SUB-MOLECULES WITH RESPECT TO ITS OWN CENTER OF
         !       MASS -----
         call ekcmr(nspm,nsp,tma,ener%cmt,xr,v,amass,1,nr)
      end if
      
      ! ----- CALCULATE THE FORCE -----

      !   ---   set irespa to get full energies calculated on step "0":
      irespa = 0
      iprint = 1

      if(ipimd==NMPIMD .or. ipimd==CMD) then
         call trans_pos_nmode_to_cart(x,cartpos)
         call force(xx,ix,ih,ipairs,cartpos,f,ener,ener%vir, &
                  xx(l96),xx(l97),xx(l98),xx(l99), qsetup, &
                  do_list_update,nstep)
#if defined(MPI) && defined(LES)
         if ( ievb == 1 .and. i_qi > 0) then
            call evb_umb ( f, cartpos, real_mass, natom, istart3, iend3 )
! 03132009            if( i_qi == 2 ) call qi_corrf_les ( cartpos, amass )
            if( i_qi == 2 ) call qi_corrf_les ( cartpos, real_mass )
            evb_nrg(1) = evb_frc%evb_nrg
            evb_nrg(2) = evb_vel0%evb_nrg
            if( nbias > 0 ) evb_nrg(3) = sum( evb_bias%nrg_bias(:) )
         endif
#endif

         call trans_frc_cart_to_nmode(f)
         i3 = 3*(istart-1)

#if defined(MPI) && defined(LES)
         if ( ievb /= 0 .and. i_qi == 0 ) then
            call evb_umb ( f, x, real_mass, natom, istart3, iend3 )
            evb_nrg(1) = evb_frc%evb_nrg
            evb_nrg(2) = evb_vel0%evb_nrg
            if( nbias > 0 ) evb_nrg(3) = sum( evb_bias%nrg_bias(:) )
         endif
#endif

      else if ( ilscivr == 1 )then
         ! prepare the Hessian Matrix of the potential for the LSC-IVR
         ! at this point, x is the position a bead at equilibrium
         ! initialize the LSC-IVR variables
         natom_lsc = natom
         ndof_lsc = natom * 3

         call lsc_init
         do ilsc = 1, natom_lsc
            mass_lsc(3*ilsc-2) = amass(ilsc)
            mass_lsc(3*ilsc-1) = amass(ilsc)
            mass_lsc(3*ilsc  ) = amass(ilsc)
         end do
         v2_lsc = 0.0d0
         do ilsc = 1, ndof_lsc
            ! ith vector of the Hesian matrix
            x_lsc = 0.0d0
            x_lsc(1:ndof_lsc) = x(1:ndof_lsc)
            x_lsc(ilsc) = x(ilsc) + dx_lsc
            call force(xx,ix,ih,ipairs,x_lsc,f_lsc,ener,ener%vir, &
                       xx(l96),xx(l97),xx(l98),xx(l99), qsetup, &
                       do_list_update,nstep)

#ifdef MPI
            call xdist( f_lsc, xx(lfrctmp), natom )
#endif
            v2_lsc(1:ndof_lsc,ilsc) = f_lsc(1:ndof_lsc)
         enddo

         call force(xx,ix,ih,ipairs,x,f,ener,ener%vir, &
               xx(l96),xx(l97),xx(l98),xx(l99), qsetup, &
               do_list_update,nstep)

#ifdef MPI
         call xdist(f, xx(lfrctmp), natom)
#endif
         ! 2nd derivative of the potential:
         do ilsc = 1, ndof_lsc
            v2_lsc(1:ndof_lsc,ilsc) = &
                ( f(1:ndof_lsc) - v2_lsc(1:ndof_lsc,ilsc) )/dx_lsc
         end do

         ! get the iniital position of the momentum:
         call lsc_xp(x,v)

      else

         ! -- ti decomp
         if(idecomp > 0) then
            decpr = .false.
            if(mod(nstep+1,ntpr) == 0) decpr = .true.
         end if
         call force(xx,ix,ih,ipairs,x,f,ener,ener%vir, &
                  xx(l96),xx(l97),xx(l98),xx(l99), qsetup, &
                  do_list_update,nstep)
#ifdef MPI
         if ( ievb /= 0 ) then
#ifdef LES
            call evb_umb_primitive ( f, x, real_mass, natom, istart, iend )
#else
            call evb_umb_primitive ( f, x, amass, natom, istart, iend )
#endif
            evb_nrg(1) = evb_frc%evb_nrg
            evb_nrg(2) = evb_vel0%evb_nrg
            if( nbias > 0 ) evb_nrg(3) = sum( evb_bias%nrg_bias(:) )
         endif
#endif

      endif

      if (sebomd_obj%do_sebomd) then
        ! computes hessian matrix if necessary
        if (sebomd_obj%ntwh /= 0) then
           ! don't output atomic charges
           sebomd_obj%iflagch_old = sebomd_obj%iflagch
           sebomd_obj%iflagch = 0
           call sebomd_gradient_write(f,3*natom)
           call sebomd_hessian_compute(xx,ix,ih,ipairs,x,f,ener, &
                             qsetup, do_list_update, nstep)
           sebomd_obj%iflagch = sebomd_obj%iflagch_old
        endif
      endif


      if (icnstph /= 0 .and. master .and. &
         ((rem /= 0 .and. mdloop > 0) .or. rem == 0)) call cnstphwrite(rem)

      for(1:nr3) = f(1:nr3)             ! lam81
#ifdef MPI
      call xdist(for,xx(lfrctmp),natom) ! lam81
#endif

      if(abfqmmm_param%abfqmmm == 1) then                                        ! lam81

        if(abfqmmm_param%system == 1) abfqmmm_param%f1(1:nr3) = for(1:nr3)       ! lam81

        if(abfqmmm_param%system == 2) then                                       ! lam81
         abfqmmm_param%f2(1:nr3) = for(1:nr3)                                    ! lam81
         call abfqmmm_combine_forces()                                           ! lam81 
         for(1:nr3) = abfqmmm_param%f(1:nr3)                                     ! lam81
         f(1:nr3) = abfqmmm_param%f(1:nr3)                                       ! lam81
        end if                                                                   ! lam81

      end if                                                                     ! lam81

      ! This FORCE call does not count as a "step". CALL NMRDCP to decrement
      ! local NMR step counter and MTMDUNSTEP to decrease the local MTMD step
      ! counter
      call nmrdcp
      call mtmdunstep

   !PLUMED force added
   plumed_stopflag=0
   if(plumed.eq.1) then
     call plumed_f_gcmd("setStep"//char(0),nstep)
     call plumed_f_gcmd("setPositions"//char(0),x)
     call plumed_f_gcmd("setMasses"//char(0),amass)
     call plumed_f_gcmd("setCharges"//char(0),xx(l15))
     call plumed_f_gcmd("setEnergy"//char(0),ener%pot)
     call plumed_f_gcmd("setForces"//char(0),f)
     call plumed_f_gcmd("setStopFlag"//char(0),plumed_stopflag)
     plumed_box=0.0
     if(ifbox==0) then
       continue
     else if(ifbox==1) then
       plumed_box(1,1)=box(1)
       plumed_box(2,2)=box(2)
       plumed_box(3,3)=box(3)
     else if(ifbox==2) then
! truncated octahedron, corresponding to a bcc lattice
! in AMBER convention, box(1) is the length of the lattice vector
! a is defined so as the bcc lattice is (a/2,a/2,a/2) (-a/2,-a/2,a/2) (a/2,-a/2,-a/2)
       plumed_box(1,1)=sqrt(1.0/3.0)*box(1)
       plumed_box(2,1)=sqrt(1.0/3.0)*box(1)
       plumed_box(3,1)=sqrt(1.0/3.0)*box(1)
       plumed_box(1,2)=-sqrt(1.0/3.0)*box(1)
       plumed_box(2,2)=-sqrt(1.0/3.0)*box(1)
       plumed_box(3,2)=sqrt(1.0/3.0)*box(1)
       plumed_box(1,3)=sqrt(1.0/3.0)*box(1)
       plumed_box(2,3)=-sqrt(1.0/3.0)*box(1)
       plumed_box(3,3)=-sqrt(1.0/3.0)*box(1)
     else
      write (6,*) "!!!!! PLUMED ERROR: Only orthorhombic and truncted octahedron cells are supported in this release."
      write (6,*) "!!!!! ABORTING RUN"
      stop
     endif
     plumed_virial=0.0
     plumed_virial(1,1)=2.0*ener%vir(1)
     plumed_virial(2,2)=2.0*ener%vir(2)
     plumed_virial(3,3)=2.0*ener%vir(3)
     call plumed_f_gcmd("setVirial"//char(0),plumed_virial)
     call plumed_f_gcmd("setBox"//char(0),plumed_box)
     call plumed_f_gcmd("calc"//char(0),0);
#ifdef MPI
! this is required since PLUMED only updates virial on master processor
#ifdef DPREC
     call mpi_bcast(plumed_virial,9,MPI_DOUBLE_PRECISION,0,commsander,ierr)
#else
     call mpi_bcast(plumed_virial,9,MPI_REAL,0,commsander,ierr)
#endif
#endif
     ener%vir(1)=0.5*plumed_virial(1,1)
     ener%vir(2)=0.5*plumed_virial(2,2)
     ener%vir(3)=0.5*plumed_virial(3,3)
   end if

   !PLUMED end


#ifdef MPI /* SOFT CORE */
      ! If softcore potentials are used, collect their dvdl contributions:
      if ( ifsc /= 0 ) then
         call mpi_reduce(sc_dvdl, sc_tot_dvdl, 1, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, 0, commsander, ierr)
         sc_dvdl=0.0d0 ! zero for next step
         call mpi_reduce(sc_dvdl_ee, sc_tot_dvdl_ee, 1, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, 0, commsander, ierr)
         sc_dvdl_ee=0.0d0 ! zero for next step
         call mpi_reduce(sc_ener, sc_ener_tmp, ti_ene_cnt, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, 0, commsander, ierr)
         sc_ener(1:ti_ene_cnt) = sc_ener_tmp(1:ti_ene_cnt)
      end if
      if ( ifsc == 2 ) then
         ! If this is a perturb to nothing run, scale forces and calculate dvdl
         call sc_nomix_frc(f,nr3,ener)
         if( numtasks>1 ) then
            call mpi_bcast(f,nr3,MPI_DOUBLE_PRECISION,0,commsander,ierr)
            call mpi_bcast(ener,state_rec_len,MPI_DOUBLE_PRECISION,0,commsander,ierr)
         end if
      end if

      if( icfe /= 0 )then
         ! ---free energies using thermodynamic integration (icfe /= 0)

         if( master ) then
            !  --- first, send the forces and energy to your partner:
            partner = ieor(masterrank,1)
            call mpi_sendrecv( f, nr3, MPI_DOUBLE_PRECISION, partner, 5, &
                               frcti, nr3+3*extra_atoms, MPI_DOUBLE_PRECISION, &
                               partner, 5, commmaster, ist, ierr )
            call mpi_sendrecv( ener, state_rec_len, MPI_DOUBLE_PRECISION, partner, 5, &
                               ecopy, state_rec_len, MPI_DOUBLE_PRECISION, partner, 5, &
                               commmaster, ist, ierr)
            ! exchange sc-dvdl contributions between masters
            call mpi_sendrecv( sc_tot_dvdl, 1, MPI_DOUBLE_PRECISION, partner, &
                               5, sc_tot_dvdl_partner, 1, &
                               MPI_DOUBLE_PRECISION, partner, 5, &
                               commmaster, ist, ierr )
            call mpi_sendrecv( sc_tot_dvdl_ee, 1, MPI_DOUBLE_PRECISION, partner, &
                               5, sc_tot_dvdl_partner_ee, 1, &
                               MPI_DOUBLE_PRECISION, partner, 5, &
                               commmaster, ist, ierr )
            if( masterrank==0 ) then
               call mix_frcti(frcti,ecopy,f,ener,nr3,clambda,klambda)
            else
               call mix_frcti(f,ener,frcti,ecopy,nr3,clambda,klambda)
            end if
         end if

         if( numtasks>1 ) then
            call mpi_bcast(f,nr3,MPI_DOUBLE_PRECISION,0,commsander,ierr)
            call mpi_bcast(ener,state_rec_len,MPI_DOUBLE_PRECISION,0,commsander,ierr)
         end if
 
      end if

#endif /* MPI SOFT CORE */

      irespa = 1
      
      ! Reset quantities depending on TEMP0 and TAUTP (which may have been
      ! changed by MODWT during FORCE call).
      ! Recalculate target kinetic energies.
      
      ekinp0 = fac(2) * temp0

#ifdef LES
      
      ! modified for LES temperature, not solvent
      
      ekins0 = 0.d0
      ekinles0 = 0.d0
      if (temp0les < 0.d0) then
         ekins0 = fac(3) * temp0
         ekin0 = fac(1) * temp0
      else
         ekinles0 = fac(3) * temp0les
         ekin0 = ekinp0 + ekinles0
      end if
#else
      ekins0 = fac(3) * temp0
      ekin0 = fac(1) * temp0
#endif

      if (ntt == 1) dttp = dt / tautp
      
      if (ntp > 0) then
         ener%volume = volume
         ener%density = tmass / (0.602204d0*volume)
         
         if( iamoeba == 0 ) then
            ener%cmt(4)  = 0.d0
            ener%vir(4)  = 0.d0
            ener%pres(4) = 0.d0
            do m = 1,3
               ener%cmt(m)  = ener%cmt(m) * 0.5d0
               ener%cmt(4)  = ener%cmt(4) + ener%cmt(m)
               ener%vir(4)  = ener%vir(4) + ener%vir(m)
               ener%pres(m) = (pconv+pconv) * (ener%cmt(m)-ener%vir(m)) / volume
               ener%pres(4) = ener%pres(4) + ener%pres(m)
            end do
            ener%pres(4) = ener%pres(4) / 3.d0
         end if
      end if
      
      ntnb = 0
      i3 = 0
      tempsu = 0.0d0
      
#ifdef LES
      ! added LES tempsu (actual LES sum of m*v**2 )
      tempsules = 0.0d0
#endif
      eke_cmd = 0.d0
      do j = 1,nrp
         winf = winv(j) * dt5
         aamass = amass(j)
         do m = 1,3
            i3 = i3+1
            rterm = v(i3)*v(i3) * aamass
#ifdef LES
            if (temp0les < 0.d0) then
               tempsu = tempsu + rterm
               if (ipimd.eq.CMD.and.(cnum(j).eq.0.or.cnum(j).eq.1)) then
                  eke_cmd = eke_cmd + aamass*v(i3)*v(i3) 
               endif
            else
               if (cnum(j) == 0) then
                  tempsu = tempsu + rterm
               else
                  tempsules = tempsules + rterm
               end if
            end if
#else
            if(ipimd.eq.CMD.and.mybeadid==1) then
               eke_cmd = eke_cmd + aamass*v(i3)*v(i3)
            end if
            tempsu = tempsu + rterm
#endif
            if(ipimd.ne.NMPIMD.and.ipimd.ne.CMD) v(i3) = v(i3) - f(i3) * winf
            if (vlim) v(i3) = sign(min(abs(v(i3)),vlimit),v(i3))
         end do
      end do

#ifdef MPI /* SOFT CORE */
      if ( ifsc /= 0 ) then
        call calc_softcore_ekin(amass,v,v,istart,iend)
        sc_ener(13) = sc_ener(6) + sc_ener(12)
      end if
#endif
      
      do im=1,iscale
         v(nr3+im) = v(nr3+im) - f(nr3+im) * dt5 / scalm
         tempsu = tempsu + scalm * v(nr3+im)*v(nr3+im)
      end do
      ener%kin%solt = tempsu * 0.5d0

#ifdef LES
      
      ! added for LES temperature using old solvent variable for ener(4)
      
      if (temp0les < 0.d0) then
         ener%kin%solv = 0.d0
         ener%kin%tot  = ener%kin%solt
         ! for CMD:
         if( ipimd > 0 ) then
            ener%kin%solv = equal_part + Epot_deriv ! "virial" estimate of KE
            ener%tot = ener%kin%solv + ener%pot%tot
         else
            ener%tot = ener%kin%tot + ener%pot%tot
         endif
         if (ipimd.eq.CMD) then
            ener%kin%tot  = eke_cmd*0.5d0
            ener%kin%solv = ener%kin%tot
         endif
      else
         ener%kin%solv = tempsules * 0.5d0
         ener%kin%tot  = ener%kin%solt + ener%kin%solv
      end if
#else
      ! for better output for parallel PIMD/NMPIM/CMD/RPMD
      if (ipimd>0) then
         ener%tot      = 0.d0
         ener%kin%tot  = 0.d0
         ener%kin%solt = 0.d0
         ener%kin%solv = 0.d0
         ener%volume   = 0.d0
      endif
      ener%kin%tot = ener%kin%solt
      ener%tot  = ener%kin%tot+ener%pot%tot

#endif

      if(ntt == 1) then
#ifdef LES
         if (temp0les >= 0.d0) then
            ekmh = max(ener%kin%solt,fac(2)*10.d0)
            ekmhles = max(ener%kin%solv,fac(3)*10.d0)
         else
            ekmh = max(ener%kin%solt,fac(1)*10.d0)
         end if
#else
         ekmh = max(ener%kin%solt,fac(1)*10.d0)
#endif
      end if
      
   end if  ! ( init == 3 )
   
   !-------------------------------------------------------------------------
   ! init = 4:  continuation of a previous trajectory
   !            this code also done for init=3
   !
   ! Note: if the last printed energy from the previous trajectory was
   !       at time "t", then the restrt file has velocities at time
   !       t + 0.5dt, and coordinates at time t + dt
   !-------------------------------------------------------------------------
   
   ! -------------------------------------------------------------------
   ekmh = 0.0d0
#ifdef LES
   ekmhles = 0.0d0
#endif

   i3 = 0
   do j = 1,nrp
      aamass = amass(j)
      do m = 1,3
         i3 = i3+1
         rterm = v(i3)*v(i3) * aamass
#  ifdef LES
         ! use copy number, not solute/solvent
         if (temp0les < 0.d0) then
            ! 1 bath
            ekmh = ekmh + rterm
         else
            if (cnum(j) == 0) then
               ekmh = ekmh + rterm
            else
               ekmhles = ekmhles + rterm
            end if
         end if
#  else
         ekmh = ekmh + rterm
#  endif
      end do
   end do

#ifdef MPI /* SOFT CORE */
   if ( ifsc /= 0 ) then
     call calc_softcore_ekin(amass,v,v,istart,iend)
     sc_ener(13) = sc_ener(6) + sc_ener(12)
   end if
#endif
   
   do im=1,iscale
      ekmh = ekmh + scalm*v(nr3+im)*v(nr3+im)
   end do
   ekmh = ekmh * 0.5d0
#ifdef LES
   ekmhles = ekmhles * 0.5d0
#endif

   do i=1,nr3+iscale
      vold(i) = v(i)
   end do

#ifdef EMIL
   !--Setup the emil calculation if required
   if ( emil_do_calc .gt. 0 ) then
       call emil_init( natom, nstep, 1.0/(temp0 * 2 * boltz2 ), &
                     mass, xx(lcrd), f, v, ener%box)
   end if
#endif
   
   if (abfqmmm_param%abfqmmm == 1) then          ! lam81
    nstep=abfqmmm_param%qmstep                   ! lam81
    if(abfqmmm_param%maxqmstep == 0) nstep = 0   ! lam81
   end if                                        ! lam81

   if (init /= 4 .or. nstlim == 0 .or. (abfqmmm_param%abfqmmm == 1 .and. abfqmmm_param%system == 1)) then  ! lam81
      
      !-------------------------------------------------------------------
      !           PRINT THE INITIAL ENERGIES AND TEMPERATURES
      !-------------------------------------------------------------------
#ifdef RISMSANDER
      if ( rismprm%irism == 1 .and. rismprm%write_thermo==1 &
           .and. nstep <= 0 .and. facc /= 'A') then
         if( rism_calc_type(0) == RISM_FULL)&
              call rism_solvdist_thermo_calc(.false.,0)
      end if
#endif /*RISMSANDER*/      

      if ( (nstep <= 0 .and. master .and. facc /= 'A') .or. &                      ! lam81
           (master .and. abfqmmm_param%abfqmmm == 1 .and. mod(abfqmmm_param%qmstep,ntpr) == 0) ) then  ! lam81

         if (isgld > 0) call sgenergy(ener)
         rewind(7)
#ifdef LES
         if (.not.ipimd.gt.0) &
           ener%tot = ener%kin%tot+ener%pot%tot
#endif /* LES */
         if(abfqmmm_param%abfqmmm /= 1 .or. abfqmmm_param%system == 1 .or. nstep == 0) &  ! lam81
         call prntmd(nstep,nitp,nits,t,ener,onefac,7,.false.)
#ifdef MPI /* SOFT CORE */
         if (ifsc /= 0) call sc_print_energies(6, sc_ener)
         if (ifsc /= 0) call sc_print_energies(7, sc_ener)
#endif
         if ( ifcr > 0 .and. crprintcharges > 0 ) then
            call cr_print_charge( xx(l15), nstep )
         end if

         !--- BEGIN DIPOLE PRINTING CODE ---
         ! See code further on for comments-explanations
         call nmlsrc('dipoles',5,prndipfind)
         if(prndipfind /= 0 ) then
            write(6,*) '------------------------------- DIPOLE &
                       &INFO ----------------------------------'
            write(6,9018) nstep,t
            read (5,'(a)') prndiptest
            call rgroup(natom,natc,nres,prndipngrp,ix(i02),ih(m02), &
                 ih(m04),ih(m06),ih(m08),ix(icnstrgp), &
                 jgroup,indx,irespw,npdec, &
                 xx(l60),xx(lcrdr),0,0,0,idecomp,5,.false.)
            rewind(5)
            if(prndipngrp > 0) then
               call printdip(prndipngrp,ix(icnstrgp),xx(lcrd), &
                    xx(l15),xx(linddip),xx(Lmass), natom)
            end if
            write(6,*) '----------------------------- END DIPOLE &
                       &INFO --------------------------------'
         end if
         !--- END DIPOLE PRINTING CODE ---

         if (nmropt > 0) then
            call nmrptx(6)
         end if
         call amflsh(7)
      end if

      if (abfqmmm_param%abfqmmm == 1 .and. abfqmmm_param%system == 1) then   ! lam81
       deallocate(for, stat=ier)                                     ! lam81
       REQUIRE(ier==0)                                               ! lam81
       return                                                        ! lam81
      end if                                                         ! lam81
      if (nstlim == 0) then                                          ! lam81
       if(abfqmmm_param%abfqmmm == 1) v(1:nr3) = abfqmmm_param%v(1:nr3) ! lam81
#ifdef MPI
       call xdist(x, xx(lfrctmp), natom)                             ! lam81
       call xdist(v, xx(lfrctmp), natom)                             ! lam81

       if(master) then                                               ! lam81 
#endif
        if(ntwr>0) call mdwrit(nstep,nrp,nr,nres,ntxo,ntr,ntb, &     ! lam81
                               x,v,xx(lcrdr),box,t,temp0)            ! lam81
        if(ntwx>0) call corpac(x,1,nrx,MDCRD_UNIT,loutfm)            ! lam81
        if(ntwv>0) call corpac(v,1,nrx,MDVEL_UNIT,loutfm)            ! lam81
        if(ntwf>0) call corpac(for,1,nrx,MDFRC_UNIT,loutfm)          ! lam81
        if(ntwe>0) call mdeng(15,nstep,t,ener,onefac,ntp,csurften)   ! lam81
#ifdef MPI
       end if                                                        ! lam81
#endif
       return                                                        ! lam81
      end if                                                         ! lam81
      init = 4
   end if


   if(ntp > 0 .and. ipimd > 0 ) then
      REQUIRE(ipimd.eq.NMPIMD)
#ifdef LES
      call part_setup_cnst_press_pimd(Nkt,tau_vol)
#else     
      call full_setup_cnst_press_pimd(Nkt,tau_vol)
#endif
      e2 = 1.0/ (2.0*3.0)
      e4 = e2 / (4.0*5.0)
      e6 = e4 / (6.0*7.0)
      e8 = e6 / (8.0*9.0)
      x_lnv = log( box(1)*box(2)*box(3) ) / 3 
   end if

   ! For CMD.
   if ( ipimd==CMD ) then

      if ( .not.eq_cmd ) then

         ! De-activate thermostat for path-centroid.
#ifdef LES
         do iatom = 1, natom
         do idim  = 1, 3
            if ( cnum(iatom)==0 .or. cnum(iatom)==1 )  then
               activate = .false.
            else
               activate = .true.
            end if
            call Thermostat_switch(thermo(idim,iatom),activate)
         enddo
         enddo
         if ( .not.restart_cmd ) then
            ! Scale path-centroid velocity and set total momentum equal to zero.
            call part_scale_vel_centroid(v,amass,istart,iend)
            nstep_cmd = 0
            t_cmd = 0.d0
         else
            t_cmd = t
            nstep_cmd = int( t / dt )
         end if
#else
         if ( mybeadid.eq.1 ) then
            activate = .false.
         else
            activate = .true.
         end if
         do iatom = 1, natom
         do idim  = 1, 3
            call Thermostat_switch(thermo(idim,iatom),activate)
         enddo
         enddo
         if ( .not.restart_cmd ) then
            ! Scale path-centroid velocity and set total momentum equal to zero.
            call full_scale_vel_centroid(v,amass,istart,iend)
            nstep_cmd = 0
            t_cmd = 0.d0
         else
            nstep_cmd = nstep
            t_cmd = t
         end if
#endif /* LES */

      else

         nstep_cmd = nstep
         t_cmd = t

      end if

   end if  ! ipimd.eq.CMD and adiab_param<1.d0

#ifdef MPI
   ! If this is a replica run and we are on exchange > 1, restore the 
   !    old ekmh value since it was reset after we left runmd last time.
   !    DAN ROE: Only for ntt==1??
   if (rem /= 0 .and. mdloop >= 1) then
      ekmh = remd_ekmh
   endif
#endif

  
   !=======================================================================
   !     ----- MAIN LOOP FOR PERFORMING THE DYNAMICS STEP -----
   !           (at this point, the coordinates are a half-step "ahead"
   !           of the velocities; the variable EKMH holds the kinetic
   !           energy at these "-1/2" velocities, which are stored in
   !           the array VOLD.)
   !=======================================================================
   
   260 continue
   onstep = mod(irespa,nrespa) == 0

   ! Constant pH setup
   if (icnstph /= 0 .and. &
       ((rem /= 0 .and. mdloop > 0) .or. rem == 0)) then

      if (ntnb == 1) then ! rebuild pairlist
            call cnstphupdatepairs(x)
      end if

      if (mod(irespa+nstlim*mdloop,ntcnstph) == 0) then
            if (icnstph .eq. 1) then
               call cnstphbeginstep(xx(l190))
            else
               call cnstph_explicitmd( xx,ix,ih,ipairs,x,winv,amass,f,v,vold, &
                                       xr,xc,conp,skip,nsp,tma,erstop,qsetup, &
                                       do_list_update,rem)
            end if
      end if

   end if
  
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  |  EVB reactive flux                                            |
!  +:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::+
!  |  Driver for coordinating backward and forward propagation as  |
!  |  well as for enforcing stopping criteria                      |
!  +---------------------------------------------------------------+

#if defined(MPI)
   if( ievb /= 0 .and. trim( adjustl( evb_dyn) ) == "react_flux" ) then
      REQUIRE( ipimd.eq.0 .or. ipimd.eq.NMPIMD )
      call react_flux ( x, v, f, winv, tempi * factt, dt5, dtx &
                      , nr, nstep, nstlim )
   endif
#endif
 
   !---------------------------------------------------------------
   !  ---Step 1a: do some setup for pressure calculations:
   !---------------------------------------------------------------

   if (ntp > 0 .and. iamoeba == 0 .and. ipimd==0) then
      ener%cmt(1:3) = 0.d0
      xr(1:nr3) = x(1:nr3)
      
      ! ----- CALCULATE THE CENTER OF MASS ENERGY AND THE COORDINATES
      !       OF THE SUB-MOLECULES WITH RESPECT TO ITS OWN CENTER OF
      !       MASS -----
      
      call timer_start(TIME_EKCMR)
      call ekcmr(nspm,nsp,tma,ener%cmt,xr,v,amass,istart,iend)
#ifdef MPI
      call trace_mpi('mpi_allreduce', &
            3,'MPI_DOUBLE_PRECISION',mpi_sum)
# ifdef USE_MPI_IN_PLACE
      call mpi_allreduce(MPI_IN_PLACE,ener%cmt,3, &
            MPI_DOUBLE_PRECISION,mpi_sum,commsander,ierr)
# else
      call mpi_allreduce(ener%cmt,mpitmp,3, &
            MPI_DOUBLE_PRECISION,mpi_sum,commsander,ierr)
      ener%cmt(1:3) = mpitmp(1:3)
# endif
#endif
      call timer_stop(TIME_EKCMR)
   end if

   ! If we're using the MC barostat, go ahead and do the trial move now
   if (ntp > 0 .and. barostat == 2 .and. mod(total_nstep+1, mcbarint) == 0) &
      call mcbar_trial(xx, ix, ih, ipairs, x, xc, f, ener%vir, xx(l96), &
               xx(l97), xx(l98), xx(l99), qsetup, do_list_update, &
               nstep, nsp, amass)

   !--------------------------------------------------------------
   !  ---Step 1b: Get the forces for the current coordinates:
   !--------------------------------------------------------------

   npbstep = nstep
   iprint = 0
   if( nstep == 0 .or. nstep+1 == nstlim ) iprint = 1

   if (sebomd_obj%do_sebomd) then
     ! write down atomic charges and density matrix if needed
     sebomd_obj%iflagch = 0
     if (sebomd_obj%ntwc /= 0) then
        if (mod(nstep+1,sebomd_obj%ntwc) == 0) sebomd_obj%iflagch = 1
     endif
!    sebomd_obj%pdmx = 0
!    if (sebomd_obj%pdump /= 0) then
!       if (mod(nstep+1,ntwr) == 0) sebomd_obj%pdmx = 1
!       if (nstep+1 == nstlim) sebomd_obj%pdmx = 1
!    endif
   endif

#ifdef MPI
   ! set do_mbar for the force contributions
   if (ifmbar /= 0) then
      do_mbar = .false.
      if ( mod(nstep+1,bar_intervall) == 0) then
         do_mbar = .true.
      end if
   end if
#endif

   if ( ipimd==NMPIMD .or. ipimd==CMD) then
      call trans_pos_nmode_to_cart(x,cartpos)
      call force(xx,ix,ih,ipairs,cartpos,f,ener,ener%vir, &
               xx(l96),xx(l97),xx(l98),xx(l99), qsetup, &
               do_list_update,nstep)

#if defined(MPI) && defined(LES)
         if ( ievb == 1 .and. i_qi > 0) then
            call evb_umb ( f, cartpos, real_mass, natom, istart3, iend3 )
! 03132009            if( i_qi == 2 ) call qi_corrf_les ( cartpos, amass )
            if( i_qi == 2 ) call qi_corrf_les ( cartpos, real_mass )
            evb_nrg(1) = evb_frc%evb_nrg
            evb_nrg(2) = evb_vel0%evb_nrg
            if( nbias > 0 ) evb_nrg(3) = sum( evb_bias%nrg_bias(:) )
         endif
#endif

      call trans_frc_cart_to_nmode(f)

#if defined(MPI) && defined(LES)
         if ( ievb /= 0 .and. i_qi == 0 ) then
            call evb_umb ( f, x, real_mass, natom, istart3, iend3 )
            evb_nrg(1) = evb_frc%evb_nrg
            evb_nrg(2) = evb_vel0%evb_nrg
            if( nbias > 0 ) evb_nrg(3) = sum( evb_bias%nrg_bias(:) )
         endif

#endif

   else
      ! -- ti decomp
      if(idecomp > 0) then
         decpr = .false.
         if(mod(nstep+1,ntpr) == 0) decpr = .true.
      end if
      call force(xx,ix,ih,ipairs,x,f,ener,ener%vir, &
               xx(l96),xx(l97),xx(l98),xx(l99), qsetup, &
               do_list_update,nstep)
#if defined(MPI)
         if ( ievb /= 0 ) then
#ifdef LES
            call evb_umb_primitive ( f, x, real_mass, natom, istart, iend )
#else
            call evb_umb_primitive ( f, x, amass, natom, istart, iend )
#endif
            evb_nrg(1) = evb_frc%evb_nrg
            evb_nrg(2) = evb_vel0%evb_nrg
            if( nbias > 0 ) evb_nrg(3) = sum( evb_bias%nrg_bias(:) )
         endif
#endif

   endif

   if (sebomd_obj%do_sebomd) then
     ! computes hessian matrix if necessary
     if (sebomd_obj%ntwh /= 0 .and. mod(nstep+1,sebomd_obj%ntwh) == 0) then
        ! don't output atomic charges
        sebomd_obj%iflagch_old = sebomd_obj%iflagch
        sebomd_obj%iflagch = 0
        call sebomd_gradient_write(f,3*natom)
        call sebomd_hessian_compute(xx,ix,ih,ipairs,x,f,ener, &
                          qsetup, do_list_update, nstep)
        sebomd_obj%iflagch = sebomd_obj%iflagch_old
     endif
   endif

   for(1:nr3) = f(1:nr3)               ! lam81
#ifdef MPI
   call xdist(for,xx(lfrctmp),natom)   ! lam81
#endif

   if(abfqmmm_param%abfqmmm == 1) then          ! lam81
     abfqmmm_param%f2(1:nr3) = for(1:nr3)       ! lam81
     call abfqmmm_combine_forces()              ! lam81
#ifdef MPI
     call mpi_bcast(abfqmmm_param%f, 3*natom, mpi_double_precision, 0, commsander, ierr)                   ! lam81
#endif 
     for(1:nr3) = abfqmmm_param%f(1:nr3)                                                                   ! lam81
     f(1:nr3) = abfqmmm_param%f(1:nr3)                                                                     ! lam81
   end if                                       ! lam81

   ! Constant pH transition evaluation for GB CpHMD (not explicit CpHMD)
   if ((icnstph == 1) .and. (mod(irespa+mdloop*nstlim,ntcnstph) == 0)) then
      call cnstphendstep(xx(l190), xx(l15), ener%pot%dvdl, temp0, solvph)
      if (master) call cnstphwrite(rem)
   end if

   !PLUMED force added
   if(plumed.eq.1) then
     plumed_stopflag=0
     call plumed_f_gcmd("setStep"//char(0),nstep)
     call plumed_f_gcmd("setPositions"//char(0),x)
     call plumed_f_gcmd("setMasses"//char(0),amass)
     call plumed_f_gcmd("setCharges"//char(0),xx(l15))
     call plumed_f_gcmd("setEnergy"//char(0),ener%pot)
     call plumed_f_gcmd("setForces"//char(0),f)
     call plumed_f_gcmd("setStopFlag"//char(0),plumed_stopflag)
     plumed_box=0.0
     if(ifbox==0) then
       continue
     else if(ifbox==1) then
       plumed_box(1,1)=box(1)
       plumed_box(2,2)=box(2)
       plumed_box(3,3)=box(3)
     else if(ifbox==2) then
! truncated octahedron, corresponding to a bcc lattice
! in AMBER convention, box(1) is the length of the lattice vector
! a is defined so as the bcc lattice is (a/2,a/2,a/2) (-a/2,-a/2,a/2) (a/2,-a/2,-a/2)
       plumed_box(1,1)=sqrt(1.0/3.0)*box(1)
       plumed_box(2,1)=sqrt(1.0/3.0)*box(1)
       plumed_box(3,1)=sqrt(1.0/3.0)*box(1)
       plumed_box(1,2)=-sqrt(1.0/3.0)*box(1)
       plumed_box(2,2)=-sqrt(1.0/3.0)*box(1)
       plumed_box(3,2)=sqrt(1.0/3.0)*box(1)
       plumed_box(1,3)=sqrt(1.0/3.0)*box(1)
       plumed_box(2,3)=-sqrt(1.0/3.0)*box(1)
       plumed_box(3,3)=-sqrt(1.0/3.0)*box(1)
     else
      write (6,*) "!!!!! PLUMED ERROR: Only orthorhombic and truncted octahedron cells are supported in this release."
      write (6,*) "!!!!! ABORTING RUN"
      stop
     endif
     plumed_virial=0.0
! It's not completely clear where the factor 2.0 comes from
! Anyway, I was able to match a change in press of 1000 bar with
! a corresponding SLOPE=66.02 added to VOLUME CV in PLUMED
! GB
     plumed_virial(1,1)=2.0*ener%vir(1)
     plumed_virial(2,2)=2.0*ener%vir(2)
     plumed_virial(3,3)=2.0*ener%vir(3)
     call plumed_f_gcmd("setVirial"//char(0),plumed_virial)
     call plumed_f_gcmd("setBox"//char(0),plumed_box)
     call plumed_f_gcmd("calc"//char(0),0);
#ifdef MPI
! this is required since PLUMED only updates virial on master processor
#ifdef DPREC
     call mpi_bcast(plumed_virial,9,MPI_DOUBLE_PRECISION,0,commsander,ierr)
#else
     call mpi_bcast(plumed_virial,9,MPI_REAL,0,commsander,ierr)
#endif
#endif
     ener%vir(1)=0.5*plumed_virial(1,1)
     ener%vir(2)=0.5*plumed_virial(2,2)
     ener%vir(3)=0.5*plumed_virial(3,3)
   end if

   !PLUMED end


#ifdef MPI
   ! If softcore potentials are used, collect their dvdl contributions:
   if ( ifsc /= 0 ) then
      call mpi_reduce(sc_dvdl, sc_tot_dvdl, 1, MPI_DOUBLE_PRECISION, &
                      MPI_SUM, 0, commsander, ierr)
      sc_dvdl=0.0d0 ! zero for next step
      call mpi_reduce(sc_dvdl_ee, sc_tot_dvdl_ee, 1, MPI_DOUBLE_PRECISION, &
                      MPI_SUM, 0, commsander, ierr)
      sc_dvdl_ee=0.0d0 ! zero for next step
      call mpi_reduce(sc_ener, sc_ener_tmp, ti_ene_cnt, MPI_DOUBLE_PRECISION, &
                      MPI_SUM, 0, commsander, ierr)
      sc_ener(1:ti_ene_cnt) = sc_ener_tmp(1:ti_ene_cnt)
   end if
   if ( ifsc == 2 ) then
      ! If this is a perturb to nothing run, scale forces and calculate dvdl
      call sc_nomix_frc(f,nr3,ener)

      if( numtasks>1 ) then
         call mpi_bcast(f,nr3,MPI_DOUBLE_PRECISION,0,commsander,ierr)
         call mpi_bcast(ener,state_rec_len,MPI_DOUBLE_PRECISION,0,commsander,ierr)
      end if
   end if

   if (ifmbar /=0 .and. do_mbar) then
      call bar_collect_cont()
   end if

   if ( icfe /= 0 )then

      ! --- free energies using thermodynamic integration (icfe /= 0)

      !  --- first, send the forces, energy, and virial to your partner:

      if( master ) then
         partner = ieor(masterrank,1)
         call mpi_sendrecv( f, nr3, MPI_DOUBLE_PRECISION, partner, 5, &
                  frcti, nr3+3*extra_atoms, MPI_DOUBLE_PRECISION, partner, 5, &
                  commmaster, ist, ierr )
         call mpi_sendrecv( ener, state_rec_len, MPI_DOUBLE_PRECISION, partner, 5, &
                            ecopy, state_rec_len, MPI_DOUBLE_PRECISION, partner, 5, &
                            commmaster, ist, ierr)

         ! exchange sc-dvdl contributions between masters:
         call mpi_sendrecv( sc_tot_dvdl, 1, MPI_DOUBLE_PRECISION, partner, 5, &
                    sc_tot_dvdl_partner, 1, MPI_DOUBLE_PRECISION, partner, 5, &
                    commmaster, ist, ierr )
         call mpi_sendrecv( sc_tot_dvdl_ee, 1, MPI_DOUBLE_PRECISION, partner, 5, &
                    sc_tot_dvdl_partner_ee, 1, MPI_DOUBLE_PRECISION, partner, 5, &
                    commmaster, ist, ierr )

         ! ---- collect statistics for free energy calculations:

         if( onstep ) then
            if( masterrank==0 ) then
               if( klambda == 1 ) then
                  edvdl = edvdl - ener + ecopy
                  edvdl_r = edvdl_r - ener + ecopy
               else
                  clfac = klambda*(1.d0 - clambda)**(klambda-1)
                  edvdl = edvdl - (ener - ecopy)*clfac
                  edvdl_r = edvdl_r - (ener - ecopy)*clfac
               end if
            else
               if( klambda == 1 ) then
                  edvdl = edvdl + ener - ecopy
                  edvdl_r = edvdl_r + ener - ecopy
               else
                  clfac = klambda*(1.d0 - clambda)**(klambda-1)
                  edvdl = edvdl + (ener - ecopy)*clfac
                  edvdl_r = edvdl_r + (ener - ecopy)*clfac
               end if
            end if
            ! This includes the sc-dvdl contribution into the vdw-part 
            !    and potential energy parts of the dvdl-statistics
            if (ifsc == 1) then
               call adj_dvdl_stat(edvdl, edvdl_r)
            end if
         end if

         ! Do energy collection for MBAR FEP runs
         if (ifmbar /= 0 .and. do_mbar) then
            call calc_mbar_energies(ener%pot%tot, ecopy%pot%tot)
         end if

         if( masterrank==0 ) then
            call mix_frcti(frcti,ecopy,f,ener,nr3,clambda,klambda)
         else
            call mix_frcti(f,ener,frcti,ecopy,nr3,clambda,klambda)
         endif
      endif

      if( numtasks>1 ) then
         call mpi_bcast(f,nr3,MPI_DOUBLE_PRECISION,0,commsander,ierr)
         call mpi_bcast(ener,state_rec_len,MPI_DOUBLE_PRECISION,0,commsander,ierr)
      end if

   end if  ! ( icfe /= 0 )

#endif /* MPI */

#ifdef EMIL
   ! Call the EMIL absolute free energy calculation.
   if ( emil_do_calc .gt. 0 ) then

        call emil_step(natom, nstep, 1.0 / (temp0 * 2 * boltz2),&
                            mass, xx(lcrd), f, v, ener%pot, ener%pot, ener%box)
   end if
#endif



   ! Reset quantities depending on TEMP0 and TAUTP (which may have been
   ! changed by MODWT during FORCE call).
   ekinp0 = fac(2)*temp0

#ifdef LES
   ! TEMP0LES may have changed too
   
   ekinles0=0.d0
   ekins0=0.d0
   if (temp0les >= 0.d0) then
      ekinles0 = fac(3)*temp0les
      ekin0 = ekinp0 + ekinles0
   else
      ekins0 = fac(3)*temp0
      ekin0 = fac(1)*temp0
   end if
#else
   ekins0 = fac(3)*temp0
   ekin0 = fac(1)*temp0
#endif

   if (ntt == 1) dttp = dt/tautp
   
   !  Pressure coupling:
   if (ntp > 0.and.ipimd>0) then
      REQUIRE(ipimd.eq.NMPIMD)
      centvir=0.0

#ifdef LES
      do iatom=istart,iend
         if(cnum(iatom).eq.0.or.cnum(iatom).eq.1) then
            centvir=centvir-x(3*iatom-2)*f(3*iatom-2)
            centvir=centvir-x(3*iatom-1)*f(3*iatom-1)
            centvir=centvir-x(3*iatom  )*f(3*iatom)
         end if
      end do
#else
      if(mybeadid.eq.1) then
         do iatom=istart,iend
            centvir=centvir-x(3*iatom-2)*f(3*iatom-2)
            centvir=centvir-x(3*iatom-1)*f(3*iatom-1)
            centvir=centvir-x(3*iatom  )*f(3*iatom)
         end do         
      end if
#endif /* LES */

      if(iamoeba.eq.1) then
         atomvir=sum(ener%vir(1:3))
#ifdef MPI
#  ifdef USE_MPI_IN_PLACE
         call mpi_allreduce(MPI_IN_PLACE,centvir,1,MPI_DOUBLE_PRECISION, &
                            mpi_sum,commworld,ierr)
         call mpi_allreduce(MPI_IN_PLACE,atomvir,1,MPI_DOUBLE_PRECISION, &
                            mpi_sum,commworld,ierr)
#  else
         call mpi_allreduce(centvir,mpitmp,1,MPI_DOUBLE_PRECISION, &
                            mpi_sum,commworld,ierr)
         centvir=mpitmp(1)
         tmp=0.0
         call mpi_allreduce(atomvir,tmp,1,MPI_DOUBLE_PRECISION, &
                            mpi_sum,commworld,ierr)
         atomvir=tmp
#  endif
#endif
      else
#ifdef MPI
#  ifdef USE_MPI_IN_PLACE
         call mpi_allreduce(MPI_IN_PLACE,centvir,1,MPI_DOUBLE_PRECISION, &
                            mpi_sum,commworld,ierr)
         call mpi_allreduce(MPI_IN_PLACE,bnd_vir,9,MPI_DOUBLE_PRECISION, &
                            mpi_sum,commworld,ierr)
         call mpi_allreduce(MPI_IN_PLACE,e14vir,9,MPI_DOUBLE_PRECISION, &
                            mpi_sum,commworld,ierr)
#    ifndef LES
         if (master) &
         call mpi_allreduce(MPI_IN_PLACE,atvir,9,MPI_DOUBLE_PRECISION, &
                            mpi_sum,commmaster,ierr)
#    endif
#  else  
         call mpi_allreduce(centvir,tmp,1,MPI_DOUBLE_PRECISION, &
                            mpi_sum,commworld,ierr)
         centvir=tmp
         tmpvir=0.0
         call mpi_allreduce(bnd_vir,tmpvir,9,MPI_DOUBLE_PRECISION, &
                            mpi_sum,commworld,ierr)
         bnd_vir=tmpvir

#    ifndef LES
         if (master) then
            tmpvir=0.0
            call mpi_allreduce(e14vir,tmpvir,9,MPI_DOUBLE_PRECISION, &
                               mpi_sum,commmaster,ierr)
            e14vir=tmpvir

            tmpvir=0.0
            call mpi_allreduce(atvir,tmpvir,9,MPI_DOUBLE_PRECISION, &
                            mpi_sum,commmaster,ierr)
            atvir=tmpvir
         endif
#    else
         tmpvir=0.0
         call mpi_allreduce(e14vir,tmpvir,9,MPI_DOUBLE_PRECISION, &
                            mpi_sum,commworld,ierr)
         e14vir=tmpvir
#    endif
#  endif
         call mpi_bcast(atvir,9,MPI_DOUBLE_PRECISION,0,commsander,ierr)
         call mpi_bcast(e14vir,9,MPI_DOUBLE_PRECISION,0,commsander,ierr)
#endif
         atomvir=0.0
         atomvir=atomvir+atvir(1,1)+bnd_vir(1,1)+e14vir(1,1)
         atomvir=atomvir+atvir(2,2)+bnd_vir(2,2)+e14vir(2,2)
         atomvir=atomvir+atvir(3,3)+bnd_vir(3,3)+e14vir(3,3)
      end if
      pressure = (Nkt*3.0-centvir-(atomvir-Eimp_virial))/(3.0*volume)
      f_lnv_p = (pressure-pres0/pconv)*volume*3.0
   end if


   if (ntp > 0) then
      ener%volume = volume
      ener%density = tmass / (0.602204d0*volume)
      if( iamoeba == 0 .and. ipimd==0 ) then
         ener%cmt(4) = 0.d0
         ener%vir(4) = 0.d0
         ener%pres(4) = 0.d0
         do m = 1,3
            ener%cmt(m)  = ener%cmt(m)*0.5d0
            ener%cmt(4)  = ener%cmt(4)+ener%cmt(m)
            ener%vir(4)  = ener%vir(4)+ener%vir(m)
            ener%pres(m) = (pconv+pconv)*(ener%cmt(m)-ener%vir(m))/volume
            ener%pres(4) = ener%pres(4)+ener%pres(m)
         end do
         ener%pres(4) = ener%pres(4)/3.d0

         ! Constant surface tension output:
 
         if (csurften > 0) then

            if (csurften == 1) then        ! Surface tension in the x direction
               ener%surface_ten = &
                    box(1) * (ener%pres(1) - 0.5d0 * &
                    (ener%pres(2) + ener%pres(3))) / (ninterface * ten_conv)

            else if (csurften .eq. 2) then ! Surface tension in the y direction
               ener%surface_ten = &
                    box(2) * (ener%pres(2) - 0.5d0 * &
                    (ener%pres(1) + ener%pres(3))) / (ninterface * ten_conv)
 
            else ! if (csurften .eq. 3)    ! Surface tension in the z direction
               ener%surface_ten = &
                    box(3) * (ener%pres(3) - 0.5d0 * &
                    (ener%pres(1) + ener%pres(2))) / (ninterface * ten_conv)

           end if

         end if

      end if
   end if

#ifdef MPI
! ------====== REMD ======------ 
! If rem /= 0 and mdloop == 0, this is the first sander call and we don't want to 
!   actually do any MD or change the initial coordinates.
! Exit here since we only wanted to get the potential energy for the first 
!  subrem exchange probability calc.
   if (rem /= 0 .and. mdloop == 0) then
#  ifdef VERBOSE_REMD
      if (master) write (6,'(a,i3)') &
         'REMD: Exiting runmd after getting initial energies for replica',repnum
#  endif
      goto 480 ! Go to the end of the runmd loop.
   endif  ! (rem /= 0 and mdloop == 0)

   !REB Do adaptive QMMM
   if ( qmmm_nml%vsolv > 1 ) then
      ! mix forces for adaptive QM/MM and
      ! calculate adaptive energy if requested
      ! note: nstep is zero during first call; this is the energy/force calculation
      !       with the starting geometry / velocities
      call adaptive_qmmm(nstep,natom,x,xold,f,ener%pot%tot, ntpr, ntwx, &
                xx, ix, ih, ipairs, qsetup, do_list_update, &
                corrected_energy, aqmmm_flag)

! ALTERNATIVE APPROACH:
!      if (ad_qmmm%calc_wbk) then
!         call ad_qmmm_check_matching_partitions()
!         if (ad_qmmm%mismatch) then
!            call force()
!         end if
!         call ad_qmmm_energy()
!      end if

! test
      i3 = 3*(istart-1)
      do j=istart,iend
         do idim = 1, 3
            xold(i3+idim)=x(i3+idim)
         enddo
         i3 = i3 + 3
      enddo
! test
   endif

#endif

   !----------------------------------------------------------------
   !  ---Step 1c: do randomization of velocities, if needed:
   !----------------------------------------------------------------
   ! ---Assign new random velocities every Vrand steps, if ntt=2

   resetvelo=.false.
   if (vrand /= 0 .and. ntt == 2) then
      if (mod((nstep+1),vrand) == 0) resetvelo=.true.
   end if

#ifdef MMTSB
   if ( mmtsb_switch == mmtsb_temp_rex .and. mmtsb_is_exchanged )  &
      resetvelo = .true.
#endif
   
   if (resetvelo) then
      ! DAN ROE: Why are only the masters doing this? Even if the velocities 
      !  are broadcast to the child processes, the wont the different # of random
      !  calls put the randomg num generators out of sync, or do we not care?
      
      if (master) then
         write (6,'(a,i8)') 'Setting new random velocities at step ', &
               nstep + 1
         call setvel(nr,v,winv,temp0*factt,init,iscale,scalm)

#ifdef MPI /* SOFT CORE */
         ! Make sure all common atoms have the same v (that of V0) in TI runs:
         if (icfe /=0 .and. ifsc /=0) call sc_sync_x(v,nr3)
#endif

#ifdef LES
     
         ! newvel call is fixed for the dual target temperatures

         if (temp0les >= 0.d0.and.temp0 /= temp0les) then
            vscalt = sqrt (temp0les/temp0) 
            do j=1,natom
              if(cnum(j) > 0) then
                i3 = 3*(j-1)
                v(i3+1) = v(i3+1) * vscalt
                v(i3+2) = v(i3+2) * vscalt
                v(i3+3) = v(i3+3) * vscalt
              endif
            end do
         end if
#endif
         if (ibelly > 0) call bellyf(nr,ix(ibellygp),v)
      end if
# ifdef MPI
      call trace_mpi('mpi_bcast',3*natom,'MPI_DOUBLE_PRECISION',0)
      call mpi_bcast(v, 3*natom, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
# endif

      ! At this point in the code, the velocities lag the positions
      ! by half a timestep.  If we intend for the velocities to be drawn 
      ! from a Maxwell distribution at the timepoint where the positions and 
      ! velocities are synchronized, we have to correct these newly 
      ! redrawn velocities by backing them up half a step using the 
      ! current force.
      ! Note that this fix only works for Newtonian dynamics.
      if( gammai==0.d0.and.(ipimd.ne.NMPIMD.or.ipimd.ne.CMD)) then
         i3 = 3*(istart-1)
         do j=istart,iend
            wfac = winv(j) * dt5
            v(i3+1) = v(i3+1) - f(i3+1)*wfac
            v(i3+2) = v(i3+2) - f(i3+2)*wfac
            v(i3+3) = v(i3+3) - f(i3+3)*wfac
            i3 = i3+3
         end do
      end if
      
   end if  ! (resetvelo)

   call timer_start(TIME_VERLET)

   !-----------------------------------------------------
   !  ---Step 2: Do the velocity update:
   !-----------------------------------------------------
   
   !step 2a: apply quenched MD if needed.  This is useful in NEB>0
   if (vv==1) call quench(f,v)
   
   !  Car-Parrinello on dipoles:  note that the (small?) kinetic energy
   !       of the dipoles is included in the epol energy
! M_WJ
!  if ( induced == 1 .and. indmeth == 3 ) call cp_dips(natom,xx(lpol),xx,dt)
   if ( induced > 0 .and. indmeth == 3 ) call cp_dips(natom,xx(lpol),xx,dt)



!   i3 = 3*(istart-1) !! Add Brownian noise for testing.        ! APJ
!   do j=istart,iend                                            ! APJ
!      do idim=1,3                                              ! APJ
!         call gauss( 0.d0, sqrt(0.1d0*boltz2*temp0)/dtx,fln )  ! APJ
!         f(i3+idim) = f(i3+idim) + fln                         ! APJ
!      enddo                                                    ! APJ
!      i3 = i3+3                                                ! APJ
!   end do                                                      ! APJ

   
   ! Nose'-Hoover thermostat (1st step).
   if ( ntt == 4 ) then

      Ekin2_tot = 0.d0
      i3 = 3*(istart-1)
      do j=istart,iend
         wfac = dtx/amass(j)
         do idim = 1, 3
#ifdef LES
            if( ntp>0.and.ipimd.eq.NMPIMD .and. &
               (cnum(j).eq.0.or.cnum(j).eq.1) ) then
#else
            if(ntp>0.and.ipimd.eq.NMPIMD.and.mybeadid.eq.1) then
#endif
               exp1 = exp(-dt5*thermo(idim,j)%v(1)-dt5*v_lnv*c2_lnv)
               Ekin2_tot = Ekin2_tot + amass(j)*v(i3+idim)*v(i3+idim)
            else
               exp1 = exp( -dt5 * thermo(idim,j)%v(1) )
            end if
            exp2 = exp1*exp1
            vold(i3+idim)=v(i3+idim)
            v(i3+idim) = v(i3+idim) * exp2 + f(i3+idim) * wfac * exp1
         end do
         i3 = i3+3
      end do

      if(ntp>0.and.ipimd.eq.NMPIMD) then
#ifdef MPI
#  ifdef USE_MPI_IN_PLACE
         call mpi_allreduce(MPI_IN_PLACE,Ekin2_tot,1,MPI_DOUBLE_PRECISION, &
                            mpi_sum,commworld,ierr)
#  else
         call mpi_allreduce(Ekin2_tot,mpitmp,1,MPI_DOUBLE_PRECISION, &
                            mpi_sum,commworld,ierr)
         Ekin2_tot=mpitmp(1)
#  endif
#endif
         f_lnv_v = Ekin2_tot*(c2_lnv-1)
         tmp = exp(-dt5*thermo_lnv%v(1))
         v_lnv = tmp*(tmp*v_lnv+dtx*(f_lnv_v+f_lnv_p)/mass_lnv)
      end if

      call Thermostat_integrate_1(nchain,thermo,nthermo,dtx,ntp)

   else if ( ntt>4 .and. ntt<=8 ) then      ! APJ

      Ekin2_tot = 0.d0
      i3 = 3*(istart-1)
      do j=istart,iend
         wfac = dtx/amass(j)
         do idim = 1, 3
#ifdef LES
            if( ntp>0.and.ipimd.eq.NMPIMD .and. &
               (cnum(j).eq.0.or.cnum(j).eq.1) ) then
#else
            if(ntp>0.and.ipimd.eq.NMPIMD.and.mybeadid.eq.1) then
#endif
               Ekin2_tot = Ekin2_tot + amass(j)*v(i3+idim)*v(i3+idim)
               !exp1 = exp(-dt5*thermo(idim,j)%v(1)-dt5*v_lnv*c2_lnv)  ! APJ
               exp1 = exp(-dt5*v_lnv*c2_lnv)                           ! APJ
            else
               !exp1 = exp( -dt5 * thermo(idim,j)%v(1) )               ! APJ
               exp1 = 1.d0                                             ! APJ
            end if
            exp2 = exp1*exp1
            vold(i3+idim)=v(i3+idim)
            !v(i3+idim) = v(i3+idim) * exp2 + f(i3+idim) * wfac * exp1 ! APJ
            v(i3+idim)=v(i3+idim)*exp2                                 ! APJ
            f(i3+idim)=f(i3+idim)*exp1                                 ! APJ
         end do
         i3 = i3+3
      end do

      if(ntp>0.and.ipimd.eq.NMPIMD) then
#ifdef MPI
#  ifdef USE_MPI_IN_PLACE
         call mpi_allreduce(MPI_IN_PLACE,Ekin2_tot,1,MPI_DOUBLE_PRECISION, &
                            mpi_sum,commworld,ierr)
#  else
         call mpi_allreduce(Ekin2_tot,mpitmp,1,MPI_DOUBLE_PRECISION, &
                            mpi_sum,commworld,ierr)
         Ekin2_tot=mpitmp(1)
#  endif
#endif
         f_lnv_v = Ekin2_tot*(c2_lnv-1)
         !tmp = exp(-dt5*thermo_lnv%v(1))                         ! APJ
         !v_lnv = tmp*(tmp*v_lnv+dtx*(f_lnv_v+f_lnv_p)/mass_lnv)  ! APJ
      end if

      if (abfqmmm_param%abfqmmm == 1) then                                           ! lam81
#ifdef MPI
       call xdist(v, xx(lfrctmp), natom)                                             ! lam81
       call xdist(f, xx(lfrctmp), natom)                                             ! lam81
#endif
       abfqmmm_param%v(1:nr3+iscale)=v(1:nr3+iscale)                                 ! lam81
       abfqmmm_param%f(1:nr3+iscale)=f(1:nr3+iscale)                                 ! lam81
      end if                                                                         ! lam81
      call Adaptive_Thermostat_integrate(nchain,thermo,nthermo,dtx,ntp,1)  ! APJ
      if (abfqmmm_param%abfqmmm == 1) then                                           ! lam81
       v(1:nr3+iscale)=abfqmmm_param%v(1:nr3+iscale)                                 ! lam81
#ifdef MPI
       call xdist(v, xx(lfrctmp), natom)                                             ! lam81
#endif
       abfqmmm_param%v(1:nr3+iscale)=v(1:nr3+iscale)                                 ! lam81
      end if                                                                         ! lam81


   else if( gammai == 0.d0 ) then

      !       ---Newtonian dynamics:
      
      ! Applying guiding force effect:
      if (isgld > 0) then
         call sgmdw(natom,istart,iend,ntp,dtx,ener,amass,winv,x,f,v)
      end if  

      i3 = 3*(istart-1)
      do j=istart,iend
         wfac = winv(j) * dtx
         v(i3+1) = v(i3+1) + f(i3+1)*wfac
         v(i3+2) = v(i3+2) + f(i3+2)*wfac
         v(i3+3) = v(i3+3) + f(i3+3)*wfac
         i3 = i3+3
      end do
      
   else if (isgld > 0) then
      !  Using SGLD algorithm:
      call sgldw(natom,istart,iend,ntp,dtx,temp0,ener,amass,winv,x,f,v)
   else  !  gamma_ln .ne. 0, which also implies ntt=3 (see mdread.f)

      !       ---simple model for Langevin dynamics, basically taken from
      !          Loncharich, Brooks and Pastor, Biopolymers 32:523-535 (1992),
      !          Eq. 11.  (Note that the first term on the rhs of Eq. 11b
      !          should not be there.)

      !  Update Langevin parameters, since temp0 might have changed:
         sdfac = sqrt( 4.d0*gammai*boltz2*temp0/dtx )
#  ifdef LES
         sdfacles = sqrt( 4.d0*gammai*boltz2*temp0les/dtx )
#  endif


#ifdef MPI /* SOFT CORE */
      if (ifsc == 1) then
         call sc_lngdyn(winv,amass,v,f,sdfac,c_explic,c_implic, &
                        istart, iend, nr, dtx)
      else
#endif


      if (no_ntt3_sync == 1) then                                               ! APJ
         !We don't worry about synchronizing the random number stream           ! APJ
         !across processors.                                                    ! APJ
         iskip_start = 0                                                        ! APJ
         iskip_end = 0                                                          ! APJ
      else                                                                      ! APJ
         ! In order to generate the same sequence of pseudorandom numbers that  ! APJ
         ! you would using a single processor you have to go through the atoms  ! APJ
         ! in order: skip those that have are being used on other processors    ! APJ
         iskip_start = 3*(istart-1)                                             ! APJ
         iskip_end = 3*(nr-iend)                                                ! APJ
#ifndef LES
         ! Always sync random number stream for PIMD          ! APJ
         ! (AWG: not sure if this is required)                ! APJ
         if (ipimd>0) then                                    ! APJ
            iskip_start = iskip_start + 3*nr*(mybeadid-1)     ! APJ
            iskip_end = iskip_end + 3*nr*(nbead-mybeadid)     ! APJ
         end if                                               ! APJ
#endif
      endif                                                   ! APJ

      do j=1,iskip_start                                      ! APJ
         ! Skip some random numbers                           ! APJ
         call gauss( 0.d0, 1.d0, fln )                        ! APJ
      end do                                                  ! APJ

      ! Do Langevin step      ! APJ
      i3 = 3*(istart-1)       ! APJ
      do j=istart,iend        ! APJ

         wfac = winv(j) * dtx  ! APJ
         aamass = amass(j)     ! APJ
#  ifdef LES
         if (temp0les >= 0 .and. temp0 /= temp0les .and. cnum(j) /= 0 ) then  ! APJ
            rsd =sdfacles*sqrt(aamass)                                        ! APJ
         else                                                                 ! APJ
            rsd = sdfac*sqrt(aamass)                                          ! APJ
         endif                                                                ! APJ
#  else
         rsd = sdfac*sqrt(aamass)                                             ! APJ
#  endif
         call gauss( 0.d0, rsd, fln )                                         ! APJ
         v(i3+1) = (v(i3+1)*c_explic + (f(i3+1)+fln)*wfac) * c_implic         ! APJ
         call gauss( 0.d0, rsd, fln )                                         ! APJ
         v(i3+2) = (v(i3+2)*c_explic + (f(i3+2)+fln)*wfac) * c_implic         ! APJ
         call gauss( 0.d0, rsd, fln )                                         ! APJ
         v(i3+3) = (v(i3+3)*c_explic + (f(i3+3)+fln)*wfac) * c_implic         ! APJ
         
         i3 = i3 + 3                                                          ! APJ
      end do                                                                  ! APJ
      
      do j=1,iskip_end                       ! APJ
         ! Skip some random numbers          ! APJ
         call gauss( 0.d0, 1.d0, fln )       ! APJ
      end do                                 ! APJ


#ifdef MPI /* SOFT CORE */
      end if ! for (ifsc==1) call sc_lngdyn 
#endif
   end if  ! ( gammai == 0.d0 )

   ! Update EMAP rigid domains 
   if(temap) call emap_move()

   !     --- consider vlimit
   
   if (vlim.and.ipimd==0) then
      vmax = 0.0d0
      do i=istart3,iend3
         vmax = max(vmax,abs(v(i)))
         v(i) = sign(min(abs(v(i)),vlimit),v(i))
      end do

      ! Only violations on the master node are actually reported
      ! to avoid both MPI communication and non-master writes.
      if (vmax > vlimit) then
         if (master) then
            write(6,'(a,i6,a,f10.4)') 'vlimit exceeded for step ',nstep, &
              '; vmax = ',vmax
         end if
      end if
   end if
   
   do im=1,iscale
      v(nr3+im) = (v(nr3+im) + f(nr3+im)*dtx/scalm)
   end do

   ! We do the force dump here if requested, since the 'old' positions are about
   ! to be dumped into the force array...

   if (master) then
      ifdump = .false.             ! Write forces this step?
      if (ntwf>0) ifdump = mod(total_nstep+1,ntwf) == 0 ! forces
      if (ntwf == -1 .and. mod(total_nstep+1,ntwx) == 0) &
         ifdump = .true. !Combined crdfrc file
      if (abfqmmm_param%abfqmmm == 1) ifdump = .false.  ! lam81
#ifdef MPI
      ! For adaptive QM/MM, only the master does a dump.
      if ( qmmm_nml%vsolv > 1 ) then
         if ( nodeid /= 0 ) then
            ifdump = .false.
         end if
      end if

      if (ifdump) then
         call xdist(f, xx(lfrctmp), natom)
      end if
#endif
      ! Force archive:
      if (ifdump) then

#ifdef MPI
         ! Write out current replica#, exchange#, step#, and mytargettemp
         ! If mdloop==0 this is a normal md run (since REMD never calls corpac
         !  when mdloop==0) and we don't want the REMD header.
         if (mdloop>0.and.loutfm) then
            if (trxsgld) then
               write (MDFRC_UNIT,'(a,4(1x,i8))') "RXSGLD ", repnum, mdloop, &
                              total_nstep+1, stagid
            else
               write (MDFRC_UNIT,'(a,3(1x,i8),1x,f8.3)') "REMD ", repnum, mdloop, &
                              total_nstep+1, my_remd_data%mytargettemp
            end if
         end if
#endif

          ! ipimd forces will probably not be right if some type of
          ! transformation is necessary. This is from the vel dump code -- keep
          ! it here as a holder in case somebody wants to fix it.
!         if(ipimd.eq.NMPIMD.or.ipimd.eq.CMD) then
!            call corpac(cartvel,1,nrx,MDVEL_UNIT,loutfm)
!         else
             call corpac(f,1,nrx,MDFRC_UNIT,loutfm)
!         endif
      end if

   else ! slaves need to participate in force distribution
      ifdump = .false.             ! Write forces this step?
      if (ntwf>0) ifdump = mod(total_nstep+1,ntwf) == 0 ! forces
      if (ntwf == -1 .and. mod(total_nstep+1,ntwx) == 0) &
         ifdump = .true. !Combined crdfrc file
      if (abfqmmm_param%abfqmmm == 1) ifdump = .false.  ! lam81
#ifdef MPI
      if (ifdump) call xdist(f, xx(lfrctmp), natom)
#endif
   end if ! master
   !-------------------------------------------------------------------
   !   Step 3: update the positions, putting the "old" positions into F:
   !-------------------------------------------------------------------

# ifdef LES
   if(ntp>0.and.ipimd.eq.NMPIMD) then
      aa = exp(dt5*v_lnv)
      arg2 = v_lnv*dt5*v_lnv*dt5
      poly = 1.0d0+arg2*(e2+arg2*(e4+arg2*(e6+arg2*e8)))
   endif

   i3 = 3*(istart-1)
   do j=istart,iend
      if(ntp>0.and.ipimd.eq.NMPIMD.and.(cnum(j).eq.0.or.cnum(j).eq.1)) then
         do idim = 1, 3
            f(i3+idim)=x(i3+idim)
            x(i3+idim)=aa*(x(i3+idim)*aa+v(i3+idim)*poly*dtx)
         enddo
      else
         do idim = 1, 3
            f(i3+idim) = x(i3+idim)
            x(i3+idim) = x(i3+idim)+v(i3+idim)*dtx
         enddo
      endif
      i3 = i3 + 3
   enddo

# else

   if(ntp>0.and.ipimd.eq.NMPIMD.and.mybeadid==1) then
      aa = exp(dt5*v_lnv)
      arg2 = v_lnv*dt5*v_lnv*dt5
      poly = 1.0d0+arg2*(e2+arg2*(e4+arg2*(e6+arg2*e8)))
      do i3=istart3,iend3
         f(i3)=x(i3)
         x(i3)=aa*(x(i3)*aa+v(i3)*poly*dtx)
      end do
   else
      do i3 = istart3, iend3
         f(i3) = x(i3)
         x(i3) = x(i3) + v(i3)*dtx
      end do
   end if

# endif /* LES */

   !Nose'-Hoover thermostat (2nd step).
   if ( ntt==4 ) then
      call Thermostat_integrate_2(nchain,thermo,nthermo,dtx,ntp)
      E_nhc = Thermostat_hamiltonian(nchain,thermo,nthermo)
   else if ( ntt>=4 .and. ntt<=8 ) then                                  ! APJ
      if(abfqmmm_param%abfqmmm == 1) then                                           ! lam81
#ifdef MPI
       call xdist(v, xx(lfrctmp), natom)                                            ! lam81
#endif
       abfqmmm_param%v(1:nr3+iscale)=v(1:nr3+iscale)                                ! lam81
      end if                                                                        ! lam81
      call Adaptive_Thermostat_integrate(nchain,thermo,nthermo,dtx,ntp,2)    ! APJ
      if (abfqmmm_param%abfqmmm == 1) then                                          ! lam81
       v(1:nr3+iscale)=abfqmmm_param%v(1:nr3+iscale)                                ! lam81
#ifdef MPI
       call xdist(v, xx(lfrctmp), natom)                                            ! lam81
#endif
       abfqmmm_param%v(1:nr3+iscale)=v(1:nr3+iscale)                                ! lam81
      end if                                                                        ! lam81
      E_nhc = Adaptive_Thermostat_hamiltonian(nchain,thermo,nthermo)
   end if

   do i = 1,iscale
      f(nr3+i) = x(nr3+i)
      x(nr3+i) = x(nr3+i)+v(nr3+i)*dtx
   end do

   call timer_stop(TIME_VERLET)
   
   if (ntc /= 1) then

      !-------------------------------------------------------------------
      !   Step 4a: if shake is being used, update the new positions to fix
      !      the bond lengths.
      !-------------------------------------------------------------------
   
      call timer_start(TIME_SHAKE)
      if (isgld > 0) call sgfshake(istart,iend,dtx,amass,x,.false.)   
      qspatial=.false.
      call shake(nrp,nbonh,nbona,0,ix(iibh),ix(ijbh),ix(ibellygp), &
      winv,conp,skip,f,x,nitp,belly,ix(iifstwt),ix(noshake), &
      shkh,qspatial)
      call quick3(f,x,ix(iifstwr),natom,nres,ix(i02))
      if(nitp == 0) then
         erstop = .true.
         goto 480
      end if
      !  Including constraint forces in self-guiding force calculation
      if (isgld > 0) call sgfshake(istart,iend,dtx,amass,x,.true.)   

      !  Need to synchronize coordinates for linearly scaled atoms after shake
#ifdef MPI
      if( icfe /= 0 ) then
         call timer_barrier( commsander )
         call timer_stop_start(TIME_SHAKE,TIME_DISTCRD)
         if ( .not. mpi_orig .and. numtasks > 1 ) then
            call xdist(x, xx(lfrctmp), natom)
         end if
         ! In dual-topology this is done within softcore.f
         if (ifsc /= 1) then
            if( master ) call mpi_bcast(x,nr3,MPI_DOUBLE_PRECISION, &
                                        0,commmaster,ierr)
         else
            if( master ) call sc_sync_x(x,nr3) 
         end if
         if( numtasks>1 ) call mpi_bcast(x,nr3,MPI_DOUBLE_PRECISION, &
                                         0,commsander,ierr)
         call timer_stop_start(TIME_DISTCRD,TIME_SHAKE)
      end if
#endif  /* MPI */
      !-----------------------------------------------------------------
      !   Step 4b:   Now fix the velocities and calculate KE
      !-----------------------------------------------------------------
      
      !  ---re-estimate the velocities from differences in positions:

      if( .not.(ipimd==NMPIMD.and.ipimd==CMD.and.mybeadid.ne.1) ) then
         v(istart3:iend3) = (x(istart3:iend3)-f(istart3:iend3)) * dtxinv
      end if

      call timer_stop(TIME_SHAKE)
   end if
   call timer_start(TIME_VERLET)

   if(ineb>0.and.(mybeadid==1.or.mybeadid==neb_nbead) ) then
      x(1:3*natom)=f(1:3*natom)
     ! CARLOS: NEB- remove velocities but ONLY for the end beads so V doesn't
     ! accumulate if high forces
      v(1:3*natom)=0.d0

   end if

   if( ntt == 1 .or. onstep ) then

      !-----------------------------------------------------------------
      !   Step 4c: get the KE, either for averaging or for Berendsen:
      !-----------------------------------------------------------------

      eke = 0.d0
      ekph = 0.d0
      ekpbs = 0.d0
#ifdef LES
      ekeles = 0.d0
      ekphles = 0.d0 
#endif
      eke_cmd = 0.d0

      if (gammai == 0.0d0) then
         i3 = 3*(istart-1)
         do j=istart,iend
            aamass = amass(j)
            do m = 1,3
               i3 = i3+1
#ifdef LES
               if (temp0les < 0.d0) then
                  eke = eke + aamass*0.25d0*(v(i3)+vold(i3))**2
                  ekph = ekph + aamass*v(i3)**2
                  if(ipimd.eq.CMD.and.(cnum(j).eq.0.or.cnum(j).eq.1)) then
                    eke_cmd = eke_cmd + aamass*0.25d0*(v(i3)+vold(i3))**2
                  endif
               else
                  if (cnum(j) == 0) then
                     eke = eke + aamass*0.25d0*(v(i3)+vold(i3))**2
                     ekph = ekph + aamass*v(i3)**2
                  else
                     ekeles = ekeles + aamass*0.25d0*(v(i3)+vold(i3))**2
                     ekphles = ekphles + aamass*v(i3)**2
                  end if
               end if

#else
               eke = eke + aamass*0.25d0*(v(i3)+vold(i3))**2

               if(mybeadid==1) then
                  eke_cmd = eke_cmd + aamass*0.25d0*(v(i3)+vold(i3))**2
               end if
               ! try pseudo KE from Eq. 4.7b of Pastor, Brooks & Szabo,
               ! Mol. Phys. 65, 1409-1419 (1988):

               ekpbs = ekpbs + aamass*v(i3)*vold(i3)
               ekph = ekph + aamass*v(i3)**2

#endif
            end do
         end do

      else

         i3 = 3*(istart-1)
         do j=istart,iend
            aamass = amass(j)
            do m = 1,3
               i3 = i3+1
#ifdef LES
               if (temp0les < 0.d0) then
                  eke = eke + aamass*0.25d0*c_ave*(v(i3)+vold(i3))**2
               else
                  if (cnum(j) == 0) then
                     eke = eke + aamass*0.25d0*c_ave*(v(i3)+vold(i3))**2
                  else
                     ekeles = ekeles + aamass*0.25d0*c_ave*(v(i3)+vold(i3))**2
                  end if
               end if
#else
               eke = eke + aamass*0.25d0*c_ave*(v(i3)+vold(i3))**2

#endif
            end do

         end do

      end if ! (if gammai == 0.0d0)

#ifdef MPI
      
      !  ---   sum up the partial kinetic energies:

      if ( ipimd.eq.CMD ) then
         call mpi_reduce(eke_cmd,tmp_eke_cmd,1,MPI_DOUBLE_PRECISION, &
              mpi_sum,0,commsander,ierr)
         eke_cmd = tmp_eke_cmd
      endif
      
#  ifdef LES
      !if ( ipimd.eq.CMD ) then
      !   call mpi_reduce(eke_cmd,tmp_eke_cmd,1,MPI_DOUBLE_PRECISION, &
      !        mpi_sum,0,commsander,ierr)
      !   eke_cmd = tmp_eke_cmd
      !endif
      if ( .not. mpi_orig .and. numtasks > 1 ) then
        if ( temp0les < 0 ) then
          mpitmp(1) = eke
          mpitmp(2) = ekph
#    ifdef USE_MPI_IN_PLACE
          call mpi_allreduce(MPI_IN_PLACE,mpitmp,2, &
               MPI_DOUBLE_PRECISION,mpi_sum,commsander,ierr)
          eke = mpitmp(1)
          ekph = mpitmp(2)
#    else
          call mpi_allreduce(mpitmp,mpitmp(3),2, &
               MPI_DOUBLE_PRECISION,mpi_sum,commsander,ierr)
          eke = mpitmp(3)
          ekph = mpitmp(4)
#    endif
        else
          mpitmp(1) = eke
          mpitmp(2) = ekph
          mpitmp(3) = ekeles
          mpitmp(4) = ekphles
#    ifdef USE_MPI_IN_PLACE
          call mpi_allreduce(MPI_IN_PLACE,mpitmp,4, &
               MPI_DOUBLE_PRECISION,mpi_sum,commsander,ierr)
          eke = mpitmp(1)
          ekph = mpitmp(2)
          ekeles = mpitmp(3)
          ekphles = mpitmp(4)
#    else
          call mpi_allreduce(mpitmp,mpitmp(5),4, &
               MPI_DOUBLE_PRECISION,mpi_sum,commsander,ierr)
          eke = mpitmp(5)
          ekph = mpitmp(6)
          ekeles = mpitmp(7)
          ekphles = mpitmp(8)
#    endif
        endif
      end if
#  else

      if ( .not. mpi_orig .and. numtasks > 1 ) then
         call trace_mpi('mpi_allreduce', &
               1,'MPI_DOUBLE_PRECISION',mpi_sum)
         mpitmp(1) = eke
         mpitmp(2) = ekph
         mpitmp(3) = ekpbs
#    ifdef USE_MPI_IN_PLACE
         call mpi_allreduce(MPI_IN_PLACE,mpitmp,3, &
               MPI_DOUBLE_PRECISION,mpi_sum,commsander,ierr)
         eke = mpitmp(1)
         ekph = mpitmp(2)
         ekpbs = mpitmp(3)

#    else

         call mpi_allreduce(mpitmp,mpitmp(4),3, &
               MPI_DOUBLE_PRECISION,mpi_sum,commsander,ierr)
         eke = mpitmp(4)
         ekph = mpitmp(5)
         ekpbs = mpitmp(6)

#    endif
      end if
#  endif
      
      ! Calculate Ekin of the softcore part of the system
      if (ifsc /= 0 ) then
        call calc_softcore_ekin(amass,v,vold,istart,iend)
        sc_ener(13) = sc_ener(6) + sc_ener(12)
      end if      
#endif
      
      !         --- all processors handle the "extra" variables:

      do im=1,iscale
         eke = eke + scalm*0.25d0*(v(nr3+im)+vold(nr3+im))**2
         ekpbs = ekpbs + scalm*v(nr3+im)*vold(nr3+im)
         ekph = ekph + scalm*v(nr3+im)**2
      end do
      
      eke = eke * 0.5d0
      ekph = ekph * 0.5d0
      ekpbs = ekpbs * 0.5d0
#ifdef LES
      ekeles = ekeles * 0.5d0
      ekphles = ekphles * 0.5d0
#endif

      if( ntt == 1 ) then
#ifdef LES
         
         if (temp0les < 0.d0) then
            scaltp = sqrt(1.d0 + 2.d0*dttp*(ekin0-eke)/(ekmh+ekph))
         else
            scaltp = sqrt(1.d0+2.d0*dttp*(ekinp0-eke)/(ekmh+ekph))
            scaltles = sqrt(1.d0+2.d0*dttp*(ekinles0-ekeles)/(ekmhles+ekphles))
         end if
#else
         
         !    --- following is from T.E. Cheatham, III and B.R. Brooks,
         !        Theor. Chem. Acc. 99:279, 1998.
         
         scaltp = sqrt(1.d0 + 2.d0*dttp*(ekin0-eke)/(ekmh+ekph))

         !    --- following is the "old" (amber7 and before) method:

         !  scaltpo = sqrt(1.d0 + dttp*(ekin0/ekph - 1.d0))
         !  write(6,*) 'scaltp: ',2.d0*dttp*(ekin0-eke)/(ekmh+ekph), &
         !            dttp*(ekin0/ekmh - 1.d0)

         !  following line reverts to the "old" behavior:
         !  scaltp = scaltpo

#endif

#ifdef MPI /* SOFT CORE */
         if (icfe /= 0) then
            if (ifsc == 1) then
               if (master) then
                  ! Linearly combine the scaling factors from both processes
                  ! the combined factor is broadcast to all nodes
                  ! the subroutine also correctly scales the softcore atom v's
                  call mix_temp_scaling(scaltp,clambda,v)
               end if
               call mpi_bcast(scaltp,1,MPI_DOUBLE_PRECISION,0,commsander,ierr)
            end if
         end if
#endif

         do j = istart,iend
            i3=(j-1)*3+1
#ifdef LES
            if (temp0les > 0.d0 .and. cnum(j) /= 0 ) then
               v(i3  ) = v(i3  )*scaltles
               v(i3+1) = v(i3+1)*scaltles
               v(i3+2) = v(i3+2)*scaltles
            else
               v(i3  ) = v(i3  ) *scaltp
               v(i3+1) = v(i3+1) *scaltp
               v(i3+2) = v(i3+2) *scaltp
            end if
#else
            v(i3  ) = v(i3  ) *scaltp
            v(i3+1) = v(i3+1) *scaltp
            v(i3+2) = v(i3+2) *scaltp
#endif
         end do
         do im=1,iscale
            v(nr3+im) = v(nr3+im)*scaltp
         end do
      end if  ! (ntt == 1 )
      
   end if  ! ( ntt == 1 .or. onstep; end of step 4c )

   !-----------------------------------------------------------------
   !   Step 5:  several tasks related to dumping of trajectory information
   !-----------------------------------------------------------------
   
   itdump = .false.             ! Write coordinates this step?
   ivdump = .false.             ! Write velocities this step?
   ifdump = .false.             ! Write forces this step?  lam81
   ixdump = .false.             ! Write restart this step?
   ifdump = .false.             ! Write forces this step?
   ivscm  = .false.             ! Do com removal this step?
#ifdef RISMSANDER
   irismdump = .false.          ! Write RISM files this step?
#endif

   !  --- Determine if trajectory, velocity, or restart
   !      writing is imminent, or if the center of mass
   !      motion will be removed.
   !      These require xdist of velocities or dipoles in parallel runs:
   !
   ! Modified so that when running REMD, writing can occur less often
   !  than exchanges (e.g. ntwx > nstlim)
   ! DAN ROE: Added two new variables, total_nstep and total_nstlim.
   !          For non-REMD runs, total_nstep=nstep+1 and total_nstlim=nstlim 
   !           just like before.
   !          For REMD runs, total_nstep=(mdloop-1)*nstlim+nstep+1, where
   !           mdloop is the current exchange - this is the current
   !           replica exchange MD step. total_nstlim=numexchg*nstlim, which is
   !           the maximum number of REMD steps.
   total_nstep=nstep+1
   total_nstlim=nstlim
   if(abfqmmm_param%abfqmmm == 1) total_nstep=abfqmmm_param%qmstep  ! lam81

#ifdef MPI
   if (rem /= 0) then
      total_nstep = (mdloop - 1) * nstlim + nstep + 1
      total_nstlim = nstlim * numexchg
   endif
#endif
   if (ntwx>0) itdump = mod(total_nstep,ntwx) == 0 ! Trajectory coords
   if (ntwv>0) ivdump = mod(total_nstep,ntwv) == 0 ! Velocity
   if (ntwf>0) ifdump = mod(total_nstep,ntwf) == 0 ! Force
   if( ntwr /= 0 ) then
      if ( mod(total_nstep, ntwr ) == 0 ) ixdump = .true. ! Restart
   endif
   if( total_nstep >= total_nstlim ) ixdump = .true. ! Final restart
   if ( nscm > 0 ) then
      if( mod(total_nstep,nscm) == 0 ) ivscm =.true. ! C.o.M. removal
   end if
   if (ntwv == -1 .and. itdump) ivdump = .true. !Combined crdvel file

#ifdef MPI
   ! adaptive QM/MM via multisander
   ! all groups have identical coords and velocities
   ! only master of first group needs to dump results
   ! We have to leave the dump values for all threads in the group, though
   ! since for dumping the coords, these are broadcast within the group
   ! (see call to xdist() below)
   if ( qmmm_nml%vsolv > 1 ) then
      if ( nodeid /= 0 ) then
         ixdump = .false.
         itdump = .false.
         ivdump = .false.
      end if
   end if
#endif

#ifdef RISMSANDER
   if(rismprm%irism ==1)then
      if(rismprm%ntwrism > 0 )then
         irismdump = mod(nstep+1,rismprm%ntwrism) == 0
         if( nstep+1 >= nstlim ) then !! do we want to do this?
            irismdump = .true.
         end if
      end if
   end if
#endif


#ifdef MPI

   !-----------------------------------------------------------------
   !  --- now distribute the coordinates, and if necessary, dipoles and vel:
   !-----------------------------------------------------------------

   call timer_barrier( commsander )
   call timer_stop_start(TIME_VERLET,TIME_DISTCRD)
   if ( .not. mpi_orig .and. numtasks > 1 ) then
      call xdist(x, xx(lfrctmp), natom)
   end if
   ! dac/knut change: force the coordinates to be the same on both masters.
   ! For certain compilers, addition may not be strictly commutative, so
   ! the forces on group 0 may be different by roundoff from the forces on 
   ! group 1.  This can lead to divergent trajectories.  The interval at
   ! which they are resynchronized is hard-wired here to 20, which seems to
   ! work fine in our tests.
   ! jwk change: coordinates are synchronized when shake is enabled above
   if( icfe /= 0 .and. mod(nstep+1,20) == 0 .and. ntc == 1 ) then

      ! In dual-topology this is done within softcore.f
      if (ifsc /= 1) then
         if( master ) call mpi_bcast(x,nr3,MPI_DOUBLE_PRECISION, &
                                     0,commmaster,ierr)
      else
         if( master ) then
            call sc_compare(x,nr3,'CRD') ! first, check if coordinates have desynced
            if (numtasks==1 ) call sc_compare(v,nr3,'VEL') ! do the same for velocities
            call sc_sync_x(x,nr3) ! then resync them
         end if
      end if
      if( numtasks>1 ) call mpi_bcast(x,nr3,MPI_DOUBLE_PRECISION, &
                                      0,commsander,ierr)
   end if
   call timer_stop(TIME_DISTCRD)

#endif  /* MPI */

   !           ----fix lone pair positions:
   if( numextra > 0 )call local_to_global(x,xx,ix)

#ifdef MPI
   if ( .not. mpi_orig .and. numtasks > 1 ) then
      call timer_start(TIME_DISTCRD)
      
      !  ---Here we provide every processor a full copy of the velocities
      !     for removal of center of mass motion, or for archiving.
      !     (Note: this is actually over-kill: for example, only the master
      !     node really needs the velocities for archiving.  But the extra
      !     overhead of doing it this way is probably small in most cases.)
      
      if( ivdump .or. ivscm .or. ixdump ) then
         call xdist(v, xx(lfrctmp), natom)
      endif
  
! M-WJ
!      if( ixdump .and. (induced == 1 .and. indmeth == 3 ) )then
      if( ixdump .and. (induced > 0 .and. indmeth == 3 ) )then
!
         call xdist(xx(ldipvel), xx(lfrctmp), natom)
         call xdist(xx(linddip), xx(lfrctmp), natom)
      end if
      call timer_stop(TIME_DISTCRD)
   end if
   call timer_start(TIME_VERLET)

   !     ========================= END AMBER/MPI =========================
#endif  /* MPI */
   
   !-------------------------------------------------------------------
   !   Step 6: zero COM velocity if requested; used for preventing
   !   ewald "block of ice flying thru space" phenomenon, or accumulation
   !   of rotational momentum in vacuum simulations
   !-------------------------------------------------------------------

   if (ivscm) then
      if (mod(nstep,nsnb) == 0) ntnb = 1
      if( ifbox == 0 ) then
         if (is_langevin) then
            ! Get current center of the system 
            call get_position(nr,x,vcmx,vcmy,vcmz,sysrange,0)
        
#ifdef MMPI /* SOFT CORE */
            if (ifsc == 1) call sc_mix_position(vcmx,vcmy,vcmz,clambda)
#endif  
            ! Center the system to the original center
            call re_position(nr,ntr,x,xc, &
                             vcmx,vcmy,vcmz,sysx,sysy,sysz,sysrange,mv_flag,0)
         else
            !  ---Non-periodic simulation: remove both translation and rotation.
            !     Back the coords up 1/2 step, so that the correspond to the
            !     velocities; temporarily store in the F() array:
            f(1:nr3) = x(1:nr3) - v(1:nr3)*dt5
            !     --- now compute the com motion, remove it, and recompute (just
            !         to check that it is really gone.....)
            call cenmas(nr,f,v,amass,ekcm,xcm,vcm,acm,ekrot,ocm,4)
            call stopcm(nr,f,v,xcm,vcm,ocm, .true.)
            call cenmas(nr,f,v,amass,ekcm,xcm,vcm,acm,ekrot,ocm,4)
         end if
      else
         if (.not. is_langevin) then        
            !    ---Periodic simulation: just remove the translational velocity:
            vcmx = 0.d0
            vcmy = 0.d0
            vcmz = 0.d0
            j = 1
            do i = 1, 3*natom,3
               aamass = amass(j)
               vcmx = vcmx + aamass * v(i)
               vcmy = vcmy + aamass * v(i+1)
               vcmz = vcmz + aamass * v(i+2)
               j = j + 1
            end do
            vcmx = vcmx * tmassinv
            vcmy = vcmy * tmassinv
            vcmz = vcmz * tmassinv
            vel2 = vcmx*vcmx + vcmy*vcmy + vcmz*vcmz
            atempdrop = 0.5d0 * tmass * vel2 * onefac(1) !onefac(1) = 1.0d0/fac(1)
            vel = sqrt(vel2)
            if ( master ) write (6,'(a,f15.6,f9.2,a)') &
               'check COM velocity, temp: ',vel,atempdrop, '(Removed)'
            do i = 1, 3*natom, 3
               v(i)   = v(i)   - vcmx
               v(i+1) = v(i+1) - vcmy
               v(i+2) = v(i+2) - vcmz
            end do
           
#ifdef MPI /* SOFT CORE */
            if (icfe==1) then
               if (ifsc==1) then
                  if (master) then
                     call sc_mix_velocities(v,nr3,clambda)
                  end if
                  call mpi_bcast(v,nr3,MPI_DOUBLE_PRECISION,0,commsander,ierr)
               end if
            end if
#endif
         end if  ! (.not. is_langevin)
      end if  ! ( ifbox == 0 )
   end if  ! (ivscm)
  
   !  Also zero out the non-moving velocities if a belly is active:
   if (belly) call bellyf(nr,ix(ibellygp),v)

   !-----------------------------------------------------------------
   !  --- put current velocities into VOLD
   !-----------------------------------------------------------------
   
   vold(istart3:iend3) = v(istart3:iend3)
   do im=1,iscale
      vold(nr3+im) = v(nr3+im)
   end do

   !-------------------------------------------------------------------
   !  Step 7: scale coordinates if NPT with Berendsen barostat:
   !-------------------------------------------------------------------
   if( ntp > 0 .and. ipimd > 0 .and. barostat == 1 ) then
      x_lnv_old = x_lnv
      x_lnv = x_lnv_old + v_lnv * dtx
      rmu(1:3) = exp( ( x_lnv - x_lnv_old ) )
      box(1:3) = box(1:3) * rmu(1:3)
      volume = box(1) * box(2) * box(3)
      ener%box(1:3) = box(1:3)
      ! only for NMPIMD in sander.LES
      ! (in sander.MPI volume, pressure and density printed in pimdout)
#ifdef LES
      ener%volume = volume
#else
      ener%volume = 0.
      totener%volume = volume
#endif
      call redo_ucell(rmu)
      call fill_tranvec()
      call ew_pscale(natom,x,amass,nspm,nsp,2)
   end if

   if( iamoeba == 0 .and. barostat == 1 ) then

      ! ntp = 1, isotropic pressure coupling

     if (ntp == 1) then
         rmu(1) = (1.d0-dtcp*(pres0-ener%pres(4)))**third
         rmu(2) = rmu(1)
         rmu(3) = rmu(1)


     ! ntp = 2, anisotropic pressure scaling

     else if (ntp == 2) then

        if (csurften > 0) then

          ! Constant surface tension adjusts the tangential pressures
          ! See Zhang, Feller, Brooks, Pastor. J. Chem. Phys. 1995

          if (csurften == 1) then        ! For surface tension in the x direction
            pres0y = pres0x - gamma_ten_int * ten_conv / box(1)
            pres0z = pres0y

          else if (csurften == 2) then   ! For surface tension in the y direction
            pres0x = pres0y - gamma_ten_int * ten_conv / box(2)
            pres0z = pres0x

          !else if (csurften == 3) then   ! For surface tension in the z !direction
          else
            pres0x = pres0z - gamma_ten_int * ten_conv / box(3)
            pres0y = pres0x

          end if

          rmu(1) = (1.d0 - dtcp * (pres0x - ener%pres(1)))**third
          rmu(2) = (1.d0 - dtcp * (pres0y - ener%pres(2)))**third
          rmu(3) = (1.d0 - dtcp * (pres0z - ener%pres(3)))**third

        else

          rmu(1) = (1.d0-dtcp*(pres0-ener%pres(1)))**third
          rmu(2) = (1.d0-dtcp*(pres0-ener%pres(2)))**third
          rmu(3) = (1.d0-dtcp*(pres0-ener%pres(3)))**third

        end if

     ! ntp = 3, semiisotropic pressure coupling
     ! (currently only for csurften>0, constant surface tension)

     !else if (ntp > 2) then
     else

        if (csurften > 0) then

          if (csurften == 1) then        ! For surface tension in the x direction
            pres0y = pres0x - gamma_ten_int * ten_conv / box(1)
            pres0z = pres0y
            press_tan_ave = (ener%pres(2) + ener%pres(3))/2
            rmu(1) = (1.d0 - dtcp * (pres0x - ener%pres(1)))**third
            rmu(2) = (1.d0 - dtcp * (pres0y - press_tan_ave))**third
            rmu(3) = (1.d0 - dtcp * (pres0z - press_tan_ave))**third

          else if (csurften == 2) then   ! For surface tension in the y direction
            pres0x = pres0y - gamma_ten_int * ten_conv / box(2)
            pres0z = pres0x
            press_tan_ave = (ener%pres(1) + ener%pres(3))/2
            rmu(1) = (1.d0 - dtcp * (pres0x - press_tan_ave))**third
            rmu(2) = (1.d0 - dtcp * (pres0y - ener%pres(2)))**third
            rmu(3) = (1.d0 - dtcp * (pres0z - press_tan_ave))**third

          !else if (csurften == 3) then   ! For surface tension in the z !direction
          else
            pres0x = pres0z - gamma_ten_int * ten_conv / box(3)
            pres0y = pres0x
            press_tan_ave = (ener%pres(1) + ener%pres(2))/2
            rmu(1) = (1.d0 - dtcp * (pres0x - press_tan_ave))**third
            rmu(2) = (1.d0 - dtcp * (pres0y - press_tan_ave))**third
            rmu(3) = (1.d0 - dtcp * (pres0z - ener%pres(3)))**third

          end if
        end if
        ! Add semiisotropic pressure scaling in any direction with no constant 
        ! surface tension here
      end if

      if (ntp > 0) then
         box(1:3) = box(1:3)*rmu(1:3)
         ener%box(1:3) = box(1:3)
         
         !    WARNING!!   This is not correct for non-orthogonal boxes if
         !    NTP > 1 (i.e. non-isotropic scaling).  Currently general cell
         !    updates which allow cell angles to change are not implemented.
         !    The viral tensor computed for ewald is the general Nose Klein,
         !    however the cell response needs a more general treatment.
      
         call redo_ucell(rmu)
         ! keep tranvec up to date, rather than recomputing each MD step.
         call fill_tranvec()  ! tranvec is dependent on only ucell

#ifdef MPI /* SOFT CORE */
         ! if softcore potentials and the dual topology approach are used
         ! C.O.M. scaling has to be changed to account for different masses 
         ! of the same molecule in V0 and V1. This is quite inefficient and is
         ! therefore done in a separate routine in softcore.f
         ! only both masters actually do the computation for ifsc==1
         ! the scaled coordinates are then broadcast to the nodes
         if (icfe /= 0 .and. ifsc == 1) then
            if (master) then
               call sc_pscale(natom,x,amass,nspm,nsp,oldrecip,ucell)
            end if
            call mpi_bcast(x,nr3,MPI_DOUBLE_PRECISION,0,commsander,ierr)
         else
#endif
            call ew_pscale(natom,x,amass,nspm,nsp,npscal)
#ifdef MPI /* SOFT CORE */
         end if
#endif
         if (ntr > 0 .and. nrc > 0) &
            call ew_pscale(natom,xc,amass,nspm,nsp,npscal)
      endif
      if (ipimd==NMPIMD.and.ntp>0) then
         ener%cmt(4) = 0.d0
         ener%vir(4) = 0.d0
         ener%pres(4) = pressure*pconv
       endif
   else if (barostat == 1) then
      if (ntp>0) then
         if (ipimd==0) then     ! for classical AMOEBA
            ener%cmt(4)  = eke      !  for printing in prntmd()
            ener%vir(4)  = ener%vir(1) + ener%vir(2) + ener%vir(3)
            ener%pres(4) = (pressure_constant/volume)*(2.d0*eke - ener%vir(4)) / 3.d0
         elseif (ipimd==NMPIMD) then     ! for NMPIMD AMOEBA
            ener%cmt(4)  = 0.d0
            ener%vir(4)  = 0.d0
            ener%pres(4) = pressure*pconv
         endif
         call AM_RUNMD_scale_cell(natom,ener%pres(4),dt,pres0,taup,x)
         call fill_tranvec()
      end if
   end if
   
#ifdef LES
   ener%kin%solt = eke
   ener%kin%solv = ekeles
   ener%kin%tot  = ener%kin%solt + ener%kin%solv
   if (ntt == 1 .and. onstep) then
      if ( temp0les < 0 ) then
        ekmh = max(ekph,fac(1)*10.d0)
      else
        ekmh = max(ekph,fac(2)*10.d0)
        ekmhles = max(ekphles,fac(3)*10.d0)
      endif
   end if

   if( ipimd > 0 ) then
      ener%kin%solv = equal_part + Epot_deriv  ! "virial" estimate of KE
      ener%tot = ener%kin%solv + ener%pot%tot
   endif
#else
   if( ipimd > 0 ) then
      ! use a "virial" estimator for the KE, rather than one derived from the 
      ! bead velocities:
      totener%kin%solv = equal_part + Epot_deriv
   else
      ener%kin%solv = ekpbs + ener%pot%tot  
                                  ! Pastor, Brooks, Szabo conserved quantity
                                  ! for harmonic oscillator: Eq. 4.7b of Mol.
                                  ! Phys. 65:1409-1419, 1988
   endif
   ener%kin%solt = eke
   ener%kin%tot  = ener%kin%solt
   if (ntt == 1 .and. onstep) then
      ekmh = max(ekph,fac(1)*10.d0)
   end if
#endif
   
   !     ---if velocities were reset, the KE is not accurate; fudge it
   !        here to keep the same total energy as on the previous step.
   !        Note that this only affects printout and averages for Etot
   !        and KE -- it has no effect on the trajectory, or on any averages
   !        of potential energy terms.
   
   if( resetvelo ) ener%kin%tot = etot_save - ener%pot%tot
   
   !     --- total energy is sum of KE + PE:
   
   if( ipimd >  0 ) then
      totener%tot = totener%kin%solv + totener%pot%tot
      etot_save   = totener%kin%tot  + totener%pot%tot
      if (ipimd==CMD) then
         etot_cmd = eke_cmd*0.5 + ener%pot%tot

         totener%tot= etot_cmd

         ener%tot   = etot_cmd
         ener%kin%tot  = eke_cmd*0.5
         ener%kin%solv = ener%kin%tot
      endif
   else
      ener%tot  = ener%kin%tot + ener%pot%tot
      etot_save = ener%tot
   end if

   !-------------------------------------------------------------------
   !  Step 8:  update the step counter and the integration time:
   !-------------------------------------------------------------------

   if(abfqmmm_param%abfqmmm /= 1) then  ! lam81   
    nstep = nstep+1
    t = t+dt
   end if                               ! lam81

   !For CMD
   if ( ipimd==CMD ) then
      nstep_cmd = nstep_cmd + 1
      t_cmd = t_cmd + dt
   end if
   
   !     ---full energies are only calculated every nrespa steps
   !     nvalid is the number of steps where all energies are calculated
   
   if (onstep .or. aqmmm_flag > 0) then
      nvalid = nvalid + 1
      ! Update all elements of these sequence types
      enert  = enert + ener
      enert2 = enert2 + (ener*ener)
#ifdef MPI
      if( ievb /= 0 ) then
         evb_nrg_ave(:) = evb_nrg_ave(:) + evb_nrg(:)
         evb_nrg_rms(:) = evb_nrg_rms(:) + evb_nrg(:)**2
      endif
      if ( ifsc /= 0 ) then
         sc_ener_ave(1:ti_ene_cnt) = sc_ener_ave(1:ti_ene_cnt) + sc_ener(1:ti_ene_cnt)
         sc_ener_rms(1:ti_ene_cnt) = sc_ener_rms(1:ti_ene_cnt) + sc_ener(1:ti_ene_cnt)**2
      end if
#endif
      if( nvalid == 1 ) etot_start = ener%tot

#ifndef LES
      if ( ipimd>0 .or. ineb>0 ) then
#  ifdef MPI
         if (master) call mpi_reduce(ener%kin%tot,totener%kin%tot,1,MPI_DOUBLE_PRECISION, &
                         mpi_sum,0,commmaster,ierr)
        
#  endif
      endif

      ! Passing of dvdl=dV/dl for TI w.r.t. mass
      ! Note that ener(39) (in runmd and mix_frcti) = 
      !       = ener(17) = ene(21) (in force). All denote dvdl.
      ! Note, ener() is now historical, MJW Feb 2010
      if (ipimd>0 .and. itimass>0) totener%pot%dvdl = ener%pot%dvdl

      if(ipimd.eq.NMPIMD.and.ntp>0) then
         totener%pres(4) = pressure * pconv
         totener%density = tmass / (0.602204d0*volume)
      endif
      if(ipimd.eq.CMD) then
         totener%kin%tot  = eke_cmd*0.5d0
         totener%kin%solv = totener%kin%tot
         totener%tot      = totener%kin%tot + totener%pot%tot
      endif
      totenert  = totenert + totener
      totenert2 = totenert2 + (totener*totener)

#endif /* LES */

      kinetic_E_save(2) = kinetic_E_save(1)
      kinetic_E_save(1) = ener%kin%tot

   end if

   ! added for rbornstat
!!FIX: TL - do we need to put in rismnrespa here?
   if (mod(irespa,nrespai) == 0 .or. irespa < 2) nvalidi = nvalidi + 1

   ntnb = 0
   if (mod(nstep,nsnb) == 0) ntnb = 1

   ! Since nstep has been incremented, total_nstep is now equal to
   !    (mdloop-1)*nstlim+nstep for REMD and nstep for MD.
   lout = mod(total_nstep,ntpr) == 0 .and. onstep

   irespa = irespa + 1
     
   ! reset pb-related flags
#ifdef MPI
   if(mytaskid == 0)then
#endif   
      if ( igb == 10 .or. ipb /= 0 ) then
         if ( mod(nstep,npbgrid) == 0 .and. nstep /= nstlim ) pbgrid = .true.
         if ( mod(nstep,ntpr) == 0 .or. nstep == nstlim ) pbprint = .true.
         if ( mod(nstep,nsnbr) == 0 .and. nstep /= nstlim ) ntnbr = 1
         if ( mod(nstep,nsnba) == 0 .and. nstep /= nstlim ) ntnba = 1
      end if
#ifdef MPI
   endif
#endif   

   !-------------------------------------------------------------------
   !  Step 9:  output from this step if required:
   !-------------------------------------------------------------------
   
#ifdef RISMSANDER
   !some 3D-RISM files require all processes to participate in output
   !due to the distributed memory
   !     RISM archive:
   if(rismprm%irism==1)then
!!$      if(irismdump)&
!!$           call rism_writeSolvDistF(rism_3d,nstep)
     ! combined thermodynamics and distribution output 
     ! Execute if we need to do either 
      if(irismdump .or. (rism_calc_type(nstep) == RISM_FULL &
           .and. rismprm%write_thermo==1 .and. lout))&
           call rism_solvdist_thermo_calc(irismdump,nstep)
   endif
#endif

   !     ...only the master needs to do the output
   if (ixdump) then
      if(ipimd.eq.NMPIMD.or.ipimd.eq.CMD) then
         call trans_pos_nmode_to_cart(x,cartpos)
         call trans_vel_nmode_to_cart(v,cartvel)
      endif
   endif

   if (itdump) then
      if(ipimd.eq.NMPIMD.or.ipimd.eq.CMD) then
         call trans_pos_nmode_to_cart(x,cartpos)
      endif
!AMD Flush amdlog file
      if(iamd.gt.0)then
#  ifdef MPI
        if (worldrank.eq.0) &
#  endif
          call write_amd_weights(ntwx,total_nstep)
      end if 
!scaledMD Flush scaledMDlog file
      if(scaledMD.gt.0)then
#  ifdef MPI
        if (worldrank.eq.0) &
#  endif
          call write_scaledMD_log(ntwx,total_nstep)
      end if 

   endif

   if (ivdump) then
      if(ipimd.eq.NMPIMD.or.ipimd.eq.CMD) then
         call trans_vel_nmode_to_cart(v,cartvel)
      endif
   endif


   if (master) then
      
      !        -- restrt:
      
      if (ixdump) then

         ! NOTE - This assumes that if numextra > 0, then velocities are
         !        found in the array v...
         if (numextra > 0) call zero_extra_pnts_vec(v,ix)

         if( iwrap == 0 ) then
            nr = nrp
#ifdef LES
            if(ipimd.eq.NMPIMD.or.ipimd.eq.CMD) then
               call mdwrit(nstep,nrp,nr,nres,ntxo,ntr,ntb, &
                    cartpos,cartvel,xx(lcrdr),box,t,temp0)
            else
               call mdwrit(nstep,nrp,nr,nres,ntxo,ntr,ntb, &
                    x,v,xx(lcrdr),box,t,temp0les)
            endif
#else
            if(ipimd.eq.NMPIMD.or.ipimd.eq.CMD) then
               call mdwrit(nstep,nrp,nr,nres,ntxo,ntr,ntb, &
                    cartpos,cartvel,xx(lcrdr),box,t,rem_val)
            else
               call mdwrit(nstep,nrp,nr,nres,ntxo,ntr,ntb, &
                    x,v,xx(lcrdr),box,t,rem_val)
            endif
#endif
         else if (iwrap == 1) then
            
            ! --- use temp. array to hold coords. so that the master's values
            !     are always identical to those on all other nodes:
            
            call get_stack(l_temp,nr3,routine)
            if(.not. rstack_ok)then
               deallocate(r_stack)
               allocate(r_stack(1:lastrst),stat=alloc_ier)
               call reassign_rstack(routine)
            endif
            REQUIRE(rstack_ok)
            if(ipimd.eq.NMPIMD.or.ipimd.eq.CMD) then
               do iatom=1,natom
               do m=1,3
                  r_stack(l_temp+3*(iatom-1)+m-1)=cartpos(m,iatom)
               end do
               end do
            else
               do m=1,nr3
                  r_stack(l_temp+m-1) = x(m)
               end do
            end if

            call wrap_molecules(nspm,nsp,r_stack(l_temp))
            if(ifbox == 2) call wrap_to(nspm,nsp,r_stack(l_temp),box)
            nr = nrp
#ifdef LES
            call mdwrit(nstep,nrp,nr,nres,ntxo,ntr,ntb, &
                  r_stack(l_temp),v,xx(lcrdr),box,t,temp0les)
#else
            call mdwrit(nstep,nrp,nr,nres,ntxo,ntr,ntb, &
                  r_stack(l_temp),v,xx(lcrdr),box,t,rem_val)
#endif
            call free_stack(l_temp,routine)
         else if (iwrap == 2) then
           ! GMS ------------------------------------------
           ! We are wrapping around a pre-determined mask
           ! Need to center it on the mask COM first, then
           ! wrap it normally as it happens on the iwrap=1 
           ! case.
           ! GMS ------------------------------------------
            call get_stack(l_temp,nr3,routine)
            if(.not. rstack_ok)then
               deallocate(r_stack)
               allocate(r_stack(1:lastrst),stat=alloc_ier)
               call reassign_rstack(routine)
            endif
            REQUIRE(rstack_ok)
            if(ipimd.eq.NMPIMD.or.ipimd.eq.CMD) then
               do iatom=1,natom
               do m=1,3
                  r_stack(l_temp+3*(iatom-1)+m-1)=cartpos(m,iatom)
               end do
               end do
            else
               do m=1,nr3
                  r_stack(l_temp+m-1) = x(m)
               end do
            end if
            nr = nrp

            ! Now, wrap the coordinates around the iwrap_mask:
            call iwrap2(n_iwrap_mask_atoms,iwrap_mask_atoms,r_stack(l_temp), &
                        box_center)
#ifdef LES
            call mdwrit(nstep,nrp,nr,nres,ntxo,ntr,ntb, &
                  r_stack(l_temp),v,xx(lcrdr),box,t,temp0les)
#else
            call mdwrit(nstep,nrp,nr,nres,ntxo,ntr,ntb, &
                  r_stack(l_temp),v,xx(lcrdr),box,t,rem_val)
#endif
            call free_stack(l_temp,routine)
         end if  ! ( iwrap == 0 )

! M-WJ
!        if( igb == 0 .and. induced == 1 .and. indmeth == 3) &
         if( igb == 0 .and. ipb == 0 .and. induced > 0 .and. indmeth == 3) &
!
            call wrt_dips(xx(linddip),xx(ldipvel),nr,t,title)

         if (icnstph /= 0 .and. ((rem /= 0 .and. mdloop > 0) .or. rem == 0)) then
            call cnstphwriterestart(chrgdat)
         end if

      end if  ! (ixdump)
      
      !     -- Coordinate archive:
      ! For formatted writes and replica exchange, write out a header line.

      if (itdump) then
#ifdef MPI
         ! Write out current replica#, exchange#, step#, and mytargettemp
         ! If mdloop==0 this is a normal md run (since REMD never calls corpac
         !  when mdloop==0) and we don't want the REMD header.
         ! total_nstep is set in step 5. 
         if (mdloop > 0 .and. loutfm) then
            if (trxsgld) then
               write (MDCRD_UNIT,'(a,4(1x,i8))') "RXSGLD ", repnum, mdloop, &
                              total_nstep, stagid
            else
               write (MDCRD_UNIT,'(a,3(1x,i8),1x,f8.3)') "REMD ", repnum, mdloop, &
                              total_nstep, my_remd_data%mytargettemp
            end if
         end if
#endif

         if( iwrap == 0 ) then
            if(ipimd.eq.NMPIMD.or.ipimd.eq.CMD) then
               call corpac(cartpos,1,nrx,MDCRD_UNIT,loutfm)
            else
               call corpac(x,1,nrx,MDCRD_UNIT,loutfm)
            endif
            if(ntb > 0)  call corpac(box,1,3,MDCRD_UNIT,loutfm)
         else if (iwrap == 1) then
            call get_stack(l_temp,nr3,routine)
            if(.not. rstack_ok)then
               deallocate(r_stack)
               allocate(r_stack(1:lastrst),stat=alloc_ier)
               call reassign_rstack(routine)
            endif
            REQUIRE(rstack_ok)
            if(ipimd.eq.NMPIMD.or.ipimd.eq.CMD) then
               do iatom=1,natom
               do m=1,3
                  r_stack(l_temp+3*(iatom-1)+m-1) = cartpos(m,iatom)
               end do
               end do
            else
               do m=1,nr3
                  r_stack(l_temp+m-1) = x(m)
               end do
            endif
            
            call wrap_molecules(nspm,nsp,r_stack(l_temp))
            if (ifbox == 2) call wrap_to(nspm,nsp,r_stack(l_temp),box)
            
            call corpac(r_stack(l_temp),1,nrx,MDCRD_UNIT,loutfm)
            call corpac(box,1,3,MDCRD_UNIT,loutfm)
            call free_stack(l_temp,routine)
         else if (iwrap == 2) then
           ! GMS ------------------------------------------
           ! We are wrapping around a pre-determined mask
           ! Need to center it on the mask COM first, then
           ! wrap it normally as it happens on the iwrap=1 
           ! case.
           ! GMS ------------------------------------------
            call get_stack(l_temp,nr3,routine)
            if(.not. rstack_ok)then
               deallocate(r_stack)
               allocate(r_stack(1:lastrst),stat=alloc_ier)
               call reassign_rstack(routine)
            endif
            REQUIRE(rstack_ok)
            if(ipimd.eq.NMPIMD.or.ipimd.eq.CMD) then
               do iatom=1,natom
               do m=1,3
                  r_stack(l_temp+3*(iatom-1)+m-1) = cartpos(m,iatom)
               end do
               end do
            else
               do m=1,nr3
                  r_stack(l_temp+m-1) = x(m)
               end do
            endif
            
            call iwrap2(n_iwrap_mask_atoms,iwrap_mask_atoms, r_stack(l_temp), &
                        box_center)

            call corpac(r_stack(l_temp),1,nrx,MDCRD_UNIT,loutfm)
            call corpac(box,1,3,MDCRD_UNIT,loutfm)
            call free_stack(l_temp,routine)
           
           
          end if ! if (iwrap == 0) ...


         !GMS: If using variable QM solvent, try to write a new pdb file
         !     with the QM coordinates for this step. This is done here
         !     to keep the PDB file in sync with the mdcrd file, which
         !     makes it easier to check later.
         if (qmmm_nml%vsolv > 0 .and. qmmm_nml%verbosity == 0) &
            call qm_print_coords(nstep,.false.)
      end if  ! (itdump)

      !     Velocity archive:
      
      if (ivdump) then

         ! NOTE - This assumes that if numextra > 0, then velocities are
         !        found in the array v...
         if (numextra > 0) call zero_extra_pnts_vec(v,ix)

#ifdef MPI
         ! Write out current replica#, exchange#, step#, and mytargettemp
         ! If mdloop==0 this is a normal md run (since REMD never calls corpac
         !  when mdloop==0) and we don't want the REMD header.
         if (mdloop>0.and.loutfm) then
            if (trxsgld) then
               write (MDVEL_UNIT,'(a,4(1x,i8))') "RXSGLD ", repnum, mdloop, &
                              total_nstep, stagid
            else
               write (MDVEL_UNIT,'(a,3(1x,i8),1x,f8.3)') "REMD ", repnum, mdloop, &
                              total_nstep, my_remd_data%mytargettemp
            end if
         end if
#endif

          if(ipimd.eq.NMPIMD.or.ipimd.eq.CMD) then
             call corpac(cartvel,1,nrx,MDVEL_UNIT,loutfm)
          else
             call corpac(v,1,nrx,MDVEL_UNIT,loutfm)
          endif
      end if
      !     Force archive lam81
      if (ifdump .and. (abfqmmm_param%abfqmmm == 1)) call corpac(for,1,nrx,MDFRC_UNIT,loutfm)        ! lam81

      !     Energy archive:
      !      (total_nstep set in Step 5.)
      if (ntwe > 0) then
         if (mod(total_nstep,ntwe) == 0.and.onstep) &
               call mdeng(15,nstep,t,ener,onefac,ntp,csurften)
      end if

      if (ioutfm > 0) then
         if (itdump) call end_binary_frame(MDCRD_UNIT)
         if (ivdump .and. ntwv>0 ) call end_binary_frame(MDVEL_UNIT)
         if (ifdump .and. ntwf>0 ) call end_binary_frame(MDFRC_UNIT)
      end if

#ifdef MPI
      if( ievb /= 0 ) call out_evb ( nstep )
#endif
      
      !     General printed output:
      
      if (lout) then
         if (facc /= 'A') rewind(7)

         ! Conserved quantity for Nose'-Hoover based thermostats.          ! APJ
         if (ipimd.eq.0 .and. ntt > 4 .and. ntt <= 8 ) then                ! APJ
            Econserved = ener%kin%tot + ener%pot%tot + E_nhc               ! APJ
            if( ntp>0 ) Econserved = Econserved + pres0 / pconv * volume   ! APJ
#  ifdef MPI
            if ( worldrank.eq.0 ) &                                        ! APJ
#  endif
            write(file_nhc,'(I10,F14.4)') nstep, Econserved                ! APJ
         endif                                                             ! APJ
#ifdef LES
         if (ipimd>0.and.ntt==4) then
            Econserved = ener%kin%tot + ener%pot%tot + E_nhc
            Econserved = Econserved   + Epot_spring
            if( ntp>0 ) Econserved = Econserved + pres0 / pconv * volume
            write(file_nhc,'(I10,F14.4)') nstep, Econserved
         endif
         if ( ipimd.eq.CMD ) then
            ener%kin%tot  = eke_cmd*0.5d0
            ener%kin%solv = ener%kin%tot
            ener%tot   = ener%kin%tot + ener%pot%tot

         end if
#else
         if ( ipimd>0 ) then
               ener%tot     = 0.d0
               ener%kin%tot = 0.d0
            ! Conserved quantity for Nose'-Hoover thermostat.
            if ( ntt==4 ) then
               Econserved = totener%kin%tot + totener%pot%tot + E_nhc
               Econserved = Econserved + Epot_spring
               if ( ntp>0 ) Econserved=Econserved+pres0/pconv*volume
#  ifdef MPI
               if ( worldrank.eq.0 ) &
#  endif
                  write(file_nhc,'(I10,F14.4)') nstep, Econserved
            endif
#  ifdef MPI
            if(worldrank.eq.0) &
#  endif
               call pimd_report(nstep,t,pimd_unit,totener,onefac)
         end if
#endif /* LES */
         call prntmd(total_nstep,nitp,nits,t,ener,onefac,7,.false.)

#  ifdef MPI
         ! print corrected energy for adaptive qm/mm runs
         ! note: nstep has already been increased here
         !       (it was not increased when adaptive_qmmm() was called above)
         if ( qmmm_nml%vsolv > 1 ) then

            if ( masterrank == 0 ) then

               if (aqmmm_flag > 0 .and. nstep > aqmmm_flag) then

                  etotcorr = corrected_energy + kinetic_E_save(aqmmm_flag)
                  nstepadc = nstep - aqmmm_flag + 1
                  tadc = t - dt * (dble( aqmmm_flag - 1) )

                  write(6,'(a)')' Adaptive QM/MM energies:'
                  write(6,'(x,a,i5,x,a,f11.4,x,2(a,f15.4,x))') &
                       'adQMMM STEP=', nstepadc, &
                       'TIME(PS)=', tadc, &
                       'ETC=', etotcorr, &
                       'EPC=', corrected_energy

                  ! print total energy for adaptive qm/mm into a separate file
                  ! when qmmm_vsolv%verbosity > 0
                  ! set reference energy to zero only for energy dumping purposes
                  if (flag_first_energy) then
                     flag_first_energy = .false.
                     adqmmm_first_energy = etotcorr
                     etotcorr = 0.0d0
                  else
                     etotcorr = etotcorr - adqmmm_first_energy
                  end if

                  if (qmmm_vsolv%verbosity > 0) then
                     open(80,file='adqmmm_tot_energy.dat',position='append')
                     write(80,'(i9,5x,f11.4,5x,f15.4)') nstepadc, tadc, etotcorr
                     close(80)
                  end if

               end if
            end if
         end if
#  endif

#ifdef MPI /* SOFT CORE */
         if (ifsc /= 0) call sc_print_energies(6, sc_ener)
         if (ifsc /= 0) call sc_print_energies(7, sc_ener)
#endif
         if ( ifcr > 0 .and. crprintcharges > 0 ) then
            call cr_print_charge( xx(l15), total_nstep )
         end if

         ! Output for CMD.
#ifdef LES
         if (ipimd.eq.CMD) then

            ncmd = 0
            do iatom = 1, natom
               if ( cnum(iatom)==0 .or. cnum(iatom)==1 )  then
                  xcmd(ncmd+1) = x(3*iatom-2)
                  xcmd(ncmd+2) = x(3*iatom-1)
                  xcmd(ncmd+3) = x(3*iatom)
                  vcmd(ncmd+1) = v(3*iatom-2)
                  vcmd(ncmd+2) = v(3*iatom-1)
                  vcmd(ncmd+3) = v(3*iatom)
                  ncmd = ncmd+3
               endif
            enddo
            write(file_pos_cmd,'(10f8.3)') xcmd(1:ncmd)
            write(file_vel_cmd,'(10f8.3)') vcmd(1:ncmd)
            write(file_pos_cmd,'(10f8.3)') box(1:3)

            eke_cmd = eke_cmd * 0.5d0
            etot_cmd = eke_cmd + ener%pot%tot

            if (eq_cmd) then
               temp_cmd = eke_cmd/boltz2/dble(3*natomCL)
            else
               temp_cmd = eke_cmd/boltz2/dble(3*(natomCL-1))
            endif

         endif
#else
         if (ipimd.eq.CMD.and.mybeadid.eq.1) then
            write(file_pos_cmd,'(10f8.3)') x(1:3*natom)
            write(file_vel_cmd,'(10f8.3)') v(1:3*natom)
            write(file_pos_cmd,'(10f8.3)') box(1:3)

            eke_cmd = eke_cmd * 0.5d0
            etot_cmd = eke_cmd + totener%pot%tot

            if (eq_cmd) then
               temp_cmd = eke_cmd/boltz2/dble(3*natom)
            else
               temp_cmd = eke_cmd/boltz2/dble(3*(natom-1))
            endif
         end if
#endif /* LES */

         !--- Print QMMM Muliken Charges if needed ---
         if (qmmm_nml%ifqnt) then
           if (qmmm_nml%printcharges .and. qmmm_mpi%commqmmm_master) then
             call qm2_print_charges(nstep,qmmm_nml%dftb_chg,qmmm_struct%nquant_nlink, &
                                    qm2_struct%scf_mchg,qmmm_struct%iqm_atomic_numbers)
           end if
         end if
         if (qmmm_nml%printdipole /= 0) then
             call qmmm_dipole(x,xx(Lmass),ix(i02),ih(m02),nres)
         end if

         !--- BEGIN DIPOLE PRINTING CODE ---

         ! RCW 2nd Dec 2003 - also output dipole information if
         ! the dipoles namelist has been specified and corresponding
         ! groups defined.

         ! Check input unit 5 for namelist dipoles
         ! We expect to find &dipoles followed by a group
         ! specification of the dipoles to output.
         call nmlsrc('dipoles',5,prndipfind)

         if(prndipfind /= 0 ) then
           !We calculate the dipoles
           write(6,*) '------------------------------- DIPOLE INFO ----------------------------------'
           write(6,9018) nstep,t
           9018 format(/1x, 'NSTEP =',i7,1x,'TIME(PS) =',f10.3)

           !Get the groups for the dipoles - Ideally we only really want
           !to call this the once but for the time being I will call it
           !every time
         
           read (5,'(a)') prndiptest

           call rgroup(natom,natc,nres,prndipngrp,ix(i02),ih(m02), &
                ih(m04),ih(m06),ih(m08),ix(icnstrgp), &
                jgroup,indx,irespw,npdec, &
                xx(l60),xx(lcrdr),0,0,0,idecomp,5,.false.)

           ! Need to rewind input file after rgroup so it is available
           ! when we next loop through

           rewind(5)

           if(prndipngrp > 0) then
             !prndipngrp - holds number of groups specified + 1
             !ix(icnstrgp) - holds map of group membership for each atom
             !x(lcrd) - X,Y,Z coords of atoms - (3,*)
             !x(l15) - Partial Charges
             !x(linddip) - induced dipoles X,Y,Z for each atom (3,*)
             !x(Lmass) - Mass of each atom
             call printdip(prndipngrp,ix(icnstrgp),xx(lcrd), &
                  xx(l15),xx(linddip),xx(Lmass), natom)
           end if
           write(6,*) '----------------------------- END DIPOLE INFO --------------------------------'
         end if
         !--- END DIPOLE PRINTING CODE ---

         if (nmropt > 0) then
            call nmrptx(6)
         end if
         if (itgtmd == 2)  then
           emtmd = 0.0d0
           call mtmdcall(emtmd,xx(lmtmd01),ix(imtmd02),x,f,ih(m04),ih(m02),ix(i02),&
                    ih(m06),xx(lmass),natom,nres,'PRNT')
         end if
         call amflsh(7)
      end if
      
      !  Output running averages:
      ! DAN ROE: total_nstep==Total nstep REMD/MD, set in step 5
      if ( ntave > 0 )then
         if ( mod(total_nstep,ntave) == 0 .and. onstep )then
            write(6,542)
#ifdef RISMSANDER
            if(rismprm%irism==1)then
               tspan = ntave/mylcm(nrespa,rismprm%rismnrespa)
            else
               tspan = ntave/nrespa
            end if
#else
            tspan = ntave/nrespa
#endif

            ! Update all elements of these sequence types
            enert_tmp  = enert - enert_old
            enert2_tmp = enert2 - enert2_old
            enert_old  = enert
            enert2_old = enert2
            enert_tmp  = enert_tmp/tspan
            enert2_tmp = enert2_tmp/tspan - &
                     enert_tmp*enert_tmp
            call zero_neg_values_state(enert2_tmp)
            enert2_tmp = sqrt(enert2_tmp)

#ifdef MPI
            if( ievb /= 0 ) then
               evb_nrg_tmp (:) = evb_nrg_ave(:) - evb_nrg_old (:)
               evb_nrg_tmp2(:) = evb_nrg_rms(:) - evb_nrg_old2(:)
               evb_nrg_old (:) = evb_nrg_ave(:)
               evb_nrg_old2(:) = evb_nrg_rms(:)
               evb_nrg_tmp (:) = evb_nrg_tmp (:) / tspan
               evb_nrg_tmp2(:) = evb_nrg_tmp2(:) / tspan - evb_nrg_tmp(:)**2
               evb_nrg_tmp2(:) = max( evb_nrg_tmp2(:), 0.0d0 )
               evb_nrg_tmp2(:) = sqrt( evb_nrg_tmp2(:) )
            endif
            if ( ifsc /= 0 ) then
               do m = 1,ti_ene_cnt
                  sc_ener_tmp(m) = sc_ener_ave(m)-sc_ener_old(m)
                  sc_ener_tmp2(m) = sc_ener_rms(m)-sc_ener_old2(m)
                  sc_ener_old(m) = sc_ener_ave(m)
                  sc_ener_old2(m) = sc_ener_rms(m)
                  sc_ener_tmp(m) = sc_ener_tmp(m)/tspan
                  sc_ener_tmp2(m) = sc_ener_tmp2(m)/tspan - sc_ener_tmp(m)**2
                  if (sc_ener_tmp2(m) < 0.0d0) sc_ener_tmp2(m) = 0.0d0
                  sc_ener_tmp2(m) = sqrt(sc_ener_tmp2(m))
               end do
            end if         
            if( ievb /= 0 ) evb_frc%evb_ave = .true.
#endif
#ifdef RISMSANDER
            if(rismprm%irism==1)then
               write(6,540) ntave/mylcm(nrespa,rismprm%rismnrespa)!nrespa
            else
               write(6,540) ntave/nrespa
            end if
#else
            write(6,540) ntave/nrespa
#endif
            call prntmd(total_nstep,izero,izero,t,enert_tmp,onefac,0,.false.)
#ifdef MPI
            if (ifsc /= 0) call sc_print_energies(6, sc_ener_tmp)
            if( ievb /= 0 ) evb_frc%evb_rms = .true.
#endif
            write(6,550)
            call prntmd(total_nstep,izero,izero,t,enert2_tmp,onefac,0,.true.)
#ifdef MPI /* SOFT CORE */
            if (ifsc /= 0) call sc_print_energies(6, sc_ener_tmp2)
#endif
            if( icfe > 0 ) then
#ifdef RISMSANDER
               if(rismprm%irism==1)then
                  write(6,541) ntave/mylcm(nrespa,rismprm%rismnrespa)!nrespa
               else
                  write(6,541) ntave/nrespa
               end if
#else
               write(6,541) ntave/nrespa
#endif
               edvdl_r = edvdl_r/tspan
               edvdl_r%pot%dvdl = enert_tmp%pot%dvdl  ! fix for DV/DL output
               edvdl_r%virvsene = 0.d0 ! virvsene should not but included here
               call prntmd(total_nstep,izero,izero,t,edvdl_r,onefac,0,.false.)
               edvdl_r = null_state_rec

            end if
            write(6,542)
         end if
      end if  ! ( ntave > 0 )
      
      !     --- end masters output ---
      
   end if  ! (master)

#ifdef MPI /* SOFT CORE */
   if (ntave > 0 .and. icfe > 0 .and. dynlmb > 0) then
      if ( mod(nstep,ntave) == 0 .and. onstep ) then
         ! For runs with dynamically changing lambda, raise lambda here
         ! and flush all buffers for the next averages
         clambda = clambda + dynlmb
         call sc_change_clambda(clambda)
         if (master) then
            sc_ener(1:ti_ene_cnt) = 0.0d0
            sc_ener_ave(1:ti_ene_cnt) = 0.0d0
            sc_ener_rms(1:ti_ene_cnt) = 0.0d0
            sc_ener_old(1:ti_ene_cnt) = 0.0d0
            sc_ener_old2(1:ti_ene_cnt) = 0.0d0
            enert      = null_state_rec
            enert2     = null_state_rec
            enert_old  = null_state_rec
            enert2_old = null_state_rec
            write (6,*)
            write (6,'(a,f12.4,a,f12.4)') &
                'Dynamically changing lambda: Increased clambda by ', &
                dynlmb, ' to ', clambda
            write (6,*)
         end if
      end if
   end if
#endif
   
   !=======================================================================
   
   !  ---major cycle back to new step unless we have reached our limit:
   
#ifdef MMTSB
   if ( mmtsb_switch /= mmtsb_off ) then
      if ( mod( nstep, mmtsb_iterations ) == 0 ) then
         write(6,'(a,i8)') &
               'MMTSB Replica Exchange iterations completed at NSTEP = ', &
               nstep
         ! apparently 23 is the magic number for potential energy.
         ! ener%pot%tot is the new 23  ;) MJW
         write(6,'(a,f12.4)') &
               'MMTSB Replica Exchange potential energy = ', ener%pot%tot
         ! write coordinates; preferred format is pdb, but can't do that
         ! so write a restart file; server will post process with ambpdb.
#  ifdef LES
         call mdwrit(nstep,nrp,nr,nres,ntxo,ntr,ntb, x,v, &
               xx(lcrdr),box,t,temp0les)
#  else
         call mdwrit(nstep,nrp,nr,nres,ntxo,ntr,ntb, x,v, &
               xx(lcrdr),box,t,rem_val)
#  endif
         ! Chop up trajectory files for later continuous temp splicing.
         call close_dump_files
         ! contact server
         if ( mmtsb_switch == mmtsb_temp_rex ) then
            call mmtsb_newtemp( ener%pot%tot, temp_mmtsb, is_done_mmtsb )
         else if ( mmtsb_switch == mmtsb_lambda_rex ) then
            ! currently temp_mmtsb is ignored, but multidimensional soon
            call mmtsb_newlambda( unpert_pe_mmtsb, pert_pe_mmtsb, &
                  lambda_mmtsb, temp_mmtsb, is_done_mmtsb )
         end if
         call open_dump_files
         if ( is_done_mmtsb ) then
            goto 480
         end if

         ! in the future we may want amber based tracking of exchanges
         ! perhaps we can use the Simmerling group's code ?
         if ( mmtsb_switch == mmtsb_temp_rex ) then
            if ( abs( temp_mmtsb - temp0 ) <= TEN_TO_MINUS3 ) then
               ! no exchange, continue at the same reference temp.
               mmtsb_is_exchanged = .false.
               write(6,'(a,i8,a,f12.4)') &
                     'MMTSB Replica Exchange temperature unchanged'
            else
               ! exchange temp via changing the reference temp.
               ! the velocities will be randomly reset at the new temp via
               ! the resetvelo variable.
               mmtsb_is_exchanged = .true.
               write(6,'(a,f8.2,a,f8.2)') &
                     'MMTSB Replica Exchange temperature change from ', &
                     temp0, ' to ', temp_mmtsb
               temp0 = temp_mmtsb
            end if
         else if ( mmtsb_switch == mmtsb_lambda_rex ) then
            if ( abs( lambda_mmtsb - clambda ) <= TEN_TO_MINUS3 ) then
               ! no exchange, continue at the same lambda
               mmtsb_is_exchanged = .false.
               write(6,'(a,i8,a,f12.4)') &
                     'MMTSB Replica Exchange lambda unchanged'
            else
               ! exchange lambda
               ! the velocities will be randomly reset via
               ! the resetvelo variable.
               mmtsb_is_exchanged = .true.
               write(6,'(a,f8.2,a,f8.2)') &
                     'MMTSB Replica Exchange lambda change from ', &
                     clambda, ' to ', lambda_mmtsb
               clambda = lambda_mmtsb
            end if
         end if  ! ( mmtsb_switch == mmtsb_temp_rex )
      else
         ! not a replica exchange update iteration.
         mmtsb_is_exchanged = .false.
      end if  ! ( mod( nstep, mmtsb_iterations ) == 0 )
   end if  ! ( mmtsb_switch /= mmtsb_off )
#endif


   call trace_integer( 'end of step', nstep )
   call trace_output_mpi_tally( )
   call timer_stop(TIME_VERLET)
#if !defined(DISABLE_NCSU) && defined(NCSU_ENABLE_BBMD)
   call ncsu_on_mdstep(ener%pot%tot, v, ekmh)
#endif /* !defined(DISABLE_NCSU) && defined(NCSU_ENABLE_BBMD) */

#if defined(RISMSANDER) && defined(RISM_DEBUG)
   if(rismprm%irism == 1) then
!!$   write(6,*) "END OF STEP",natom
!     call calc_cm(x,cm,amass,natom)
      angvel=0
      do m=1,natom
         r = x((m-1)*3+1:(m-1)*3+3)-cm
!!$      write(6,*) m,v((m-1)*3+1:(m-1)*3+3)
         call cross(r,v((m-1)*3+1:(m-1)*3+3),rxv)
         angvel = angvel + rxv/sum(r**2)
      end do
      moi=0
      erot=0
      do m=1,natom
         r = x((m-1)*3+1:(m-1)*3+3)-cm
         call cross(r,v((m-1)*3+1:(m-1)*3+3),rxv)
         proj = sum(r*angvel)/sum(angvel**2)*angvel
!!$      write(6,*) "angvel   ",angvel
!!$      write(6,*) "r   ",r,sum((r)**2)
!!$      write(6,*) "proj",proj
!!$      write(6,*) "r-proj",r-proj,sum((r-proj)**2)
         moi=moi+amass(m)*sum((r-proj)**2)
         erot = erot + .5*amass(m)*sum((r-proj)**2)*sum((rxv/sum(r**2))**2)
      end do
!!$   write(6,*) moi
!!$   do m=1,3
!!$      write(6,*) m,sum(v(m:3*natom:3))
!!$      write(6,*) m,sum(amass(1:natom)*v(m:3*natom:3))
!!$      write(6,*) m,angvel(m),sum(angvel**2)
!!$   end do
!!$   write(6,*) "EROT", 0.5*moi*sum(angvel**2), erot
!!$   write(6,*) "EROT", erot
!!$   call mexit(6,1)
   end if
#endif /*RISMSANDER && RISM_DEBUG*/

   if(abfqmmm_param%abfqmmm == 1) then                 ! lam81
#ifdef MPI
    call xdist(v, xx(lfrctmp), natom)                  ! lam81
#endif
    abfqmmm_param%v(1:nr3+iscale) = v(1:nr3+iscale)    ! lam81
    deallocate(for, stat=ier)                          ! lam81 
    return                                             ! lam81
   end if                                              ! lam81

   if (plumed.eq.1 .and. plumed_stopflag/=0) goto 480

   if (nstep < nstlim) goto 260
   480 continue

#ifdef MPI     
! ------====== REMD Post-Dynamics ======------
   if(next_rem_method == 1) then
      remd_ekmh=ekmh
               
      ! ---=== HYBRID REMD ===---
      if (numwatkeep>=0) then 
         ! This is a hybrid REMD run. Get energy of stripped system for next
         !  exchange.
         call hybrid_remd_ene(xx,ix,ih,ipairs,qsetup,      &
                              numwatkeep,hybridgb,igb,ntr,nspm,t,temp0, &
                              ntb,cut,                              &
                              ener,ener%vir,do_list_update,nstep,   &
                              nitp,nits,onefac,loutfm )
      else ! numwatkeep>=0
         ! The positions are currently one step ahead of the energy ener%pot%tot,
         ! since force was called prior to the position propagation. Thus, call 
         ! force one more time to update ener%pot%tot to reflect the current 
         ! coordinates.
         call force(xx,ix,ih,ipairs,x,f,ener,ener%vir, &
                  xx(l96),xx(l97),xx(l98),xx(l99), qsetup, &
                  do_list_update,nstep)
      endif ! numwatkeep>=0

      ! Set myeptot, mytemp, and mytargettemp
!     if (mdloop>0) mytemp = ener%kin%tot * onefac(1)
      my_remd_data%mytemp = ener%kin%tot * onefac(1)
      my_remd_data%myeptot = ener%pot%tot


      my_remd_data%mytargettemp = temp0
#  ifdef VERBOSE_REMD
      if (master) write(6,'(a,f15.4,2(a,f6.2))') &
         "REMD: myEptot= ",my_remd_data%myeptot," myTargetTemp= ", &
         my_remd_data%mytargettemp," mytemp= ",my_remd_data%mytemp
#  endif
#  ifdef LES
   else if(next_rem_method == 2 ) then
      my_remd_data%mytemp       = ener%kin%solv * onefac(3)
      my_remd_data%myeptot      = ener%eptot
      my_remd_data%mytargettemp = temp0les
#  endif
   else if (next_rem_method == 3) then
      remd_ekmh = ekmh
      if (mdloop > 0) my_remd_data%mytemp = ener%kin%tot * onefac(1)
      my_remd_data%mytargettemp = temp0
!    Call force here to bring all energies up-to-date
      call force(xx,ix,ih,ipairs,x,f,ener,ener%vir,xx(l96),xx(l97),xx(l98), &
                 xx(l99),qsetup,do_list_update)

      my_remd_data%myeptot = ener%pot%tot

!   Call nmrdcp to decrement the NMR counter, since this should not count as
!   a real step (JMS 2/12). This is OK, since the counter got incremented at
!   the _very_ end of nmrcal, so we haven't already printed an unwanted value
      if (nmropt /= 0) call nmrdcp
!   Call xdist such that  master has all the  velocities(DSD 09/12)
      call xdist(v, xx(lfrctmp), natom)   
     
   else if (next_rem_method == 4) then
      remd_ekmh = ekmh

   endif ! rem == 1 
! ------====== END REMD Post-Dynamics ======------
#endif /* MPI */

   !=======================================================================
   !     ----- PRINT AVERAGES -----
   !=======================================================================
   
#  ifdef MPI
   ! -- ti decomp
   if (icfe /= 0 .and. idecomp /= 0) then
      if( idecomp == 1 .or. idecomp == 2 ) then
         call collect_dec(nrs)
      !else if( idecomp == 3 .or. idecomp == 4 ) then
      !   call collect_dec(npdec*npdec)
      end if
   end if

   ! Turn off avg. for REMD. and explicit solvent CpHMD, since it's not
   ! accumulated correctly in that case for each compiler
   if (master.and.rem == 0) then
#  else
   if (master) then
#  endif /*MPI*/
      tspan = nvalid
      if (nvalid > 0) then

         ! Update all elements of these sequence types
         enert  = enert/tspan
         enert2 = enert2/tspan - enert*enert
         call zero_neg_values_state(enert2)
         enert2 =  sqrt(enert2)
         edvdl = edvdl/tspan

         ! for PIMD/NMPIMD/CMD/RPMD averages
         if (ipimd>0) then
           totenert = totenert/tspan
           totenert2 = totenert2/tspan - (totenert*totenert)
           call zero_neg_values_state(totenert2)
           totenert2 =  sqrt(totenert2)
         endif

#ifdef MPI
         if( ievb /= 0 ) then
            evb_nrg_ave(:) = evb_nrg_ave(:) / tspan
            evb_nrg_rms(:) = evb_nrg_rms(:) / tspan - evb_nrg_ave(:)**2
            evb_nrg_rms(:) = max( evb_nrg_rms(:), 0.0d0 ) 
            evb_nrg_rms(:) = sqrt( evb_nrg_rms(:) ) 
         endif
         if ( ifsc /= 0 ) then
            do m = 1,ti_ene_cnt
               sc_ener_ave(m) = sc_ener_ave(m)/tspan
               sc_ener_rms(m) = sc_ener_rms(m)/tspan - sc_ener_ave(m)**2
               if(sc_ener_rms(m) < 0.0d0) sc_ener_rms(m) = 0.0d0
               sc_ener_rms(m) = sqrt(sc_ener_rms(m))
            end do
         end if         
         if( ievb /= 0 ) evb_frc%evb_ave = .true.
#endif
         write(6,540) nvalid
         call prntmd(total_nstep,izero,izero,t,enert,onefac,0,.false.)
#ifdef MPI /* SOFT CORE */
         if (ifsc /= 0) call sc_print_energies(6, sc_ener_ave)
         if( ievb /= 0 ) evb_frc%evb_rms = .true.
         if ( ipimd > 0 .and. worldrank==0 ) then
            write(pimd_unit,540) nvalid
            call pimd_report(nstep,t,pimd_unit,totenert,onefac)
            write(pimd_unit,550)
            call pimd_report(nstep,t,pimd_unit,totenert2,onefac)
         endif
#endif
         if (nmropt > 0) call nmrptx(6)
         write(6,550)
         call prntmd(total_nstep,izero,izero,t,enert2,onefac,0,.true.)

#ifdef MPI
         if (ifsc /= 0) call sc_print_energies(6, sc_ener_rms)
         if (ifsc /= 0) call sc_print_dvdl_values()
         
         if( icfe > 0 ) then
            write(6,541) nvalid
            edvdl%pot%dvdl = enert%pot%dvdl  ! fix for DV/DL output
            edvdl%virvsene = 0.d0 ! virvsene should not but included here
            call prntmd(total_nstep,izero,izero,t,edvdl,onefac,0,.false.)
            ! -- ti decomp
            if(worldrank == 0 .and. idecomp /= 0) then
               call checkdec(idecomp)
               if(idecomp == 1 .or. idecomp == 2) call printdec(ix)
            end if
         end if
#endif

         if (nmropt >= 1) then
            write(6,500)
            if (iredir(7) /= 0) call pcshift(-1,x,f)
            call ndvptx(x,f,ih(m04),ih(m02),ix(i02),nres,xx(l95), &
                  natom, xx(lwinv),ntb,xx(lnmr01),ix(inmr02),6)
         end if
         
         ! Print Born radii statistics
         
         if ((rbornstat == 1).and.(igb /= 0 .or. ipb /= 0)) then

            ! Born radii stats collected every nrespai step not nrespa step
            tspan = nvalidi

            write(6,580) nstep
            write(6,590)
            do m = 1,natom
               xx(l188-1+m) = xx(l188-1+m)/tspan
               xx(l189-1+m) = xx(l189-1+m)/tspan - &
                     xx(l188-1+m)*xx(l188-1+m)
               xx(l189-1+m) = sqrt(xx(l189-1+m))
               write(6,600) m, xx(l186-1+m), xx(l187-1+m), &
                     xx(l188-1+m), xx(l189-1+m)
            end do
         end if
         
         enert%kin%tot   = enert%kin%tot*onefac(1)
         enert2%kin%tot  = enert2%kin%tot*onefac(1)
         enert%kin%solt  = enert%kin%solt*onefac(2)
         enert2%kin%solt = enert2%kin%solt*onefac(2)
         enert%kin%solv  = enert%kin%solv*onefac(3)
         enert2%kin%solv = enert2%kin%solv*onefac(3)

         temp = enert%kin%tot
      end if  ! (nvalid > 0)
      
      if (ntp > 0 .and. barostat == 2) call mcbar_summary

   end if  ! (master)

#ifdef MPI
   if( ievb /= 0 ) then
      call evb_dealloc
#if defined(LES)
      if( master ) call evb_pimd_dealloc
#endif
   endif
#endif
   
   if( icfe /= 0 ) then
      deallocate( frcti, stat = ier )
      REQUIRE( ier == 0 )
   end if

   if (plumed.eq.1) then
     call plumed_f_gfinalize()
   end if

   500 format(/,' NMR restraints on final step:'/)
   540 format(/5x,' A V E R A G E S   O V E R ',i7,' S T E P S',/)
   541 format(/5x,' DV/DL, AVERAGES OVER ',i7,' STEPS',/)
   542 format('|',79('='))
   550 format(/5x,' R M S  F L U C T U A T I O N S',/)
   580 format('STATISTICS OF EFFECTIVE BORN RADII OVER ',i7,' STEPS')
   590 format('ATOMNUM     MAX RAD     MIN RAD     AVE RAD     FLUCT')
   600 format(i4,2x,4f12.4)
   call trace_exit( 'runmd' )
   return
end subroutine runmd 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Stripped-down runmd routine for running relaxation dynamics on a given mask
subroutine relaxmd(xx,ix,ih,ipairs,x,winv,amass,f, &
      v,vold,xr,xc,conp,skip,nsp,tma,erstop, qsetup, &
      relax_nstlim, mobile_atoms, increment_nmropt)

   !  Runmd operates in kcal/mol units for energy, amu for masses,
   !     and angstroms for distances.  To convert the input time parameters
   !     from picoseconds to internal units, multiply by 20.455
   !     (which is 10.0*sqrt(4.184)).

   use bintraj, only: end_binary_frame
   use barostats, only : mcbar_trial
   use constants, only : third, ten_to_minus3
   use crg_reloc, only: ifcr, crprintcharges, cr_print_charge
   use fastwt
   use file_io_dat
   use molecule, only: n_iwrap_mask_atoms, iwrap_mask_atoms
   use nblist,only: fill_tranvec,volume,oldrecip,ucell
   use qmmm_module, only : qmmm_nml,qmmm_struct, qmmm_mpi, qm2_struct, &
                           qmmm_vsolv
   use stack
   use state
   use trace

!     Variable Descriptions
!
! Passed variables
!  xx          : global real array. See locmem.f for structure/pointers
!  ix          : global integer array. See locmem.f for structure/pointers
!  ih          : global hollerith array. See locmem.f for structure/pointers
!  ipairs      : ?? Global pairlist ?? --add description (JMS 11/2010)
!  x           : global position array *
!  winv        : array with inverse masses *
!  amass       : mass array *
!  f           : force array, used to hold old coordinates temporarily, too
!  v           : velocity array
!  vold        : old velocity array, from the previous step
!  xr          : coordinates with respect to COM of molecule
!  conp        : bond parameters for SHAKE
!  skip        : logical skip array for SHAKE (and QM/MM too, I think)
!  nsp         : submolecule index array (?)
!  tma         : submolecular weight array (?)
!  erstop      : should we stop in error (?)
!  qsetup      : Not quite sure what this does, if anything anymore.
!  mobile_atoms: bellymask-style array with 1s for moving atoms and 0s for
!                frozen atoms
!  relax_nstlim: Number of relaxation dynamics steps to run
!  increment_nmropt: Do we allow the nmropt counter to increment?
!
! Local variables
!  factt       : degree-of-freedom correction factor for temperature scaling
!  nr          : local copy of nrp, number of atoms
!  nr3         : 3 * nr, used for runtime efficiency
!
! Common memory variables
!  nrp         : number of atoms, adjusted for LES copies

   implicit none
   character(kind=1,len=7) :: routine="relaxmd"
   integer   ipairs(*), ix(*), relax_nstlim
   integer, intent(in) :: mobile_atoms(*)
   logical, intent(in) :: increment_nmropt
   _REAL_ xx(*)
   character(len=4) ih(*)
   _REAL_ combination

#ifdef MPI
#  include "parallel.h"
   include 'mpif.h'
   _REAL_ mpitmp(8) !Use for temporary packing of mpi messages.
   integer ist(MPI_STATUS_SIZE), partner, ierr
#else
   ! mdloop and REM is always 0 in serial
   integer, parameter :: mdloop = 0, rem = 0
#endif

#include "../include/md.h"
#include "box.h"
#include "nmr.h"
#include "../include/memory.h"
#include "extra.h"
#include "ew_frc.h"
#include "ew_cntrl.h"
#include "ew_mpole.h"
#include "def_time.h"
#include "extra_pts.h"
#include "../lib/random.h"   
      
   _REAL_ sysx,sysy,sysz,sysrange(3,2)
   logical mv_flag

   _REAL_ , dimension(1) :: shkh
   integer, dimension(1) :: ifstwr2
   integer :: nshkh

   integer idx, iatom, iatomCL,m
   _REAL_  Ekin2_tot,tmp,f_lnv
   integer :: idim, ithermo
   _REAL_ :: E_nhc, exp1, exp2, v_sum
   
   logical ivscm
   logical qspatial
   character(len=6)fnam

   logical resetvelo
   integer nshak
   _REAL_ ekgs,eold3,eold4,etot_save,ekpbs
   
   logical do_list_update
   logical skip(*),lout,loutfm,erstop,vlim,onstep
   _REAL_ x(*),winv(*),amass(*),f(*),v(*),vold(*), &
         xr(*),xc(*),conp(*)
   type(state_rec) :: ener
   type(state_rec) :: ecopy, edvdl
   type(state_rec) :: edvdl_r
   _REAL_ rmu(3),fac(3),onefac(3),clfac, etot_start
   _REAL_ tma(*)

   _REAL_ tspan,atempdrop,fln,scaltp,scaltpo
   _REAL_ vel,vel2,vcmx,vcmy,vcmz,vmax,vx,vy,vz
   _REAL_ winf,aamass,rterm,ekmh,ekph,ekpht,wfac,rsd,ekav
   _REAL_ fit,fiti,fit2,vscalt
   logical is_langevin  ! Is this a Langevin dynamics simulation
   _REAL_ gammai,c_implic,c_explic,c_ave,sdfac,ekins0
   _REAL_ dtx,dtxinv,dt5,factt,ekin0,ekinp0,dtcp,dttp
   _REAL_ rndf,rndfs,rndfp,boltz2,pconv,tempsu
   _REAL_ xcm(3),acm(3),ocm(3),vcm(3),ekcm,ekrot
   _REAL_ emtmd

! Variables and parameters for constant surface tension:
   _REAL_, parameter :: ten_conv = 100.0d0 !ten_conv - converts
                                                    !dyne/cm to bar angstroms
   _REAL_  :: pres0x
   _REAL_  :: pres0y
   _REAL_  :: pres0z
   _REAL_  :: gamma_ten_int
   _REAL_  :: press_tan_ave

   integer nsp(*)
   integer idumar(4)
   integer l_temp
   integer i,j,im,i3,nitp,nits
   integer nstep,nrep,nrek,nren,iend,istart3,iend3
   integer nrx,nr,nr3,ntcmt,izero,istart
   logical qsetup
   
   integer nvalid, nvalidi
   _REAL_ eke,eket
   _REAL_ extent

   _REAL_ xcen,ycen,zcen,extents(3,2)
   _REAL_, allocatable, dimension(:) :: frcti
   integer ier

   _REAL_ small
   data small/1.0d-7/
   data nren/51/

   !--- VARIABLES FOR DIPOLE PRINTING ---
   integer prndipngrp
   integer prndipfind
   character(len=4) prndiptest

   _REAL_,parameter :: pressure_constant = 6.85695d+4
   ! variables used in constant pressure PIMD
   _REAL_ :: Nkt,centvir,pressure, aa, arg2, poly, e2, e4, e6, e8 
   ! variable used in CMD
   real(8) :: tmp_eke_cmd !Use for temporary packing of mpi messages.

   _REAL_ :: box_center(3)

   !==========================================================================
  
   call trace_enter( 'relaxmd' )

   !     ----- INITIALIZE SOME VARIABLES -----
   
   vlim = vlimit > small
   ntcmt = 0
   izero = 0
   lout = .true.
   loutfm = ioutfm <= 0
   nr = natom
   nr3 = 3*nr
   ekmh = 0.d0
   onstep = .true.

   do_list_update=.false.
#ifdef MPI
   istart = iparpt(mytaskid) + 1
   iend = iparpt(mytaskid+1)
#else
   istart = 1
   iend = nr
#endif
   istart3 = 3*istart -2
   iend3 = 3*iend

   ! If NTWPRT.NE.0, only print the atoms up to this value
   nrx  = nr3
   if (ntwprt > 0) nrx = ntwprt*3
   
   !=======================================================================
   ! Determine system degrees of freedom (for T scaling, reporting)
   
   ! Call DEGCNT to get the actual number of degrees of freedom for the
   ! solute and solvent. The 'belly' atoms are just the mobile ones
   
   call degcnt(1,nr,mobile_atoms,nsolut,nbonh,nbona,0, &
         ix(iibh),ix(ijbh),ix(iiba),ix(ijba),idumar, &
         idumar,ntc,idumar,0,0,0, &
         idumar,rndfp,rndfs)
   
   ! RNDFP = # degrees of freedom for solute
   ! RNDFS = # degrees of freedom for solvent
   ! RNDF = total number of degrees of freedom.

   ! qtw - substract the number of overlapping noshake QM atoms in noshakemask
   rndfp = rndfp - qmmm_struct%noshake_overlap
   !    modify RNDFP to reflect NDFMIN (set in mdread) and num_noshake
   rndfp = rndfp - ndfmin + num_noshake
   rndf = rndfp+rndfs

   call fix_degree_count(rndf) ! correct for extra points
   
   ! End of degrees of freedom stuff
   !=======================================================================
   
   boltz2 = 8.31441d-3 * 0.5d0
   pconv = 1.6604345d+04  ! factor to convert the pressure kcal/mole to bar
   
   !     ---convert to kcal/mol units
   
   boltz2 = boltz2/4.184d0   ! k-sub-B/2
   dtx = dt*20.455d+00
   dtxinv = 1.0d0 / dtx
   dt5 = dtx * 0.5d0
   pconv = pconv*4.184d0
   
   ! FAC() are #deg freedom * kboltz / 2
   ! multiply by T to get expected kinetic energy
   ! FAC(1) is for total system
   
   fac(1) = boltz2*rndf
   fac(2) = boltz2*rndfp

   if(rndfp < 0.1d0) fac(2) = 1.d-6

   fac(3) = boltz2*rndfs
   if(rndfs < 0.1d0) fac(3) = 1.d-6
   onefac(1) = 1.0d0/fac(1)
   onefac(2) = 1.0d0/fac(2)
   onefac(3) = 1.0d0/fac(3)
   factt = rndf/(rndf+ndfmin)
   
   ! these are "desired" kinetic energies based on
   ! # degrees freedom and target temperature
   ! they will be used for calculating the velocity scaling factor
   
   ekinp0 = fac(2)*temp0
   ekins0 = fac(3)*temp0
   ekin0  = fac(1)*temp0

   ! Langevin dynamics setup:

   is_langevin = gamma_ln > 0.0d0
   gammai = gamma_ln/20.455d0
   c_implic = 1.d0/(1.d0+gammai*dt5)
   c_explic = 1.d0 - gammai*dt5
   c_ave    = 1.d0+gammai*dt5
   sdfac = sqrt( 4.d0*gammai*boltz2*temp0/dtx )
   if (is_langevin .and. ifbox==0) then
      call get_position(nr,x,sysx,sysy,sysz,sysrange,0)
   end if
   if (ntt == 1) dttp = dt/tautp
   if (ntp > 0) dtcp = comp * 1.0d-06 * dt / taup

   ! Constant surface tension setup:

   if (csurften > 0) then

      ! Set pres0 in direction of surface tension.
      ! The reference pressure is held constant in on direction dependent
      ! on what the surface tension direction is set to.
      if (csurften .eq. 1) then           ! pres0 in the x direction
        pres0x = pres0

      else if (csurften .eq. 2) then      ! pres0 in the y direction
        pres0y = pres0

      !else if (csurften .eq. 3) then      ! pres0 in the z direction
      else
        pres0z = pres0

      end if

      ! Multiply surface tension by the number of interfaces
      gamma_ten_int = dble(ninterface) * gamma_ten
 
   end if

   nrek = 4
   nrep = 15
   
   nvalid = 0
   nvalidi = 0
   nstep = 0
   fit = 0.d0
   fiti = 0.d0
   fit2 = 0.d0

   ener       = null_state_rec ! Zeros all elements
   ener%kin%pres_scale_solt = 1.d0
   ener%kin%pres_scale_solv = 1.d0
   ener%box(1:3) = box(1:3)
   ener%cmt(1:4) = 0.d0
   nitp = 0
   nits = 0
   
   ekmh = 0.0d0

   i3 = 0
   do j = 1,nrp
      aamass = amass(j)
      do m = 1,3
         i3 = i3+1
         rterm = v(i3)*v(i3) * aamass
         ekmh = ekmh + rterm
      end do
   end do

   do im=1,iscale
      ekmh = ekmh + scalm*v(nr3+im)*v(nr3+im)
   end do
   ekmh = ekmh * 0.5d0
   do i=1,nr3+iscale
      vold(i) = v(i)
   end do
   
   !=======================================================================
   !     ----- MAIN LOOP FOR PERFORMING THE DYNAMICS STEP -----
   !           (at this point, the coordinates are a half-step "ahead"
   !           of the velocities; the variable EKMH holds the kinetic
   !           energy at these "-1/2" velocities, which are stored in
   !           the array VOLD.)
   !=======================================================================
   
   260 continue

   !---------------------------------------------------------------
   !  ---Step 1a: do some setup for pressure calculations:
   !---------------------------------------------------------------

   if (ntp > 0) then
      ener%cmt(1:3) = 0.d0
      xr(1:nr3) = x(1:nr3)
      
      ! ----- CALCULATE THE CENTER OF MASS ENERGY AND THE COORDINATES
      !       OF THE SUB-MOLECULES WITH RESPECT TO ITS OWN CENTER OF
      !       MASS -----
      
      call timer_start(TIME_EKCMR)
      call ekcmr(nspm,nsp,tma,ener%cmt,xr,v,amass,istart,iend)
#ifdef MPI
      call trace_mpi('mpi_allreduce', &
            3,'MPI_DOUBLE_PRECISION',mpi_sum)
# ifdef USE_MPI_IN_PLACE
      call mpi_allreduce(MPI_IN_PLACE,ener%cmt,3, &
            MPI_DOUBLE_PRECISION,mpi_sum,commsander,ierr)
# else
      call mpi_allreduce(ener%cmt,mpitmp,3, &
            MPI_DOUBLE_PRECISION,mpi_sum,commsander,ierr)
      ener%cmt(1:3) = mpitmp(1:3)
# endif
#endif
      call timer_stop(TIME_EKCMR)
   end if

   ! If we're using the MC barostat, go ahead and do the trial move now
   if (ntp > 0 .and. barostat == 2 .and. mod(nstep+1, mcbarint) == 0) &
      call mcbar_trial(xx, ix, ih, ipairs, x, xc, f, ener%vir, xx(l96), &
               xx(l97), xx(l98), xx(l99), qsetup, do_list_update, nstep, nsp, &
               amass)

   !--------------------------------------------------------------
   !  ---Step 1b: Get the forces for the current coordinates:
   !--------------------------------------------------------------

   iprint = 0
   if( nstep == 0 .or. nstep+1 == relax_nstlim ) iprint = 1

   call force(xx,ix,ih,ipairs,x,f,ener,ener%vir, &
            xx(l96),xx(l97),xx(l98),xx(l99), qsetup, &
            do_list_update,nstep)

   ! If we don't want to increment the NMROPT counter, decrement it here.
   if (.not. increment_nmropt) &
      call nmrdcp

   ! Reset quantities depending on TEMP0 and TAUTP (which may have been
   ! changed by MODWT during FORCE call).
   ekinp0 = fac(2)*temp0
   ekins0 = fac(3)*temp0
   ekin0 = fac(1)*temp0

   if (ntt == 1) dttp = dt/tautp
   
   if (ntp > 0) then
      ener%volume = volume
      ener%density = tmass / (0.602204d0*volume)
      ener%cmt(4) = 0.d0
      ener%vir(4) = 0.d0
      ener%pres(4) = 0.d0
      do m = 1,3
         ener%cmt(m)  = ener%cmt(m)*0.5d0
         ener%cmt(4)  = ener%cmt(4)+ener%cmt(m)
         ener%vir(4)  = ener%vir(4)+ener%vir(m)
         ener%pres(m) = (pconv+pconv)*(ener%cmt(m)-ener%vir(m))/volume
         ener%pres(4) = ener%pres(4)+ener%pres(m)
      end do
      ener%pres(4) = ener%pres(4)/3.d0

      ! Constant surface tension output:

      if (csurften > 0) then

        if (csurften == 1) then          ! Surface tension in the x direction
          ener%surface_ten = &
                box(1) * (ener%pres(1) - 0.5d0 * &
                (ener%pres(2) + ener%pres(3))) / (ninterface * ten_conv)

        else if (csurften .eq. 2) then   ! Surface tension in the y direction
          ener%surface_ten = &
                box(2) * (ener%pres(2) - 0.5d0 * &
                (ener%pres(1) + ener%pres(3))) / (ninterface * ten_conv)

        else ! if (csurften .eq. 3) then ! Surface tension in the z direction
          ener%surface_ten = &
                box(3) * (ener%pres(3) - 0.5d0 * &
                (ener%pres(1) + ener%pres(2))) / (ninterface * ten_conv )

        end if

      end if

   end if

   !----------------------------------------------------------------
   !  ---Step 1c: do randomization of velocities, if needed:
   !----------------------------------------------------------------
   ! ---Assign new random velocities every Vrand steps, if ntt=2

   resetvelo=.false.
   if (vrand /= 0 .and. ntt == 2) then
      if (mod((nstep+1),vrand) == 0) resetvelo=.true.
   end if

   if (resetvelo) then
      ! DAN ROE: Why are only the masters doing this? Even if the velocities 
      !  are broadcast to the child processes, the wont the different # of random
      !  calls put the randomg num generators out of sync, or do we not care?
      
      if (master) then
!        write (6,'(a,i8)') 'Setting new random velocities at step ', &
!              nstep + 1
         call setvel(nr,v,winv,temp0*factt,init,iscale,scalm)
      end if

# ifdef MPI
      call trace_mpi('mpi_bcast',3*natom,'MPI_DOUBLE_PRECISION',0)
      call mpi_bcast(v, 3*natom, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
# endif

      ! At this point in the code, the velocities lag the positions
      ! by half a timestep.  If we intend for the velocities to be drawn 
      ! from a Maxwell distribution at the timepoint where the positions and 
      ! velocities are synchronized, we have to correct these newly 
      ! redrawn velocities by backing them up half a step using the 
      ! current force.
      ! Note that this fix only works for Newtonian dynamics.
      if( gammai==0.d0 ) then                         
         i3 = 3*(istart-1)
         do j=istart,iend
            wfac = winv(j) * dt5
            v(i3+1) = v(i3+1) - f(i3+1)*wfac
            v(i3+2) = v(i3+2) - f(i3+2)*wfac
            v(i3+3) = v(i3+3) - f(i3+3)*wfac
            i3 = i3+3
         end do
      end if
      
   end if  ! (resetvelo)

   call timer_start(TIME_VERLET)

   !-----------------------------------------------------
   !  ---Step 2: Do the velocity update:
   !-----------------------------------------------------
   
   !step 2a: apply quenched MD if needed.  This is useful in NEB>0
   if (vv==1) call quench(f,v)
   
   if( gammai == 0.d0 ) then

      !       ---Newtonian dynamics:
      
      i3 = 3*(istart-1)
      do j=istart,iend
         wfac = winv(j) * dtx
         v(i3+1) = v(i3+1) + f(i3+1)*wfac
         v(i3+2) = v(i3+2) + f(i3+2)*wfac
         v(i3+3) = v(i3+3) + f(i3+3)*wfac
         i3 = i3+3
      end do
      
   else  !  gamma_ln .ne. 0, which also implies ntt=3 (see mdread.f)

      !       ---simple model for Langevin dynamics, basically taken from
      !          Loncharich, Brooks and Pastor, Biopolymers 32:523-535 (1992),
      !          Eq. 11.  (Note that the first term on the rhs of Eq. 11b
      !          should not be there.)

      !  Update Langevin parameters, since temp0 might have changed:
         sdfac = sqrt( 4.d0*gammai*boltz2*temp0/dtx )

      i3 = 3*(istart-1)

      if (no_ntt3_sync == 1) then
        !We don't worry about synchronizing the random number stream
        !across processors.
        do j=istart,iend

           wfac = winv(j) * dtx
           aamass = amass(j)         
           rsd = sdfac*sqrt(aamass)
           call gauss( 0.d0, rsd, fln )
           v(i3+1) = (v(i3+1)*c_explic + (f(i3+1)+fln)*wfac) * c_implic
           call gauss( 0.d0, rsd, fln )
           v(i3+2) = (v(i3+2)*c_explic + (f(i3+2)+fln)*wfac) * c_implic
           call gauss( 0.d0, rsd, fln )
           v(i3+3) = (v(i3+3)*c_explic + (f(i3+3)+fln)*wfac) * c_implic
           i3 = i3+3
        end do

      else

        do j=1,nr
           if( j<istart .or. j>iend ) then
           ! In order to generate the same sequence of pseudorandom numbers that
           ! you would using a single processor you have to go through the atoms 
           ! in order.  The unused results are thrown away
              call gauss( 0.d0, 1.d0, fln )
              call gauss( 0.d0, 1.d0, fln )
              call gauss( 0.d0, 1.d0, fln )
              cycle
           end if

           wfac = winv(j) * dtx
           aamass = amass(j)         
           rsd = sdfac*sqrt(aamass)
           call gauss( 0.d0, rsd, fln )
           v(i3+1) = (v(i3+1)*c_explic + (f(i3+1)+fln)*wfac) * c_implic
           call gauss( 0.d0, rsd, fln )
           v(i3+2) = (v(i3+2)*c_explic + (f(i3+2)+fln)*wfac) * c_implic
           call gauss( 0.d0, rsd, fln )
           v(i3+3) = (v(i3+3)*c_explic + (f(i3+3)+fln)*wfac) * c_implic
           i3 = i3+3
        end do
      end if ! no_ntt3_sync

   end if  ! ( gammai == 0.d0 )

   if (vlim) then
      vmax = 0.0d0
      do i=istart3,iend3
         vmax = max(vmax,abs(v(i)))
         v(i) = sign(min(abs(v(i)),vlimit),v(i))
      end do

      ! Only violations on the master node are actually reported
      ! to avoid both MPI communication and non-master writes.
      if (vmax > vlimit) then
         if (master) then
            write(6,'(a,i6,a,f10.4)') 'vlimit exceeded for step ',nstep, &
              '; vmax = ',vmax
         end if
      end if
   end if
   
   do im=1,iscale
      v(nr3+im) = (v(nr3+im) + f(nr3+im)*dtx/scalm)
   end do

   !-------------------------------------------------------------------
   !   Step 3: update the positions, putting the "old" positions into F:
   !-------------------------------------------------------------------

   i = istart - 1
   do i3 = istart3, iend3, 3
      f(i3  ) = x(i3  )
      f(i3+1) = x(i3+1)
      f(i3+2) = x(i3+2)
      if (mobile_atoms(i) == 1) then
         x(i3  ) = x(i3  ) + v(i3  )*dtx
         x(i3+1) = x(i3+1) + v(i3+1)*dtx
         x(i3+2) = x(i3+2) + v(i3+2)*dtx
      end if
      i = i + 1
   end do

   do i = 1,iscale
      f(nr3+i) = x(nr3+i)
      x(nr3+i) = x(nr3+i)+v(nr3+i)*dtx
   end do

   call timer_stop(TIME_VERLET)
   
   if (ntc /= 1) then

      !-------------------------------------------------------------------
      !   Step 4a: if shake is being used, update the new positions to fix
      !      the bond lengths.
      !-------------------------------------------------------------------
   
      call timer_start(TIME_SHAKE)
      qspatial=.false.
      call shake(nrp,nbonh,nbona,0,ix(iibh),ix(ijbh),ix(ibellygp), &
                 winv,conp,skip,f,x,nitp,.false.,ix(iifstwt),ix(noshake), &
                 shkh,qspatial)
      call quick3(f,x,ix(iifstwr),natom,nres,ix(i02))
      if(nitp == 0) then
         erstop = .true.
         goto 480
      end if

      !-----------------------------------------------------------------
      !   Step 4b:   Now fix the velocities and calculate KE
      !-----------------------------------------------------------------
      
      !  ---re-estimate the velocities from differences in positions:

      v(istart3:iend3) = (x(istart3:iend3)-f(istart3:iend3)) * dtxinv

      call timer_stop(TIME_SHAKE)
   end if
   call timer_start(TIME_VERLET)

   if( ntt == 1 .or. onstep ) then

      !-----------------------------------------------------------------
      !   Step 4c: get the KE, either for averaging or for Berendsen:
      !-----------------------------------------------------------------

      eke = 0.d0
      ekph = 0.d0
      ekpbs = 0.d0

      if (gammai == 0.0d0) then
         i3 = 3*(istart-1)
         do j=istart,iend
            aamass = amass(j)
            do m = 1,3
               i3 = i3+1
               eke = eke + aamass*0.25d0*(v(i3)+vold(i3))**2

               ! try pseudo KE from Eq. 4.7b of Pastor, Brooks & Szabo,
               ! Mol. Phys. 65, 1409-1419 (1988):

               ekpbs = ekpbs + aamass*v(i3)*vold(i3)
               ekph = ekph + aamass*v(i3)**2

            end do
         end do

      else

         i3 = 3*(istart-1)
         do j=istart,iend
            aamass = amass(j)
            do m = 1,3
               i3 = i3+1
               eke = eke + aamass*0.25d0*c_ave*(v(i3)+vold(i3))**2
            end do

         end do

      end if ! (if gammai == 0.0d0)

#ifdef MPI
      !  ---   sum up the partial kinetic energies:

      if ( numtasks > 1 ) then
         call trace_mpi('mpi_allreduce', &
               1,'MPI_DOUBLE_PRECISION',mpi_sum)
         mpitmp(1) = eke
         mpitmp(2) = ekph
         mpitmp(3) = ekpbs
#    ifdef USE_MPI_IN_PLACE
         call mpi_allreduce(MPI_IN_PLACE,mpitmp,3, &
               MPI_DOUBLE_PRECISION,mpi_sum,commsander,ierr)
         eke = mpitmp(1)
         ekph = mpitmp(2)
         ekpbs = mpitmp(3)

#    else /* USE_MPI_IN_PLACE */

         call mpi_allreduce(mpitmp,mpitmp(4),3, &
               MPI_DOUBLE_PRECISION,mpi_sum,commsander,ierr)
         eke = mpitmp(4)
         ekph = mpitmp(5)
         ekpbs = mpitmp(6)

#    endif /* USE_MPI_IN_PLACE */
      end if
#endif /* MPI */
      
      !         --- all processors handle the "extra" variables:

      do im=1,iscale
         eke = eke + scalm*0.25d0*(v(nr3+im)+vold(nr3+im))**2
         ekpbs = ekpbs + scalm*v(nr3+im)*vold(nr3+im)
         ekph = ekph + scalm*v(nr3+im)**2
      end do
      
      eke = eke * 0.5d0
      ekph = ekph * 0.5d0
      ekpbs = ekpbs * 0.5d0
      if( ntt == 1 ) then
         
         !    --- following is from T.E. Cheatham, III and B.R. Brooks,
         !        Theor. Chem. Acc. 99:279, 1998.
         
         scaltp = sqrt(1.d0 + 2.d0*dttp*(ekin0-eke)/(ekmh+ekph))

         !    --- following is the "old" (amber7 and before) method:

         !  scaltpo = sqrt(1.d0 + dttp*(ekin0/ekph - 1.d0))
         !  write(6,*) 'scaltp: ',2.d0*dttp*(ekin0-eke)/(ekmh+ekph), &
         !            dttp*(ekin0/ekmh - 1.d0)

         !  following line reverts to the "old" behavior:
         !  scaltp = scaltpo

         do j = istart,iend
            i3=(j-1)*3+1
            v(i3  ) = v(i3  ) *scaltp
            v(i3+1) = v(i3+1) *scaltp
            v(i3+2) = v(i3+2) *scaltp
         end do
         do im=1,iscale
            v(nr3+im) = v(nr3+im)*scaltp
         end do
      end if  ! (ntt == 1 )
      
   end if  ! ( ntt == 1 .or. onstep; end of step 4c )

   !-----------------------------------------------------------------
   !   Step 5:  several tasks related to dumping of trajectory information
   !-----------------------------------------------------------------
   
   !  --- Determine if trajectory, velocity, or restart
   !      writing is imminent, or if the center of mass
   !      motion will be removed.
   !      These require xdist of velocities or dipoles in parallel runs:
   !
   ! Modified so that when running REMD, writing can occur less often
   !  than exchanges (e.g. ntwx > nstlim)
   ! DAN ROE: Added two new variables, total_nstep and total_nstlim.
   !          For non-REMD runs, total_nstep=nstep+1 and total_nstlim=nstlim 
   !           just like before.
   !          For REMD runs, total_nstep=(mdloop-1)*nstlim+nstep+1, where
   !           mdloop is the current exchange - this is the current
   !           replica exchange MD step. total_nstlim=numexchg*nstlim, which is
   !           the maximum number of REMD steps.

#ifdef MPI

   !-----------------------------------------------------------------
   !  --- now distribute the coordinates, and if necessary, dipoles and vel:
   !-----------------------------------------------------------------

   call timer_barrier( commsander )
   call timer_stop_start(TIME_VERLET,TIME_DISTCRD)
   if ( numtasks > 1 ) then
      call xdist(x, xx(lfrctmp), natom)
   end if
   call timer_stop(TIME_DISTCRD)

#endif  /* MPI */

   !           ----fix lone pair positions:
   if( numextra > 0 )call local_to_global(x,xx,ix)

#ifdef MPI
   call timer_start(TIME_VERLET)
   !     ========================= END AMBER/MPI =========================
#endif  /* MPI */
   
   !-------------------------------------------------------------------
   !   Step 6: zero COM velocity if requested; used for preventing
   !   ewald "block of ice flying thru space" phenomenon, or accumulation
   !   of rotational momentum in vacuum simulations
   !-------------------------------------------------------------------

   !-----------------------------------------------------------------
   !  --- put current velocities into VOLD
   !-----------------------------------------------------------------
   
   vold(istart3:iend3) = v(istart3:iend3)
   do im=1,iscale
      vold(nr3+im) = v(nr3+im)
   end do

   !-------------------------------------------------------------------
   !  Step 7: scale coordinates if constant pressure run:
   !-------------------------------------------------------------------

      ! ntp = 1, isotropic pressure coupling

   if (ntp == 1) then
      rmu(1) = (1.d0-dtcp*(pres0-ener%pres(4)))**third
      rmu(2) = rmu(1)
      rmu(3) = rmu(1)


     ! ntp = 2, anisotropic pressure scaling

   else if (ntp == 2) then

      if (csurften > 0) then

         ! Constant surface tension adjusts the tangential pressures
         ! See Zhang, Feller, Brooks, Pastor. J. Chem. Phys. 1995

         if (csurften == 1) then        ! For surface tension in the x direction
            pres0y = pres0x - gamma_ten_int * ten_conv / box(1)
            pres0z = pres0y

         else if (csurften == 2) then   ! For surface tension in the y direction
            pres0x = pres0y - gamma_ten_int * ten_conv / box(2)
            pres0z = pres0x

         !else if (csurften == 3) then   ! For surface tension in the z !direction
         else
            pres0x = pres0z - gamma_ten_int * ten_conv / box(3)
            pres0y = pres0x

         end if

         rmu(1) = (1.d0 - dtcp * (pres0x - ener%pres(1)))**third
         rmu(2) = (1.d0 - dtcp * (pres0y - ener%pres(2)))**third
         rmu(3) = (1.d0 - dtcp * (pres0z - ener%pres(3)))**third

      else

         rmu(1) = (1.d0-dtcp*(pres0-ener%pres(1)))**third
         rmu(2) = (1.d0-dtcp*(pres0-ener%pres(2)))**third
         rmu(3) = (1.d0-dtcp*(pres0-ener%pres(3)))**third

      end if

     ! ntp = 3, semiisotropic pressure coupling
     ! (currently only for csurften>0, constant surface tension)

     !else if (ntp > 2) then
     else

        if (csurften > 0) then

           if (csurften == 1) then        ! For surface tension in the x direction
              pres0y = pres0x - gamma_ten_int * ten_conv / box(1)
              pres0z = pres0y
              press_tan_ave = (ener%pres(2) + ener%pres(3))/2
              rmu(1) = (1.d0 - dtcp * (pres0x - ener%pres(1)))**third
              rmu(2) = (1.d0 - dtcp * (pres0y - press_tan_ave))**third
              rmu(3) = (1.d0 - dtcp * (pres0z - press_tan_ave))**third

           else if (csurften == 2) then   ! For surface tension in the y direction
              pres0x = pres0y - gamma_ten_int * ten_conv / box(2)
              pres0z = pres0x
              press_tan_ave = (ener%pres(1) + ener%pres(3))/2
              rmu(1) = (1.d0 - dtcp * (pres0x - press_tan_ave))**third
              rmu(2) = (1.d0 - dtcp * (pres0y - ener%pres(2)))**third
              rmu(3) = (1.d0 - dtcp * (pres0z - press_tan_ave))**third

            !else if (csurften == 3) then   ! For surface tension in the z !direction
           else
              pres0x = pres0z - gamma_ten_int * ten_conv / box(3)
              pres0y = pres0x
              press_tan_ave = (ener%pres(1) + ener%pres(2))/2
              rmu(1) = (1.d0 - dtcp * (pres0x - press_tan_ave))**third
              rmu(2) = (1.d0 - dtcp * (pres0y - press_tan_ave))**third
              rmu(3) = (1.d0 - dtcp * (pres0z - ener%pres(3)))**third

           end if ! csurften == 1
        end if ! csurften > 0
        ! Add semiisotropic pressure scaling in any direction with no constant 
        ! surface tension here
      end if

      if (ntp > 0) then
         box(1:3) = box(1:3)*rmu(1:3)
         ener%box(1:3) = box(1:3)
         
         !    WARNING!!   This is not correct for non-orthogonal boxes if
         !    NTP > 1 (i.e. non-isotropic scaling).  Currently general cell
         !    updates which allow cell angles to change are not implemented.
         !    The viral tensor computed for ewald is the general Nose Klein,
         !    however the cell response needs a more general treatment.
      
         call redo_ucell(rmu)
         ! keep tranvec up to date, rather than recomputing each MD step.
         call fill_tranvec()  ! tranvec is dependent on only ucell

         call ew_pscale(natom,x,amass,nspm,nsp,npscal)
         if (ntr > 0 .and. nrc > 0) &
            call ew_pscale(natom,xc,amass,nspm,nsp,npscal)
      endif
   
   ener%kin%solv = ekpbs + ener%pot%tot  
                               ! Pastor, Brooks, Szabo conserved quantity
                               ! for harmonic oscillator: Eq. 4.7b of Mol.
                               ! Phys. 65:1409-1419, 1988
   ener%kin%solt = eke
   ener%kin%tot  = ener%kin%solt
   if (ntt == 1 .and. onstep) then
      ekmh = max(ekph,fac(1)*10.d0)
   end if
   
   !     ---if velocities were reset, the KE is not accurate; fudge it
   !        here to keep the same total energy as on the previous step.
   !        Note that this only affects printout and averages for Etot
   !        and KE -- it has no effect on the trajectory, or on any averages
   !        of potential energy terms.
   
   if( resetvelo ) ener%kin%tot = etot_save - ener%pot%tot
   
   !     --- total energy is sum of KE + PE:
   
   ener%tot  = ener%kin%tot + ener%pot%tot
   etot_save = ener%tot

   !-------------------------------------------------------------------
   !  Step 8:  update the step counter and the integration time:
   !-------------------------------------------------------------------
   
   nstep = nstep+1

   !     ---full energies are only calculated every nrespa steps
   !     nvalid is the number of steps where all energies are calculated
   
   ntnb = 0
   if (mod(nstep,nsnb) == 0) ntnb = 1

#if 0
! DEBUG code -- this will print out every frame of the relaxation dynamics to
! the trajectory if uncommented

   if (master) then
      
      !     -- Coordinate archive:
      if (.true.) then
         if( iwrap == 0 ) then
            call corpac(x,1,nrx,MDCRD_UNIT,loutfm)
            if(ntb > 0)  call corpac(box,1,3,MDCRD_UNIT,loutfm)
         else if (iwrap == 1) then
            call get_stack(l_temp,nr3,routine)
            if(.not. rstack_ok)then
               deallocate(r_stack)
               allocate(r_stack(1:lastrst),stat=alloc_ier)
               call reassign_rstack(routine)
            endif
            REQUIRE(rstack_ok)
            do m=1,nr3
               r_stack(l_temp+m-1) = x(m)
            end do
            
            call wrap_molecules(nspm,nsp,r_stack(l_temp))
            if (ifbox == 2) call wrap_to(nspm,nsp,r_stack(l_temp),box)
            
            call corpac(r_stack(l_temp),1,nrx,MDCRD_UNIT,loutfm)
            call corpac(box,1,3,MDCRD_UNIT,loutfm)
            call free_stack(l_temp,routine)
          end if ! if (iwrap == 0) ...


         !GMS: If using variable QM solvent, try to write a new pdb file
         !     with the QM coordinates for this step. This is done here
         !     to keep the PDB file in sync with the mdcrd file, which
         !     makes it easier to check later.
         if (qmmm_nml%vsolv > 0 .and. qmmm_nml%verbosity == 0) &
            call qm_print_coords(nstep,.false.)
      end if  ! (itdump)

      if (ioutfm > 0) then
         if (.true.) call end_binary_frame(MDCRD_UNIT)
      end if

   end if  ! (master)
#endif /* 0 */

   !=======================================================================
   
   !  ---major cycle back to new step unless we have reached our limit:

   call trace_integer( 'end of step', nstep )
   call trace_output_mpi_tally( )
   call timer_stop(TIME_VERLET)

   if (nstep < relax_nstlim) goto 260
   480 continue

end subroutine relaxmd

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Enter description for quench here.
subroutine quench(f,v) 
   implicit none
   
#include "../include/md.h" 
!need access to vv - temp verlet scaling
#include "../include/memory.h" 
!need access to natom

   _REAL_ f(*),v(*),dotproduct,force
   !f is the forces and v is the velocity
 
   integer index
   dotproduct = 0.d0
   force = 0.d0
 
   do index=1,3*natom
      force = force + f(index)**2
      dotproduct = dotproduct + v(index)*f(index)
   enddo
  
   if (force/=0.0d0) then 
     force = 1.0d0/sqrt(force)
     dotproduct = dotproduct*force
   end if
   
   if (dotproduct>0.0d0) then
      v(1:3*natom) = dotproduct*f(1:3*natom)*force
   else 
      !v(1:3*natom) = 0.0d0
      v(1:3*natom) = vfac*dotproduct*f(1:3*natom)*force
   end if
   
end subroutine quench
