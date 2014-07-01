!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE plugin_forces()
  !----------------------------------------------------------------------------
  !
  !
  USE mp_global,        ONLY : intra_image_comm
  USE mp,               ONLY : mp_bcast
  USE io_global,        ONLY : stdout, ionode, ionode_id
  USE kinds,            ONLY : DP
  USE io_files,         ONLY : outdir
  !
  USE plugin_flags
  !
  USE cell_base,        ONLY : alat, at
  USE ions_base,        ONLY : tau, nat,amass
  USE force_mod,        ONLY : force,sigma
  USE control_flags,    ONLY : istep
  USE ener,             ONLY : etot 
  !
  IMPLICIT NONE
  !
  INTEGER:: i,j
  REAL(DP) :: at_plumed(3,3)
  REAL(DP) :: virial(3,3)
  REAL(DP) :: volume
  REAL(DP), ALLOCATABLE :: tau_plumed(:,:)
  !
  IF(use_plumed) then
    IF(ionode)THEN
      at_plumed=alat*at;  ! the cell, rescaled properly
      allocate(tau_plumed(3,nat))
      tau_plumed=alat*tau
      volume=+at_plumed(1,1)*at_plumed(2,2)*at_plumed(3,3) &
             +at_plumed(1,2)*at_plumed(2,3)*at_plumed(3,1) &
             +at_plumed(1,3)*at_plumed(2,1)*at_plumed(3,2) &
             -at_plumed(1,1)*at_plumed(3,2)*at_plumed(2,3) &
             -at_plumed(1,2)*at_plumed(3,3)*at_plumed(2,1) &
             -at_plumed(1,3)*at_plumed(3,1)*at_plumed(2,2) 
      virial=-sigma*volume

      CALL plumed_f_gcmd("setStep"//char(0),istep)
      CALL plumed_f_gcmd("setMasses"//char(0),amass)
      CALL plumed_f_gcmd("setForces"//char(0),force)
      CALL plumed_f_gcmd("setPositions"//char(0),tau_plumed)
      CALL plumed_f_gcmd("setBox"//char(0),at_plumed)
      CALL plumed_f_gcmd("setVirial"//char(0),virial)
      CALL plumed_f_gcmd("setEnergy"//char(0),etot)
      CALL plumed_f_gcmd("calc"//char(0),0)

      sigma=-virial/volume

      deallocate(tau_plumed)
    ENDIF
    CALL mp_bcast(force, ionode_id, intra_image_comm)
    CALL mp_bcast(sigma, ionode_id, intra_image_comm)
  ENDIF
  !
  !
END SUBROUTINE plugin_forces
