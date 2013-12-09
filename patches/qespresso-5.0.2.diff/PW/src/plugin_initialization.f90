!
! Copyright (C) 2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE plugin_initialization()
  !----------------------------------------------------------------------------
  !
  USE io_global,        ONLY : stdout, ionode
  USE kinds,            ONLY : DP
  USE io_files,         ONLY : tmp_dir
  !
  USE plugin_flags
  !
  USE ions_base,        ONLY : amass, ityp, nat
  !
  USE dynamics_module,  ONLY : dt
  USE constants,        ONLY : au_ps
  !
  !
  IMPLICIT NONE
  !
  INTEGER  :: na
  INTEGER  :: plumedavailable
  REAL*8   :: energyUnits,lengthUnits,timeUnits 
  !
  IF(use_plumed) then

     CALL plumed_f_installed(plumedavailable)
   
     IF(plumedavailable<=0)THEN
        write(stdout,*)"YOU ARE LOOKING FOR PLUMED BUT LOOKS LIKE IT IS NOT AVAILABLE: DO YOU HAVE IT IN YOUR LD_LIBRARY_PATH?" 
        STOP 
     ELSE 
        IF (ionode) THEN

            write(stdout,*)"  CREATING PLUMED FROM THE PROGRAM" 
            call plumed_f_gcreate()
            CALL plumed_f_gcmd("setRealPrecision"//char(0),8)
            energyUnits=1312.75  ! Ry to kjoule mol 
            lengthUnits=0.0529177249 ! bohr to nm
            timeUnits=2*au_ps ! internal time to ps
            call plumed_f_gcmd("setMDEnergyUnits"//char(0),energyUnits)
            call plumed_f_gcmd("setMDLengthUnits"//char(0),lengthUnits)
            call plumed_f_gcmd("setMDTimeUnits"//char(0),timeUnits)
            call plumed_f_gcmd("setPlumedDat"//char(0),"plumed.dat"//char(0))
            call plumed_f_gcmd("setLogFile"//char(0),"PLUMED.OUT"//char(0))
            call plumed_f_gcmd("setNatoms"//char(0),nat)
            call plumed_f_gcmd("setMDEngine"//char(0),"qespresso");
            call plumed_f_gcmd("setTimestep"//char(0),dt);
            call plumed_f_gcmd("init"//char(0),0);


        ENDIF
     ENDIF
  ENDIF
  !
  !
END SUBROUTINE plugin_initialization
