PROGRAM MAIN
  IMPLICIT NONE
  INTEGER :: check,natoms
  natoms=5
  CALL plumed_f_installed(check)
  IF(check==1) THEN
    CALL plumed_f_gcreate()
    CALL plumed_f_gcmd("setMDEngine"//char(0),"AFortranCode"//char(0))
    CALL plumed_f_gcmd("setNatoms"//char(0),natoms)
    CALL plumed_f_gcmd("setPlumedDat"//char(0),"plumed.dat"//char(0))
    CALL plumed_f_gcmd("init"//char(0),0)
    CALL plumed_f_gfinalize()
  ENDIF
END PROGRAM MAIN
