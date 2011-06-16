PROGRAM MAIN
  IMPLICIT NONE
  INTEGER :: check,natoms
  natoms=5
  CALL plumed_installed(check)
  IF(check==1) THEN
    CALL plumed_g_create()
    CALL plumed_g_cmd("setMDEngine"//char(0),"AFortranCode"//char(0))
    CALL plumed_g_cmd("setNatoms"//char(0),natoms)
    CALL plumed_g_cmd("setPlumedDat"//char(0),"plumed.dat"//char(0))
    CALL plumed_g_cmd("init"//char(0),0)
    CALL plumed_g_finalize()
  ENDIF
END PROGRAM MAIN
