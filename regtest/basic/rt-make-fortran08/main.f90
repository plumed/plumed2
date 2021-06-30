PROGRAM main
  USE PLUMED_MODULE_F08
  IMPLICIT NONE
  TYPE(PLUMED) :: p
  CALL TEST3()
  CALL TEST4()
  call p%cmd("init")
  call plumed_finalize(p)
END PROGRAM main
