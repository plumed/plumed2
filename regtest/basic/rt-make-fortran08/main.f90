PROGRAM main
  USE PLUMED_MODULE_F08
  IMPLICIT NONE
  TYPE(PLUMED) :: p
  CALL TEST3()
  call p%cmd("init")
  call p%finalize()
END PROGRAM main
