PROGRAM main
  USE PLUMED_F08_MODULE
  IMPLICIT NONE
  TYPE(PLUMED) :: p
  if(.not. plumed_installed()) then
    stop "plumed not installed"
  endif
  CALL TEST3()
  CALL TEST4()
  CALL TEST5()
  CALL TEST6()
  call p%cmd("init")
  call p%finalize()
  call TESTFULL()
END PROGRAM main
