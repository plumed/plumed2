PROGRAM main
  USE PLUMED_MODULE
  IMPLICIT NONE
  CHARACTER(len=32) :: p
  INTEGER :: i
  INTEGER :: natoms
  REAL(8), ALLOCATABLE :: positions(:,:)
  REAL(8), ALLOCATABLE :: forces(:,:)
  REAL(8), ALLOCATABLE :: masses(:)
  REAL(8) :: box(3,3)
  REAL(8) :: virial(3,3)
  call plumed_f_create(p)
  call plumed_f_cmd(p,"CLTool setArgvLine"//char(0),"plumed --help"//char(0))
  call plumed_f_cmd(p,"CLTool run"//char(0),i)
  call plumed_f_finalize(p)

  ALLOCATE(positions(3,10))
  ALLOCATE(forces(3,10))
  ALLOCATE(masses(10))
  positions=0.0
  box=0.0
  virial=0.0
  forces=0.0
  do i=1,10
    positions(1,i)=1.0+10*i
    positions(2,i)=2.0+10*i
    positions(3,i)=3.0+10*i
  enddo

  ! global interface
    call plumed_f_gcreate()
    call plumed_f_gcmd("setNatoms"//char(0),10)
    call plumed_f_gcmd("init"//char(0),0)
    call plumed_f_gcmd("readInputLine"//char(0),"p: POSITION ATOM=2"//char(0))
    call plumed_f_gcmd("readInputLine"//char(0),"PRINT ARG=p.*"//char(0))
    call plumed_f_gcmd("setStep"//char(0),1)
    call plumed_f_gcmd("setPositions"//char(0),positions)
    call plumed_f_gcmd("setForces"//char(0),forces)
    call plumed_f_gcmd("setMasses"//char(0),masses)
    call plumed_f_gcmd("setBox"//char(0),box)
    call plumed_f_gcmd("setVirial"//char(0),virial)
    call plumed_f_gcmd("calc"//char(0),0)
    call plumed_f_gfinalize()
  
  ! object interface
    call plumed_f_create(p)
    call plumed_f_cmd(p,"setNatoms"//char(0),10)
    call plumed_f_cmd(p,"init"//char(0),0)
    call plumed_f_cmd(p,"readInputLine"//char(0),"p: POSITION ATOM=2"//char(0))
    call plumed_f_cmd(p,"readInputLine"//char(0),"PRINT ARG=p.*"//char(0))
    call plumed_f_cmd(p,"setStep"//char(0),1)
    call plumed_f_cmd(p,"setPositions"//char(0),positions)
    call plumed_f_cmd(p,"setForces"//char(0),forces)
    call plumed_f_cmd(p,"setMasses"//char(0),masses)
    call plumed_f_cmd(p,"setBox"//char(0),box)
    call plumed_f_cmd(p,"setVirial"//char(0),virial)
    call plumed_f_cmd(p,"calc"//char(0),0)
    call plumed_f_finalize(p)

END PROGRAM main
