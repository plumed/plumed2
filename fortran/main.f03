PROGRAM main
  USE PLUMED_MODULE
  IMPLICIT NONE
  TYPE(PLUMED) :: p
  TYPE(PLUMED_ERROR) :: error
  INTEGER :: i
  INTEGER :: natoms
  REAL(8), ALLOCATABLE :: positions(:,:)
  REAL(8), ALLOCATABLE :: forces(:,:)
  REAL(8), ALLOCATABLE :: masses(:)
  REAL(8) :: box(3,3)
  REAL(8) :: virial(3,3)
  p=plumed_f03_create()
  call plumed_f03_cmd(p,"CLTool setArgvLine","plumed --help")
  call plumed_f03_cmd(p,"CLTool run",i)
  call plumed_f03_finalize(p)

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
  call plumed_f03_gcreate()
  call plumed_f03_gcmd("setNatoms",10)
  call plumed_f03_gcmd("init")
  call plumed_f03_gcmd("readInputLine","p: POSITION ATOM=2")
  call plumed_f03_gcmd("readInputLine","PRINT ARG=p.*")
  call plumed_f03_gcmd("setStep",1)
  call plumed_f03_gcmd("setPositions",positions)
  call plumed_f03_gcmd("setForces",forces)
  call plumed_f03_gcmd("setMasses",masses)
  call plumed_f03_gcmd("setBox",box)
  call plumed_f03_gcmd("setVirial",virial)
  call plumed_f03_gcmd("calc")
  call plumed_f03_gfinalize()

! object interface
  p=plumed_f03_create()
  call plumed_f03_cmd(p,"setNatoms",10)
  call plumed_f03_cmd(p,"init")
  call plumed_f03_cmd(p,"readInputLine","p: POSITION ATOM=2")
  call plumed_f03_cmd(p,"readInputLine","PRINT ARG=p.*")
  call plumed_f03_cmd(p,"setStep",1.0,error=error)
  if(error%code /= 0) then
    print * , "ERROR:",error%code,error%what
    call plumed_f03_cmd(p,"setStep",1,error=error)
    if(error%code /= 0) then
      stop "should never arrive here"
    endif
  endif
  call plumed_f03_cmd(p,"setPositions",positions)
  call plumed_f03_cmd(p,"setForces",forces)
  call plumed_f03_cmd(p,"setMasses",masses)
  call plumed_f03_cmd(p,"setBox",box)
  call plumed_f03_cmd(p,"setVirial",virial)
  call plumed_f03_cmd(p,"calc")
  call plumed_f03_finalize(p)

END PROGRAM main

