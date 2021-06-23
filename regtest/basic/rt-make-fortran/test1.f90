SUBROUTINE TEST1(testdat,testme)
  USE PLUMED_MODULE
  IMPLICIT NONE
  CHARACTER(*) :: testdat,testme
  CHARACTER(len=32) :: p
  INTEGER :: i
  INTEGER :: natoms
  REAL(8), ALLOCATABLE :: positions(:,:)
  REAL(8), ALLOCATABLE :: forces(:,:)
  REAL(8), ALLOCATABLE :: masses(:)
  REAL(8) :: box(3,3)
  REAL(8) :: virial(3,3)
  REAL(8) :: testconvert
  INTEGER(4) :: stopflag,ene
  REAL(8) :: bias
  call plumed_f_create(p)
  call plumed_f_cmd(p,"CLTool setArgvLine"//char(0),"plumed --help"//char(0))
  call plumed_f_cmd(p,"CLTool run"//char(0),i)
  call plumed_f_finalize(p)

  open(unit=10,file=testdat,status='new')

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
    stopflag=0
    ene=1
    call plumed_f_gcreate()
    call plumed_f_gcmd("setNatoms"//char(0),8+2) ! pass a temporary
    call plumed_f_gcmd("init"//char(0),0)
    call plumed_f_gcmd("setStopFlag"//char(0),stopflag)
    call plumed_f_gcmd("readInputLine"//char(0),"p: POSITION ATOM=2"//char(0))
    call plumed_f_gcmd("readInputLine"//char(0),"g: GYRATION ATOMS=@allatoms"//char(0))
    call plumed_f_gcmd("readInputLine"//char(0),"r: RESTRAINT ARG=g AT=0 KAPPA=3"//char(0))
    call plumed_f_gcmd("readInputLine"//char(0),"COMMITTOR ARG=p.x STRIDE=1 BASIN_LL1=0 BASIN_UL1=30"//char(0))
    call plumed_f_gcmd("readInputLine"//char(0),"PRINT ARG=p.*,r.* FILE="//testme//char(0))
    call plumed_f_gcmd("setStep"//char(0),1)
    call plumed_f_gcmd("setPositions"//char(0),positions)
    call plumed_f_gcmd("setForces"//char(0),forces)
    call plumed_f_gcmd("setMasses"//char(0),masses)
    call plumed_f_gcmd("setBox"//char(0),box)
    call plumed_f_gcmd("setVirial"//char(0),virial)
    write(10,"(A,I5)") "stopflag should be 0",stopflag
    write(10,"(A,I5)") "isEnergyNeeded should be 1",ene
    call plumed_f_gcmd("isEnergyNeeded"//char(0),ene)
    write(10,"(A,I5)") "isEnergyNeeded should be 0",ene
    write(10,"(A,I5)") "stopflag should be 0",stopflag
    call plumed_f_gcmd("calc"//char(0),0)
    write(10,"(A,I5)") "stopflag should be 1",stopflag
    call plumed_f_gcmd("getBias"//char(0),bias)
    write(10,"(A,F10.4)") "bias",bias

  ! test convert
    testconvert=0.0
    call plumed_f_gcmd("convert 3*2-1"//char(0),testconvert)
    write(10,"(A,F10.4)") "convert should be 5",testconvert
    call plumed_f_gfinalize()
  
  ! object interface
    call plumed_f_create(p)
    call plumed_f_cmd(p,"setNatoms"//char(0),10)
    call plumed_f_cmd(p,"init"//char(0),0)
    call plumed_f_cmd(p,"readInputLine"//char(0),"p: POSITION ATOM=2"//char(0))
    call plumed_f_cmd(p,"readInputLine"//char(0),"g: GYRATION ATOMS=@allatoms"//char(0))
    call plumed_f_cmd(p,"readInputLine"//char(0),"r: RESTRAINT ARG=g AT=0 KAPPA=3"//char(0))
    call plumed_f_cmd(p,"readInputLine"//char(0),"PRINT ARG=p.*,r.* FILE="//testme//"2"//char(0))
    call plumed_f_cmd(p,"setStep"//char(0),1)
    call plumed_f_cmd(p,"setPositions"//char(0),positions)
    call plumed_f_cmd(p,"setForces"//char(0),forces)
    call plumed_f_cmd(p,"setMasses"//char(0),masses)
    call plumed_f_cmd(p,"setBox"//char(0),box)
    call plumed_f_cmd(p,"setVirial"//char(0),virial)
    call plumed_f_cmd(p,"calc"//char(0),0)
    call plumed_f_finalize(p)

    close(unit=10)

END SUBROUTINE
