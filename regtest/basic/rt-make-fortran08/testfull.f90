subroutine testfull()
use plumed_module_f08
use iso_c_binding
implicit none

type(plumed)      :: pl
type(c_ptr)       :: ptr
type(plumed_error) :: error
integer :: s,i
integer, target :: natoms
real(c_double), allocatable :: positions(:,:)
real(c_double), allocatable :: masses(:)
real(c_double), allocatable :: forces(:,:)
real(c_double) :: box(3,3)
real(c_double) :: virial(3,3)
integer(c_int) :: ene,stopflag
real(c_double) :: bias
open(10,file="error_codes_full")

call pl%cmd("initx","a",error=error)
write(10,*)error%code
call pl%cmd("initx","a",error=error)
write(10,*)error%code
call pl%cmd("setNatoms",3.0,error=error)
write(10,*)error%code
call pl%finalize()

call plumed_create(pl)

natoms=10
ALLOCATE(positions(3,10))
ALLOCATE(forces(3,10))
ALLOCATE(masses(10))
positions=0.0
box=0.0
virial=0.0
forces=0.0
do i=1,natoms
  positions(1,i)=1.0+10*i
  positions(2,i)=2.0+10*i
  positions(3,i)=3.0+10*i
enddo

stopflag=0
ene=1
bias=-10
call pl%cmd("setNatoms",c_loc(natoms)) ! pass a c pointer
call pl%cmd("init",0)
call pl%cmd("setStopFlag",stopflag)
call pl%cmd("readInputLine","p: POSITION ATOM=2")
call pl%cmd("readInputLine","g: GYRATION ATOMS=@allatoms")
call pl%cmd("readInputLine","r: RESTRAINT ARG=g AT=0 KAPPA=3")
call pl%cmd("readInputLine","COMMITTOR ARG=p.x STRIDE=1 BASIN_LL1=0 BASIN_UL1=30")
call pl%cmd("readInputLine","PRINT ARG=p.*,r.* FILE=testme2")
call pl%cmd("setStep",1)
call pl%cmd("setPositions",positions)
call pl%cmd("setForces",forces)
call pl%cmd("setMasses",masses)
call pl%cmd("setBox",box)
call pl%cmd("setVirial",virial)
write(10,"(A,I5)") "stopflag should be 0",stopflag
write(10,"(A,I5)") "isEnergyNeeded should be 1",ene
call pl%cmd("isEnergyNeeded",ene)
write(10,"(A,I5)") "isEnergyNeeded should be 0",ene
write(10,"(A,I5)") "stopflag should be 0",stopflag
call pl%cmd("calc",0)
write(10,"(A,I5)") "stopflag should be 1",stopflag
call pl%cmd("getBias",bias)
write(10,"(A,F10.4)") "bias",bias

open(11,file="forces")
do i=1,natoms
  write(11,"(3F10.2)") forces(:,i)
enddo

close(10)
close(11)
end subroutine testfull
