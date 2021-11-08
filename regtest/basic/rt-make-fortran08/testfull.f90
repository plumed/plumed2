subroutine testfull()
use plumed_f08_module
use iso_c_binding
implicit none

type(plumed)      :: pl
type(c_ptr)       :: ptr
type(plumed_error) :: error
integer :: s,i
integer, target :: natoms
real(c_double), target, allocatable :: positions(:,:)
real(c_double), target, allocatable :: masses(:)
real(c_double), target, allocatable :: forces(:,:)
real(c_double), target :: box(3,3)
real(c_double), target :: virial(3,3)

! Required due to INTEL bug
real(c_double), pointer :: positions_ptr(:,:)
real(c_double), pointer :: masses_ptr(:)
real(c_double), pointer :: forces_ptr(:,:)
real(c_double), pointer :: box_ptr(:,:)
real(c_double), pointer :: virial_ptr(:,:)

integer(c_int), target :: ene,stopflag
real(c_double) :: bias
open(10,file="error_codes_full")

call pl%cmd_val("initx","a",error=error)
write(10,*)error%code
call pl%cmd_val("initx","a",error=error)
write(10,*)error%code
call pl%cmd_val("setNatoms",3.0,error=error)
write(10,*)error%code
call pl%finalize()

call plumed_create(pl)

natoms=10
ALLOCATE(positions(3,10))
ALLOCATE(forces(3,10))
ALLOCATE(masses(10))

! Required due to INTEL bug
positions_ptr => positions
masses_ptr => masses
forces_ptr => forces
box_ptr => box
virial_ptr => virial

positions=0.0
box=0.0
virial=0.0
forces=0.0
masses=1.0
do i=1,natoms
  positions(1,i)=1.0+10*i
  positions(2,i)=2.0+10*i
  positions(3,i)=3.0+10*i
enddo

stopflag=0
ene=1
bias=-10
call pl%cmd_const_ptr("setNatoms",c_loc(natoms)) ! dummy test
call pl%cmd_val("setNatoms",1*natoms)
call pl%cmd("init")
call pl%cmd_ptr("setStopFlag",stopflag)
call pl%cmd_val("readInputLine","p: POSITION ATOM=2")
call pl%cmd_val("readInputLine","g: GYRATION ATOMS=@allatoms")
call pl%cmd_val("readInputLine","r: RESTRAINT ARG=g AT=0 KAPPA=3")
call pl%cmd_val("readInputLine","COMMITTOR ARG=p.x STRIDE=1 BASIN_LL1=0 BASIN_UL1=30")
call pl%cmd_val("readInputLine","PRINT ARG=p.*,r.* FILE=testme2")
call pl%cmd_val("setStep",1)
!! check 
! NB intel cannot even compile this so I have to remove it:
! call pl%cmd_const_ptr("setPositions",positions_ptr(:,::2),error=error)
! write(10,*) "code",error%code

call pl%cmd_const_ptr("setPositions",positions_ptr)
call pl%cmd_ptr("setForces",forces_ptr)
call pl%cmd_const_ptr("setMasses",masses_ptr)
call pl%cmd_const_ptr("setBox",box_ptr)
call pl%cmd_ptr("setVirial",virial_ptr)
write(10,"(A,I5)") "stopflag should be 0",stopflag
write(10,"(A,I5)") "isEnergyNeeded should be 1",ene
call pl%cmd_ref("isEnergyNeeded",ene)
write(10,"(A,I5)") "isEnergyNeeded should be 0",ene
write(10,"(A,I5)") "stopflag should be 0",stopflag
call pl%cmd("calc")
write(10,"(A,I5)") "stopflag should be 1",stopflag
call pl%cmd_ref("getBias",bias)
write(10,"(A,F10.4)") "bias",bias

open(11,file="forces")
do i=1,natoms
  write(11,"(3F10.2)") forces(:,i)
enddo

close(10)
close(11)
end subroutine testfull
