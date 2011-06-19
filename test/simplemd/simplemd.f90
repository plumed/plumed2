
program simplemd
use routines
implicit none

integer              :: natoms          ! number of atoms
real,    allocatable :: positions(:,:)  ! atomic positions
real,    allocatable :: velocities(:,:) ! velocities
                                        ! was calculated last time
real,    allocatable :: masses(:)       ! masses
real,    allocatable :: forces(:,:)     ! forces   
real                 :: cell(3)         ! cell size
real                 :: cell9(9)         ! cell size

! neighbour list variables
! see Allen and Tildesey book for details
integer              :: listsize        ! size of the list array
integer, allocatable :: list(:)         ! neighbour list
integer, allocatable :: point(:)        ! pointer to neighbour list
real,    allocatable :: positions0(:,:) ! reference atomic positions, i.e. positions when the neighbour list

! input parameters
! all of them have a reasonable default value, set in read_input()
real           :: tstep          ! simulation timestep
real           :: temperature    ! temperature
real           :: friction       ! friction for Langevin dynamics (for NVE, use 0)
real           :: listcutoff     ! cutoff for neighbour list
real           :: forcecutoff    ! cutoff for forces
integer        :: nstep          ! number of steps
integer        :: nconfig        ! stride for output of configurations
integer        :: nstat          ! stride for output of statistics
integer        :: maxneighbour   ! maximum average number of neighbours per atom
integer        :: idum           ! seed
logical        :: wrapatoms      ! if true, atomic coordinates are written wrapped in minimal cell
character(256) :: inputfile      ! name of file with starting configuration (xyz)
character(256) :: outputfile     ! name of file with final configuration (xyz)
character(256) :: trajfile       ! name of the trajectory file (xyz)
character(256) :: statfile       ! name of the file with statistics
character(256) :: parfile        ! name of the file with parameters (optional) 
character(256) :: string         ! a string for parsing 


real    :: engkin         ! kinetic energy
real    :: engconf        ! configurational energy
real    :: engint         ! integral for conserved energy in Langevin dynamics

logical :: recompute_list ! control if the neighbour list have to be recomputed
integer :: istep          ! step counter
integer :: iatom

integer :: i  ! an integer for loops 
integer :: argcount  ! a counter for the arguments
logical :: has_parfile  ! a flag for the parameter file 

logical :: plumed
integer :: plumedavailable

CALL plumed_installed(plumedavailable)

plumed=.false.
if(plumedavailable>0)plumed=.true.

IF(plumed) THEN
  CALL plumed_g_create()
END IF

argcount = IARGC()

has_parfile=.false.
do i=1,argcount
        call getarg(i,string)
        if (INDEX(string,'-in').NE.0)then
          call getarg(i+1,string)
          read(string,*) parfile 
          has_parfile=.true.
        endif
enddo

call read_input(temperature,tstep,friction,forcecutoff, &
                listcutoff,nstep,nconfig,nstat, &
                wrapatoms, &
                inputfile,outputfile,trajfile,statfile, &
                maxneighbour,idum,has_parfile,parfile)

! number of atoms is read from file inputfile
call read_natoms(inputfile,natoms)

! write the parameters in output so they can be checked
write(*,*) "Starting configuration           : ",trim(inputfile)
write(*,*) "Final configuration              : ",trim(outputfile)
write(*,*) "Number of atoms                  : ",natoms
write(*,*) "Temperature                      : ",temperature
write(*,*) "Time step                        : ",tstep
write(*,*) "Friction                         : ",friction
write(*,*) "Cutoff for forces                : ",forcecutoff
write(*,*) "Cutoff for neighbour list        : ",listcutoff
write(*,*) "Number of steps                  : ",nstep
write(*,*) "Stride for trajectory            : ",nconfig
write(*,*) "Trajectory file                  : ",trim(trajfile)
write(*,*) "Stride for statistics            : ",nstat
write(*,*) "Statistics file                  : ",trim(statfile)
write(*,*) "Max average number of neighbours : ",maxneighbour
write(*,*) "Seed                             : ",idum
write(*,*) "Are atoms wrapped on output?     : ",wrapatoms

! Since each atom pair is counted once, the total number of pairs
! will be half of the number of neighbours times the number of atoms
listsize=maxneighbour*natoms/2

! allocation of dynamical arrays
allocate(positions(3,natoms))
allocate(positions0(3,natoms))
allocate(velocities(3,natoms))
allocate(forces(3,natoms))
allocate(masses(natoms))
allocate(point(natoms))
allocate(list(listsize))


! masses are hard-coded to 1
masses=1.0

! energy integral initialized to 0
engint=0.0

! positions are read from file inputfile
call read_positions(inputfile,natoms,positions,cell)

! velocities are randomized according to temperature
call randomize_velocities(natoms,temperature,masses,velocities,idum)

!CALL init_metadyn(natoms,tstep,masses,masses,1,1.0D0,"plumed.dat"//char(0));
IF(plumed) THEN
  CALL plumed_g_cmd("setNatoms"//char(0),natoms)
  CALL plumed_g_cmd("setMDEngine"//char(0),"simplemd"//char(0))
  CALL plumed_g_cmd("setTimestep"//char(0),tstep)
  CALL plumed_g_cmd("setPlumedDat"//char(0),"plumed.dat"//char(0))
  CALL plumed_g_cmd("init"//char(0))
ENDIF


! neighbour list are computed, and reference positions are saved
call compute_list(natoms,listsize,positions,cell,listcutoff,point,list)
write(*,*) "List size: ",point(natoms)-1
positions0=positions

! forces are computed before starting md
call compute_forces(natoms,listsize,positions,cell,forcecutoff,point,list,forces,engconf)

! here is the main md loop
! Langevin thermostat is applied before and after a velocity-Verlet integrator
! the overall structure is:
!   thermostat
!   update velocities
!   update positions
!   (eventually recompute neighbour list)
!   compute forces
!   update velocities
!   thermostat
!   (eventually dump output informations)
do istep=1,nstep
  call thermostat(natoms,masses,0.5*tstep,friction,temperature,velocities,engint,idum)

  do iatom=1,natoms
    velocities(:,iatom)=velocities(:,iatom)+forces(:,iatom)*0.5*tstep/masses(iatom)
  end do

  do iatom=1,natoms
    positions(:,iatom)=positions(:,iatom)+velocities(:,iatom)*tstep
  end do

! a check is performed to decide whether to recalculate the neighbour list
  call check_list(natoms,positions,positions0,listcutoff,forcecutoff,recompute_list)
  if(recompute_list) then
    call compute_list(natoms,listsize,positions,cell,listcutoff,point,list)
    positions0=positions
    write(*,*) "Neighbour list recomputed at step ",istep
    write(*,*) "List size: ",point(natoms)-1
  end if

  call compute_forces(natoms,listsize,positions,cell,forcecutoff,point,list,forces,engconf)
!  forces=0.0
IF(plumed) THEN
  CALL plumed_g_cmd("setMasses"//char(0),masses)
  CALL plumed_g_cmd("setForces"//char(0),forces)
  CALL plumed_g_cmd("setStep"//char(0),istep)
  CALL plumed_g_cmd("setPositions"//char(0),positions)
  CALL plumed_g_cmd("calc"//char(0),0)
ENDIF

  do iatom=1,natoms
    velocities(:,iatom)=velocities(:,iatom)+forces(:,iatom)*0.5*tstep/masses(iatom)
  end do

  call thermostat(natoms,masses,0.5*tstep,friction,temperature,velocities,engint,idum)

! kinetic energy is calculated
  call compute_engkin(natoms,masses,velocities,engkin)

! eventually, write positions and statistics
  if(modulo(istep,nconfig)==0) call write_positions(trajfile,natoms,positions,cell,wrapatoms)
  if(modulo(istep,nstat)==0)   call write_statistics(statfile,istep,tstep,natoms,engkin,engconf,engint)

end do

call write_final_positions(outputfile,natoms,positions,cell,wrapatoms)
write(*,*) "Execution completed"

! deallocation of all allocatable array
deallocate(positions)
deallocate(velocities)
deallocate(forces)
deallocate(masses)
deallocate(positions0)
deallocate(point)
deallocate(list)

IF(plumed) THEN
  CALL plumed_g_finalize()
END IF

end program simplemd
