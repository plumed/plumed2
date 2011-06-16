module routines
implicit none
contains
subroutine read_input(temperature,tstep,friction, &
                      forcecutoff,listcutoff,nstep,&
                      nconfig,nstat, &
                      wrapatoms, &
                      inputfile,outputfile,trajfile,statfile, &
                      maxneighbours,idum)
  implicit none
  real,           intent(out) :: temperature
  real,           intent(out) :: tstep
  real,           intent(out) :: friction
  real,           intent(out) :: forcecutoff
  real,           intent(out) :: listcutoff
  integer,        intent(out) :: nstep
  integer,        intent(out) :: nconfig
  integer,        intent(out) :: nstat
  logical,        intent(out) :: wrapatoms
  integer,        intent(out) :: maxneighbours
  character(256), intent(out) :: inputfile
  character(256), intent(out) :: outputfile
  character(256), intent(out) :: trajfile
  character(256), intent(out) :: statfile
  integer, intent(out) :: idum
  integer :: iostat
  character(256) :: line,keyword,keyword1
  integer :: i
  logical :: foundsharp
! default values
  temperature=1.0
  tstep=0.005
  friction=0.0
  forcecutoff=2.5
  listcutoff=3.0
  nstep=1
  nconfig=10
  nstat=1
  maxneighbours=1000
  idum=0
  wrapatoms=.false.
  statfile=""
  trajfile=""
  outputfile=""
  inputfile=""
  
  do
    read(*,"(a)",iostat=iostat) line
! when the file finishes, exit the loop
    if(iostat/=0) exit
! delete everything past an eventual "#" comment
    foundsharp=.false.
    do i=1,len(line)
      if(line(i:i)=="#") foundsharp=.true.
      if(foundsharp) line(i:i)=" "
    end do
! if the remaining line is empty, skip it
    if(len_trim(line)==0) cycle
! read the first word from line
    read(line,*) keyword
! the second word is then read to the proper variable
    select case(keyword)
    case("temperature")
      read(line,*) keyword1,temperature
    case("tstep")
      read(line,*) keyword1,tstep
    case("friction")
      read(line,*) keyword1,friction
    case("forcecutoff")
      read(line,*) keyword1,forcecutoff
    case("listcutoff")
      read(line,*) keyword1,listcutoff
    case("nstep")
      read(line,*) keyword1,nstep
    case("nconfig")
      read(line,*) keyword1,nconfig,trajfile
    case("nstat")
      read(line,*) keyword1,nstat,statfile
    case("maxneighbours")
      read(line,*) keyword1,maxneighbours
    case("wrapatoms")
      read(line,*) keyword1,wrapatoms
    case("inputfile")
      read(line,*) keyword1,inputfile
    case("outputfile")
      read(line,*) keyword1,outputfile
    case("seed")
      read(line,*) keyword1,idum
      idum=-idum ! idum for ran1() needs to be negative
    case default
! an unknown word will stop the execution
      write(0,*) "Unknown keyword :",trim(keyword)
      stop
    end select
  end do
  if(inputfile=="") then
    write(0,*) "Specify input file"
    stop
  end if
  if(outputfile=="") then
    write(0,*) "Specify output file"
    stop
  end if
  if(trajfile=="") then
    write(0,*) "Specify traj file"
    stop
  end if
  if(statfile=="") then
    write(0,*) "Specify stat file"
    stop
  end if
end subroutine read_input

subroutine read_positions(inputfile,natoms,positions,cell)
! read positions and cell from a file called inputfile
! natoms (input variable) and number of atoms in the file should be consistent
  implicit none
  character(*), intent(in)  :: inputfile
  integer,      intent(in)  :: natoms
  real,         intent(out) :: positions(3,natoms)
  real,         intent(out) :: cell(3)
  integer :: iatom
  character(100) :: atomname
  open(10,file=inputfile)
  read(10,*)
  read(10,*) cell
  do iatom=1,natoms
    read(10,*) atomname,positions(:,iatom)
  end do
! note: atomname is read but not used
  close(10)
end subroutine read_positions

subroutine read_natoms(inputfile,natoms)
! read the number of atoms in file "input.xyz"
  implicit none
  character(*),     intent(in) :: inputfile
  integer,          intent(out):: natoms
  open(10,file=inputfile)
  read(10,*) natoms
  close(10)
end subroutine read_natoms

subroutine write_positions(trajfile,natoms,positions,cell,wrapatoms)
! write positions on file trajfile
! positions are appended at the end of the file
  implicit none
  character(*), intent(in) :: trajfile
  integer,      intent(in) :: natoms
  real,         intent(in) :: positions(3,natoms)
  real,         intent(in) :: cell(3)
  logical,      intent(in) :: wrapatoms
  integer :: iatom
  real :: pos(3)
  logical, save :: first=.true.
  if(first) then
    open(10,file=trajfile)
    first=.false.
  else
    open(10,file=trajfile,position="append")
  end if
  write(10,*) natoms
  write(10,*) cell
  do iatom=1,natoms
! usually, it is better not to apply pbc here, so that diffusion
! is more easily calculated from a trajectory file:
    if(wrapatoms) then
      call pbc(cell,positions(:,iatom),pos)
    else
      pos=positions(:,iatom)
    end if
    write(10,*) "Ar",pos
  end do
  close(10)
end subroutine write_positions

subroutine write_final_positions(outputfile,natoms,positions,cell,wrapatoms)
! write positions on file outputfile
  character(*), intent(in) :: outputfile
  integer,      intent(in) :: natoms
  real,         intent(in) :: positions(3,natoms)
  real,         intent(in) :: cell(3)
  logical,      intent(in) :: wrapatoms
  integer :: iatom
  real :: pos(3)
  open(10,file=outputfile)
  write(10,*) natoms
  write(10,*) cell
  do iatom=1,natoms
! usually, it is better not to apply pbc here, so that diffusion
! is more easily calculated from a trajectory file:
    if(wrapatoms) then
      call pbc(cell,positions(:,iatom),pos)
    else
      pos=positions(:,iatom)
    end if
    write(10,*) "Ar",pos
  end do
  close(10)
end subroutine write_final_positions

subroutine write_statistics(statfile,istep,tstep,natoms,engkin,engconf,engint)
! write statistics on file statfile
  character(*), intent(in) :: statfile
  integer,      intent(in) :: istep
  real,         intent(in) :: tstep
  integer,      intent(in) :: natoms
  real,         intent(in) :: engkin
  real,         intent(in) :: engconf
  real,         intent(in) :: engint
  logical, save :: first=.true.
  integer, save :: last_time_reopened=0
  if(first) then
! first time this routine is called, open the file
    open(666,file=statfile)
    first=.false.
  end if
  if(istep-last_time_reopened>100) then
! every 100 steps, reopen the file to flush the buffer
    close(666)
    open(666,file=statfile,position="append")
    last_time_reopened=istep
  end if
  write(666,"(i10,99g13.6)") istep,istep*tstep,2.0*engkin/(3.0*natoms),engconf,engkin+engconf,engkin+engconf+engint
end subroutine write_statistics

subroutine randomize_velocities(natoms,temperature,masses,velocities,idum)
! randomize the velocities according to the temperature
  implicit none
  integer, intent(in)    :: natoms
  real,    intent(in)    :: temperature
  real,    intent(in)    :: masses(natoms)
  real,    intent(out)   :: velocities(3,natoms)
  integer, intent(inout) :: idum
  real, external :: gasdev
  integer :: iatom,i
  do iatom=1,natoms
    do i=1,3
      velocities(i,iatom)=sqrt(temperature/masses(iatom))*gasdev(idum)
    end do
  end do
end subroutine randomize_velocities

subroutine compute_engkin(natoms,masses,velocities,engkin)
! calculate the kinetic energy from the velocities
  implicit none
  integer, intent(in) :: natoms
  real,    intent(in) :: masses(natoms)
  real,    intent(in) :: velocities(3,natoms)
  real,    intent(out):: engkin
  integer :: iatom
  engkin=0.0
  do iatom=1,natoms
    engkin=engkin+0.5*masses(iatom)*sum(velocities(:,iatom)**2)
  end do
end subroutine compute_engkin

subroutine pbc(cell,vin,vout)
! apply periodic boundary condition to a vector
  implicit none
  real, intent(in) :: cell(3)
  real, intent(in) :: vin(3)
  real, intent(out) :: vout(3)
  integer :: i
  do i=1,3
    vout(i)=vin(i)-nint(vin(i)/cell(i))*cell(i)
  end do
end subroutine pbc

subroutine check_list(natoms,positions,positions0,listcutoff,forcecutoff,recompute)
! check if the neighbour list have to be recomputed
  implicit none
  integer, intent(in)    :: natoms
  real,    intent(in)    :: positions(3,natoms)
  real,    intent(in)    :: positions0(3,natoms)
  real,    intent(in)    :: listcutoff
  real,    intent(in)    :: forcecutoff
  logical, intent(out)   :: recompute
  real :: displacement(3) ! displacement from positions0 to positions
  real :: delta2          ! square of the 'skin' thickness
  integer :: iatom
  recompute=.false.
  delta2=(0.5*(listcutoff-forcecutoff))**2
! if ANY atom moved more than half of the skin thickness, recompute is set to .true.
  do iatom=1,natoms
    displacement=positions(:,iatom)-positions0(:,iatom)
    if(sum(displacement**2)>delta2) recompute=.true.
  end do
end subroutine check_list

subroutine compute_list(natoms,listsize,positions,cell,listcutoff,point,list)
  implicit none
  integer, intent(in)  :: natoms
  integer, intent(in)  :: listsize
  real,    intent(in)  :: positions(3,natoms)
  real,    intent(in)  :: cell(3)
  real,    intent(in)  :: listcutoff
! see Allen-Tildesey for a definition of point and list
  integer, intent(out) :: point(natoms)
  integer, intent(out) :: list(listsize)
  integer :: iatom,jatom  ! indexes of the two involved atoms
  real :: distance(3)     ! distance of the two atoms
  real :: distance_pbc(3) ! minimum-image distance of the two atoms
  real    :: listcutoff2  ! squared list cutoff
  listcutoff2=listcutoff**2
  point(1)=1
  do iatom=1,natoms-1
    point(iatom+1)=point(iatom)
    do jatom=iatom+1,natoms
      distance=positions(:,iatom)-positions(:,jatom)
      call pbc(cell,distance,distance_pbc)
! if the interparticle distance is larger than the cutoff, skip
      if(sum(distance_pbc**2)>listcutoff2) cycle
      if(point(iatom+1)>listsize) then
! too many neighbours
        write(0,*) "Verlet list size exceeded"
        write(0,*) "Increase maxneighbours"
        stop
      end if
      list(point(iatom+1))=jatom
      point(iatom+1)=point(iatom+1)+1
    end do
  end do
end subroutine compute_list

subroutine compute_forces(natoms,listsize,positions,cell,forcecutoff,point,list,forces,engconf)
  implicit none
  integer, intent(in)  :: natoms
  integer, intent(in)  :: listsize
  real,    intent(in)  :: positions(3,natoms)
  real,    intent(in)  :: cell(3)
  real,    intent(in)  :: forcecutoff
  integer, intent(in)  :: point(natoms)
  integer, intent(in)  :: list(listsize)
  real,    intent(out) :: forces(3,natoms)
  real,    intent(out) :: engconf
  integer :: iatom,jatom  ! indexes of the two involved atoms
  integer :: jlist        ! counter for the neighbours of iatom
  real :: distance(3)     ! distance of the two atoms
  real :: distance_pbc(3) ! minimum-image distance of the two atoms
  real :: distance_pbc2   ! squared minimum-image distance
  real :: forcecutoff2    ! squared force cutoff
  real :: f(3)            ! force
  real :: engcorrection   ! energy necessary shift the potential avoiding discontinuities
  forcecutoff2=forcecutoff**2
  engconf=0.0
  forces=0.0
  engcorrection=4.0*(1.0/forcecutoff2**6-1.0/forcecutoff2**3)
  do iatom=1,natoms-1
    do jlist=point(iatom),point(iatom+1)-1
      jatom=list(jlist)
      distance=positions(:,iatom)-positions(:,jatom)
      call pbc(cell,distance,distance_pbc)
      distance_pbc2=sum(distance_pbc**2)
! if the interparticle distance is larger than the cutoff, skip
      if(distance_pbc2>forcecutoff2) cycle
      engconf=engconf+4.0*(1.0/distance_pbc2**6-1.0/distance_pbc2**3)-engcorrection
      f=2.0*distance_pbc*4.0*(6.0/distance_pbc2**7-3.0/distance_pbc2**4)
! same force on the two atoms, with opposite sign:
      forces(:,iatom)=forces(:,iatom)+f
      forces(:,jatom)=forces(:,jatom)-f
    end do
  end do
end subroutine compute_forces

subroutine thermostat(natoms,masses,dt,friction,temperature,velocities,engint,idum)
! Langevin thermostat, implemented as decribed in Bussi and Parrinello, Phys. Rev. E (2007)
! it is a linear combination of old velocities and new, randomly chosen, velocity,
! with proper coefficients
  implicit none
  integer, intent(in)    :: natoms
  real,    intent(in)    :: masses(natoms)
  real,    intent(in)    :: dt
  real,    intent(in)    :: friction
  real,    intent(in)    :: temperature
  real,    intent(inout) :: velocities(3,natoms)
  real,    intent(inout) :: engint ! contribution of phase-space compression
                                   ! its increment is equal to minus the kinetic-energy increment due to the
                                   ! thermostat
  integer, intent(inout) :: idum
  real :: c1 ! coefficient for the old velocity
  real :: c2 ! coefficient for the new velocity
  real, external :: gasdev
  integer :: i,iatom
  c1=exp(-friction*dt)
  do iatom=1,natoms
    c2=sqrt((1.0-c1**2)*temperature/masses(iatom))
    do i=1,3
      engint=engint+0.5*masses(iatom)*velocities(i,iatom)**2
      velocities(i,iatom)=c1*velocities(i,iatom)+c2*gasdev(idum)
      engint=engint-0.5*masses(iatom)*velocities(i,iatom)**2
    end do
  end do
end subroutine thermostat
end module routines

