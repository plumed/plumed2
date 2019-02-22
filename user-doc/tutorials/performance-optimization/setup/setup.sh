#! /bin/bash

git clone https://github.com/srnas/ff

sed "s/\*/\'/" GG.pdb > GGp.pdb

export GMXLIB=$PWD/ff/
{ echo 1 ; echo 1 ; } | gmx_mpi pdb2gmx -f GGp.pdb -o conf_start.gro

gmx_mpi editconf -f conf_start.gro -o conf_box.gro -bt dodecahedron -d 3.0
gmx_mpi solvate -cp conf_box.gro -cs -o conf_water.gro -p topol.top
cat > min2.mdp << EOF
title           = Minimization     ; Title of run
cpp             = /lib/cpp      ; Preprocessor
integrator      = steep         ; Algorithm (steep = steepest descent minimization)
nsteps          = 50000         ; Maximum number of (minimization) steps to perform
nstenergy       = 50            ; Write energies to disk every nstenergy steps
nstxtcout       = 50            ; Write coordinates to disk every nstxtcout steps
nstvout         = 50
nstlog          = 50
xtc_grps        = System        ; Which coordinate group(s) to write to disk
energygrps      = System        ; Which energy group(s) to write to disk
cutoff-scheme=Verlet
nstlist         = 5             ; Frequency to update the neighbor list and long range forces
ns_type         = grid          ; Method to determine neighbor list (simple, grid)
constraints     = none          ; Bond types to replace by constraints
coulombtype = PME 
rvdw = 1.0 
rlist = 1.0 
rcoulomb = 1.0 
pbc                      = xyz
EOF

gmx_mpi grompp -f min2.mdp -c conf_water.gro -p topol.top -o genion.tpr -maxwarn 1
echo SOL | gmx_mpi genion -s genion.tpr -o conf_ions1.gro -p topol.top -neutral -pname K

gmx_mpi grompp -f min2.mdp -c conf_ions1.gro -p topol.top -o genion.tpr -maxwarn 1
echo SOL | gmx_mpi genion -s genion.tpr -o conf_ions.gro -p topol.top -conc 0.1 -pname K -nname CL

awk '{
  line[NR]=$0
  nr=NR
}END{
  for(i=1;i<=nr;i++) {
    if(i==2) {
      print line[i] -1
    } else if(match(line[i],"K") && !done) {
      gsub("K ","MG",line[i])
      gsub(" K","MG",line[i])
      print(line[i])
      done=1
      i++
    } else print line[i]
  }
}' < conf_ions.gro > conf.gro


awk '{
  if($2=="molecules") inside=1
  if(inside && $1=="K") {
     print "MG",1
     print "K", $2-2
     inside=0
  } else print
}' < topol.top > topol_new.top
mv topol_new.top topol.top

cat > run-equil.mdp << EOF

title                    = MD simulation
integrator               = md
dt                       = 0.0002
nsteps                   = 25000               ;number of steps
nstxout                  =
nstvout                  =
nstlog                   = 200
nstenergy                = 200
nstxtcout                = 2500
xtc_grps                 = System     ;group(s) to write to xtc trajectory
energygrps               = System       ;group(s) to write to energy file
nstlist                  = 10                   ;Frequency to update the neighbor list (and the long-range forces,
                                                ;when using twin-range cut-off's).
ns_type                  = grid                 ;Make a grid in the box and only check atoms in neighboring grid cells
                                                ;when constructing a new neighbor list every nstlist steps.
cutoff-scheme=Verlet
coulombtype = PME
rvdw =  1.0
rlist =  1.0
rcoulomb =  1.0
pbc                      = xyz                  ; Periodic boudary conditions in all the directions
tcoupl                   = v-rescale            ;Temperature coupling
tc-grps                  = RNA Water_and_ions
tau_t                    = 0.1 0.1
ref_t                    = 300 300
Pcoupl                   = Parrinello-Rahman    ;Pressure coupling
Pcoupltype               = isotropic
tau_p                    = 2.0
compressibility          = 4.5e-5
ref_p                    = 1.0
refcoord_scaling         = com
;Constrain all bonds
constraints              = h-bonds

EOF

gmx_mpi grompp -f run-equil.mdp -c conf.gro -o topol.tpr
gmx_mpi mdrun -v -c conf.gro

cat > run.mdp << EOF

title                    = MD simulation
integrator               = md
dt                       = 0.002
nsteps                   = 25000               ;number of steps
nstxout                  =
nstvout                  =
nstlog                   = 200
nstenergy                = 200
nstxtcout                = 2500
xtc_grps                 = System     ;group(s) to write to xtc trajectory
energygrps               = System       ;group(s) to write to energy file
nstlist                  = 10                   ;Frequency to update the neighbor list (and the long-range forces,
                                                ;when using twin-range cut-off's).
ns_type                  = grid                 ;Make a grid in the box and only check atoms in neighboring grid cells
                                                ;when constructing a new neighbor list every nstlist steps.
cutoff-scheme=Verlet
coulombtype = PME
rvdw =  1.0
rlist =  1.0
rcoulomb =  1.0
pbc                      = xyz                  ; Periodic boudary conditions in all the directions
tcoupl                   = v-rescale            ;Temperature coupling
tc-grps                  = RNA Water_and_ions
tau_t                    = 0.1 0.1
ref_t                    = 300 300
Pcoupl                   = Parrinello-Rahman    ;Pressure coupling
Pcoupltype               = isotropic
tau_p                    = 2.0
compressibility          = 4.5e-5
ref_p                    = 1.0
refcoord_scaling         = com
;Constrain all bonds
constraints              = h-bonds

EOF
gmx_mpi grompp -f run.mdp -c conf.gro -o topol.tpr
echo 0 | gmx_mpi trjconv -f conf.gro -o conf.pdb
