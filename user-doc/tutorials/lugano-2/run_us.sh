for AT in -3.00  -2.75 -2.50 -2.25 -2.00 -1.75 -1.50 -1.25 -1.00 -0.75 -0.50 -0.25 +0.00 +0.25 +0.50 +0.75 +1.00 +1.25 +1.50 +1.75 +2.00 +2.25 +2.50 +2.75 3.00
do

cat >plumed.dat << EOF
# vim:ft=plumed
MOLINFO STRUCTURE=diala.pdb
phi: TORSION ATOMS=@phi-2
psi: TORSION ATOMS=@psi-2
#
# Impose an umbrella potential on CV 1 and CV 2
# with a spring constant of 500 kjoule/mol
# at fixed points on the Ramachandran plot
#
restraint-phi: RESTRAINT ARG=phi KAPPA=250.0 AT=$AT
# monitor the two variables and the bias potential from the two restraints
PRINT STRIDE=10 ARG=phi,psi,restraint-phi.bias FILE=COLVAR$AT
EOF

gmx mdrun -plumed plumed.dat -nsteps 100000 -x traj$AT.xtc -nb cpu

done
