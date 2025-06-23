# GPU parallelism in PLUMED

For certain actions in PLUMED you can use the USEGPU flag. This flag turns on an experimental GPU parallized version of the 
command. GPU parallelism in PLUMED has been implemented using [openACC](https://www.openacc.org) and is currently experimental. We are actively working
on these features at the moment. __There is thus no guarantee that the GPU accelerated versions of actions are any faster than 
the CPU versions.__ If you have experimented with these features on your own calculations we would love to hear from you (even 
if your experience was negative.)

## List of actions that can be called with the USEGPU option:

 - module:
   - ACTION

 - colvar:
   - [ANGLE](ANGLE.md)
   - [DIPOLE](DIPOLE.md)
   - [DISTANCE](DISTANCE.md)
   - [PLANE](PLANE.md)
   - [POSITION](POSITION.md)
   - [TORSION](TORSION.md)
 - crystdistrib:
   - ~~[QUATERNION_BOND_PRODUCT_MATRIX](QUATERNION_BOND_PRODUCT_MATRIX.md)~~ setup, but deactivated
 - secondarystructure:
   - [SECONDARY_STRUCTURE_DRMSD](SECONDARY_STRUCTURE_DRMSD.md), and in particular:
     - [ALPHARMSD](ALPHARMSD.md) only with **TYPE=DRMSD**
     - [ANTIBETARMSD](ANTIBETARMSD.md) only with **TYPE=DRMSD**
     - [PARABETARMSD](PARABETARMSD.md) only with **TYPE=DRMSD**
 - volumes:
   - [AROUND](AROUND.md)
   - [INCYLINDER](INCYLINDER.md)
   - [INSPHERE](INSPHERE.md)
