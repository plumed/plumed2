# GPU parallelism in PLUMED

For certain actions in PLUMED you can use the USEGPU flag. This flag turns on an experimental GPU parallized version of the 
command. GPU parallelism in PLUMED has been implemented using [openACC](https://www.openacc.org) and is currently experimental. We are actively working
on these features at the moment. __There is thus no guarantee that the GPU accelerated versions of actions are any faster than 
the CPU versions.__ If you have experimented with these features on your own calculations we would love to hear from you (even 
if your experience was negative.)

## [Experimental] Compiling plumed with openacc

_This section will likely be moved  in the proper installation page_

To compile PLUMED with openacc enabled you will need to have the [NVIDIA HPC SDK](https://developer.nvidia.com/hpc-sdk) avaiable on your path.

Plumed is tested with the [24.3](https://developer.nvidia.com/nvidia-hpc-sdk-243-downloads) version.

To prepare the compilation add `--enable-opeacc` to the `./configure` options.
It is also possible to pass some extra options by exporting or specifying the the following variables:
 - **PLUMED_ACC_TYPE**: if omitted defaults to `gpu`, can be changed to `host` or `multicore` to try a compilation that targets the CPU also for the openacc accellerate part of the code.
 - **PLUMED_ACC_GPU**: if `PLUMED_ACC_TYPE` is set to `gpu`, it can be used to specify a range of `-gpu` optios to pass to the nvhpc compiler (for example the target gpu, see the [compiler manual](https://docs.nvidia.com/hpc-sdk/compilers/hpc-compilers-user-guide/index.html) for the options) options can be comma separated or space separated


!!! warning
    (Currently) modules that use a custom openmp reduction can be compiled with nvhpc.
    Currently `membranefusion` is not compatible and should be excluded from the compilation


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
