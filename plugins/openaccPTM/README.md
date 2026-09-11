This plugin allows you to use the USEGPU flag with certain actions in PLUMED. This flag turns on an experimental GPU parallized version of the 
command. GPU parallelism in PLUMED has been implemented using [openACC](https://www.openacc.org) and is currently experimental. We are actively working
on these features at the moment. __There is thus no guarantee that the GPU accelerated versions of actions are any faster than 
the CPU versions.__ If you have experimented with these features on your own calculations we would love to hear from you (even
if your experience was negative.)

## [Experimental] Compiling the openaccPTM plugin

The openaccPTM plugin is compiled separately to the main PLUMED code. To compile the plugin you need to have the 
[NVIDIA HPC SDK](https://developer.nvidia.com/hpc-sdk) and the plumed library avaiable in your path.  You then compile 
the plugin by issuing the following commands:

```bash
cd plugin/openaccPTM
./configure 
make clean
NVCXX=mpic++ make
```  

It is worth running these commands more than once if the compilation doesn't work the first time. If they complete correctly a dynamic libarary called `plumedOpenACC.so` 
will be created in the plumed2/plugin directory.

We have tested this module with the [24.3](https://developer.nvidia.com/nvidia-hpc-sdk-243-downloads) version of the NVIDIA HPC SDK compiler and found that it works.
It currently does not work with NVHPC 26.5 compiler.

Note that we recommend running the test suite after compilation.  You can do this by issuing the following commands:

```bash
cd plugin/openaccPTM/regtest
make
```

## [Experimental] Using PLUMED with the openaccPTM plugin

To use openACC in PLUMED you need to load the dynamic libary that you compiled in the previous section.  You can do this as follows:

```plumed
LOAD FILE=/path/to/plumedOpenACC.so
```

where `/path/to` is the directory that contains the compiled plugins. If you choose not to move the library after compilation this will 
be the `plumed2/plugin` directory. To use the plugin to calculate an action using the GPU you add the `USEGPU` keyword to the list of keywords
for that action.  For example, to calculate a contact matrix using the GPU you would use the following input:

```plumed
LOAD FILE=/path/to/plumedOpenACC.so
CONTACT_MATRIX GROUP=1-1000 SWITCH={EXP D_0=0.2 R_0=0.1 D_MAX=0.66} USEGPU
``` 

The next section contains a list of the actions that can take the USEGPU command in input and that can thus be run on the GPU using this plugin.

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
- adjmat
   - [CONTACT_MATRIX](CONTACT_MATRIX.md)
- function
   - [LESS_THAN](LESS_THAN.md)
   - [MORE_THAN](MORE_THAN.md)
   - [BETWEEN](BETWEEN.md)
   - [COMBINE](COMBINE.md)
- matrixtools
   - [MATRIX_PRODUCT](MATRIX_PRODUCT.md)
   - [MATRIX_VECTOR_PRODUCT](MATRIX_VECTOR_PRODUCT.md)
- symfunc
   - [Q1](Q1.md)
   - [Q3](Q3.md)
   - [Q4](Q4.md)
   - [Q6](Q6.md)
   - [LOCAL_AVERAGE](LOCAL_AVERAGE.md)
   - [FCCUBIC](FCCUBIC.md)
   - [TETRAHEDRAL](TETRAHEDRAL.md)
   - [SIMPLECUBIC](SIMPLECUBIC.md)
   - [COORDINATION_SHELL_FUNCTION](COORDINATION_SHELL_FUNCTION.md)
   - [COORDINATION_SHELL_AVERAGE](COORDINATION_SHELL_AVERAGE.md)
   - [LOCAL_Q1](LOCAL_qQ1.md)
   - [LOCAL_Q3](LOCAL_qQ3.md)
   - [LOCAL_Q4](LOCAL_qQ4.md)
   - [LOCAL_Q6](LOCAL_qQ6.md)
