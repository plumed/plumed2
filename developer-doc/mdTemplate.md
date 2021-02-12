\page HowToPlumedYourMD How to add plumed to an MD code

\brief Learn how to use plumed in a not yet supported MD code

Plumed ships with scripts that can be used to add it to many of the standard MD packages.  Obviously though, if no patch is provided for the MD code you use then you will have to write one yourself.  Plumed has been designed so that it can be added to an MD code
either statically (i.e. as a collection of objects) or as a dynamic library.
For technical reasons it
is NOT possible to pack all the plumed objects as a static library and link the library
(if you really want to know it: when static librearies are linked the linker
typically discards all the objects which are not used; plumed is using an
automatic registration process which is not compatible with this choice).
Furthermore, unlike the previous version of plumed, the plumed source code is not compiled at the same time as your MD code is compiled.  Instead plumed now has its own makefile and is compiled separately to the MD engines that rely on it.  This makes plumed patches considerably simpler as now they only do two things:

- Modify the makefile so that the plumed is linked to the MD code
- Modify the source code to add all the required calls to plumed

\section makefile Modifying your makefile

Once the plumed source code has been compiled a file src/Plumed.inc is
generated. This file should be included in your MD code's makefile as it
informs the code where to find all the plumed source. There are three possible
ways to link PLUMED: static (linking the .o files directly), shared (linking
the libplumed.so library) or runtime (which links only a wrapper to plumed,
whereas the location of the remaining part on plumed, i.e. the libplumedKernel.so,
can be specified at runtime). The relevant variables are
- \$(PLUMED_STATIC_LOAD) the options to be added to the linker for static
  linking, i.e. all the plumed objects with their full path,
   plus special linking options, plus al the libraries used by plumed (e.g.
  matheval)
- \$(PLUMED_SHARED_LOAD) the options to be added to the linker for shared
  linking, i.e. full plumed shared library (libplumed.so) with its full path
  plus special linking options.
- \$(PLUMED_RUNTIME_LOAD) the options to be added to the linker for runtime
  linking, i.e. the wrapper object with its full path plus special linking options.

The libplumedKernel.so basically contains the same objects as the
libplumed.so except for the wrapper.

The simplest approach is to take advantage of the three variants
src/Plumed.inc.static , src/Plumed.inc.shared and src/Plumed.inc.runtime .
Including one of those, it is sufficient to add to the linker line the macro
\$(PLUMED_LOAD).

The best way to patch your MD file is:
- go into the root directory of the MD code and type
\verbatim
> plumed patch --new name-of-the-code
\endverbatim
  this will create a new (empty) patch file for your MD code
- patch the code with the empty patch
\verbatim
> plumed patch --patch --engine name-of-the-code
\endverbatim
  this will produce a symbolic link Plumed.inc in the root directory of the MD
  code
- make a backup copy of the Makefile of your MD code, naming it
  Makefile.preplumed
\verbatim 
> cp src/Makefile src/Makefile.preplumed
\endverbatim
- edit the Makefile including the Plumed.inc file
  (e.g. "include ../Plumed.inc") and adding "\$(PLUMED_LOAD)" to the linking
command
- save your modification into the plumed patch
\verbatim
> plumed patch --save
\endverbatim

Now you can patch/unpatch your MD code, automatically linking PLUMED, as if it
was an officially released interface (even though at this stage you do not
have any functionality of plumed yet in your code - see below how to add it).

Note that after you have this done, PLUMED creates a file in
$PLUMED_ROOT/patches/ that is generally named 
name-of-the-code.diff ( if your code is named name-of-the-code in the "save"
phase). You can also manually create, in the same directory,  a file  
called name-of-the-code.config so that you can use to perform five actions (totally optional).
The first one can be used to do some heuristic check to verify that the patch
is for the correct MD engine and that the MD source tree is in the correct state
(this typically means that it hass been already configured).
Then there are two functions that can
do some modifications to some
files that might vary considerable (as the Makefiles produced by the config
procedure that change with the machine) and do some actions after the
patching procedure is performed (a sanity check or a makedependence script to
be rerun).
The last two ones should revert the effect of the previous ones in such
a way that it is possible to smoothly revert the patch. Try to keep
them consistent so as the revert feature works correctly.
This file is in bash and typically looks like this 
\verbatim
function plumed_preliminary_test(){
# Just check of the code is the right one: no error means it is ok 
  grep -q  name-of-the-code Somefile 
}

function plumed_before_patch(){
# typically do something on the Makefile
}

function plumed_after_patch(){
# typically call some makedepends if necessary
}

function plumed_after_revert(){
# revert exacty the effect of plumed_before_patch
# typically undoes what was done on the Makefile
# additionally, one could add a makedepends
}

function plumed_before_revert(){
# revert exactly the effect of plumed_after_patch
# however, a makedepends script should always go in plumed_after_revert
}

\endverbatim

\attention Be careful with these scripts. You have access to a lot of flexibility
(adding files, modifying existing ones arbitrarily etc) but you should always 
take care of the revertibility with "patch -r".




\section language Language dependence

The next step is to add proper calls to plumed into the MD source code.

The interface between plumed and the MD codes is written is plain C so it is compatible with codes that are written in plain C, C++ and fortran.
To use this interface one should first create a plumed object (plumed_create()), then sent to it
messages (using plumed_cmd()) and deallocate it (plumed_finalize()) at the end
of the simulation. All of these routines have C++ and FORTRAN equivalents.
Notice that in C the plumed object is represented using a struct (plumed), in
C++ it is represented using a class (PLMD::Plumed) and in FORTRAN using a
string of 32 characters (CHARACTER(LEN=32)).
Also notice that whenever passing strings from FORTRAN you should add to them
"char(0)" as
a string terminator to be sure that the string is properly interpreted
(e.g. "pass-this-string"//char(0)).

As people might get confused by
the idea of "plumed being an object", we also added the possibility of using
a "generic global plumed instance", which can be accessed
using routines with a _g_ in the name (plumed_g_create(),
plumed_g_cmd() and plumed_g_finalize()).

This interface is very general, will not change in future releases, and is
fully contained in Plumed.h and Plumed.c files. For a reference to it, see
\ref ReferencePlumedH


\section mdImplementation A practical example

We now describe how to pass data to plumed in from an MD code that is written in C/C++.
First save all the files that you intend to modify into a .preplumed file (for
example if you want to modify myfile.c save it first in its original version
into myfile.c.preplumed): this will tell the patching procedure that this  is
a file that will enter the set of the files to be modified. 

In C or C++ files containing calls to plumed you will have to include the Plumed.h file.
In addition, you will also have to define a plumed obect that is visible in all these routines
and you will probably like to define some sort of plumedswitch,
which can be read in from the input to your code and used to tell the code whether or not this is to be a run with plumed.
Finally, you might like to include something in input so that you specify the name of the plumed input file.
How these things are best done will depend on your code and so we leave it to your discretion.  

Plumed must perform three actions inside your MD code:

- It must be initialized before the main MD loop starts so that the plumed input files are read in.
- It must be called every MD step so that the forces from the bias can be computed.
- It must be finalized when the simulation is completed

The only routine which can be used to send message to plumed is plumed_cmd
(or, equivalently, Plumed::cmd in C++ and plumed_f_cmd in FORTRAN).

Notice that as we PLUMED evolves new commands could be available.
Thus, if you want your interface to be compatible with multiple PLUMED versions, you should first
check that the PLUMED version is new enough to accept a given command.
To know which is the
currently linked PLUMED version see \ref apiversion .

The various calls that can be used during initialization are as follows:

\verbatim
plumed plumedmain; plumedmain=plumed_create();                 // Create the plumed object

// Calls to pass data to plumed
plumed_cmd(plumedmain,"setRealPrecision",&real_precision);     // Pass a pointer to an integer containing the size of a real number (4 or 8)
plumed_cmd(plumedmain,"setMDEnergyUnits",&energyUnits);        // Pass a pointer to the conversion factor between the energy unit used in your code and kJ mol-1
plumed_cmd(plumedmain,"setMDLengthUnits",&lengthUnits);        // Pass a pointer to the conversion factor between the length unit used in your code and nm 
plumed_cmd(plumedmain,"setMDTimeUnits",&timeUnits);            // Pass a pointer to the conversion factor between the time unit used in your code and ps

// This is valid only if API VERSION > 3
plumed_cmd(plumedmain,"setMDChargeUnits",&chargeUnits);        // Pass a pointer to the conversion factor between the charge unit used in your code and e

// This is valid only if API VERSION > 3
plumed_cmd(plumedmain,"setMDMassUnits",&massUnits);            // Pass a pointer to the conversion factor between the mass unit used in your code and amu

plumed_cmd(plumedmain,"setPlumedDat",&plumedInput);            // Pass the name of the plumed input file from the md code to plumed
plumed_cmd(plumedmain,"setMPIComm",&MPI_COMM_WORLD);           // Pass a pointer to the MPI communicator to plumed
// notice that from fortran the command "setMPIFComm" should be used instead
plumed_cmd(plumedmain,"setNatoms",&natoms);                    // Pass a pointer to the number of atoms in the system to plumed
plumed_cmd(plumedmain,"setMDEngine","gromacs");                // Pass the name of your md engine to plumed (now it is just a label) 
plumed_cmd(plumedmain,"setLog",fplog);                         // Pass the file on which to write out the plumed log (if the file is already open)
plumed_cmd(plumedmain,"setLogFile",fplog);		       // Pass the file  on which to write out the plumed log (to be created)
plumed_cmd(plumedmain,"setTimestep",&delta_t);                 // Pass a pointer to the molecular dynamics timestep to plumed

// This is valid only if API VERSION > 1
plumed_cmd(plumedmain,"setKbT",&kbT);                          // Pointer to a real containing the value of kbT

// This is valid only if API VERSION > 2
plumed_cmd(plumedmain,"setRestart",&res);                      // Pointer to an integer saying if we are restarting (zero means no, one means yes)

// This is valid only if API VERSION > 3
plumed_cmd(plumedmain,"readInputLine","d: DISTANCE ATOMS=1,2");// Read a single input line directly from a string

// This is valid only if API VERSION > 7
plumed_cmd(plumedmain,"readInputLines","d: DISTANCE ATOMS=1,2\n"
                                       "PRINT ARG=d");         // Read a multiple lines directly from a string. Allows comments and continuation lines.

// Calls to do the actual initialization (all the above commands must appear before this call)
plumed_cmd(plumedmain,"init",NULL);                            // Do all the initialization of plumed
plumed_cmd(plumedmain,"read",read);                            // Read the plumed input.  N.B. This is called during init and so this call is only required in special cases. 
\endverbatim

Please note that if your code is in FORTRAN you should append a "char(0)"
token to every string. Also, remember that FORTRAN is by default passing
arguments by reference, so that the "&" symbols which are required in C are
not necessary in FORTRAN.


The various calls that can be used pass data and calculate the forces due to the bias are as follows:

\verbatim
// Calls to pass data to plumed
plumed_cmd(plumedmain,"setStep",&step);                      // Pass a pointer to the current timestep to plumed
/ *** The way that you pass positions will depend on how they are stored in your code.  If the x, y and z position are all stored in a single array you may use:
plumed_cmd(plumedmain,"setPositions",&pos[0][0]);            // Pass a pointer to the first element in the atomic positions array to plumed  
                                                             // assuming they
                                                             // are stored in
                                                             // a
                                                             // x1,y1,z1,x2,y2,z2 ...
                                                             // kind of ordering
/ *** Othersize if you pass the three separate vectors of x, y and z positions using:
plumed_cmd(plumedmain,"setPositionX",&x[0]);                 // Pass a pointer to the first element in the array of x component of the atomic positions to plumed
plumed_cmd(plumedmain,"setPositionY",&y[0]);                 // Pass a pointer to the first element in the array of y component of the atomic positions to plumed
plumed_cmd(plumedmain,"setPositionZ",&z[0]);                 // Pass a pointer to the first element in the array of z component of the atomic positions to plumed
plumed_cmd(plumedmain,"setMasses",&mass[0]);                 // Pass a pointer to the first element in the masses array to plumed
plumed_cmd(plumedmain,"setCharges",&charge[0]);              // Pass a pointer to the first element in the charges array to plumed
plumed_cmd(plumedmain,"setBox",&box[0][0]);                  // Pass a pointer to the first element in the box share array to plumed
plumed_cmd(plumedmain,"setEnergy",&poteng);                  // Pass a pointer to the current value of the potential energy to plumed?
/ *** The way that you pass forces will depend on how they are stored in your code.  If the x, y and z force are all stored in a single array you may use:
plumed_cmd(plumedmain,"setForces",&f[0][0]);                 // Pass a pointer to the first element in the foces array to plumed
/ *** Othersize if you pass the three separate vectors of x, y and z forces using:
plumed_cmd(plumedmain,"setForcesX",&fx[0]);                  // Pass a pointer to the first element in the array of the x components of the atomic forces to plumed
plumed_cmd(plumedmain,"setForcesY",&fy[0]);                  // Pass a pointer to the first element in the array of the y components of the atomic forces to plumed
plumed_cmd(plumedmain,"setForcesZ",&fz[0]);                  // Pass a pointer to the first element in the array of the z components of the atomic forces to plumed
plumed_cmd(plumedmain,"setVirial",&force_vir[0][0]);         // Pass a pointer to the first element in the virial array to plumed

// Calls to do actual calculations
plumed_cmd(plumedmain,"calc",NULL);                          // Calculate and apply forces from the biases defined in the plumed input

// One can break up the "calc" command in two parts:
plumed_cmd(plumedmain,"prepareCalc",NULL);                   // Prepare to do a calculation by requesting all the atomic positions from the MD code
plumed_cmd(plumedmain,"performCalc",NULL);                   // Use the atomic positions collected during prepareCalc phase to calculate colvars and biases.
// The "performCalc" step can be further split into:
// This is valid only if API VERSION > 3
plumed_cmd(plumedmain,"performCalcNoUpdate",NULL);           // Same as "performCalc", skipping the update phase. Could be called multiple time per step
// This is valid only if API VERSION > 3
plumed_cmd(plumedmain,"update",NULL);                        // Only performs the update phase. Should be called once per step

// After the first part it will be possible to ask PLUMED e.g. if the energy is required with
plumed_cmd(plumedmain,"isEnergyNeeded,&flag);                // assuming flag is an int, that will be set to 0 if energy is not needed and 1 if it is needed

// The "prepareCalc" can be further split into:
plumed_cmd(plumedmain,"prepareDependencies",NULL);           // Work out what we are calculating during this MD step (this is the first step of prepareCalc)
plumed_cmd(plumedmain,"shareData",NULL);                     // Request all the atomic positions from the MD code (this is the second step of prepareCalc)
// This will allow to overlap sharing of atoms across multiple processors (sent in shareData) and calculation

// Some extra calls that might come in handy
plumed_cmd(plumedmain,"createFullList",&n);                  // Create a list containing of all the atoms plumed is using to do calculations (return the number of atoms in n)
plumed_cmd(plumedmain,"getFullList",&list);                  // Return a list (in list) containing all the indices plumed is using to do calculations
plumed_cmd(plumedmain,"clearFullList",NULL);                 // Clear the list of all the atoms that plumed is using to do calculations
plumed_cmd(plumedmain,"clear",clear);                        // Clear and delete all the pointers inside plumed.
\endverbatim

The plumed calls for the finalization tasks is as follows:

\verbatim
plumed_finalize(plumedmain);          // Call the plumed destructor
\endverbatim

\section mpicodes Dealing with parallelism

Plumed has functionality to deal with parallel MD codes.  The particular form of the functionality used to do this (and the frequency with which you will have to call these routines) will depend on whether your code is parallelism using a domain decomposition or particle decomposition strategy.  The calls for required for using this functionality are as follows:

\verbatim
plumed_cmd(plumedmain,"setAtomsNlocal",&nlocal);            // Pass a pointer to the number of atoms on this node
plumed_cmd(plumedmain,"setAtomsGatindex",gatindex);         // Pass an array (from a code in c) containing the indices of all the atoms on this node (used for domain decomposition)
plumed_cmd(plumedmain,"setAtomsFGatindex",gatindex);        // Pass an array (from a code in fortran) containing the indices of all the atoms on this node (used for domain decomposition)
plumed_cmd(plumedmain,"setAtomsContiguous",&start);         // Number the atoms on this node from start to start+nlocal   (used for particle decomposition)
\endverbatim

\section apiversion Inquiring for the plumed version

New functionalities might be added in the future to plumed. The description of
all the possible "commands" sent to plumed is called its API (application
programming interface). If you want to know which API is supported by plumed
you can use this call:
\verbatim
plumed_cmd(plumedmain,"getApiVersion",&api);                 // Pass the api version that plumed is using
\endverbatim 
With the current version, this will set the api variable (an integer) to 2. As
we add new features, this number will be increased.

\section Saving the diffs

This is similar to plumed 1. All the files that you want to modify should be
first copied to .preplumed files. Then use "plumed patch -s" to save the diff.

\section debugging-patch Debugging your PLUMED+MD code

If you want to be sure that your MD code is going to work properly with PLUMED we strongly suggest to go through the checks below.

\subsection debugging-patch-positions Check positions

The first thing to do is to checkout if positions are properly passed to PLUMED.
Let's assume your system has 100 atoms.
Just run PLUMED with an input such as the following one
\verbatim
DUMPATOMS ATOMS=1-100 FILE=test.xyz
\endverbatim
File test.xyz should look like:
\verbatim
100
 1000.000000 1000.000000 1000.000000
X 0.2 0.9 1.2
X 0.7 2.9 3.6
....... other 98 lines ......
100
 1000.000000 1000.000000 1000.000000
X 0.3 0.8 1.1
X 0.8 2.8 3.5
....... and so on ......
\endverbatim
That is, foreach snapsnot: first line, number of atoms; second line, box; then for each atom name/x/y/z coordinates.
Notice that coordinates are expected to be in in nanometers (PLUMED units).
For orthorhombic boxes the box line will contain three numbers, namely box size in x, y, and z directions.
For non orthorombic boxes, the box line will contain nine numbers, namely (ax,ay,az,bx,by,bz,cx,cy,cz),
where ax is the x component of the first lattice vector.

The produced test.xyz should match the trajectory file produced by your MD code.
Check all the atoms (not just the first ones).
Also, check that results are consistent also when you run in parallel (with particle or domain decomposition).

\subsection debugging-timestep Check timestep

PLUMED relies on the timestep from the MD code for several calculations. To be sure it has been passed properly,
run PLUMED with
\verbatim
d: DISTANCE ATOMS=1,10
PRINT ARG=d FILE=colvar
\endverbatim
Check the first column of `colvar` file (time). It should increase by one timestep at every line. Notice
that time here should appear in picoseconds (PLUMED units).

\subsection debugging-pass-energy Check energy

(This only applies if your MD codes actually passes energy to PLUMED)

Just try the following
\verbatim
e: ENERGY
PRINT ARG=e FILE=colvar
\endverbatim
Check the second column of `colvar` file (energy). It should be equivalent to the energy computed
in the MD code. Notice that PLUMED will write energy in kj/mol.

\subsection debugging-mass-charges Check masses and charges

The best way to debug the masses and charges is to print their values from PLUMED and
compare it with the MD code. This will be possible with PLUMED 2.2.

Meanwhile, you can try to compute center of masses or dipole moments and compare with what you expect.

\subsection debugging-patch-forces Check passed forces

Then you should double check whether PLUMED is able to return proper forces to the MD engine.
A first check could be the following. Run a short MD with the following PLUMED input
\verbatim
d: DISTANCE ATOMS=1,10
PRINT ARG=d FILE=colvar
\endverbatim
The run again with
\verbatim
d: DISTANCE ATOMS=1,10
RESTRAINT ARG=d AT=0.0 SLOPE=10
PRINT ARG=d FILE=colvar-biased
\endverbatim
Now plot both files with
\verbatim
gnuplot> plot "colvar" u 1:2 , "colvar-bias" u 1:2
\endverbatim

The two lines should start at the same value, but the second one should be systematically below the first one.

Notice that this test is not quantitative. Since an error in the force units would just give a qualitative
change, it is better to do a more rigorous check.

The best way to do a quantitative check is to alternatively add a restraint with the MD code and
with PLUMED and checkout that the obtained results are equivalent.

\todo better explain here

If your code passes these tests, you can likely start to do biased MD simulations in the NVT ensemble.
If you need NPT ensemble or if you want to bias the total energy you should continue with further tests.

\subsection debugging-patch-virial Check virial contribution

Most of the people use plumed to bias a small number of coordinates. This makes the contribution of plumed forces to the total virial
very small and makes it particularly difficult to find errors in this part. In case you combine plumed with a new MD code we highly
suggest to use the following procedure to debug the virial.

First run a short simulation (1000 steps should be enough) at constant pressure with pressure=1bar and the following plumed input:
\verbatim
v: VOLUME
PRINT ARG=v FILE=volume
\endverbatim

Then run another short simulation starting from identical conditions (i.e. same positions, velocities, random seed, etc) at constant pressure
with pressure=1001bar and the following plumed input:
\verbatim
v: VOLUME 
# slope should be just 10 times the Avogadro constant:
RESTRAINT AT=0.0 ARG=v SLOPE=-60.2214129
PRINT ARG=v FILE=volume2
\endverbatim
In this way, the negative pressure provided by plumed should exactly compensate for the extra 1000bar set in the barostat.
Thus, the two files `volume` and `volume2` should be virtually identical. Notice that small differences
due to different precision in the storage of avogadro number and other issues related to difficulties in exactly reproducing
a simulation could make the two trajectory slightly different.

If you pass this test, you can safely run biased MD in the NPT ensemble. Otherwise, there could be some
issue in the way the virial is passed to the MD code.

\subsection debugging-patch-energy Check forces on energy

(This only applies if your MD codes actually passes energy to PLUMED)

First run a short simulation (1000 steps should be enough) at constant temperature
\f$T=300\f$K and a given integration timestep \f$\Delta t\f$.
and the following PLUMED input:
\verbatim
e: ENERGY
PRINT ARG=e FILE=energy1
\endverbatim

Then run another short simulation starting from identical conditions (i.e. same positions, velocities, random seed, etc) at constant temperature
\f$T'=\alpha^2 T\f$K where  \f$\alpha=1.1\f$ (that is \f$T=363\f$) and a shorter timestep \f$\Delta t'=\frac{\Delta t}{\alpha}\f$.
Notice that all absolute times provided in the MD input (e.g. relaxation time for the thermostat) should be consistently divided by \f$\alpha\f$.
Use the following PLUMED input:
\verbatim
e: ENERGY
# slope is such that 
PRINT ARG=e FILE=energy2
# slope should be (alpha-1)=0.1
RESTRAINT AT=0.0 ARG=e SLOPE=0.1
\endverbatim

The two files `energy1` and `energy2` should be virtually identical.

In case you were able to have the virial properly working (see previous section), then you can try the same with a constant temperarure-constant pressure
simulation. In this case, please also monitor the volume of the resulting trajectory.
