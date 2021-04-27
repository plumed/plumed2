\page regtests Adding regressions tests

When you write new functionality for PLUMED it is important that you document these features AND that you provide suitable
regression tests so as to ensures that developers working on the code in the future do not break any existing features.
The regression tests that are currently within PLUMED are all in the regtests directory separated.  The tests are contained 
in a number of separate directories that (roughly) correspond to the modules of PLUMED.  These tests are run whenever the 
user executes the

\verbatim
make tests
\endverbatim

command.  Furthermore, they are also executed every time one of the core developers pushes commits to the main plumed2 fork
on github.  In fact these tests are executed across multiple architectures whenever a push is performed.  To do this we use 
a web-app called travis-ci.  The current status of the code on travis-ci can be viewed here: https://travis-ci.org/plumed/plumed2
<b> We take code testing various seriously and will not consider including contributions from developers unless appropriate 
regression tests are included. </b> 

\section addingregtests Creating a regression test

The first step in creating a regression test is to create a directory for it somewhere in the regtest directory of PLUMED.  If you
are implementing your own module for PLUMED then it is best to keep your regression tests in a subdirectory named after your module.
If you are not doing this then you can create a directory in one of the already available module directories.  Your regression test
module <b> must </b> be named:

\verbatim
rt-something
\endverbatim

where the something should be replaced with a suitable description of the feature that is being tested within the module.  The directory 
name must begin with rt-*, however, as this is how the various Makefiles than run the regression tests identify the directories that 
contain tests that must be run.

\subsection regconfig Creating the config file and the Makefile

Once you have created the directory for your regression test you will need to create two files: a Makefile and a config file.  These files
can be copied from one of the regtest directories that are already there.  The Makefile does not need to be changed as it essentially 
just includes a script from elsewhere in PLUMED.  You may have to edit the config file, however, so in order to understand what the various 
lines in this file do here is an example config file: 

\verbatim
mpiprocs=2
type=driver
plumed_modules=adjmat
arg="--plumed plumed.dat --trajectory-stride 50 --timestep 0.005 --ixyz diala_traj_nm.xyz --dump-forces forces --dump-forces-fmt=%10.6f"
extra_files="../../trajectories/diala_traj_nm.xyz  ../trajectories/path_msd/all.pdb"
\endverbatim

- The first line in this file - the one containing the string mpiprocs=2 - tells PLUMED that when this regtest should be run using MPI on two nodes.
If mpiprocs is not specified then PLUMED will, by default, run the test one node.

- The second line tells PLUMED that the test is to be performed using the command line tool <a href="../../user-doc/html/driver.html">driver</a>.  
There are four options that you can use here:
type=driver tells PLUMED to run the test by reading in a trajectory and analysing it using <a href="../../user-doc/html/driver.html">driver</a>, 
type=simplemd tells PLUMED that the test will involve
running some Lennard Jones MD using <a href="../../user-doc/html/simplemd.html">simplemd</a>, 
type=sum_hills tells PLUMED that the test will involve summing gaussians using the 
<a href="../../user-doc/html/sum_hills.html">sum_hills</a> untility
 and type=make tells PLUMED that a main function (called main.cpp) is included in the regtest directory and that the test should compile this function, link in
the PLUMED library and run the resulting executible. The vast majority of features are best tested using type=driver here although the type=make function
can be useful for testing tools that are used in many places in the code (see the code in regtest/basic/rt-make-4 for an example of how use type=make to test
PLUMED's internal matrix class).  If you need to test PLUMED in a way that is not included in these various types please contact the PLUMED developers before
making any changes.

- The third line - the one containing the string plumed_modules=adjmat - tells PLUMED that this regtest needs the optional module adjmat to be installed in 
order for it to be run.  Before trying to run this test PLUMED will thus check whether or not the module is installed.  If it is not installed then it will 
abandon the attempt without running the test.  <b> You must incorporate a line like this if you are testing some feature that is included in an optional module </b>

- The fourth line specifies the command line arguments that need to be passed to the PLUMED tool.  In the above example the final command the test will 
run is thus as follows:

\verbatim
plumed driver --plumed plumed.dat --trajectory-stride 50 --timestep 0.005 --ixyz diala_traj_nm.xyz --dump-forces forces --dump-forces-fmt=%10.6f
\endverbatim    

- The final line tells PLUMED that a few files need to be copied from elsewhere in the PLUMED directory hierarchy in order to run the test.  This copying of files
is useful as it ensures that our repository does not grow large because the same trajectory is duplicated in many regtest directory.  However, if all the files you 
need to run the test  are contained in the test directory this line is not necessary.  Having said that, however, it may be useful to reuse the trajectories contained 
in the directory regtest/trajectories in creating your own regression tests.  

\subsection reffiles Creating the reference files

The next step in creating your regression test will be to copy suitable trajectories to run the calculation into the test directory and to write a plumed input file.
Obviously, we cannot tell you what to do here (although there are some suggestions in the next section) as the tests you will need to write will depend on the feature
you are incorporating.  <b>In general though it should be possible it should be possible to run the test in a few seconds and the output files should be small.</b>
Once you have created this input you can run the test by simply executing 

\verbatim
make
\endverbatim

inside the directory you have created for the regression test.  When this process completes a new directory called tmp is created, which contains the output from your
calculation.  Regression tests work by comparing this output to a reference output and you should thus create this reference now.  This is simple to do you simply copy 
the output file named <i>name</i> from the tmp directory to the regtest directory and rename it <i>name</i>.reference.  As an example suppose we were testing the following
PLUMED input by analysing a trajectory using driver:

\verbatim
d1: DISTANCE ATOMS=1,2
PRINT ARG=d1 FMT=%8.4f FILE=colvar
\endverbatim 

This input creates an output file called colvar so in future we can check the code is running in the same way by comparing with a reference version of the colvar file.  To
create this reference file we execute the following command:

\verbatim
cp tmp/colvar colvar.reference
\endverbatim

Immediately after our first successful run of the test.

\subsection regtrepo Adding the regtest to PLUMED's git repository

Suppose you use the instructions outlined in the previous section to create a regression test in a directory called rt-mynewtest.  You can add this test to the git 
repository by running the commands:

\verbatim
git add rt-mynewtest
git commit
\endverbatim

in the directory containing rt-mynewtest. By adding the test in this way you will ensure that your feature is not broken in future versions of the code.

\section reg-suggestions-1 Some suggestions on how to test new collective variables

If the new feature you are implementing is a new collective variable or a new function of a collective variable and you inherited from either 
\ref PLMD::Colvar, \ref PLMD::function::Function or \ref PLMD::multicolvar::MultiColvarBase then there is a reasonably well established way of testing
that your implementing is correct.  In essence you have to ensure that the numerical derivatives and the analytical derivatives of the variable 
are equal.  This sort of test is easy to setup with PLUMED.  To explain how to setup the test let's suppose for a moment that you want to test 
\ref DISTANCE.  You would do this using an input file something like the one below:

\verbatim
d1: DISTANCE ATOMS=1,2
d1n: DISTANCE ATOMS=1,2 NUMERICAL_DERIVATIVES
PRINT ARG=d1 FILE=colvar FMT=%8.4f 
DUMPDERIVATIVES ARG=d1,d1n FILE=deriv FMT=%8.4f
\endverbatim 

The input above outputs the analytical and numerical derivatives of the distance between atom 1 and 2 to a file called deriv.  The analytical 
derivatives will be in the third column of this file while the numerical derivatives are in the fourth.  If you had run a regtest with this 
input you can thus ensure that the two sets of derivatives were equal using the command:

\verbatim
awk '{print $3-$4}' tmp/deriv
\endverbatim

Lets suppose the derivatives are equal. You might thus conclude the creation of your regtest by copying the deriv and colvar file to files called 
deriv.reference and colvar.reference and checking in the regtest directory.  We would urge you not to do the following, however, as we have found 
that numerical derivatives calculated on different architectures vary, which makes it difficult to determine if the regtests are failing for real 
or if it is just small numerical errors in the calculated numerical derivatives.  On top of this calculating numerical derivatives is expensive.
We would thus ask that once you have ensured that the numerical and analytical derivatives match that you <b> check in deriv files that contain the 
analytical derivatives only. </b>  Also notice the use of the FMT keywords in the above input.  This is important - there will be small numerical 
differences when the code is run by different architectures.  These small differences are not a problem, however, if the format to use 
for real numbers in the output is stated explicitly in the input using FMT.

Notice that this method of testing will also work if you are testing a \ref PLMD::Function.  As an example we can test \ref COMBINE
as follows:

\verbatim
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=3,4
c1: COMBINE ARG=d1,d2 PERIODIC=NO 
c1n: COMBINE ARG=d1,d2 PERIODIC=NO NUMERICAL_DERIVATIVES
PRINT ARG=c1 FILE=colvar FMT=%8.4f
DUMPDERIVATIVES ARG=c1,c1n FILE=deriv FMT=%8.4f 
\endverbatim

Now, however, the numerical and analytical derivatives in deriv are the derivatives of the value c1 with respect to the distances d1 and d2.  In
other words the file deriv does not contain derivatives with respect to atomic positions.

\section reg-suggestions-2 More complicated tests

If your feature is more complicated than those covered in the previous section - if for example you are using multiple lines in the PLUMED input
to calculate a final CV or if you have had to write a new apply method in your method - then you will need to do more to test your implementation.
It is probably still a good idea to test analytical derivatives against numerical derivatives in these complicated cases, however.  You can do this
by adding a bias to the system and by exploiting the --debug-forces feature that is contained in <a href="../../user-doc/html/driver.html">driver</a>.  
To understand how this works suppose you have a PLUMED input that looks like this:

\verbatim
d1: DISTANCE ATOMS=1,2
RESTRAINT ARG=d1 AT=1 KAPPA=10
\endverbatim

This input applies a harmonic restraint on the value of the distance between atoms 1 and 2.  If you use the above as input to a driver command such as the one
below:

\verbatim
plumed driver --ixyz trajectory.xyz --debug-forces forces.num
\endverbatim

Then a file called forces.num will be output that contains three columns. The first column is a numerical indexr. The second contains the forces on the atoms
in your system calculated using the analytical equations for the forces that are contained within PLUMED.  The final column also contains the forces on the atoms
in your system but this time these forces are calculated numerically using finite differences.  Calculating the numerical derivatives using --debug-forces is 
very expensive so you definitely should not include regtests that exploit this option in PLUMED.   Instead you can just output the analytical forces using a command
such as:

\verbatim
plumed driver --ixyz trajectory.xyz --dump-forces forces --dump-forces-fmt %8.4f 
\endverbatim

Notice, however, that with both these commands, even thought there are only forces on the 6 coordinates that are used to specify the positions of 1 and 2 
and the cell vectors, the forces with respect to all the atomic positions are output.  In other words, the output files produced usig --debug-forces and 
-dump-foces will contain a lot of zeros. 

