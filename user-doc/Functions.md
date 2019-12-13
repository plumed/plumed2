\page Function Functions

When performing biased dynamics or analyzing a trajectory you may wish to analyze/bias the value of
some function of a set of collective variables rather than the values of the collective variables
directly.  You can do this with PLUMED by using any one of the following list of functions.

Notice that in many functions you should explicitly say to PLUMED whether the result
is a periodic variable or not using the keyword `PERIODIC`.
This is crucial to allow a variable to be properly based.
To know if a function is periodic
of not you should answer to the following question:

- Can my function change with a discontinuity when I move my atoms in a continuous manner?

In case the answer is no, than you should use `PERIODIC=NO`. In case the answer is yes, then you should
consider the following question:

- Are the values of the function at the discontinuity always the same or do they change?

In case the answer is that they are the same, you should use `PERIODIC=A,B` where `A`
is the smallest value and `B` is the largest value. In case the answer is that the
values at the discontinuity are not always the same, then you cannot construct a variable that
can be biased with PLUMED. Consider the following examples:

\plumedfile
t: TORSION ATOMS=1,2,3,4
# When atoms are moved, t could jump suddenly from -pi to +pi

c: MATHEVAL ARG=t FUNC=x*x*x PERIODIC=-31.0062766802998,31.0062766802998
# When atoms are moved, c could jump suddenly from -pi**3 to +pi**3

# equivalently, we could have used:
# c: COMBINE ARG=t POWERS=3 PERIODIC=-31.0062766802998,31.0062766802998

# compute x/y/z components of the distance between atoms 1 and 10
d: DISTANCE ATOMS=1,10 COMPONENTS

# make a new variable equal to d.z but with the correct periodicity
dz: COMBINE ARG=d.z PERIODIC=-10,10
# here we assumed the system is in a orthorhombic box with z side = 20
\endplumedfile

@FUNCTION@

