\page AddingAMultiColvar How to add a new MultiColvar

As you are no doubt aware within plumed 2 you can calculate multiple 
instances of a collective coorinate from a single line in the input file.
One can then calculate functions such as the minimum, number less than,...
from the resulting distribution of collective variables.  To create these kinds
of collective variables we use the functionality implemented in PLMD::multicolvar::MultiColvar.
In fact in writing a single PLMD::multicolvar::MultiColvar you actually write many CVs at once
as the minimum of a distribution, the number less than and so on come with no 
additional effort on your part.

To better understand how to go about writing a new MultiColvar examine the interior
of one of the existing MultiColvars, e.g. PLMD::multicolvar::Distances or PLMD::multicolvar::CoordinationNumbers.
In fact a good way to start is to copy one of these files as certain features (e.g. the fact that you
have to include MultiColvar.h and ActionRegister.h and that all your code should be inside the namespace
PLMD) never change. 

The essential idea with multicolvar is that you can reuse the same functionality within many collective variables.
In particular, you can exploit the same parallelisation strategy in a wide variety of different collective variables.
Similarly one can use the same implementation of link cells in a wide variety of different cvs.  Lastly, one can
create functions of multicolvars and thus use the calculated cvs in a wide variety of different contexts.  If you are 
reading this you are probably aware of this and are (hopefully) impressed by the flexibility of this objects. 

The reason for this flexibility is that many CVs can be calculated using a variant on the following formula:

\f[
s = \sum_i g[f(\{X\}_i)]
\f]

where \f$g\f$ and \f$f\f$ are functions and the sum runs over a sets of atomic positions \f$\{X\}_i\f$.
In the above formula the sum and a set of likely \f$g\f$ functions are implemented within the base classes
base classes of multicolvar PLMD::multicolvar::MultiColvarBase, PLMD::multicolvar::MultiColvar, PLMD::multicolvar::MultiColvarFunction and 
PLMD::multicolvar::BridgedMultiColvarFunction.  When you implement a new multicolvar what you thus need to do is 
write code to calculate the function \f$f\f$ from a set of atomic positions.  One way to understand how to do this would be to look
through these classes and attempt to understand the code contained within them.  The following is an alternative.  It
provides a rather prescriptive description as to how to implement a new multicolvar.  I hope this helps but I am keenly
aware that any description is bound to be incomplete so if you have questions please do email the user list to ask.

\section AddingAMultColvarDocs Creating documentation  

The first thing you will have to change is the documentation.  As discussed on the \ref usingDoxygen page
of this manual the documentation is created using Doxygen.  You are implementing a cv so your PLMEDOC line
should read:

\verbatim
//+PLUMEDOC MCOLVAR MYCVKEYWORD 
\endverbatim 

Your documentation should contain a description of what your CV calculates and some examples.  You do not
need to write a description of the input syntax as that will be generated automatically.  For more information
on how to write documentation go to \ref usingDoxygen.

\section class Creating the class

The first step in writing the executable code for your CV is to create a class that calculates your CV.  The
declaration for your class should appear after the documentation and will look something like:

\verbatim
class MyNewMultiColvar : public MultiColvar {
private:
   // Declare all the variables you need here  
public:
  static void registerKeywords( Keywords& keys );
  MyNewMultiColvar(const ActionOptions&);
  virtual double compute( const unsigned& tindex, AtomValuePack& myatoms );
  bool isPeriodic(){ return false; }
};

PLUMED_REGISTER_ACTION(MyNewMultiColvar,"MYCVKEYWORD")
\endverbatim

This new class (MyNewMultiColvar) inherits from MultiColvar and so contains much of the functionality we 
require already.  Furthermore, by calling PLUMED_REGISTER_ACTION we have ensured that whenever the keyword
MYCVKEYWORD is found in the input file an object of type MyNewMultiColvar will be generated and hence that your 
new CV will be calculated wherever it is required.

The four functions that are defined in the above class (registerKeywords, the constructor, isPeriodic and compute) are mandatory.  
Without these functions the code will not compile and run.  Writing your new CV is simply a matter of writing these
four subroutines.

\section register RegisterKeywords

RegisterKeywords is the routine that is used by plumed to create the remainder of the documentation.  As much of this
documentation is created inside the MultiColvar class itself the first line of your new registerKeywords routine must
read:

\verbatim
MultiColvar::registerKeywords( keys )
\endverbatim

as this creates all the documentation for the features that are part of all PLMD::MultiColvar objects.  To see how to create 
documentation that describes the keywords that are specific to your particular CV please read \ref usingDoxygen.

\par Reading the atoms

Creating the lists of atoms involved in each of the colvars in a PLMD::MultiColvar is quite involved as you have to generate all the sets of 
atoms for which you want to calculate your \f$f\f$ function.  What is more this is
an area where we feel it is important to maintain some consistency in the input. Many of the multicolvars that are currently
implemented in the code thus use one or multiple of the following keywords: 

<table align=center frame=void width=95%% cellpadding=5%%>
<tr>
<td width=5%> ATOMS </td> <td> The atoms keyword specifies that one collective coordinate is to be calculated for each set of atoms specified.  
                               Hence, for MultiColvarDistance the command DISTANCES ATOMS1=1,2 ATOMS2=2,3 ATOMS3=3,4 specifies 
                               that three distances should be calculated. </td>
</tr> <tr>
<td width=5%> GROUP GROUPA GROUPB </td> <td> The GROUP keyword is used for quantities such as distances and angles.  A single GROUP specifies that
                                             a CV should be calculated for each distinct set of atoms that can be made from the group.  Hence,
                                             MUTIDISTANCE GROUP=1,2,3 specifies that three distances should be calculated (1,2), (1,3) and (2,3).
                                             If there is a GROUPA and GROUPB one CV is calculated for each set of atoms that includes at least one
                                             member from each group.  Thus MULTIDISTANCE GROUPA=1 GROUPB=2,3 calculates two distance (1,2) and (1,3). </td>
</tr> <tr>
<td width=5%> SPECIES SPECIESA SPECIESB </td> <td> The SPECIES keywords is used for quantities like coordination numbers.  The way this works is
                                                   best explained using an example.  Imagine a user working on NaCl wishes to calculate the average
                                                   coordination number of the sodium ions in his/her system.  To do this he/she uses COORDINATIONNUMBER SPECIES=1-100
                                                   which tells plumed that atoms 1-100 are the sodium ions.  To calculate the average coordination number
                                                   plumed calculates 100 coordination numbers (one for each sodium) and averages.  Obviously, each of these coordination
                                                   numbers involves the full set of 100 sodium atoms.  By contrast if the user wanted to calculate the coordination number
                                                   of Na with Cl he/she would do COORDINATIONNUMBER SPECIESA=1-100 SPECIESB=101-200, where obviously 101-200 are the 
                                                   chlorine ions.  Each of these 100 heteronuclear coordination numbers involves the full set of atoms speciefied using
                                                   SPECIESB and one of the ions speciefied by SPECIESA. </td>
</tr>
</table>

If you wish to use the ATOMS or SPECIES options above you can do so by adding:

\verbatim
keys.use("ATOMS");
\endverbatim

or 

\verbatim
keys.use("SPECIES"); keys.use("SPECIESA"); keys.use("SPECIESB")
\endverbatim

The documentation for these keywords will then be generated automatically.  If you wish to use the 
something like the GROUP keyword above then you need to register GROUP.  Here is how this is done within 
PLMD::multicolvar::Distances

\verbatim
keys.add("atoms-1","GROUP","Calculate the distance between each distinct pair of atoms in the group");
keys.add("atoms-2","GROUPA","Calculate the distances between all the atoms in GROUPA and all "
                            "the atoms in GROUPB. This must be used in conjunction with GROUPB.");
keys.add("atoms-2","GROUPB","Calculate the distances between all the atoms in GROUPA and all the atoms "
                            "in GROUPB. This must be used in conjunction with GROUPA.");
\endverbatim

To be clear, if these keywords are registered the actual readin of the atoms is looked after by the method
PLMD::multicolvar::MultiColvar::readAtoms.  However, there are a number of things I would like to point out 
about the registering syntax above.  Firstly, writing atoms-1 for GROUP and atoms-2 for GROUPA and GROUPB 
ensures nice formatting in the manual.  In particular it is made clear to the user that these are two 
distinct options for specificying the atoms in the input.   Notice also that we use similar commands to 
register keywords in PLMD::multicolvar::Angles:

\verbatim
keys.add("atoms-1","GROUP","Calculate angles for each distinct set of three atoms in the group");
keys.add("atoms-2","GROUPA","A group of central atoms about which angles should be calculated");
keys.add("atoms-2","GROUPB","When used in conjunction with GROUPA this keyword instructs plumed "
                            "to calculate all distinct angles involving one atom from GROUPA "
                            "and two atoms from GROUPB. The atom from GROUPA is the central atom.");
keys.add("atoms-3","GROUPC","This must be used in conjunction with GROUPA and GROUPB.  All angles "
                            "involving one atom from GROUPA, one atom from GROUPB and one atom from "
                            "GROUPC are calculated. The GROUPA atoms are assumed to be the central "
                            "atoms");
\endverbatim

Again though the actual reading of the GROUP, GROUPA, GROUPB and GROUPC keywords are looked after within
PLMD::multicolvar::MultiColvar::readAtoms

In some cases you may wish to use keywords other than ATOMS, GROUP, SPECIES etc.  In these cases the following
sets of routines may be useful:

<table align=center frame=void width=95%% cellpadding=5%%>
<tr>
<td width=10%> PLMD::multicolvar::MultiColvar::readAtomsLikeKeyword </td> <td> This is used to read in the ATOMS keyword.  By using this function you can change the ATOMS keyword to whatever you wish. </td>
</tr> <tr>
<td width=10%> PLMD::multicolvar::MultiColvar::readTwoGroups </td> <td> This is used to read in the GROUPA and GROUPB keywords.  When you use this funciton you can change GROUPA and GROUPB to whatever you desire </td>
</tr> <tr>
<td width=10%> PLMD::multicolvar::MultiColvar::readThreeGroups </td> <td> This is used to read in the GROUPA, GROUPB and GROUPC keywords. 
                                                                          When you use this function you can change GROUPA, GROUPB and GROUPC to whatever you desire.
                                                                          For an example see PLMD::multicolvar::Bridge </td>
</tr> <tr>
<td width=10%> PLMD::multicolvar::MultiColvar::readSpeciesKeyword </td> <td> This is used to read in SPECIESA and SPECIESB keywords.  Once again you can use this funciton and use different keywords. </td>
</tr>
</table>

\section label4 The constructor

Here is the constructor for the COORDINATIONNUMBERS multicolvar:

\verbatim
CoordinationNumbers::CoordinationNumbers(const ActionOptions&ao):
PLUMED_MULTICOLVAR_INIT(ao)
{
  // Read in the switching function
  std::string sw, errors; parse("SWITCH",sw);
  if(sw.length()>0){
     switchingFunction.set(sw,errors);
     if( errors.length()!=0 ) error("problem reading SWITCH keyword : " + errors );
  } else {
     double r_0=-1.0, d_0; int nn, mm;
     parse("NN",nn); parse("MM",mm);
     parse("R_0",r_0); parse("D_0",d_0);
     if( r_0<0.0 ) error("you must set a value for R_0");
     switchingFunction.set(nn,mm,r_0,d_0);
  }
  log.printf("  coordination of central atom and those within %s\n",( switchingFunction.description() ).c_str() );
  // Set the link cell cutoff
  setLinkCellCutoff( switchingFunction.get_dmax() );
  rcut2 = switchingFunction.get_dmax()*switchingFunction.get_dmax();

  // Read in the atoms
  int natoms=2; readAtoms( natoms );
  // And setup the ActionWithVessel
  checkRead();
}
\endverbatim

This code does three things:

- It reads in all the keywords that are specific to your new CV.  In this case it is a switching function that is read by parse("SWITCH",sw).
- It reads in the atoms involved in the CV.  In this case this is all looked after by the base class method PLMD::multicolvar::MultiColvar::readAtoms
- It sets the cutoff for link cells by exploiting the method PLMD::multicolvar::MultiColvarBase::setLinkCellCutoff.  This is really important - your code will run much faster if you exploit link cells. 
- It checks that readin was successful.  This is the final call to checkRead

Within the constructor of a multicolvar there must be a call to PLMD::multicolvar::MultiColvar::readAtoms even if you choose to read in the atoms using
one of the methods described in the table above.  

As you can see from the above example the majority of the constructor will be concerned with reading keywords from the input directive. More information
and guidance on PLUMED's parsing functionality can be found in \ref parsing.

\section per Implementing isPeriodic

Some \f$f\f$ functions (e.g. TORSIONS) have a co-domain that is periodic the majority, however, do not.  If the \f$f\f$ function of your 
new multicolvar falls into this category your isPeriodic function will be as follows:

\verbatim
bool isPeriodic(){ return false; }
\endverbatim

If your \f$f\f$ function does have a periodic co-domain then you will need to have methods like the following:

\verbatim
bool isPeriodic(){ return true; }
void retrieveDomain( std::string& min, std::string& max ){ min="-pi"; max="pi"; }
\endverbatim

As is hopefully obvious the retrieveDomain should return the minimum and maximum values (as strings) that your function, \f$f\f$, can take. 
The above example is taken from the code in PLMD::multicolvar::Torsions. 

\section label5 Compute

Compute is the routine that actually calculates the function \f$f\f$ and this function's derivatives from the set of atomic positions \f$\{X\}\f$.
This function should have as its return value the value of \f$f\f$.  The atomic positions can be extracted by using the PLMD::multicolvar::AtomValuePack::getPosition method 
that belongs to the PLMD::multicolvar::AtomValuePack class that is passed to the compute method.  One can then calculate the distance vector connecting two atoms by using
PLMD::multicolvar::MultiColvarBase::getSeparation.  Once you have calculated the derivatives of your \f$f\f$ function you can store these within the 
PLMD::multicolvar::AtomValuePack by using the methods PLMD::multicolvar::AtomValuePack::addAtomsDerivatives and PLMD::multicolvar::AtomValuePack::addBoxDerivatives.
 
\section mregtests Writing regression tests

Once you have written all the above function you will need to write regression tests for your new feature if you would like
us to incorporate it in the release version of PLUMED.  Details on how to write regression tests are provided here: \ref regtests  

