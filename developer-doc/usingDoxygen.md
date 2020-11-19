\page usingDoxygen Creating plumed documentation

To create the plumed manual you should go to the <b> user-doc </b> directory and type <b> make </b>. 
This command works because user documentation for all the PLMD::Action is inside the source code.  If
you look at the documentation page for any of the actions that are implemented in plumed you will
see that it is composed of three pars:

- A short introduction which describes what the method does.
- A description of the various keywords for the calculation.
- An example/some examples of how the PLMD::Action can be used.

Furthermore, you will also have noticed that if you make an error in the input for any PLMD::Action the 
descriptions of all the keywords from the manual appears in the log file.  This is possible because manual
pages for PLMD::Action are inside the code and because the manual is created using the 
following packages:

- Doxygen:   http://www.doxygen.org
- Graphviz:  http://www.graphviz.org/ 

In addition a special class, PLMD::Keywords, is used to store the descriptions of the syntax for any given
action so that this data can be produced in the output when the user makes a mistake in input.  In the following
a step-by-step explanaition as to how to use the documentation prodcuing functionality of plumed is provided.  When
you have finished creating your new documentation they should be incorporated automatically in the manual the next 
time you make the manual.  Furthermore, a number of checks of your manual are performed when you make the manual
so any errors should be straightforward to find. 

The plumed manual does not only contain descriptions of the various actions that have been implemented and their 
keywords. There are also pages that place these actions in context and which attempt to explain the interplay 
between colvars, functions, biases, analysis methods and so on.  More importantly, the plumed manual contains
instructive examples that explain to users how to use particular features of the code.  We would encourage all
users of the code to submit these sort of things to the plumed repository particularly if they have found it difficult to 
start using a particular method.  Instructions as to how to go about writing a material that will appear in the How-tos 
section of the manual can be found \ref tutorials here

\section registerkeys Registering Keywords

When you implement any new PLMD::Action in plumed you must first create the documentation for the keywords that you
will use to read the input for your new method.  In fact you cannot read in undocumented keywords using plumed.  The
documentation for keywords is created in a static method of the action called registerKeywords.  This method
should be declared in the definition of the class as follows:

\verbatim
static void registerKeywords(Keywords& keys);
\endverbatim

The static attribute allows one to use the registerKeywords routine when no instances of the class have been created.
This is essential as keywordRegistration in plumed is done before the list of PLMD::Action is created.  This means
that when the keywords are created plumed has no understanding of the hierarchy of inherited Actions.  Hence, before
adding your own Keywords you must ensure that the keywords for the class from which your new class inherits have been
added.  In pracise this is done by calling PLMD::Colvar::registerKeywords, PLMD::function::Function::registerKeywords or 
PLMD::Bias::registerKeywords for a new PLMD::Colvar, PLMD::function::Function or PLMD::bias::Bias respectively.  To be clear these
functions will ensure that generic keywords such as LABEL, NUMERICAL_DERIVATIVES or ARG are registered and explained in the
manual. If your method requires the derivatives of some value and you have no way of implementing the analytical derivatives
you should also call PLMD::ActionWithValue::noAnalyticalDerivatives.  This routine will ensure that plumed's numerical
derivatives routines are used to calculation your derivatives automatically and will ensure that a message is put in the plumed
manual so that other users are aware that numerical derivatives are being used. 

Once you have called the reigsterKeywords routine for the PLMD::Action above yours in the hierarchy you can begin to add
the keywords for your new method.  These keywords will have one of 5 attributes:

<table align=center frame=void width=95%% cellpadding=5%%>
<tr>
<td width=5%%> <b> compulsory </b> </td> <td> These are the quantities that must be defined in order to perform your action </td> 
</tr> <tr>
<td> <b> optional </b> </td> <td> If there is some alternate way of performing your calculation that requires numerical input you should declare your keyword as optional </td>
</tr> <tr>
<td> <b> flag </b> </td> <td> This is used to declare keywords such as NOPBC that tell plumed to turn on/off some feature of the calculation </td>
</tr> <tr>
<td> <b> atoms </b> </td> <td> If you are reading a list of atoms after the keyword then you should use this keyword. You can easily specify if there are multiple ways of defining the atoms involved in the action in the manual.  To register the keywords for this first method of specifying the atoms using atoms-1. Then register the keywords for the second way of specifying the atoms using atoms-2 and so on.  A manual that states that these keywords can be used in an either or fashion, much that for <a href="../../user-doc/html/_t_o_r_s_i_o_n.html"> TORSION </a>, will then be generated.  </td>
</tr> <tr>
<td> <b> numbered </b> </td> <td> If you need to read in a list of similar keywords such as keyword0, keyword1, keyword2... then you must use this option.  These keywords will be assumed to be optional.  However, you can set them to be atoms or whatever by using reset_style(keyword,newstyle).  </td>
</table>

All keywords (other than flags) are added using the add method of PLMD::Keywords.  This command has the following syntax:

\verbatim
keys.add( attribute, keyword, explanation );
\endverbatim

where <i> attribute </i> is one of the options from the above table, <i> keyword </i> is the word that appears on the input line and <i> explanation </i> is an explantion
of what the keyword does.  If your keyword is compulsory it can also be added using:

\verbatim
keys.add( attribute, keyword, default, explanation );
\endverbatim

where <i> default </i> is a string containing the default value to use for the quantity.

Flags are added using the add flag method, this has syntax:

\verbatim
keys.addFlag( keyword, default, explantion );   
\endverbatim

where default is a bool that tells plumed if by default this option is/is not in use.  

\section reading Reading the input keywords

Keywords are read in using either PLMD::Action::parse, PLMD::Action::parseVector, PLMD::Action::parseNumberedVector or PLMD::Action::parseFlag.  
These routines will use the information provided during keyword registration to check the sanity of any input.  For instance if you declare a 
compulsory keyword and do not specify a default value then the code will automatically complain if the particular keyword is missing from input.  
In addition, if the vector you pass to PLMD::Action::parseVector and PLMD::Action::parseNumberedVector has a
size greater than 0 plumed will assume that the input should contain a vector of this size and will complain if it finds a different sized vector.

\section components Registering components

In plumed 2.1 we will also begin registering all the components that are calculated by your action.
In plumed 2.2 registering components will become compulsory and your features will not work if this is not done.
This registering of components means that in the registerKeywords method for commands such as:

\verbatim
d1: DISTANCE COMPONENTS
\endverbatim

which calculates quantities that can be referenced in the input as d1.x, d1.y and d1.z, you will have to provide
documentation for the manual that describes what information is stored in the x, y and z components. As an example
for the distances components this documentation takes the following form:

\verbatim
keys.addOutputComponent("x","COMPONENTS","the x-component of the vector connecting the two atoms");
keys.addOutputComponent("y","COMPONENTS","the y-component of the vector connecting the two atoms");
keys.addOutputComponent("z","COMPONENTS","the z-component of the vector connecting the two atoms");
\endverbatim

As you can see this new feature works in a similar manner to the feature by which the keywords are registered.  Both
these features serve to keep the manual up to date with only a relatively small amount of effort on the part of the
developer.

All components registered by plumed should be added using a variant of the addOutputComponent method of PLMD::Keywords.
This command has the following syntax:

\verbatim
keys.addOutputComponent( name, keyword, explanation )
\endverbatim

where <i> name </i> is the name the component will be given i.e. the quantity will be referencable in input as <i> label.name </i>.
<i> keyword </i> is the name of any keyword that must be added to the input for the action in order to calculate this quantity
and <i> explanation </i> is an explanation of the quantity that will appear in the manual.

If your Action always generates a particular set of components then the form of this command changes slightly. That is to say
if it makes no sense to use reference the isolated label for your command then this command should read:

\verbatim
componentsAreNotOptional(keys);
keys.addOutputComponent( name, "default", explanation )
\endverbatim

So for example in RESTRAINT, which always generates a component called bias and a component called force2, the components are registered
using the following code:

\verbatim
componentsAreNotOptional(keys);
keys.addOutputComponent("bias","default","the instantaneous value of the bias potential");
keys.addOutputComponent("force2","default","the instantaneous value of the squared force due to this bias potential");
\endverbatim

Lastly, if you have some method that takes in arguments and gives back one component per argument then these components should be
labelled <i>argument-name</i>_<i>description</i>.  The MOVINGRESTRAINT command below gives an example of how this is done in practise.

\verbatim
DISTANCE ATOMS=1,2 LABEL=d1
DISTANCE ATOMS=3,4 LABEL=d2
MOVINGRESTRAINT ARG=d1,d2 AT0=2,2 AT1=6,6 STEP0=0 STEP1=100 KAPPA=1
\endverbatim

This command has components called d1_steer, d2_steer, d1_cntr and d2_cntr.  These components describe the work done in moving the 
system along the d1 and d2 axis and the instantaneous positions of the harmonic potential on the d1 and d2 axis respectively.

The various components in MOVINGRESTRAINT are registered using:

\verbatim
componentsAreNotOptional(keys);
keys.addOutputComponent("bias","default","the instantaneous value of the bias potential");
keys.addOutputComponent("force2","default","the instantaneous value of the squared force due to this bias potential");
keys.addOutputComponent("_cntr","default","one or multiple instances of this quantity will be refereceable elsewhere in the input file. "
                                          "these quantities will named with  the arguments of the bias followed by "
                                          "the character string _cntr. These quantities give the instantaneous position "
                                          "of the center of the harmonic potential.");
keys.addOutputComponent("_work","default","one or multiple instances of this quantity will be refereceable elsewhere in the input file. "
                                          "These quantities will named with the arguments of the bias followed by "
                                          "the character string _work. These quantities tell the user how much work has "
                                          "been done by the potential in dragging the system along the various colvar axis.");
\endverbatim

Appropriately labelled components can be created using the following command:

\verbatim
comp=getPntrToArgument(i)->getName()+"_cntr";
addComponent(comp); componentIsNotPeriodic(comp);
\endverbatim

They can then be set by using something like:

\verbatim
getPntrToComponent(getPntrToArgument(i)->getName()+"_work")->set(val[i]);
\endverbatim

\section reserved Reserved Keywords

To maintain some consistency for end users of the code certain keywords (e.g. ARG, PERIODIC) are reserved.
The reserved keywords for PLMD::Colvar, PLMD::function::Function and PLMD::bias::Bias are explained inside the documentation for
these actions.  To use one of the registered keywords you shold insert the following command into the registerKeywords
method of your new function.

\verbatim
keys.use( keyword );
\endverbatim

where <i> keyword </i> is a string that tells the method which reserved keyword you wish to use.  To be clear
when you use a reserved keyword all the parsing and registering for it is looked after automatically. In addition,
if new components are generated in the output of the action when the keyword is present the registering of those
components is also looked after elsewhere.  So, for example, if you do something like:

\verbatim
keys.use("MIN")
\endverbatim

In a new PLMD::multicolvar::MultiColvar then there is no need to add the command:

\verbatim
keys.addOutputComponent("min","MIN","the minimum value. This is calculated using the formula described in the description of the "
                                    "keyword so as to make it continuous.");
\endverbatim 

As this documentation will already have been created for you elsewhere in the code.

\section errors Generating errors

You may need to check for other mistakes in input.  When you find these mistakes you can report them to users using PLMD::Action::error.  This routine will 
output a description of how the input to your Action should appear.  It takes as input a string that describes the particular nature of the error that the user has made.

\section manual Creating the rest of the manual

The remainder of the manual - the detailed description of your action and some examples of how the PLMD::Action can be used - is created 
from the comments at the top of the cpp file that contains the various subroutines that your PLMD::Action performs.  This is converted
to manual using Doxygen as this allows one to incorporate equations, bibliographic information and bits and pieces of HTML.  At the
start of the block of manual information the following lines should appear:

\verbatim
//+PLUMEDOC TYPE ACTIONNAME 
/*
\endverbatim

ACTIONAME is the first word that appears on the input line - i.e. it is the command that a user would use in input in order to make
use of your particular PLMD::Action.  TYPE, meanwhile, tells Doxygen where in the manual the Docuementation should be placed.  TYPE
should be one of the following:

<table align=center frame=void width=95%% cellpadding=5%%>
<tr> 
<td width=5%%> <b> COLVAR </b> </td> <td> This is used if your PLMD::Action is the calculation of a CV </td>
</tr> <tr>
<td width=5%%> <b> DCOLVAR </b> </td> <td> This is used if your PLMD::Action is a colvar that measures the distance from a reference frame in some metric </td>
</tr> <tr>
<td width=5%%> <b> MCOLVAR </b> </td> <td> This is used if your PLMD::Action calculates some Function of a distribution of CVs. (Action inherits from PLMD::multicolvar::MultiColvar)  </td>
</tr> <tr>
<td width=5%%> <b> MCOLVARF </b> </td> <td> This is used if your PLMD::Action calculates a paricularly complicated function of a distribution of CVs. (Action inherits from PLMD::multicolvar::MultiColvarFunction) </td>
</tr> <tr>
<td width=5%%> <b> FUNCTION </b> </td> <td> This is used if your PLMD::Action calculates some Function of a set of CVs </td>
</tr> <tr>
<td width=5%%> <b> VATOM </b> </td> <td> This is used if your PLMD::Action calculates the position of a new atom e.g. for COM </td>
</tr> <tr>
<td width=5%%> <b> ANALYSIS </b> </td> <td> This is used if your PLMD::Action does some analysis of the trajectory </td>
</tr> <tr>
<td width=5%%> <b> BIAS </b> </td> <td> This is used if your PLMD::Action is a bias that adds supplemental forces to the potential in order to enhance sampling </td> 
</tr> <tr>
<td width=5%%> <b> TOPOLOGY </b> </td> <td> PLMD::setup::MolInfo has documentation of this type.  This command is used to provide information about the chemistry of the system under study.  For MolInfo this is what constitute the backbone atoms of the protin, what the residues are etc. </td>
</tr> <tr>
<td width=5%%> <b> GENERIC </b> </td> <td> This should be used if you want to specify manually where in the manual your documentation should appear.  If you feel this really is the correct way to incorporate your new feature please contact the core developers so as to discuss it. </td>
</tr> 
</table>

Immediately after the start comment symbol you should place a single line that describes in a sentence or two what it is your PLMD::Action does.  This information will
appear beside the link to your more detailed manual page in the general pages of the user manual.  The code will use everything up to the first blank
line in input to create this brief description.  You can then write a longer description of your PLMD::Action to appear at the start of its
particular page in the manual.  As described below this description can incorporate equations and bibliographic information.

\subsection Equations

You can add formulae in latex using:

\verbatim
This is an inline equation \f$s=y+x\f$ but this is an equation:

\f[
r = \sqrt{ \mathbf{s}^T \mathbf{C}^{-1} \mathbf{s} }
\f]

And this is an equation array:

\f{eqnarray*}{
 f &=& \frac{1}{2} \\
 g &=& \frac{2}{3}
\f}
\endverbatim

In the manual this will be translated into:

This is an inline equation \f$s=y+x\f$ but this is an equation:
 
\f[
r = \sqrt{ \mathbf{s}^T \mathbf{C}^{-1} \mathbf{s} }
\f]

And this is an equation array:

\f{eqnarray*}{
 f &=& \frac{1}{2} \\
 g &=& \frac{2}{3}
\f} 

\subsection Lists

You can create lists of data using: 

\verbatim
- First item in list
- Second item in list
\endverbatim

which becomes:

- First item in list
- Second item in list

\subsection Formatting 

You can create a new section in your documentation using:

\verbatim
\section manual Creating the rest of the manual
\endverbatim

In fact I used this very command earlier in writing this page.  I can therefore reference it here (\ref manual) by using:

\verbatim
\ref manual 
\endverbatim

You can also reference external webpages by typing web addresses directly in the documentation.

\subsection Citations

You can create citations using:

\verbatim
\cite bibtex-tag
\endverbatim

This command uses an interface between Doxygen and bibtex to create bibliographic data.  Inside
the user-doc directory you will find a bibtex file called bibliography.bib that contains all
the references that are included in the user documentation for plumed.  To add your reference
you should add bibliographic data for the article you want to cite in this file.

\section Examples

Manual entries for actions and tutorials <b>must</b> contain some examples.  The most basic way to include these is as follows:

\verbatim
\par Example

The following input tells plumed to print the distance between atoms 3 and 5,
the distance between atoms 2 and 4 and the x component of the distance between atoms 2 and 4.
\plumedfile
DISTANCE ATOMS=3,5             LABEL=d1
DISTANCE ATOMS=2,4 COMPONENTS  LABEL=d2
PRINT ARG=d1,d2,d2.x
\ endplumedfile /*** But with no space between the \ and the endplumedfile
\endverbatim 

In the manual this will be converted to:

\par Example

The following input tells plumed to print the distance between atoms 3 and 5,
the distance between atoms 2 and 4 and the x component of the distance between atoms 2 and 4.
<pre class="fragment">
<a href="../../user-doc/html/_d_i_s_t_a_n_c_e.html" style="color:green">DISTANCE</a> ATOMS=3,5             LABEL=d1
<a href="../../user-doc/html/_d_i_s_t_a_n_c_e.html" style="color:green">DISTANCE</a> ATOMS=2,4 COMPONENTS  LABEL=d2
<a href="../../user-doc/html/_p_r_i_n_t.html" style="color:green">PRINT</a> ARG=d1,d2,d2.x
</pre>

Please be aware of the blank line between after the title of the paragraph.  If this line is not present your manual will look ugly.  
Also be aware that your Examples section <b> must </b> be called Examples and not Example because of a perculiarity in the 
script that generates the manual.

By including the example input in a plumedfile environment you ensure two things:

- That the action names are converted to links to the relevant pages in the manual when the manual is constructed.
- That the code to construct the user manual will test to see if your example input can be parsed by PLUMED whenever the user manual is built.

To achieve the second of these objectives with the input shown above it is sufficient to include the example input in a plumedfile environment.
As detailed in the following sections, however, there are some cases where things are a little more complicated.

\subsection multirepeg Including example inputs for multiple replica simulations

If you have an input for a simulation that is to be run with three replicas such as the one below:

<pre class="fragment">
<span style="color:blue"># Compute a distance</span>
d: <a href="../../user-doc/html/_d_i_s_t_a_n_c_e.html" style="color:green">DISTANCE</a> ATOMS=1,2
<span style="color:blue"># Apply a restraint.</span>
<a href="../../user-doc/html/_r_e_s_t_r_a_i_n_t.html" style="color:green">RESTRAINT</a> ARG=d AT=@replicas:1.0,1.1,1.2 KAPPA=1.0
</pre>

Then you must specify that the input is to be run on three replicas in the first (SETTINGS) line of the input file as shown below: 

\verbatim
\plumedfile{3}
#SETTINGS NREPLICAS=3
# Compute a distance
d: DISTANCE ATOMS=1,2
# Apply a restraint.
RESTRAINT ARG=d AT=@replicas:1.0,1.1,1.2 KAPPA=1.0
\ endplumedfile /*** But with no space between the \ and the endplumedfile
\endverbatim

Notice that there should not be a space between the hash sign at the start of this line and word settings. 

\subsection auxfileeg Including example inputs that require an auxiliary file

Suppose that you have an input such as the one below:

<pre class="fragment">
<a href="../../user-doc/html/_r_m_s_d.html" style="color:green">RMSD</a> REFERENCE=file.pdb TYPE=OPTIMAL
</pre>

As RMSD has been used here you are also required to provide an input file which in this case would be called file.pdb.  You can include 
this input in an auxfile environment as shown below:

\verbatim
\auxfile{file.pdb}
ATOM      1  CL  ALA     1      -3.171   0.295   2.045  1.00  1.00
ATOM      5  CLP ALA     1      -1.819  -0.143   1.679  1.00  1.00
ATOM      6  OL  ALA     1      -1.177  -0.889   2.401  1.00  1.00
ATOM      7  NL  ALA     1      -1.313   0.341   0.529  1.00  1.00
ATOM      8  HL  ALA     1      -1.845   0.961  -0.011  1.00  1.00
END
\ endauxfile /*** But with no space between the \ and the endauxfile
\endverbatim

Obviously, the file.pdb inside the curly braces in the top line here indicates that the auxiliary file to be constructed from this data should be named 
file.pdb.  Files input in this way can be given any name but:

- If two auxfiles are used on the same page they must be given different names (if they are on different pages it does not matter)
- auxfiles should not be named *.dat as the script that builds the user manual assumes that all *.dat files are plumed input files. 

\subsection incfileeg Using INCLUDE in your example input files

Suppose that you have split your input by using an INCLUDE file as shown below:

<pre class="fragment">
<a href="../../user-doc/html/_d_i_s_t_a_n_c_e.html" style="color:green">DISTANCE</a> ATOMS=1,2 LABEL=dist
<a href="../../user-doc/html/_i_n_c_l_u_d_e.html" style="color:green">INCLUDE</a> FILE=toBeIncluded.inc
</pre>

<pre class="fragment">
<span style="color:blue"># this is toBeIncluded.inc</span>
<a href="../../user-doc/html/_r_e_s_t_r_a_i_n_t.html" style="color:green">RESTRAINT</a> ARG=dist AT=2.0 KAPPA=1.0
</pre>

To include an input like this in the manul you would write the following:

\verbatim
\plumedfile
DISTANCE ATOMS=1,2 LABEL=dist
INCLUDE FILE=toBeIncluded.inc
\ endplumedfile    /*** But with no space between the \ and the endplumedfile

\plumedfile
#SETTINGS FILENAME=toBeIncluded.inc  
RESTRAINT ARG=dist AT=2.0 KAPPA=1.0
\ endplumedfile   /*** But with no space between the \ and the endplumedincludefile
\endverbatim

By including the FILENAME attribute on the SETTINGS line you can set the name of the plumed input file that is generated when the input is tested.
Also notice that if, as in the example above, the included file is not (by itself) a valid plumed input it CANNOT be called *.dat as the script that 
checks the input will complain.  

\subsection molfileeg Using MOLFILE in your example input files

If you use have used a \ref MOLINFO command in the example input that you specified as has been done here:

<pre class="fragment">
<a href="./_m_o_l_i_n_f_o.html" style="color:green">MOLINFO</a> STRUCTURE=helix.pdb
<a href="./_w_h_o_l_e_m_o_l_e_c_u_l_e_s.html" style="color:green">WHOLEMOLECULES</a> ENTITY0=1-100
alpha: <a href="./_a_l_p_h_a_r_m_s_d.html" style="color:green">ALPHARMSD</a> RESIDUES=all TYPE=OPTIMAL R_0=0.1
</pre> 

Then you must provide information on the location from whence PLUMED can the reference input so that the example checking script can copy the input
for the MOLINFO.   The above input would thus be included in the manual as shown below:

\verbatim
\plumedfile
#SETTINGS MOLFILE=regtest/basic/rt32/helix.pdb
MOLINFO STRUCTURE=helix.pdb
WHOLEMOLECULES ENTITY0=1-100
alpha: ALPHARMSD RESIDUES=all TYPE=OPTIMAL R_0=0.1
\ endplumedfile    /*** But with no space between the \ and the endplumedfile
\endverbatim

\subsection otherfiles Other actions requiring external files/folder

Other actions in plumed may require reading input files, examples include reading gromacs .ndx files in \ref GROUP, reading chemical shifts in \ref CS2BACKBONE, etc.
To make these example work correctly in the manual you can use the keywords AUXFILE and AUXFOLDER as in the following:

\verbatim
\plumedfile
#SETTINGS MOLFILE=regtest/basic/rt77/peptide.pdb
MOLINFO MOLTYPE=protein STRUCTURE=peptide.pdb
WHOLEMOLECULES ENTITY0=1-111

# This allows us to select only non-hydrogen atoms
#SETTINGS AUXFILE=regtest/basic/rt77/index.ndx
protein-h: GROUP NDX_FILE=index.ndx NDX_GROUP=Protein-H

# We extend the cutoff by 0.1 nm and update the neighbor list every 40 steps
solv: EEFSOLV ATOMS=protein-h

# Here we actually add our calculated energy back to the potential
bias: BIASVALUE ARG=solv

PRINT ARG=solv FILE=SOLV

#SETTINGS AUXFOLDER=regtest/isdb/rt-cs2backbone/data NREPLICAS=2
cs: CS2BACKBONE ATOMS=1-174 DATADIR=data/
encs: ENSEMBLE ARG=(cs\.hn-.*),(cs\.nh-.*)
stcs: STATS ARG=encs.* SQDEVSUM PARARG=(cs\.exphn-.*),(cs\.expnh-.*)
RESTRAINT ARG=stcs.sqdevsum AT=0 KAPPA=0 SLOPE=24

PRINT ARG=(cs\.hn-.*),(cs\.nh-.*) FILE=RESTRAINT STRIDE=100
\ endplumedfile    /*** But with no space between the \ and the endplumedfile
\endverbatim

\section tutorials Writing how-to instructions

On every page of the plumed user manaul there are three tabs: Main-page, Glossary and How-tos.  Here we are going to describe how to
go about writing something that will appear on the How-tos page.  On that page you will find that there are two kinds of resources.  
The first set of resources are a set of links to shortist descriptions of how to particular things with plumed.  The second set of 
resources are links to a set of external websites that we think might be of interest to users of plumed.

\subsection websites Adding a link to your website 

If you have a website that you think it would be useful for us to link from the plumed How-tos page please send us a file that contains
the following information in the following format:

\verbatim
link: http://en.wikipedia.org/wiki/Metadynamics

description: A wikipedia article on metadynamics
\endverbatim 

In case it isn't abundantly clear the line that starts "link:" contains the hyperlink to your page, while the second line (the one starting "description:")
contains the description of the page that will appear beside the link on the how-tos page.  If this file is placed in the user-doc/tutorials directory of the
plumed source and if it is given a name that ends in .site then the link will appear on the how-tos page.

\subsection tute Writing a how-to

Lets say you now want to write a set of how-to instructions.  You will first need to create an additional file in the user-doc/tutorials directory of the plumed
source.  Much like the rest of plumed manual this file is written in Doxygen so you should read the instructions about how to go about using \ref Equations, \ref Lists, 
\ref Formatting, \ref Citations and \ref Examples.  You will also need to ensure that the text that you want to appear on the your page is contained between
the special doxygen comment symbols as shown below:

\verbatim
This will not appear on the final how-to page.

/**
This will appear on the final how-to page.
*/

This will also not appear on the final how-to page.
\endverbatim 

To ensure that a link is created on the main How-tos page you need to include the following instructions after the closing doxygen comment as shown below:

\verbatim
/**
\page filename My explanation of something in plumed

Text of how-to page

*/
link: @subpage filename

description: The description of what is described on your page
\endverbatim

For this particular example the above should be contained in a file inside the user-doc/tutorials directory called filename.txt.  If the user comes to the 
How-tos page they will see a link with the text "My explanation of something in plumed," which they can click on to get to your page.  Beside that link the 
description "The description of what is described on your page" will appear.

One final thing, if your How-to is a tutorial - that is to say if you have a set of exercises for users to work through - then it may be useful to 
provide them with some example files that they can download.  This is possible you add the files to our repository and you can put a link to download them
somewhere on your page.
To keep things manageable there should only be one (or a few) tar-ball file per tutorial.
As such if you need to provide users with multiple files please put them in a single directory. The script that 
build the documentation will make one or more tar-balls from there.
please put them in a zipped tar ball.
The simplest way is thus to add a directory
in the user-doc/tutorials directory and ideally given a name so that it can be identified
with your corresponding *.txt file.
For example if the name of the tutorial is `mytute` you should make a directory named
`user-doc/tutorials/mytute`. 
You can then download it by including the following into your *.txt file.

\verbatim
/**
\page mytute A tutorial that helps users learn to do something with plumed

Please download the following <a href="tutorial-resources/mytute.tar.gz" download="mytute.tar.gz"> file </a>.

*/

link: @subpage mytute

description: A tutorial for users to work through

additional-files: mytute
\endverbatim 

In this case the tar ball you add is called mytute.tar.gz.  The user can download this file by clicking on the word file.

\section updating-web-manuals Updating web manuals

Precompiled versions of PLUMED manuals can be found on github at an address such as http://www.plumed.org/doc-v2.1/user-doc/html/index.html
(replace v2.1 with the actual version number). These manuals take advantage of a nice github feature: any branch named gh-pages
is shown as a webpage. In this example, the repository is located at http://github.com/plumed/doc-v2.1 .
Before version 2.1.1 it was necessary to upload the precompiled manual by hand. Since version 2.1.1, this is done
from Travis CI automatically whenever a commit containing in its log the tag [makedoc] is pushed into the plumed2 github repository.
Since version 2.3.3, manual is always updated, and tag [makedoc] is ignored.
Notice that Travis CI will try to push the manual on a repository named http://github.com/plumed/doc-NAMEOFTHEBRANCH , so that 
this should work for all the future release branches as long as an appropriate repository is created on the github.com/plumed
organization.
We could even easily create repositories to host the documentation of temporary branches.
Also notice that these repositories (plumed/doc-xxx) need to give write access to a dummy github account (PlumedBot). A token
for that user enabling html access is stored in the environment variable GIT_TOKEN which is saved (not visible) on travis-ci.org.
In this way, any commit made in the plumed repository by one of the developers will have access to the variable and will trigger
manual build and push. Conversely, pull requests by external users should not be able to
access the token and won't update manual changes.
Starting with PLUMED 2.3 the web manual also contains a coverage scan that shows which parts of the code are actually
covered by the regtests. As of PLUMED 2.5, the coverage scan ends up on a separate repository named 
http://github.com/plumed/coverage-NAMEOFTHEBRANCH.

Notice that to solve [this issue](https://github.com/plumed/plumed2/issues/239) as of PLUMED 2.3.2 the script that
pushes the documentation to travis-ci adds special information to remove from search engine results pages from
unofficial or unsupported branch (see .travis/pushdoc script).

Bottom line: manual will always be updated after a commit that can pass the tests.
Twenty minutes or so after your push the manual should be up to date, remember to double check on the web
and to revert the commit if there are errors!

It is possible to generate PLUMED manuals for your own personal forks 
using a similar procedure as described above. 
For this to work you need to enable Travis CI for your forked repository 
and define appropriately the environment variables on Travis CI. 
The github account used to automatically push the generated manuals 
should be defined using the `GIT_BOT` variable, 
preferably this should be a dummy account. A github token
enabling html access for that account should be defined using the `GIT_TOKEN` variable. 
Furthermore, you need to define an email address associated to the account using the `GIT_BOT_EMAIL` variable. 
It is better to make all these environment variable hidden such that they are 
not shown in the public logs on travis-ci.org. 
To generate a manual for a specific branch you need to create a repository 
`USERNAME/doc-NAMEOFTHEBRANCH` and give write access to the account given in 
`GIT_BOT`. The generated manuals will be accessible on 
https://USERNAME.github.io/doc-NAMEOFTHEBRANCH. Note that manuals generated in 
this way will always be labeled as unofficial and not shown in search engine results.
Starting with PLUMED 2.5, if you want to show the results of the coverage scan you should
similarly create arepository `USERNAME/coverage-NAMEOFTHEBRANCH`.


