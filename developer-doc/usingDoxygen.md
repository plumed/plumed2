\page usingDoxygen Creating plumed documentation

Our main documentation effort for plumed is the plumed tutorials site (www.plumed-tutorials.org) and to a lesser 
extent the nest (www.plumed-nest.org). The manual is written in a way that supports these two efforts. A key 
part of our phillosophy is that <b> there is plenty of documentation for PLUMED in the papers that use PLUMED
and there should be links to these papers in the manual.</b>.  If you are writing documentation for an action 
the most useful thing you can do is add DOIs to your papers (see \ref dois and \ref moduledois). 

To create the plumed manual you should go to the root directory and type <b> make docs </b>. 
This command works because user documentation for all the PLMD::Action is inside the source code.  
Much of the text in the manual is generated automatically, which is possible as manual
pages for PLMD::Action are inside the code and because the manual is created using the 
following package:

- mkdocs:   https://www.mkdocs.org   (for user documentation)
- Doxygen:   http://www.doxygen.org  (for developer documentation)

In addition a special class, PLMD::Keywords, is used to store the descriptions of the syntax for any given
action so that this data can be produced in the output when the user makes a mistake in input.  In the following
a step-by-step explanaition as to how to use the documentation producing functionality of plumed is provided.  When
you have finished creating your new documentation they should be incorporated automatically in the manual the next 
time you make the manual.  Furthermore, a number of checks of your manual are performed when you make the manual
so any errors should be straightforward to find. 

\section structure The structure of the manual

When you visit the manual pages at https://www.plumed.org/doc-master/user-doc/html/ you should notice that the manual is structured
like the code.  There are pages for each of the actions and command line tools that are implemented with PLUMED and pages for each 
of the modules.  Every module page contains a table of the actions and command line tools that are implemented within that particular module
and a tag menu that allows you to search the actions.  There is then a page where all the actions are listed and a page where all the modules 
are listed.

The division of the code into modules allows us to acknowledge the many authors who have contributed functionality to PLUMED as you 
see in the documentation for the modules in the manual.  

In earlier versions of the manual we imposed a more complicated taxonomy on the actions in plumed. We chose to discard this as it is difficult 
to maintain when you have multiple contributors. Furthermore, for any person who wishes to impose some taxnonmy in PLUMED you can do by using
tags (see \ref tags) or by writing a tutorial. 

<b>Ultimately, our view is that tutorials that explain exactly how to perform a particular calculation
are the most useful documentation we can provide to users.</b> We learn to use codes by finding examples that are similar to what we want to 
do and then modifying them to our own purposes. It is rare to find anyone who learns by carefully reading and internalising the user manual. In short, 
<b>our aim with this documentation is to support people who want to show how the calculations that formed the basis of the argument in their papers were performed.</b>

\section incode In code documentation

A lot of the documentation for PLUMED is inside the PLUMED executible itself. This documentation is used when we construct the annotated examples on 
www.plumed-nest.org and www.plumed-tutorials.org.  The following sections explain how to construct this in code documentation.

\subsection registerkeys Registering Keywords

When you implement any new PLMD::Action in plumed you must first create the documentation for the keywords that you
will use to read the input for your new method.  In fact you cannot read in undocumented keywords using plumed.  You will also 
need to document the values that can be assed out of the action.  Again you cannot pass values unless they are documented.

The documentation for keywords and values is created in a static method of the action called registerKeywords.  This method
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

\subsection reading Reading the input keywords

Keywords are read in using either PLMD::Action::parse, PLMD::Action::parseVector, PLMD::Action::parseNumberedVector or PLMD::Action::parseFlag.  
These routines will use the information provided during keyword registration to check the sanity of any input.  For instance if you declare a 
compulsory keyword and do not specify a default value then the code will automatically complain if the particular keyword is missing from input.  
In addition, if the vector you pass to PLMD::Action::parseVector and PLMD::Action::parseNumberedVector has a
size greater than 0 plumed will assume that the input should contain a vector of this size and will complain if it finds a different sized vector.

\subsection argument-keywords Keywords for reading arguments

If you have an action that takes the values output by other actions in input.  In other words, if your action reads argument you need to use the 
special method for documenting these keywords that is shown below:

\verbatim
keys.addInputKeyword( attribute , keyword , type , explanation );
\endverbtim

The addInputKeyword method above is part of PLMD::Keywords.  This function takes the <i> attribute </i>, <i> keyword </i> and <i> explanation </i> arguments that were
discussed in the previous section as well as a new <i> type </i> attribute.  This <i> type </i> tells PLUMED what kind of input value is expected of the four basic 
types that PLUMED can pass about:

<table align=center frame=void width=95%% cellpadding=5%%>
<tr>
<td width=5%%> <b> scalar </b> </td> <td> a scalar </td> 
</tr> <tr>
<td width=5%%> <b> vector </b> </td> <td> a vector </td> 
</tr> <tr>
<td width=5%%> <b> matrix </b> </td> <td> a matrix </td> 
</tr> <tr>
<td width=5%%> <b> grid </b> </td> <td> a function whose values are stored on a grid </td> 
</tr>
</table>

You do not have to use a single <i> type </i> when specifying input keywords.  You can specity any combination of the four types in the table above by using types specified 
by slashes as indicated in the example below:

\verbatim
keys.addInputKeyword("compulsory","ARG","scalar/vector/matrix","the values input to this function");
\endverbatim

\subsection switching Documenting switching functions

If you have a keyword that takes in the input to a switching function you can add the following command in the registerKeywords function:

\verbatim
keys.linkActionInDocs("SWITCH","LESS_THAN");
\endverbatim

This command automatically adds the sentence, "Options for this keyword are explained in the documentation for LESS_THAN" into the documentation for the SWITCH keyword.
Furthermore, the LESS_THAN in this sentence is a link to the documentation for the LESS_THAN command.  You can obviously add a link to any action in the descriptions for keywords
by using this same mechanism.

\subsection shortcuts Documenting shortcut actions

If you have written a shortcut action that works by creating a more complicated PLUMED input from a simple one line command there are more things that you need to do in order
to construct the documentation.  I doubt anyone will do this so I have not documented this here.  If I am wrong and you think this documentation would be useful then email
gareth.tribello@gmail.com.

\subsection components Registering values and components

In plumed 2.1 we will also begin registering all the components that are calculated by your action.
In plumed 2.2 registering components will become compulsory and your features will not work if this is not done.

The registering of values means that in the registerKeywords method for commands such as:

\verbatim
d1: DISTANCE
\endverbatim

which calculates a quantity that can be reference in the input as d1 you have to provide documentation 
fr the manual that describes hwat information is stored in the d1 value.  As an example, for the
distances command above this documentation takes the following form:

\verbatim
keys.setValueDescription(scalarOrVector,"the DISTANCE between this pair of atoms");
\endverbatim

All values registered by plumed should be added using a variant of the setValueDescription of PLMD::Keywords.
This command has the following syntax:

\verbatim
keys.setValueDescription( type, explanation ) 
\endverbatim 

where <i> type </i> documents whether the output is a scalar, vector, matrix or grid (see \ref argument-keywords) 
and <i> explanation </i> is an explanation of the quantity that will appear in the manual.

The registering of components means that in the registerKeywords method for commands such as:

\verbatim
d1: DISTANCE COMPONENTS
\endverbatim

which calculates quantities that can be referenced in the input as d1.x, d1.y and d1.z, you will have to provide
documentation for the manual that describes what information is stored in the x, y and z components. As an example
for the distances components this documentation takes the following form:

\verbatim
keys.addOutputComponent("x","COMPONENTS",scalarOrVector,"the x-components of the vectors connecting the two atoms");
keys.addOutputComponent("y","COMPONENTS",scalarOrVector,"the y-components of the vectors connecting the two atoms");
keys.addOutputComponent("z","COMPONENTS",scalarOrVector,"the z-components of the vectors connecting the two atoms");
\endverbatim

As you can see this new feature works in a similar manner to the feature by which the keywords are registered.  Both
these features serve to keep the manual up to date with only a relatively small amount of effort on the part of the
developer.

All components registered by plumed should be added using a variant of the addOutputComponent method of PLMD::Keywords.
This command has the following syntax:

\verbatim
keys.addOutputComponent( name, keyword, type, explanation )
\endverbatim

where <i> name </i> is the name the component will be given i.e. the quantity will be referencable in input as <i> label.name </i>.
<i> keyword </i> is the name of any keyword that must be added to the input for the action in order to calculate this quantity, <i> type </i>
documents whether the output is a scalar, vector, matrix or grid (see \ref argument-keywords) and <i> explanation </i> is an explanation of the quantity that will appear in the manual.

If your Action always generates a particular set of components then the form of this command changes slightly. That is to say
if it makes no sense to use reference the isolated label for your command then this command should read:

\verbatim
componentsAreNotOptional(keys);
keys.addOutputComponent( name, "default", "scalar", explanation )
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

\subsection dois Referening papers

<b>If you have contributed a method to PLUMED you will have generally written extensive documentation about it in the articles 
you have written. You can add links to these papers on the action's manual page by using the addDOI method from PLMD::Keywords
as shown below:</b>

\verbatim
keys.addDOI(doi);
\endverbatim

This is the easist way to quickly document your action.

\subsection reserved Reserved Keywords

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

\subsection errors Generating errors

You may need to check for other mistakes in input.  When you find these mistakes you can report them to users using PLMD::Action::error.  This routine will 
output a description of how the input to your Action should appear.  It takes as input a string that describes the particular nature of the error that the user has made.

\section actiondocs Documenting actions

The remainder of the manual - the detailed description of your action and some examples of how the PLMD::Action can be used - is created 
from the comments at the top of the cpp file that contains the various subroutines that your PLMD::Action performs.  This documentation is written in 
github markdown as this allows one to incorporate equations, bibliographic information and bits and pieces of HTML. At the
start of the block of manual information the following lines should appear:

\verbatim
//+PLUMEDOC TAG1 ACTIONNAME TAG2 TAG3 
/*
\endverbatim

ACTIONAME is the first word that appears on the input line - i.e. it is the command that a user would use in input in order to make
use of your particular PLMD::Action.  TAG1, TAG2 and TAG3, are words that you can use for organising the actions in the manual (see \ref tags).

Immediately after the start comment symbol you should place a single line that describes in a sentence or two what it is your PLMD::Action does.  This information will
appear beside the link to your more detailed manual page in the general pages of the user manual.  The code will use everything up to the first blank
line in input to create this brief description.  You can then write a longer description of your PLMD::Action to appear on the 
particular page in the manual.  

You write your documentation in the github markdown syntax with additional plumed features that is explained here:

https://www.plumed-tutorials.org/instructions.html

\section moduledocs Documenting a module

The documentation for modules is included in two files that appear in the module directory that are called <i> module.md </i> and <i> module.yml </i>.  The <i> module.md </i>
file is a markdown file that contains a description of the module.  This is written in the github markdown syntax that is explained here:

https://www.plumed-tutorials.org/instructions.html

\subsection moduledois The module.json file

The <i> module.yml </i> file is a Yaml file with the following structure:

\verbatim
name: name of module
authors: authors of module
description: description of module
dois: ["doi1", "doi2", "doi2"]
tags:
  TAG1: description of TAG1
  TAG2: description of TAG2
\endverbatim

As you can see you use this file to provide information on the name, authors and description of the module.  You can also include a list of DOIs for papers that should be cited on the 
module's documentation page.

\subsection tags Using action tags

As you can see in the previous section you can define tags the <i>module.yml</i> file and give them a short description.  You can then attach these tags to actions by including them on
the line that includes the string PLUMEDOC.  These tags provide a community-sourced way of grouping actions that can serve similar purposes in the manual.  

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
unofficial or unsupported branch (see .ci/push script).

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


