\mainpage Introduction

This is the developer manual. Please first have a look at the <a href="../../user-doc/html/index.html"> user manual </a>.

Plumed 2 is written in C++ and uses many of the advanced, object-oriented features of this language.  This structure makes the implementation of collective coordinates and free energy methods straightforward.  In fact, it should be possible to implement methods and collective coordinates (CV) by creating a single file and without touching any other part of the code. Futhermore, to implement new methodology does not require one to be some sort of C++ wizzard. Rather, the code has been specifically redisigned to make the implementation of new CVs and new free energy methods straightforward so as to encourage people to implement whatever new functionality they require.  This document serves then to provide an introduction as to how to go about implementing new functionality in plumed. A good starting point is \ref INHERIT as this page contains links to parts of the manual where you can find information on how to go about implementing CV, functions and biases. Another useful page is the \subpage TOOLBOX page, which contains information on the many reusable objects that have been implemented in plumed.  

If you want to understand a little more about the code and the way that we use the various features of C++ before you start then we describe this briefly here:

\ref ABriefIntroduction 

And finally, for the developers of MD codes, we provide information as to how to incorperate plumed into your codes here:

\ref HowToPlumedYourMD

If you would like to contribute new functionalities to PLUMED please read the following guidance:

\ref HowToContributeToPlumed

We ask that contributors endeavor to maintain the portability of plumed by, as much as possible, by only using the STL library and lapack in modifications.  
If you need to use any less standard library (e.g. Boost, Sockets) please ensure that your functionality is not installed during a default compilation.  
However, do feel free to provide alternative compilation options that incorperate your functionality.

Information about C++
http://www.parashift.com/c++-faq-lite/

\par Code Coverage

This manual might  also contain a detailed analysis of which parts of the PLUMED code have been tested when compiling it.
In case so, you will find it at <a href="../coverage/index.html"> this link </a>.
If this manual was compiled on Travis-CI, notice that as of PLUMED 2.5 the coverage scan ends up on a separate 
repository named `github.com/plumed/coverage-branchname`.

\defgroup TOOLBOX Tool Box
@{
Classes providing basic tools in plumed.

Classes of this group are designed to be reusable and to incorporate all sorts of functionality in plumed.
We try to keep their documentation as complete and clear as possible so as to increase the
chance that they will be reused.

If you implement a new class that you think might be useful to others please add it to the 
list by including the following inside the header file.
\verbatim
\ingroup TOOLBOX
\endverbatim
@}

\defgroup MULTIINHERIT Classes for multiple inheritance
@{
Classes for multiple inheritance.

Each of these classes implements some special feature which can
be then used to compose complex Actions.
All of them are "public virtual" derivatives of PLMD::Action,
so that it is possible to build ad Action which is based on multiple
classes from this group. This is the only place in the Action hierarchy
where multiple inheritance should be used.

Multiple inheritance allows for immediate combination of these features,
but add some C++ subtleties. If you do not fully understand them don't worry
and directly inherits from classes in the \ref INHERIT group.

To add a class to this group, just put a
\verbatim
\ingroup MULTIINHERIT 
\endverbatim
statement somewhere inside the header.
@}

\defgroup INHERIT Base classes for CVs, functions, biases, etc.
@{
Classes which can be used to create CVs, functions, biases and so on.

The typical way to add a new feature to plumed is to create a new class
which just inherits from one of the classes of this group. For example,
a new collective variable can be created by inheriting from PLMD::Colvar.
Most of the \ref INPUTDIRECTIVES are based on classes from this group.

To add a class to this group, just put a
\verbatim
\ingroup INHERIT
\endverbatim
statement somewhere inside the header.
@}

\defgroup INPUTDIRECTIVES Classes providing input directives
@{
Classes which implement directive that we used in the plumed input file.

Most of these classes are only used to provide a new feature which will be
available from the plumed input file. As such, they are typically not reused in
other places of the code. For this reason, almost all of them are directly
provided into an implementation file (.cpp), and have no accociated header file (.h).
A notable exceptions is PLMD::SetupMolInfo, which needs to be accessed directly
from other classes.




Each of these classes provides one directive for the plumed input file.
This list is built automatically based on the PLUMED_REGISTER_ACTION macro.
@}


