\page HowToContributeToPlumed How to contribute new functionality to PLUMED

We welcome researchers to contribute new functionality to the PLUMED code.  In fact, we would 
argue that new biasing methods, cvs and analysis tools should be implemented in common software 
libaries and that the community should adopt shared standards when it comes to the format of input files
and output files in order to make collaboration and communication more straightforward.  There are 
some (not particularly onerous) caveats, however. We would thus ask any person who is 
considering contributing some functionality to PLUMED to read the following page carefully before 
commencing. 

\section maintanence Our phillosophy on sharing code

Writing programs is not so difficult.  Writing software that other people can use and that can
be maintaned is a lot of work, however.   In fact much of the work that we (the core developers)
do to maintain PLUMED does not involve coding.  It involves ensuring that the code has a sufficiently
large test suite, ensuring that old features are not broken when new methods are added to the code,
updating the manual, answering users questions on the mail list and organising PLUMED
meetings and tutorials.  None of the core developers are employed to do these things full-time and as 
such we would like to minimise the amount of time we spend on this sort of maintence.  As such it is
very important to us that we <b>ensure that we understand who wrote every file within PLUMED</b>.  There are
two reasons why we believe this is important:

- It is important that credit is given where credit is due.  In other words, we want users to know who
wrote the functionalities they are using and not to simply assume that we (the core developers) wrote 
everything in PLUMED.  
- The author of any feature is responsible for ensuring their feature is maintained in perpetuity.  
That is to say we (the core developers) will not take responsibility for maintaining features that we
did not ourselves write.  Furthermore, this maintenance involves answering users questions on the mail
lists as well as the occasional fix to the relevant cpp files.

In addition, to these two requirements we do not want the core-developer team to expand as we believe that
if it does expand our meetings will become unweildy.  

We decided on these priorities early during the 
development of PLUMED 2 as we learnt a lot of hard lessons about software management based on our experience 
with PLUMED 1.  We have thus always thought of PLUMED 2 as consisting a small set of <b>core</b> functionalities
together with various extensions.  This led us to write the code so that developers <b>can implement new 
functionalities without changing any other files within PLUMED</b>. 

We believe this distinction between core code and extensions divides the community of PLUMED developers into two
distinct groups.   There is a small group of core-developers who work to ensure that the core of the code is maintained
and a second larger group of contributors who work on new features.  We should be clear, however, that we do not
make this distinction because we believe that what we (the core developers) do is more valuable than the work 
of contributors.  The core developers are only called the core developers because they maintain the core parts of the code
that are used in every PLUMED calculation.  This maintence work is not even software development per say as
<b>the core of the code has not changed much since the publication of the original paper.</b>  In fact we consider almost every
addition we have made since publication and even some of the functionalities that were described in the original papers to 
be extensions upon the original core code.

With all the above in mind we ask that developers who wish to contribute features to PLUMED put any new cpp files in a 
separate module (see what follows for an explanation of how to create a module).  Authors of these modules can use 
different copyright information at the top of their source code files to make it clear that they (and not the core 
developers of PLUMED) wrote these functions.  Furthermore, information about the modules that have been contributed to 
PLUMED will be put here.  On this page information on the authors of the module, the relevant papers and a graphical
image to illustrate the purpose of the module will be provided for the PLUMED users to see.  

Ultimately we would like PLUMED to be a community code that serves the requirements of all the users and developers in this field.
We feel that the best model for achieving this is to have a code that is composed of a number of semi-autonomous modules with their 
own individual identities that are developed in separate research groups.  These modules should (as much as possible) make use of 
a common input/output syntax and a common interface with the various large MD codes.  Furthermore, it should be possible to use 
functionalities from different modules concurrently.  Dividing the code into a core and extensions is what will allow us to achieve 
this aim.       

\section cmodule Creating a module 

In the following sections a set of step-by-step instructions for creating a new module for incorporating some new functionality 
within PLUMED is provided.

\subsection forking Fork the PLUMED repository

Git is amazing!  If you learn to use it properly it will make your life much easier.  It is perhaps hard to get started but it is
well worth the effort and there are some excellent tutorials online e.g. https://try.github.io/levels/1/challenges/1.  The first 
things you should do when you start working on developing PLUMED is to learn a bit about git, create an account on github 
https://github.com and create your own fork of plumed on github (see \ref https://help.github.com/articles/fork-a-repo/).  
By forking plumed on github you are creating your own independent repository on a remote server.  You can change this repository
without breaking the main plumed2 repository as it is your own personal repository.  Furthermore, because it is on a remote server
it is easy for you transfer your code between your laptop and desktop computers.  

The main advantages of having your own fork are that you can merge changes to the main
plumed2 repository into your own repository so that you have all the new features added by other developers using the following
instructions (https://help.github.com/articles/syncing-a-fork/).  Furthermore, as we shall see in later section working in a fork makes it 
straightforward to merge your changes (once they are ready) into the main plumed2 repository and the release version of the code. 

\subsection cdirectory Create a directory for the module source code

Once you have your own fork of PLUMED you can begin to add your new features.  Ideally when you do so you should not need to modify
any of the cpp files that are already part of PLUMED.  In other words, all your features should be implemented in new cpp files and new
header files.  It is a good idea to bundle all these files together into a single directory as has been done for other PLUMED modules 
such as crystallisation, adjmat and metainference.  Within your cpp files it is also a good idea to put all your new code into its own
separate namespace within the main PLMD namespace.  By doing so you prevent naming conflicts with other developers.  

Notice that conflicts may still happen if you pick a name for your collective variable that collides with something
else existing in the code. In this respect, before merging your contribution, the core developers may ask you to change
the name of some of the keywords that you added.

The first step in writing your new feature will thus be to create a sub-directory within src in which to hold your new feature.  Before you 
write any c++ code you will need to create two files within this directory. The first of these files will be called module.type and will 
just contain the following text:

\verbatim
default-off
\endverbatim

By setting this file up this way you ensure that your module is not compiled unless the user specifically asks for it to be compiled at
configure time.  Explanations on how to configure PLUMED to include your module will appear in the documentation for the module automatically.
Furthermore, the procedure for compiling with your code enabled is relatively straightforward.

The second file you will need to create is the Makefile.  An example module Makefile is shown below:

\verbatim
USE=core tools vesselbase multicolvar
# generic makefile
include ../maketools/make.module
\endverbatim

You should only need to modify the first line of this file - the line starting with USE=.  This line is used to tell PLUMED at compile time
what other modules are required in order for this module to function.  This module relies on functionality that is contained in the core, tools, 
vesselbase and multicolvar modules and so these modules are all required during compilation of this particular module.  For your new module you 
will most likely always need to use core and tools.  The other modules that are required will depend on what you are implementing.

The last thing you will need to do before you start programming is that you will need to modify the .gitignore file in the src directory in order 
to stop git from ignoring your new module directory.   If you look in the .gitignore file you will see that it reads something like this:

\verbatim
/*

# Only track modules that are part of the plumed 
!/Makefile
!/README
!/adjmat
!/analysis
\endverbatim

If your new module is called newmodule then you need to add a !/newmodule in the .gitignore file in order to prevent git from ignoring the directory.

\subsection cwrite Writing your source code

Obviously, the source code you will write in your modules directory will depend on the particular feature you are implementing
and there is thus little generic advice that we can give at this stage.  We would ask that you observe a few rules, however.  
In particular:

- Please ensure that you fully document all the new features that you add to the code and that you include examples in your documentation.
There is information on how to write documentation for PLUMED on this page: \ref usingDoxygen.  Note that whenever you make a commit that 
changes the documentation for the code you need to add the string [makedoc] somewhere in the commit message to update the website.  
- We ask that you maintain the portability of plumed by only using the STL library and lapack in modifications.
If you need to use any less standard library (e.g. Boost, Sockets) please ensure that your functionality is not 
installed during a default compilation.  Instead add new options to the configure script in order to make it search for 
libraries so as to make compilation straightforward for users who are using/not using your new feature.  There is information 
on how to add complilation options on this page: \ref UsingExternalLibs.  N.B. you should only need to modify the configure
scripts if you are using modules.  Flags for activating your module at configure time will be generated automatically.
- We ask you not to include C++ features that are too new and not supported by the vast majority of compilers.
Until PLUMED v2.3, we were not using any C++11 feature. Since PLUMED v2.4, a C++11-compliant compiler is explicitly
requested and so you can use C++11 features.

\subsection ctesting Writing regression tests 

It is really important for us (the core developers) to be able to tell if a changes to the code affect the way PLUMED behaves and the answers that 
it produces.  It is thus really important that you write suitable regression tests for new features.  In fact you should be writing regtests as
you code - every time you implement a new feature write a regression tests immediately after you are done.  Testing is an important part of 
developing code and you shouldn't leave it to the end as an afterthought - doing so is horrendously bad software development practise.  

There are instructions as to how to create regression tests on this 
page: \ref regtests.  If you put all your new rt* directories in a single subdirectory in the regtest directory and if you give the directory you
create the same name as your module then that makes our life easier.  Just remember, if you do it this way, that you need to copy the Makefile 
from one of the other module folders to your new directory and that you need to modify the .gitignore file in the regtest folder to make it not
ignore the files in your new folder.  If you look in the .gitignore file you will see it reads soemthing like this:

\verbatim
/*

# Only track modules that are part of the plumed 
!/Makefile
!/README
!/adjmat
!/analysis
!/bias
\endverbatim 

If your module is called newmodule then you need to add a !/newmodule in the .gitignore file in order to prevent git from ignoring the directory.

\subsection pull-request Do a git pull request

Once you have finished coding and once you have written your regression tests submit a git pull request.  There is an explanation of how 
to do a pull request here: https://help.github.com/articles/using-pull-requests/  By doing the pull request you will make the core developers aware
that you want to add a new feature to the release version of PLUMED.  We will be able to see all the changes that you would like to make 
and any new files that you would like to add.  We may come back with a few questions or suggestions at this stage but we will eventually merge your
code into the main PLUMED repository.  Once we have merged your code it is very important that you <b> do not delete your fork of the plumed repository</b>.
If you need to make changes to your source code files you will first need to change them in your own forked repository.  Once you have made these changes
you will then need to <b> submit another pull request</b> in order to change the release version of the code.  Remember also that you are responsible for
maintaining your source code in perpetuity and that we may therefore contact you and ask you to fix something in response to some change we have made.  
In case you will not have time to do it, we will be forced to delete your module from the main PLUMED repository. People
will still be able to download it from your fork, but obviosuly your module won't benefit anymore from the further
enhancement we will add to the main PLUMED repository.
As always this fix should be done on your fork of the repository first and then merged into the release version of PLUMED using a pull request. 

