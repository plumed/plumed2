\page ABriefIntroduction A brief introduction to the plumed core

Plumed 2, unlike its predecessor plumed 1, which was written in plain C, is written in C++.  
C++, unlike C, fortran and many of the other languages that are commonly used in the 
atomistic simulation community, is an object oriented programming language.  As such the way things
are done in plumed may feel unfamiliar to developers in the scientific community who,
if our experience is anything to go by, are more used to programming in non-object
oriented languages.  For this reason we have tried in what follows to explain how
we have used the features of object oriented programming in plumed 2.  We
hope that this guide is helpful and appologize in advance to any developers
who feel patronized.     

\section intro Object oriented programming

The main objective in object oriented programming is to write code that is more 
resilient to bugs.  There are two ways that object oriented programing allows
us to acchieve these aims:

- In object oriented programs one generally needs fewer lines of code
- Object oriented programming allows us to use the compiler to do many more checks of the code 

To be clear though object oriented programming does not allow us to do things that 
would be impossible with other programming languages.  All programs perform some set of 
mathematical operations and any programming language can be used to implement these mathematical 
operations. The only advantage of C++ is that the advanced, object-oriented features of the 
language make implementing things more straighforward. 

As you are no doubt aware, in C one can create structures to order the variables in your code.  
A naive way of thinking about the objects, or more correctly classes, that one uses in C++ is
that these are structures that also contain functions.  This is useful for making neat header files
as the parameters are kept near the functions.  However, at this level of thinking the C++ way 
of doing things:

\verbatim
class fclass {
bool param0;
int param1,param2;
double param3;
double f( double x );
}; 
\endverbatim

is not much better than the C way of doing things:

\verbatim
struct fparams {
bool param0;
int param1,param2;
double param3;
}; 
double f( double x, struct fparams myparams ); 
\endverbatim

<b>
Nevertheless for reasons that will hopefully become clear as you read this document every bias, 
colvar and function that is implemented in plumed 2 is inside its own separate class.
</b>

\section publicprivate Public and private members

The variables in a C struct can be accessed anywhere in the code - any function in the code
can copy the information from a structure or change the values of the variables in the structure.  This
was a particularly fun feature of plumed 1.0, every function in the old code could change any of the variables
in the code!  Obviously this causes problems as new functions can accidentally change the 
values of some variable in a widely used structure that should never have been changed.  As one can imagine this can cause 
drastic problems.  To prevent errors like this C++ provides a set of functionalities to 
allow one to specify what any given function can do to the members of a class.  This is not possible in 
C and it was the ability to use this functionality to create flexible, easily-extendable code that motivated
our choice of C++ for plumed 2.  The example class below shows how this is done in practise:

\verbatim
class myclass{
private:
  bool b;  //--This can only be accessed when you are in one of the methods in the class
public:
  int i;  //--This can be acessed by anywhere in the code
};

\section constructors Constructors

As someone who learnt to program in fortran it was this aspect of C++, more than any other, that confused me 
the most. In actually though it is rather simple and I don't really know why I was confused. In essence every 
class must contain some method for creating it. This class should set the initial values of all the variables
inside the class. Obviously the functions (or more correctly methods) that the class
contains cannot be used until an instance of the class has been created using the constructor.

An example of how all this works in practise is given below:

\verbatim
class myclass{
private:
  bool b;
  int i;
  double d;
public:
  myclass( bool bb, int ii, double dd ) : b(bb), i(ii), d(dd) {}
  static double g(int j);
  double f( double x );
};

// Here there are currently no instances of myclass and so I cannot run f
// I can run g however as it is a static member of the class - I run it using
double d=myclass::g(j);

// Now I create an instance of the class using the constructor
myclass thisIsTheInstance(false, 3, 6.5);
// so that I can run the method called f 
double s=thisIsTheInstance.f(4.0); 
\endverbatim         

<b>
In plumed 2 all the lines in the input file are read in inside the constructors.  This ensures
that the parameters inside any given method are set correctly from the outset.
</b> 

\section operators Operators

Addition, subtraction, multiplication, division and so on are all functions (they are obviously
not variables).  We usually don't think of them as functions however because we use 
these operations all the time.  C++ recognizes that the short cuts of +, -, *, / and so on
 are very useful.  It thus allows one to define operators in our new classes that explain
to the compiler what a given symbol means for a given class.  Among other things we can define:

- How to perform the logical operators !=, ==, etc
- How to perform arithmatic for a class: +, -, /, *, +=, -= etc
- What brackets mean i.e. the meanings of (), [] 

We do not use this extensively in plumed 2 but it does occasionally appear.

\section inclusion Including the functionality of one class in a new class 1: Inclusion

There are various ways that one can include the functionality of one class inside a second class.
By far the simplest is to create an instance of class 1 inside class 2 as shown below:

\verbatim
class class1 {
private:
  double d1,d2,d3;
public:
  class1();
  void f(double x);
};

class class2 {
private:
  class1 myclass;
public:
  class2();
  // The methods of class 2 here
};
\endverbatim

This is really simple one includes a class in this way in exactly the same way that one includes a 
double, int or whatever variable.

<b>
This kind of inclusion is used extensively in plumed 1.0 and there are a huge number of classes that you
can re-use in this way to create new colvars, functions or biases.  For a full list of the classes that are 
available see \subpage TOOLBOX.
</b>

\section inheritance Including the functionality of one class in a second class 2: Inheritance

There is an alternate way of reusing the functionality from one class in a second class that is available
in C++.  This method is called inheritance and it offers some advantages over simply including class A
inside class B as described above.  To create a class called B that inherits from A one writes:

\verbatim
class B : public A {
// Contents of class B
};  
\endverbatim

One advantage of this method over inclusion is that I can use protected members to more closely control
what members of A class B is allowed to access.  Hence, rather than simply having private and public members
I now have:

<table align=center frame=void width=95%% cellpadding=5%%>
<tr>
<td> <b> public </b> </td> <td> These members can be accessed by anyone </td>
</tr> <tr> 
<td> <b> protected </b> </td> <td> These members can only be accessed by the methods of class A and class B </td>
</tr> <tr>
<td> <b> private </b> </td> <td> These members can only by accessed by the methods of class A (not by class B) </td>
</tr>
</table>

In addition, I can use inheritance to treat pointers to objects of class B as if they were pointers to objects of 
class A.  In other words, if I create an object of type B I can convert it to an object of type A using dynamic_cast
as shown below:

\verbatim
B* mynewB=new B();   // This is a special way of calling the constructor so you get a pointer
A* BpretendingToBeA=dynamic_cast<A*>(mynewB); 
\endverbatim

All the colvars and free energy methods of plumed use inheritance. In fact all these methods are
built on a single base class called PLMD::Action. This class contains all the functionality for
reading stuff from input, the stuff for controlling the dependencies Actions and a set of controls
that decide which actions are performed when. All the functionality for the different methods is
then built on this root. As you can see (PLMD::Action) the inheritance tree for the code is
quite complicated.  However, in practise if you are implementing a CV, function, bias or virtual atom
the correct start point is with one of the classes listed on this page \ref INHERIT all of which contain
detailed descriptions of how to use them. 

\section minheritance Including the functionality of one class in a second class 3: Multiple inheritance

Immediately above the PLMD::Action root of the inheritance tree in plumed there is a very complicated looking
layer in the inheritance structure of the code.  This layer looks ugly because in this layer
we are using multiple inheritance - the classes in the layer above inherit from multiple classes simultaneously.
This way of incorporating functionality from classes is unique to C++ and brings with it a special set of
 difficulties in programming.  Its great advantage is though that one can 
create classes that incorporate bring a set of particular attributes.  This will perhaps be most clear
if you look at what each of the classes in the multiple inheritance layer is doing (see \ref MULTIINHERIT) and see how these
functionalities are used in Colvars, Functions and Biases.  Please be aware that, unless you are doing something really wacky, 
you should be able to implement whatever you need to implement without writing classes that take advantage of multiple inheritance. Furthermore,
you should not need to touch the classes in this region of the code.  The information here is there for
completeness only.  If you feel you really must change something in this part of the code please contact the
developers before doing anything.     

\section static-poly Static Polymorphism

Polymorhpism is a way of using the same code to do many different things.  As an example consider a 
Matrix.  The elements of a Matrix can be ints, doubles, floats or even some fancy new class
but we would still want the operator (i,j) to return the element in row i and column j.  That is to
say the operator (const int i, const int j) of a matrix is independent of what is actually inside the
matrix.  Using C++ we can use so called template classes to implement thse kinds of things and can then 
re-use them to do an enormous variety of different operations.  To see how this works in practise take
a look at PLMD::Matrix, which is a working version of our Matrix example.  Be aware that all the 
routines in a template class must be inside the header file.  To use a template within the code
you declare it as follows:

\verbatim
Matrix<double> mat;   // This is a matrix of doubles
\endverbatim

<b>
The most common way we use this kind of functionality in plumed 2 is when we take advantage of the
features that are available in the C++ standard library.  For more details on the standard library
visit http://www.cplusplus.com/reference/
</b>  

\section dynamic-poly Dynamic Polymorhpism

When you run a calculation with plumed the code calculates a number of CVs.  The bias and the forces
due to the bias are then calculated and in the final step these forces are propegated back onto the 
atoms using the chain rule.  For example PLMD::colvar::Distance
contains the function that calculates a distance between atoms, while PLMD::bias::MetaD contains
the function for doing metadynamics.  What may thus seem remarkable to the programmer unfamiliar with
C++ is that the class that calls the functions that calculate the CVs, biases and so on only uses 
PLMD::Action.  To make that clear it looks like the code can calculate the distances between
atoms without ever calling any of the routines from PLMD::colvar::Distance!      

We can program in this way because we take advantage of dynamic polymorhpism.  If you look at the
documenation for PLMD::Action you will see that the method PLMD::Action::calculate is declare inside
PLMD::Action as:

\verbatim
virtual void calculate()=0;
\endverbatim

This kind of declaration promises two things to a class:

- That the class will only ever be used in derived classes.  No PLMD::Action class is ever 
constructed in the code. The functionality in PLMD::Action is only ever used in the derived classes 
that inherit PLMD::Action. Classes like PLMD::Action are called abstract base classes.
- That in one of the classes that inherits from PLMD::Action a method called calculate will be defined.

The great advantage of declaring calculate() inside PLMD::Action in this way is that the calculate 
routine that we declare in the derived class is a member of PLMD::Action.  We thus can thus write
a class for doing all the business of plumed in the manner described previously.

\section ForwardDeclaration Forward declaration

One problem of including classes inside other classes in C++ is that this enforces one
to include one .h file into another one, thus leading to a large set of objects
needing to be recompiled just because a single .h file was touched. In some cases
this is not avoidable, e.g. when classes inherits from other classes. However,
when only a pointer (or a reference) to another class is used, it might be better
to just use a forward declaration as in this example:
\verbatim
/////////////////////////////////////////////
// This is file A.h
namespace PLMD{

class A{
  int pippo;
};

}

/////////////////////////////////////////////
// This is file B-bad.h
// it has to include A.h
#include "A.h"
namespace PLMD{

class B{
public:
// notice that here we only use a reference to class A
  int do_something(A&a);
};

}

/////////////////////////////////////////////
// This is file B-good.h
namespace PLMD{

// this command just instructs the compiler that A is a class:
class A;
// no inclusion of A.h is required!

class B{
public:
// notice that here we only use a reference to class A
  int do_something(A&a);
};

}

\endverbatim

This trick however does not work is a class is including an instance of another class.
E.g., if B _contains_ an instance of A one should know exactly the A declaration
to build a B object.
In this case, a similar effect can be obtained at the price of adding some
more lines of code in constructor and destructor of B as in the following example

\verbatim

/////////////////////////////////////////////
// This is file B-bad.h
// it has to include A.h
#include "A.h"
namespace PLMD{

class B{
  A content;
};

}

/////////////////////////////////////////////
// This is file B-good.h
namespace PLMD{

// this command just instructs the compiler that A is a class:
class A;
// no inclusion of A.h is required!

class B{
// As "content" is a reference, it can be used exactly as if it was a normal object
// However, it is represents by the compiler as a pointer.
  A& content;
public:
  B();
  ~B();
};

}

/////////////////////////////////////////////
// Using B-good.h enforces to add something in B-good.cpp

#include "A.h"
#include "B-good.h"

using namespace PLMD;

B::B():
// now "content" needs to be explicitly allocated ...
  content(*new A){
}

B::~B(){
// ... and deallocated
  delete &content;
}

\endverbatim

Notice that this trick cannot be always be applied, e.g., if the constructor to be A needs parameter,
or if object "content" is to be accessed by inline methods of B for efficiency. Another example where
this does not work is when inline methods are used because of template expressions.

This trick is extensively used in plumed so as to avoid too many indirect dependencies
among .h files.

\section cxx11features C++11 Features

Since PLUMED 2.4 we systematically use C++11 features. Some of the most important ones are discussed here.

\subsection cxx11features-auto Using auto

Auto allows to implicitly declare the type of a variable.
This is particularly handy when using iterators from STL containers.
For instance, you can replace the following:
\verbatim
   map<string,vector<AtomNumber> >::const_iterator m=atoms.groups.find(strings[i]);
\endverbatim
with:
\verbatim
   const auto m=atoms.groups.find(strings[i]);
\endverbatim
Notice that the syntax is significantly simpler, especially if you do not
remember which exact type of map is the variable `groups`.

Iterators are often used in loops. Thus, you can now replace
\verbatim
  for(std::map<AtomNumber,Tensor>::const_iterator p=gradients.begin();p!=gradients.end();++p){
    a+=(*p).second;
  }
\endverbatim
with
\verbatim
  for(auto p=gradients.begin();p!=gradients.end();++p){
    a+=(*p).second;
  }
\endverbatim
However, in cases where you do not need to explicitly use `p` as an iterator, you might find
even more convenient to use range-based loops:
\verbatim
  for(const auto & p : gradients){
    a+=p.second;
  }
\endverbatim
Notice that now `p` is a constant reference, so it is not anymore necessary to use the `*` operator.

\subsection cxx11features-smart-pointers Using smart pointers

There are many resources on the web about this topic. Have a look at <a href="https://mbevin.wordpress.com/2012/11/18/smart-pointers/"> this link </a> for a
concise introduction.

Smart pointers can be most of the time used in place of regular pointers so as to better manage memory.
In PLUMED you can find many times sections of code such as
\verbatim
  object* obj;
  if(...) obj=new type1;
  else obj=new type2;

  ...

  obj->method();

  ..

  delete obj;
\endverbatim
Here we use a pointer to allow dynamic polymorphism.

In this case, the object pointed by `obj` is not transferred anywhere else.
In other words, you can think that the `obj` pointer owns the object itself.
You can replace it with a `std::unique_ptr` as follows:
\verbatim
  std::unique_ptr<object> obj;
  if(...) obj.reset(new type1);
  else obj.reset(new type2);

  ...
  
  obj->method();

  ..
\endverbatim

Notice that instead of assigning it with `=` you should assign it with `reset()`. This is because
the `std::unique_ptr` cannot be copied and so does not understand the assignment operator.
More importantly, notice that the delete command has disappeared. Indeed, when `obj` goes
out of scope, the pointee is automatically deleted.

You can also use vectors of pointers. Consider the following example
\verbatim
  std::vector<object*> objs;
  for(unsigned i=0;i<10;i++) objs.push_back(new object);

  ...

  for(unsigned i=0;i<10;i++) delete objs[i];
\endverbatim
This can be replaced with
\verbatim
  std::vector<std::unique_ptr<object>> objs;
  for(unsigned i=0;i<10;i++) objs.emplace_back(new object);

  ...
\endverbatim

Notice that instead of using `push_back()` we used `emplace_back()`. The reason is that the latter move
the pointer instead of copying it. More importantly, notice that the delete command has disappeared.
When the vector is cleared, also the contained objects are deleted.

Notice that `emplace_back` needs to be fed with a so-called rvalue. In case you created the
unique_ptr in advance, you should insert it with the following syntax
\verbatim
  std::vector<std::unique_ptr<object>> objs;
  std::unique_ptr<object> to_insert;
  to_insert.reset(new object);
  objs.emplace_back(std::move(to_insert));
\endverbatim

\subsection cxx11features-forward Forward declarations using C++11

Notice that also forward declarations discussed above are a bit simpler
to implement using C++11 syntax. This can be done using a std::unique_ptr
or, even better, using the class ForwardDecl, which is a small utility class that
only implement two methods:
- A constructor, which takes an arbitrary number of parameters and use them
  to construct an internally stored `std::unique_ptr`.
- A `operator *`, which returns a pointer to the object.

An example usage is below:


\verbatim
/////////////////////////////////////////////
// This is file B-good.h
#include "tools/ForwardDecl.h"

namespace PLMD{

// this command just instructs the compiler that A is a class:
class A;
// no inclusion of A.h is required!

class B{
  ForwardDecl<A> content_fwd;
  A& content1=*content1_fwd;
  ForwardDecl<A> content_fwd;
  A& content2=*content2_fwd;
public:
  B();
  ~B();
};

}

/////////////////////////////////////////////
// Using B-good.h enforces to add something in B-good.cpp

#include "A.h"
#include "B-good.h"

using namespace PLMD;

B::B():
// constructors that need no argument can be omitted
  content2_fwd(argument)
{
}

B::~B(){
// empty destructor
}

\endverbatim

Notice that it is necessary to add a destructor, even though it is empty.
The reason is that if the compiler tries to construct an inline destructor for this class
it will not be able to create it (the class is not completely defined in `B.h`.
However, the advantage is that objects are deallocated in the correct order as if they were
normal members of class B, that is the inverse of the initialization order.

\section conc Conclusion

The above is meant to give you some feel as to how plumed works.  If there is stuff you do not
understand it is not necessarily that important.  The great advantage of the code as it is currently
written is that you can implement new methods without touching the core of the code.  So in conclusion
give it a go and have fun!

\section Notes Notes

More information about C++
http://www.parashift.com/c++-faq-lite/

\page UsingExternalLibs Using external libraries in PLUMED

When implementing new CVs and methods you may at time find it useful to make use functionality that has
been implemented in some external library (eg <a href="http://www.boost.org"> boost </a>).  This is strongly 
encouraged and good practise - you introduce fewer bugs if you write less code.  Having said that we would 
prefer it if any functionality that is reliant on external libraries is <b> not </b> enabled by default.  The 
reason for this is that we would like to ensure that the code remains easy to compile.  We would rather not 
have users struggling to resolve lots of dependencies on external libraries just so that they can compile 
PLUMED to run metadynamics with an distance and a torsion.

The first step in ensuring that the code will compile even if your fancy library is not available is to thus 
put all the calls to the functions in these libraries inside an ifdef block.  So for example here is a block
of code that uses the boost graph library:

\verbatim
#ifdef __PLUMED_HAS_BOOST_GRAPH
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/graph_utility.hpp>
#endif
\endverbatim

Here the `#ifdef __PLUMED_HAS_BOOST_GRAPH` ensures that this part of the code is only compiled when the 
`-D__PLUMED_HAS_BOOST_GRAPH` flag is set at compile time.  

This example is a bit of a special case as the particular PLMD::Action I took this example from both when the 
library is available and when it is not.  Obviously, if your PLMD::Action does not work without the library
then you should place everything in the file inside the file (the definition of the class, the definitions of 
the methods and the PLUMED_REGISTER_ACTION command) an ifdef block like the one shown above. 

Once you have written the code and surrounded all the calls to the library by an ifdef block like the one above
you then need to adjust the configure scripts so that the user is provided with the option to compile the code
with your new functionality and the link to the external library available.  To do this you need to edit the
file configure.ac in the plumed2 directory.  The first edit you need to make to this file should read as follows:

\verbatim
PLUMED_CONFIG_ENABLE([library_name],[search for library_name],[no])
\endverbatim

Add this in the part of the file where there are similar names.  This command allows users to enable the package
using `--enable-library_name` when they run the configure command.  Obviously, library_name here should be replaced
with the name of the particular library you are using. A shell variable named `library_name` will be set
to `true` if the library has been requested. The last argument says that, by default, the library is not requested.
If you replace it with a `[yes]`, then you will be able to use `--disable--library_name` to disable the library.

The second edit you need to make to the configure.ac file is as follows:

\verbatim
if test $library_name == true ; then
  PLUMED_CHECK_PACKAGE([header_file.hpp],[function],[__PLUMED_HAS_LIBRARY_NAME],[library_name])
fi
\endverbatim 

This command checks if the library is available on the machine that PLUMED is being compiled on.  If it is the 
Makefiles are generated with the appropriate precompiler directives and lists of libraries to link.  If it is not
available a warning is issued to the user that tells him/her that PLUMED will not be compiled with the requested
features.  The PLUMED_CHECK_PACKAGE command that we use to do this check here takes in four arguments.  The first
is the name of one of the header files from the library that you are using.  You specify the location of this header
file in the same way as you specified its location within the code.  So for example if you had #include <boost/graph/adjacency_list.hpp> 
in your cpp code you would replace header_file.hpp in the above with boost/graph/adjacency_list.hpp.  The next 
argument to PLUMED_CHECK_PACKAGE is one of the functions that is within the library that you are trying to link.  If in doubt
you can use the exit function here, which should work even if you are using template functions.  The third argument is the precompiler
directive that appears around your library calling code and that we discussed earlier.  The last argument meanwhile is the name of the library - 
the part that appears after the -l.  In the example above the code would thus try and link -llibrary_name. 

Notice that in the for C++ libraries sometimes one has to check something more complicated than the presence of a single function.
In this case you can use a more advanced version of the command which allows you to write a short test. Look at this example:
\verbatim
if test $boost_serialization == true ; then
  PLUMED_CHECK_CXX_PACKAGE([boost serialization],[
#include <fstream>
#include <boost/archive/text_oarchive.hpp>
int main() {
    std::ofstream ofs("filename");
    boost::archive::text_oarchive oa(ofs);
    return 0;
}
  ], [__PLUMED_HAS_BOOST_SERIALIZATION],[boost_serialization boost_serialization-mt])
fi
\endverbatim

The first argument here (`[boost serialization]`) is just used in the configure log. The second argument is a
complete C++ program that should compile and link. The last arguments are similar to those of the `PLUMED_CHECK_PACKAGE` macro.

Notice that for both macros (`PLUMED_CHECK_PACKAGE` and `PLUMED_CHECK_CXX_PACKAGE`) you can omit the final library, which
means that the function should be found without adding any extra `-l` option. In addition, you can put multiple libraries.
In this case, autoconf will scan them to find the appropriate one.

Once you have made these edits issue the command:

\verbatim
autoconf
\endverbatim

and the necessary changes will be made to the configure script.  <b>You should never edit the configure script directly</b>

