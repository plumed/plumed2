/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2019 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_tools_Exception_h
#define __PLUMED_tools_Exception_h

#include <string>
#include <stdexcept>
#include <sstream>

namespace PLMD {

/**
\ingroup TOOLBOX
Class to deal with Plumed runtime errors.

This class and the related macros can be used to detect programming
errors. Typical cases are internal inconsistencies or errors in the plumed<->MD
interface. Mistakes made by final users (i.e. in the `plumed.dat` file)
should probably be documented in some better way (e.g. printing parts of the manual in the output).
However, also this class allows for significant information to be attached.
Let's try to make error messages as informative as possible!

\note This class has been rewritten in PLUMED 2.5. It works in a backward compatible manner,
but is much more flexible. The main novelty is that we can use insertion operators to
add arbitrary messages, as in `plumed_error()<<"check this vector "<<v;`
See below for more details.

To throw an error, just throw a c++ exception
\verbatim
  if(something_bad) throw Exception();
\endverbatim
or better add an error message to that
\verbatim
  if(something_bad) throw Exception("describe the error here");
\endverbatim

As of PLUMED 2.5 you can add multiple messages, they will just be concatenated,
but to do se you should use the insertion operator. Notice that anything that
can be formatted with an insertion operator can go to the exception, even a \ref Vector
\verbatim
  Vector v;
  if(something_bad) throw Exception()<<"problem with this "<<v;
\endverbatim
In principle you can mix the two syntax (add a message as an argument and insert others with `<<`),
however it is not very clear and should be avoided.
We only allow using arguments in parenthesis in order to keep backward compatibility.

\par Using macros

In order to provide more context, especially for debugging, it might be useful to know where the exception
originated from. The macros below add information about the exact location of the error in the file (filename, line
and, when available, function name). Macros ending in "error" unconditionally throw
the exception, whereas macros ending in "assert" first perform a conditional check
(similarly to standard assert()).
An extra `m` in the name (e.g. `plumed_merror`) indicates a macro that provides a message as its argument.
However, as of PLUMED 2.5 we should prefer adding messages using insertion operators.
\verbatim
// this is correct but not recommended. add a message please!
  plumed_assert(a>0);

// this is the old syntax (with argument).
// this syntax is basically available for backward compatibility.
  plumed_massert(a>0,"a should be larger than zero);

// this is the recommended syntax, with insertion operators.
// it allows to easily insert multiple objects
  plumed_assert(a>0)<<"a should be larger than zero. a="<<a;

// same as above, but the test is made explicitly:
  if(a<=0) plumed_error();
  if(a<=0) plumed_error("a should be larger than zero);
  if(a<=0) plumed_error()<<"a should be larger than zero. a="<<a;
\endverbatim

The additional macros
plumed_dbg_assert() and plumed_dbg_massert() are similar
to plumed_assert() and plumed_massert() respectively, but the corresponding
check is only performed when NDEBUG macro is not defined. They should
be used when the check is expensive and should be skipped in production
code. So, for instance, in the following case:
\verbatim
  plumed_dbg_assert(expensive_function(i)>0)<<"message";
\endverbatim
`expensive_function()` is not called in the production code.
Notice that the compiler should be able to completely optimize away the
whole statement including functions used to produce the message as in this example:
\verbatim
  plumed_dbg_assert(expensive_function(i)>0)<<"I did this check "<<other_expensive_function(i);
\endverbatim

Finally, notice that there is another macro available, \ref plumed_here.
In can be used in order to create an exception with information about the
line/file coordinates without trowing it. That is, the two following syntaxes
are equivalent
\verbatim
// First way, all at once
plumed_error()<<"some message";
/////////////////////////////////
// Second way, one step at a time
// Create exception
Exception e;
// Append information about line and file
e<<plumed_here;
// Append some other message
e<<"some message";
// Throw the resulting exception
throw e;
\endverbatim

Exceptions can be caught within plumed or outside of it.
E.g., in an external c++ code using PLUMED as a library, one can type
\verbatim
  try{
    plumed.cmd("setPrecision",n);
  } catch (std::exception & e) {
    printf("ee %s",e.what());
    exit(1);
  }
\endverbatim
This can be useful if an external code wants to exit in a controlled manner
(e.g. flushing files, printing the error message in a specific file, etc.)
but is anyway limited to c++ codes. Moreover,
since these errors are expected to be unrecoverable, the MD code will
usually not be able to do something more clever than exiting.

\note
We store message and stack trace in growing strings. This is in
principle not recommended, since copying the exception might fail if
copying the string throw another exception. However, this has been like
this in all previous PLUMED versions. In case it is necessary, we can replace
it later with a fixed size array placed on the stack.

*/
class Exception : public std::exception
{
/// Reported message
  std::string msg;
/// Stack trace at exception
  std::string stackString;
/// Flag to remembed if we have to write the `+++ message follows +++` string.
/// Needed so that the string appears only at the beginning of the message.
  bool note;
/// Stream used to insert objects.
/// It is not copied when the Exception is copied.
  std::stringstream stream;

public:

/// Auxiliary containing the location of the exception in the file.
/// Typically used from the macros below.
  class Location {
  public:
    const char*file;
    const unsigned line;
    const char* pretty;
    explicit Location(const char*file,unsigned line,const char* pretty=nullptr):
      file(file),
      line(line),
      pretty(pretty)
    {}
  };

/// Auxiliary containing the failed assertion.
/// Typically used from the macros below.
  class Assertion {
  public:
    const char*assertion;
    explicit Assertion(const char*assertion=nullptr):
      assertion(assertion)
    {}
  };

/// Default constructor with no message.
/// Only records the stack trace.
  Exception();

/// Constructor compatible with PLUMED <=2.4.
  explicit Exception(const std::string & msg):
    Exception()
  {
    *this << msg;
  }

/// Copy constructor.
/// Needed to make sure stream is not copied
  Exception(const Exception & e):
    msg(e.msg),
    stackString(e.stackString),
    note(e.note)
  {
  }

/// Assignment.
/// Needed to make sure stream is not copied
  Exception & operator=(const Exception & e) {
    msg=e.msg;
    stackString=e.stackString;
    note=e.note;
    stream.str("");
    return *this;
  }

/// Returns the error message.
/// In case the environment variable PLUMED_STACK_TRACE was defined
/// and equal to `yes` when the exception was raised,
/// the error message will contain the stack trace as well.
  virtual const char* what() const noexcept {return msg.c_str();}

/// Returns the stack trace.
/// Stack trace stored only if the required functions were found at configure time.
  virtual const char* stack() const noexcept {return stackString.c_str();}

/// Destructor should be defined and should not throw other exceptions
  virtual ~Exception() noexcept {}

/// Insert location.
/// Format the location properly.
  Exception& operator<<(const Location&);

/// Insert assertion.
/// Format the assertion properly
  Exception& operator<<(const Assertion&);

/// Insert string.
/// Append this string to the message.
  Exception& operator<<(const std::string&);

/// Insert anything else.
/// This allows to dump also other types (e.g. double, or even Vector).
/// Anything that can be written on a stream can go here.
  template<typename T>
  Exception& operator<<(const T & x) {
    stream<<x;
    (*this)<<stream.str();
    stream.str("");
    return *this;
  }
};

/// Class representing a generic error
class ExceptionError :
  public Exception {
public:
  using Exception::Exception;
  template<typename T>
  ExceptionError& operator<<(const T & x) {
    *((Exception*) this) <<x;
    return *this;
  }
};

/// Class representing a debug error (can only be thrown when using debug options)
class ExceptionDebug :
  public Exception {
public:
  using Exception::Exception;
  template<typename T>
  ExceptionDebug& operator<<(const T & x) {
    *((Exception*) this) <<x;
    return *this;
  }
};

#ifdef __GNUG__
// With GNU compiler, we can use __PRETTY_FUNCTION__ to get the function name
#define __PLUMED_FUNCNAME __PRETTY_FUNCTION__
#else
// Otherwise, we use the standard C++11 variable
#define __PLUMED_FUNCNAME __func__
#endif

/// \relates PLMD::Exception
/// Auxiliary macro that generates a PLMD::Exception::Location object.
/// Might be useful if we want to use derived exceptions that could
/// be thrown using `throw DerivedException()<<plumed_here<<" "<<other stuff"`.
/// It is used in the macros below to throw PLMD::Exception.
#define plumed_here PLMD::Exception::Location(__FILE__,__LINE__,__PLUMED_FUNCNAME)

/// \relates PLMD::Exception
/// Throw an exception with information about the position in the file.
/// Messages can be inserted with `plumed_error()<<"message"`.
#define plumed_error() throw PLMD::ExceptionError() << plumed_here

/// \relates PLMD::Exception
/// Throw an exception with information about the position in the file
/// and a message. Mostly available for backward compatibility
#define plumed_merror(msg) plumed_error() << msg

/// \relates PLMD::Exception
/// Launches plumed_merror only if test evaluates to false.
/// The string describing the test is also reported.
/// Further messages can be inserted with `<<`.
#define plumed_assert(test) if(!(test)) plumed_error() << PLMD::Exception::Assertion(#test)

/// \relates PLMD::Exception
/// Launches plumed_merror only if test evaluates to false.
/// The string describing the test is also reported, in addition to
/// messages reported in the extra argument. Mostly available for backward compatibility.
#define plumed_massert(test,msg) plumed_assert(test) << msg

#ifdef NDEBUG

// These are the versions used when compiling with NDEBUG flag.
// The condition is always true, so that the rest of the statement
// should be optimized away.
#define plumed_dbg_assert(test) plumed_assert(true)
#define plumed_dbg_massert(test,msg) plumed_massert(true,msg)

#else

/// \relates PLMD::Exception
/// Same as \ref plumed_assert, but only evaluates the condition if NDEBUG is not defined.
#define plumed_dbg_assert(test) if(!(test)) throw PLMD::ExceptionDebug() << plumed_here << PLMD::Exception::Assertion(#test)

/// \relates PLMD::Exception
/// Same as \ref plumed_massert, but only evaluates the condition if NDEBUG is not defined.
#define plumed_dbg_massert(test,msg) plumed_dbg_assert(test) << msg

#endif

}
#endif
