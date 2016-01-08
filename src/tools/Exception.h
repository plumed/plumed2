/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2016 The plumed team
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

namespace PLMD{

/**
\ingroup TOOLBOX
Class to deal with Plumed runtime errors.

This class, or better the related macros, can be used to detect programming
errors. Typical cases are internal inconsistencies or errors in the plumed<->MD
interface. Mistakes made by final users (i.e. in the plumed.dat file)
should be documented in some better way (e.g. printing parts of the manual in the output).

To throw an error, just throw a c++ exception
\verbatim
  if(something_bad) throw Exception();
\endverbatim
or better add an error message to that
\verbatim
  if(something_bad) throw Exception("describe the error here);
\endverbatim

Even better, you can use the predefined macros 
plumed_error(), plumed_assert(), plumed_merror() and plumed_massert(),
which add information about the exact location of the error in the file (filename, line
and, for g++, function name). Macros ending in "error" unconditionally throw
the exception, whereas macros ending in "assert" first perform a conditional check
(similarly to standard assert()). The extra "m" in the name means that an
extensive error message can be added.
\verbatim
// examples:
  plumed_assert(a>0);
  plumed_massert(a>0,"a should be larger than zero");
  if(a<=0) plumed_error();
  if(a<=0) plumed_merror("a should be larger than zero");
\endverbatim
The additional macros
plumed_dbg_assert() and plumed_dbg_massert() are similar
to plumed_assert() and plumed_massert() respectively, but the corresponding
check is only performed when NDEBUG macro is not defined. They should
be used when the check is expensive and should be skipped in production
code.

By default, execution is terminated imediately and a message is printed on stderr.

If PLUMED is compiled with -D__PLUMED_EXCEPTIONS execution will continue
and the exception will be passed to c++, so that it will be possible
to intercepted it at a higher level, even outside plumed.
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
In the future we might decide to extend the usage of exceptions to
detect recoverable errors (e.g. optional arguments not found, etc),
even if I am not fully convinced that this is a good idea.
Notice that sometime people claim that code compiled with exception enabled
is slower (GB)
*/
class Exception : public std::exception
{
  std::string msg;
/// Common tool, invoked by all the constructor to build the message string
  static std::string format(const std::string&,const std::string&,unsigned,const std::string&);
/// Method which aborts in case exceptions are disabled
  void abortIfExceptionsAreDisabled();
public:
/// Without message
  Exception();
/// With message
  explicit Exception(const std::string&);
/// With message plus file, line and function (meant to be used through a preprocessor macro)
  Exception(const std::string&,const std::string&,unsigned,const std::string&);
/// Returns the error message
  virtual const char* what() const throw(){return msg.c_str();}
/// Destructor should be defined and should not throw other exceptions
  virtual ~Exception() throw(){}
};

// With GNU compiler, we can use __PRETTY_FUNCTION__ to get the function name
#if !defined(__PRETTY_FUNCTION__)
#define __PRETTY_FUNCTION__ ""
#endif

/// \relates PLMD::Exception
/// Just print file/line/function information and exit
#define plumed_error() throw PLMD::Exception("",__FILE__,__LINE__,__PRETTY_FUNCTION__)
/// \relates PLMD::Exception
/// Print file/line/function information plus msg and exit
#define plumed_merror(msg) throw PLMD::Exception(msg,__FILE__,__LINE__,__PRETTY_FUNCTION__)
/// \relates PLMD::Exception
/// Conditionally print file/line/function information and exit
#define plumed_assert(test) if(!(test)) throw PLMD::Exception("assertion failed " #test,__FILE__,__LINE__,__PRETTY_FUNCTION__)
/// \relates PLMD::Exception
/// Conditionally print file/line/function information plus msg and exit
#define plumed_massert(test,msg) if(!(test)) throw PLMD::Exception("assertion failed " #test ", " msg,__FILE__,__LINE__,__PRETTY_FUNCTION__)

#ifdef NDEBUG
#define plumed_dbg_assert(test) 
#define plumed_dbg_massert(test,msg)
#else
/// \relates PLMD::Exception
/// Conditionally print file/line/function information and exit when NDEBUG flag is not present
#define plumed_dbg_assert(test) if(!(test)) throw PLMD::Exception("assertion failed " #test,__FILE__,__LINE__,__PRETTY_FUNCTION__)
/// \relates PLMD::Exception
/// Conditionally print file/line/function information plus msg and exit when NDEBUG flag is not present
#define plumed_dbg_massert(test,msg) if(!(test)) throw PLMD::Exception("assertion failed " #test ", " msg,__FILE__,__LINE__,__PRETTY_FUNCTION__)
#endif

}
#endif
