#ifndef __PLUMED_Plumed_h
#define __PLUMED_Plumed_h

/*
  Plumed.h and Plumed.c contain the external plumed interface, which is used to
  integrate it with MD engines. This interface is very general, and is expected
  not to change across plumed versions. Plumed.c also implements a dummy version
  of the interface. These files could be directly included in the official
  host MD distribution. In this manner, it will be sufficient to link the plumed
  library at link time (on all systems) or directly at runtime (on some system)
  to include plumed features.
  The C++ interface is more or less equivalent to the PlumedMain class, but it is
  based on the plumed_{create,cmd,finalize} C interface. In this manner it
  can take advantage of dynamic binding.

  The available actions are:
  create   --- creates a plumed object (in C++ this is done by the constructor)
  cmd      --- sends a command to the plumed object
  finalize --- destroys a plumed object (in C++ this is done by the destructor)

  The dummy interface does:
  create   --- don't do anything
  cmd      --- report an error and die
  finalize --- don't do anything

  There is also a "global" interface, where a static plumed object is used. It can be called
  also from FORTRAN.

  C++ interface:
    bool installed();
    class PLMD::Plumed{
      cmd(const char* key,const void* val=NULL);
    }
  C interface (can be used in C/C++)
    void   plumed_installed(int*)
    struct plumed;
    plumed plumed_create();
    void   plumed_cmd(plumed p,const char* key,const void* val);
    void   plumed_finalize(plumed p);
  Global interface (can be used in C/C++/FORTRAN)
    plumed_installed(int*)
    plumed_g_create();
    plumed_g_command(const char* key,const void* val);
    plumed_g_finalize();

  The main routine is "cmd", which accepts two arguments:
  key is a string containing the name of the command
  val is the argument. it is declared const so as to use allow passing const objects, but in practice plumed
      is going to modify val in several cases (using a const_cast).
  In some cases val can be omitted: just pass a NULL pointer (in C++, val is optional and can be omitted).
  The set of possible keys is the real API of the plumed library, and will be expanded with time.
  New commands will be added, but backward compatibility will be retained as long as possible.

  The plumed_install routine set its argument to a positive integer if plumed is installed,
  to zero otherwise. Similarly, the C++ version installed() returns true in the first case and
  false in the second. This is useful for the host MD code to know if
  plumed has been linked or not, and is the only command which can be executed in the dummy interface.

  To pass plumed a callback function use the following syntax:
    plumed_function_holder ff;
    ff.p=your_function;
    plumed_cmd(plumed,"xxxx",&ff);
  (this is passing the your_function() function to the "xxxx" command)
*/

#ifdef __cplusplus
 extern "C" {
#endif

/**
  \brief Main plumed object

  This is an object containing a Plumed instance, which should be used in
  the MD engine. It should first be initialized with plumed_create(),
  then it communicates with the MD engine using plumed_cmd(). Finally,
  before the termination, it should be deallocated with plumed_finalize().
  Its interface is very simple and general, and is expected
  not to change across plumed versions.
*/
typedef struct {
/**
  \private
  \brief Void pointer holding the real PlumedMain structure
*/
  void*p;
} plumed;

/* Generic function pointer */
typedef void (*plumed_function_pointer)(void);

/**
  \brief Holder for function pointer.

  To pass plumed a callback function use the following syntax:

    plumed_function_holder ff;

    ff.p=your_function;

    plumed_cmd(plumed,"xxxx",&ff);

  (this is passing the your_function() function to the "xxxx" command)
*/

typedef struct {
  plumed_function_pointer p;
} plumed_function_holder;

/* C interface: */

/** \relates plumed
    \brief Constructor

    \return The constructed plumed object
*/
plumed plumed_create(void);

/** \relates plumed
    \brief Tells p to execute a command

    \param p The plumed object on which command is acting
    \param key The name of the command to be executed
    \param val The argument. It is declared as const to allow calls like plumed_cmd("A","B"),
               but for some choice of key it can change the content
*/
void plumed_cmd(plumed p,const char*key,const void*val);

/** \relates plumed
    \brief Destructor

    \param p The plumed object to be deallocated
*/
void plumed_finalize(plumed p);

/** \relates plumed
    \brief Check if plumed is installed (for runtime binding)

    \param flag Is set to 1 if plumed is installed, to 0 otherwise
*/
void plumed_installed(int*flag);

/* global C interface, working on a global object */

/** \relates plumed
    \brief Constructor for the global interface.

    Equivalent to plumed_create(), but initialize a static global plumed object
*/
void plumed_g_create(void);

/** \relates plumed
    \brief Tells to the global interface to execute a command.

    Equivalent to plumed_cmd(), but skipping the plumed argument
*/
void plumed_g_cmd(const char*,const void*);

/** \relates plumed
    \brief Destructor for the global interface.

    Equivalent to plumed_finalize(), but skipping the plumed argument
*/
void plumed_g_finalize(void);

/* fortran (underscored) wrappers, only to the global interface */


/** \relates plumed
    \brief Same as plumed_g_create().

*/
void plumed_g_create_(void);

/** \relates plumed
    \brief Same as plumed_g_cmd().

*/

void plumed_g_cmd_(const char*,const void*);
/** \relates plumed
    \brief Same as plumed_g_finalize().

*/
void plumed_g_finalize_(void);

/** \relates plumed
    \brief Same as plumed_installed().

*/
void plumed_installed_(int*);

#ifdef __cplusplus
 }
#endif

#ifdef __cplusplus

/* this is to include the NULL pointer */
#include <cstdlib>

/* C++ interface is hidden in PLMD namespace (same as plumed library) */
namespace PLMD {

/**
  C++ wrapper for \link plumed
*/

class Plumed{
  plumed main;
public:
/**
  Constructor - equivalent to plumed_create()
*/
  Plumed(){main=plumed_create();};
/**
  Cmd - equivalent to plumed_cmd()
*/
  void cmd(const char*key,const void*val=NULL){plumed_cmd(main,key,val);};
/**
  Destructor - equivalent to plumed_finalize()
*/
  ~Plumed(){plumed_finalize(main);};
};

/** \relates Plumed
    \brief Check if plumed is installed (for runtime binding)

    \return true if plumed is installed, false otherwise
*/
inline
bool installed(){int i;plumed_installed(&i);if(i>0) return true;else return false;}

}

#endif


#endif
