/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
University of Illinois Open Source License
Copyright 2003 Theoretical and Computational Biophysics Group, 
All rights reserved.

Developed by:		Theoretical and Computational Biophysics Group
			University of Illinois at Urbana-Champaign
			http://www.ks.uiuc.edu/

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the Software), to deal with 
the Software without restriction, including without limitation the rights to 
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
of the Software, and to permit persons to whom the Software is furnished to 
do so, subject to the following conditions:

Redistributions of source code must retain the above copyright notice, 
this list of conditions and the following disclaimers.

Redistributions in binary form must reproduce the above copyright notice, 
this list of conditions and the following disclaimers in the documentation 
and/or other materials provided with the distribution.

Neither the names of Theoretical and Computational Biophysics Group, 
University of Illinois at Urbana-Champaign, nor the names of its contributors 
may be used to endorse or promote products derived from this Software without 
specific prior written permission.

THE SOFTWARE IS PROVIDED AS IS, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL 
THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR 
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
OTHER DEALINGS WITH THE SOFTWARE.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_molfile_largefiles_h
#define __PLUMED_molfile_largefiles_h
namespace PLMD{
namespace molfile{
/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2016 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: largefiles.h,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.4 $       $Date: 2016/11/28 05:01:54 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   Platform dependent defines for enabling 64-bit file I/O on 32-bit machines
 *
 ***************************************************************************/
 
#if defined(_AIX)
/* Define to enable large file extensions on AIX */
#define _LARGE_FILE
#define _LARGE_FILES
#else
/* Defines which enable LFS I/O interfaces for large (>2GB) file support
 * on 32-bit machines.  These must be defined before inclusion of any
 * system headers.
 */
#ifndef _LARGEFILE_SOURCE
#define _LARGEFILE_SOURCE
#endif
#define _FILE_OFFSET_BITS 64
#endif
}
}
#endif
