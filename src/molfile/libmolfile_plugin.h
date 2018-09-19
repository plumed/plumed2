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
#ifndef __PLUMED_molfile_libmolfile_plugin_h
#define __PLUMED_molfile_libmolfile_plugin_h
#include "vmdplugin.h"
namespace PLMD{
namespace molfile{


 int molfile_dcdplugin_init(void);
 int molfile_dcdplugin_register(void *, vmdplugin_register_cb);
 int molfile_dcdplugin_fini(void);
 int molfile_crdplugin_init(void);
 int molfile_crdplugin_register(void *, vmdplugin_register_cb);
 int molfile_crdplugin_fini(void);
 int molfile_gromacsplugin_init(void);
 int molfile_gromacsplugin_register(void *, vmdplugin_register_cb);
 int molfile_gromacsplugin_fini(void);
 int molfile_pdbplugin_init(void);
 int molfile_pdbplugin_register(void *, vmdplugin_register_cb);
 int molfile_pdbplugin_fini(void);

#define MOLFILE_INIT_ALL \
    molfile_dcdplugin_init(); \
    molfile_crdplugin_init(); \
    molfile_gromacsplugin_init(); \
    molfile_pdbplugin_init(); \

#define MOLFILE_REGISTER_ALL(v, cb) \
    molfile_dcdplugin_register(v, cb); \
    molfile_crdplugin_register(v, cb); \
    molfile_gromacsplugin_register(v, cb); \
    molfile_pdbplugin_register(v, cb); \

#define MOLFILE_FINI_ALL \
    molfile_dcdplugin_fini(); \
    molfile_crdplugin_fini(); \
    molfile_gromacsplugin_fini(); \
    molfile_pdbplugin_fini(); \

}
}
#endif
