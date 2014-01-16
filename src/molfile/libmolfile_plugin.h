#ifndef LIBMOLFILE_PLUGIN_H
#define LIBMOLFILE_PLUGIN_H
#include "vmdplugin.h"

#ifdef __cplusplus
extern "C" {
#endif

extern int molfile_dcdplugin_init(void);
extern int molfile_dcdplugin_register(void *, vmdplugin_register_cb);
extern int molfile_dcdplugin_fini(void);
extern int molfile_gromacsplugin_init(void);
extern int molfile_gromacsplugin_register(void *, vmdplugin_register_cb);
extern int molfile_gromacsplugin_fini(void);

#define MOLFILE_INIT_ALL \
    molfile_dcdplugin_init(); \
    molfile_gromacsplugin_init(); \

#define MOLFILE_REGISTER_ALL(v, cb) \
    molfile_dcdplugin_register(v, cb); \
    molfile_gromacsplugin_register(v, cb); \

#define MOLFILE_FINI_ALL \
    molfile_dcdplugin_fini(); \
    molfile_gromacsplugin_fini(); \

#ifdef __cplusplus
}
#endif
#endif
