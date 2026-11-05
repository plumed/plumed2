/* Default storage capacities; override consistently with compiler -D flags. */
#ifndef PLUMED_MPIW_CAPACITY_H
#define PLUMED_MPIW_CAPACITY_H

#ifndef PLUMED_MPIW_HANDLE_CAPACITY
#define PLUMED_MPIW_HANDLE_CAPACITY ((sizeof(void *) > 8U) ? sizeof(void *) : 8U)
#endif

#ifndef PLUMED_MPIW_STATUS_CAPACITY
#define PLUMED_MPIW_STATUS_CAPACITY 128
#endif

#endif
