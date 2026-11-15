#ifndef GMX_MDRUN_PLUMEDEXCHANGE_H
#define GMX_MDRUN_PLUMEDEXCHANGE_H

#include "config.h"

#include "gromacs/applied_forces/plumed/plumedforceprovider.h"

namespace gmx
{
//! Find this simulation/rank's optional provider; ownership stays with its MD module.
inline PlumedForceProvider* findPlumedForceProvider(const ForceProviders* providers)
{
    return providers ? providers->plumedForceProvider() : nullptr;
}
} // namespace gmx

#if GMX_LIB_MPI
#    include <cstdint>
#    include <cstdio>

#    include "gromacs/utility/basedefinitions.h"
#    include "gromacs/utility/real.h"

struct t_commrec;
struct gmx_multisim_t;
struct gmx_repl_ex;
class t_state;

namespace gmx
{
//! Select this step's HREX partner and broadcast its replica index within each replica.
//! Called by every rank of every replica; returns the local replica index when unpaired.
int plumed_hrex_partner(FILE* fplog, const t_commrec* cr, const gmx_multisim_t* ms, gmx_repl_ex* re, int64_t step);

//! Cache this rank's physical potential contribution at the foreign configuration.
//! Uses the local Hamiltonian and GROMACS units; GREX sums collectively within the replica.
void plumed_hrex_cache_foreign(PlumedForceProvider* provider, real* energy);

//! Exchange global coordinates, velocities, box and thermostat/barostat state with a partner.
//! Called on replica main ranks; domain redistribution is the caller's responsibility.
void exchange_state(const gmx_multisim_t* ms, int partner, t_state* state);

//! Native serial-state handoff: src and dest currently must alias; otherwise this asserts.
void copy_state_serial(const t_state* src, t_state* dest);
} // namespace gmx
#endif // GMX_LIB_MPI
#endif
