/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Declares Plumed force provider class
 *
 * \author Daniele Rapetti <drapetti@sissa.it>
 * \ingroup module_applied_forces
 */
#ifndef GMX_APPLIED_FORCES_PLUMEDFORCEPROVIDER_H
#define GMX_APPLIED_FORCES_PLUMEDFORCEPROVIDER_H

#include "config.h"

#include <memory>
#include <optional>
#include <vector>

#include "gromacs/mdtypes/iforceprovider.h"

struct gmx_multisim_t;

namespace PLMD
{
class Plumed;
}
namespace gmx
{
template<typename>
class ArrayRef;
struct PlumedOptions;
// Owned separately by each simulation/rank, including thread-MPI ranks.
// Registered buffers stay valid for GREX and lambda trials.
struct PlumedContext
{
    //! Whether this MD step attempts replica exchange.
    bool exchangeStep = false;
    //! Engine eligibility for pre-update trials; shells and native MTS retain legacy ordering.
    bool preUpdateExchange = true;
    //! Whether general Hamiltonian replica exchange was requested by -hrex.
    bool hrex = false;
    //! Foreign HREX pass: suppress normal callbacks, or evaluate only its bias with API14.
    bool hrexTrial = false;
    //! Current or trial alchemical state index exposed as EXTRACV NAME=lambda.
    real lambdaState = 0;
    //! Force-output buffer registered for lambda; discrete-state MCMC uses bias energies instead.
    real lambdaForce = 0;
    //! Physical potential for an MCMC candidate, in GROMACS units: global on main rank, zero on workers.
    real trialEnergy = 0;
    //! Physical potential differences already included in native lambda acceptance, indexed by replica.
    std::vector<real> nativePhysicalDeltas;
    //! Selected partner replica, shared by the early trials and later acceptance.
    int partner = 0;
    //! Partner's native lambda index, or -1 for ordinary/Hamiltonian exchange.
    int partnerLambda = -1;
    //! Whether the pre-update foreign bias has been evaluated on this step.
    bool foreignReady = false;
    //! Whether the complete pre-update correction is ready for acceptance.
    bool exchangeReady = false;
    //! Foreign physical potential, summed over the evaluating replica, in MD units.
    real foreignEnergy = 0;
    //! Foreign bias before the ordinary local update, in MD units.
    real foreignBias = 0;
    //! Bias difference, plus physical potential difference only for -hrex.
    real exchangeDelta = 0;
    //! Whether this is the replica's main MPI rank; used to supply global ENERGY only once.
    bool mainRank = true;
};

/*! \internal \brief
 * Implements IForceProvider for PLUMED.
 */
class PlumedForceProvider final : public IForceProvider
{
public:
    /*! \brief Initialize the PLUMED interface with the given options
     *
     * \param options PLUMED options
     */
    PlumedForceProvider(const PlumedOptions& options);
    ~PlumedForceProvider();
    //! Typed discovery through this simulation/rank's force-provider collection.
    PlumedForceProvider* plumedForceProvider() override { return this; }
    //! Mutable exchange/MCMC data, valid for this provider's lifetime.
    PlumedContext& context() { return context_; }
    //! Non-owning access to the provider-owned PLUMED object from the MD engine.
    PLMD::Plumed* plumedHandle() const;
    //! Kernel capability for side-effect-protected calculation-only trials.
    bool hasCalculationOnlyProtection() const { return plumedAPIversion_ >= 14; }
    //! Use early cached exchange energies only when both kernel and engine allow it.
    bool usesPreUpdateExchange() const { return hasCalculationOnlyProtection() && context_.preUpdateExchange; }
#if GMX_LIB_MPI
    //! Evaluate the foreign bias from existing cross potentials, before the local update.
    //! Called collectively on PP ranks after the physical potentials are accumulated.
    void prepareExchange(const gmx_enerdata_t& energy);
#endif
    /*! @brief Sets the global atom indices list
     * param[in] globalAtomIndices  the, optional, list of global atom indices
     */
    void setGlobalAtomIndices(const std::optional<ArrayRef<const int>>& globalAtomIndices);
    /*! \brief Tells PLUMED to output the checkpoint data
     *
     * This is called on the main rank only, when GROMACS writes a checkpoint. It is only a
     * fallback for the simulators that do not call prepareStep(), such as the modular
     * simulator, which writes its checkpoints before the force calculation of the step.
     * calculateForces() then tells PLUMED on all ranks.
     * If the PLUMED API version is not greater than 3 it will do nothing.
     */
    void writeCheckpointData();
    /*! \brief Calculate the forces with PLUMED
     * \param[in] forceProviderInput input for force provider
     * \param[out] forceProviderOutput output for force provider
     */
    void calculateForces(const ForceProviderInput& forceProviderInput,
                         ForceProviderOutput*      forceProviderOutput) override;
    /*! \brief Prepares PLUMED for \p step, and returns whether it needs the potential energy
     * \param[in] step                 the step that is about to be computed
     * \param[in] isCheckpointingStep  whether GROMACS writes a checkpoint at \p step
     */
    bool prepareStep(int64_t step, bool isCheckpointingStep) override;
    //! Whether PLUMED asked to stop the simulation during the last calculation
    bool requestsStop() const override;
    /*! \brief Passes the potential energy to PLUMED and applies energy biases
     * \param[in]     energyWasComputed  whether the potential energy is valid this step
     * \param[in,out] potentialEnergy    the MD potential energy
     * \param[in,out] force              the complete force array
     * \param[in,out] virial             the force virial
     * \returns whether the atoms that PLUMED used include virtual sites
     */
    bool applyAfterPotentialEnergy(bool           energyWasComputed,
                                   real*          potentialEnergy,
                                   ArrayRef<RVec> force,
                                   tensor         virial) override;

private:
    //! Register local or trial inputs, retaining the provider-owned box binding.
    void bindInput(const ForceProviderInput& input);
#if GMX_LIB_MPI
    //! Current force call's input: owns its box and refers to MD-owned atom arrays.
    std::optional<ForceProviderInput> stepInput_;
    //! Simulation-owned replica communication and metadata.
    const MpiComm* mpiComm_ = nullptr;
    const gmx_multisim_t* multiSim_ = nullptr;
    //! Global atom count for temporary coordinate communication.
    int numAtoms_ = 0;
    //! Snapshot/validate complete local inputs and read the already calculated bias.
    //! Performs no additional forward calculation under the strict STRIDE policy.
    real localExchangeBias(real physicalEnergy);
#endif
    //! Whether the atoms that PLUMED used on this step include virtual sites
    bool usesVirtualSites();

    PlumedContext context_;
    std::unique_ptr<PLMD::Plumed> plumed_;
    int                           plumedAPIversion_;
    bool                          replex_;
    //! Whether lambda is supplied for expanded-ensemble MCMC.
    bool                          expanded_ = false;
    //! Whether each atom is a virtual site, empty if there are none
    std::vector<bool>             isVirtualSite_;

    std::optional<ArrayRef<const int>> globalAtomIndices_;
    int                                plumedNeedsEnergy_ = 0;
    int                                plumedWantsToStop_ = 0;
    matrix                             plumedVir_ = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
    //! PLUMED reuses this binding after the short-lived ForceProviderInput is destroyed.
    matrix                             plumedBox_ = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
    //! Step for which prepareStep() already prepared the dependencies, if any
    std::optional<int64_t> preparedStep_;
    //! Whether prepareStep() is called, so that PLUMED is told about checkpoints before the calculation
    bool prepareStepIsCalled_ = false;
    //! Set on the main rank by writeCheckpointData(), passed to PLUMED by calculateForces()
    int checkpointRequested_ = 0;
    //! Step for which prepareStep() was last called
    int64_t currentStep_ = 0;
    //! Last step that PLUMED calculated
    std::optional<int64_t> calculatedStep_;
    //! Whether PLUMED already calculated this MD step (shell relaxation): it is then not updated again
    bool repeatedCall_ = false;
};

} // namespace gmx

#endif // GMX_APPLIED_FORCES_PLUMEDFORCEPROVIDER_H
