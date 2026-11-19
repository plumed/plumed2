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
 * Implementation of the Plumed force provider class
 *
 * \author Daniele Rapetti <drapetti@sissa.it>
 * \ingroup module_applied_forces
 */

#include "plumedforceprovider.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <limits>

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/math/units.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdrunutility/handlerestart.h"
#include "gromacs/mdrunutility/multisim.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/mpicomm.h"
#include "gromacs/utility/vec.h"
#include "gromacs/version.h"

#include "plumedOptions.h"

#define __PLUMED_WRAPPER_FORTRAN 0 // NOLINT(bugprone-reserved-identifier)

#define __PLUMED_WRAPPER_LINK_RUNTIME 1 // NOLINT(bugprone-reserved-identifier)
#define __PLUMED_WRAPPER_EXTERN 0       // NOLINT(bugprone-reserved-identifier)

#define __PLUMED_WRAPPER_CXX 1            // NOLINT(bugprone-reserved-identifier)
#define __PLUMED_WRAPPER_LIBCXX11 1       // NOLINT(bugprone-reserved-identifier)
#define __PLUMED_WRAPPER_LIBCXX17 1       // NOLINT(bugprone-reserved-identifier)
#define __PLUMED_WRAPPER_IMPLEMENTATION 1 // NOLINT(bugprone-reserved-identifier)
#define __PLUMED_HAS_DLOPEN               // NOLINT(bugprone-reserved-identifier)

// Embed the adapter using this host build's MPI types and functions.
#if (GMX_LIB_MPI || GMX_THREAD_MPI) && !defined(__PLUMED_WRAPPER_MPI_HOST_ADAPTER)
#define __PLUMED_WRAPPER_MPI_HOST_ADAPTER 1
#endif
#include "../../../../Plumed.h"
#include "../../mdrun/plumedexchange.h"

namespace gmx
{

PLMD::Plumed* PlumedForceProvider::plumedHandle() const
{
    return plumed_.get();
}

PlumedForceProvider::~PlumedForceProvider() = default;
PlumedForceProvider::PlumedForceProvider(const PlumedOptions& options)
try : plumed_(std::make_unique<PLMD::Plumed>()),
      replex_(options.replex_),
      isVirtualSite_(options.isVirtualSite_)
{
    // I prefer to pass a struct with data because it stops the coupling
    // at the implementation and not at the function signature:
    // less code to edit when adding new options :)
    expanded_ = options.inputrec_ != nullptr && options.inputrec_->bExpanded;
    if (expanded_ && replex_)
        GMX_THROW(InvalidInputError("Expanded-ensemble lambda MCMC combined with replica exchange is unsupported"));
    context_.mainRank = options.mpiComm_->isMainRank();
    context_.hrex     = options.hrex_;
#if GMX_LIB_MPI
    mpiComm_          = options.mpiComm_;
    multiSim_         = options.ms_;
    numAtoms_         = options.natoms_;
#endif
    plumed_->cmd("getApiVersion", &plumedAPIversion_);
#if __PLUMED_WRAPPER_MPI_HOST_ADAPTER
    if (plumedAPIversion_ >= 13)
    {
        // Inject before the first command that can use MPI.
        plumed_->cmd("injectMPIImplementation", plumed_mpi_implementation());
    }
    else if (options.mpiComm_->isMainRank())
    {
        std::fprintf(stderr,
                     "WARNING: This PLUMED kernel (API %d) does not support MPI injection. "
                     "For multi-rank runs, use a PLUMED kernel built with the same MPI library "
                     "as GROMACS, or upgrade PLUMED.\n",
                     plumedAPIversion_);
    }
#endif
#if GMX_THREAD_MPI
    if (plumedAPIversion_ < 13 && options.mpiComm_->isParallel())
    {
        GMX_THROW(InvalidInputError(
                "plumed MPI interface is not compatible with THREAD_MPI when uses more than one "
                "rank"));
    }
#endif
    int real_precision = sizeof(real);
    plumed_->cmd("setRealPrecision", &real_precision);
    real energyUnits = 1.0;
    plumed_->cmd("setMDEnergyUnits", &energyUnits);
    real lengthUnits = 1.0;
    plumed_->cmd("setMDLengthUnits", &lengthUnits);
    real timeUnits = 1.0;
    plumed_->cmd("setMDTimeUnits", &timeUnits);

    plumed_->cmd("setPlumedDat", options.plumedFile_.c_str());

    if (plumedAPIversion_ > 1)
    {
        /* setting kbT is only implemented with api>1) */
        if (options.ensembleTemperature_.has_value())
        {
            real kbT = options.ensembleTemperature_.value() * gmx::c_boltz;
            plumed_->cmd("setKbT", &kbT);
        }
    }

    if (plumedAPIversion_ > 2)
    {
        if ((options.startingBehavior_ != StartingBehavior::NewSimulation))
        {
            int res = 1;
            plumed_->cmd("setRestart", &res);
        }
    }

    if (plumedAPIversion_ > 11)
    {
        /* The step from which the simulation continues after -cpi: PLUMED removes from its
           output files the lines written at or after it, as GROMACS does with its own files */
        if (options.startingBehavior_ != StartingBehavior::NewSimulation && options.inputrec_ != nullptr)
        {
            long long int restartStep = options.inputrec_->init_step;
            plumed_->cmd("setRestartStep", &restartStep);
        }
    }

    if (plumedAPIversion_ > 5)
    {
        /* Tell PLUMED how many OpenMP threads GROMACS is using, so that it can
           parallelise its own work. Without this PLMD::OpenMP::getNumThreads()
           always reports 1 unless PLUMED_NUM_THREADS is set by hand. */
        int numOmpThreads = gmx_omp_nthreads_get(ModuleMultiThread::Default);
        plumed_->cmd("setNumOMPthreads", &numOmpThreads);
    }

    if (isMultiSim(options.ms_))
    {
        /* Hand PLUMED the multi-replica communicators. This is what populates
           PLMD::PlumedMain::multi_sim_comm, without which every multi-replica
           action silently degrades to a single replica. */
        if (options.mpiComm_->isMainRank())
        {
            plumed_->cmd("GREX setMPIIntercomm", &options.ms_->mainRanksComm_);
        }
        MPI_Comm intracomm = options.mpiComm_->comm();
        plumed_->cmd("GREX setMPIIntracomm", &intracomm);
        plumed_->cmd("GREX init", nullptr);
    }

    if (options.mpiComm_->isParallel())
    {
        /* PLUMED reads this argument as an MPI_Comm*, so it must be handed the
           address of a communicator rather than the communicator itself. */
        MPI_Comm plumedComm = options.mpiComm_->comm();
        plumed_->cmd("setMPIComm", &plumedComm);
    }

    plumed_->cmd("setNatoms", options.natoms_);
    plumed_->cmd("setMDEngine", "gromacs");
    {
        const char*       logFileEnv = std::getenv("PLUMED_LOG_FILE");
        const std::string logFile =
                (logFileEnv && logFileEnv[0] != '\0') ? logFileEnv : "PLUMED.OUT";
        plumed_->cmd("setLogFile", logFile.c_str());
    }

    plumed_->cmd("setTimestep", &options.simulationTimeStep_);
    plumed_->cmd("init", nullptr);
}
catch (const GromacsException&)
{
    // Keep the type and message of the errors raised on purpose, e.g. NotImplementedError
    throw;
}
catch (const std::exception& ex)
{
    GMX_THROW(InternalError(
            std::string("An error occurred while initializing the PLUMED force provider:\n") + ex.what()));
}

void PlumedForceProvider::setGlobalAtomIndices(const std::optional<ArrayRef<const int>>& globalAtomIndices)
{
    globalAtomIndices_ = globalAtomIndices;
}

void PlumedForceProvider::writeCheckpointData()
try
{
#if GMX_LIB_MPI
    if (context_.hrexTrial)
        return;
#endif
    if (!prepareStepIsCalled_)
    {
        checkpointRequested_ = 1;
    }
}
catch (const GromacsException&)
{
    // Keep the type and message of the errors raised on purpose, e.g. NotImplementedError
    throw;
}
catch (const std::exception& ex)
{
    GMX_THROW(InternalError(
            std::string("An error occurred while PLUMED was writing the checkpoint data\n:") + ex.what()));
}

void PlumedForceProvider::bindInput(const ForceProviderInput& forceProviderInput)
{
    long int lstep = prepareStepIsCalled_ ? currentStep_ : forceProviderInput.step_;
    // prepareStep() already opened this input epoch. Repeating setStep would
    // clear older kernels' prepared energy requirements before sharing data.
    if (preparedStep_ != lstep)
    {
        plumed_->cmd("setStepLong", &lstep);
    }
    if (expanded_)
    {
        context_.lambdaForce = 0;
        plumed_->cmd("setExtraCV lambda", &context_.lambdaState);
        plumed_->cmd("setExtraCVForce lambda", &context_.lambdaForce);
    }
    if (globalAtomIndices_.has_value())
    {
        int nat_home = forceProviderInput.homenr_;
        plumed_->cmd("setAtomsNlocal", &nat_home);
        plumed_->cmd("setAtomsGatindex", globalAtomIndices_.value().data());
    }

    plumed_->cmd("setPositions", &(forceProviderInput.x_.data()->as_vec()[0]));
    plumed_->cmd("setMasses", forceProviderInput.massT_.data());
    plumed_->cmd("setCharges", forceProviderInput.chargeA_.data());
    // GREX snapshots and lambda trials read the box after this callback's
    // short-lived input container has gone out of scope. Keep the binding alive.
    copy_mat(forceProviderInput.box_, plumedBox_);
    plumed_->cmd("setBox", plumedBox_);

}

void PlumedForceProvider::calculateForces(const ForceProviderInput& forceProviderInput,
                                          ForceProviderOutput*      forceProviderOutput)
try
{
#if GMX_LIB_MPI
    if (context_.hrexTrial)
    {
        if (usesPreUpdateExchange())
            stepInput_.emplace(forceProviderInput);
        return;
    }
    if (context_.exchangeStep && usesPreUpdateExchange())
    {
        // Delay the local bias until the foreign trial and physical energy are complete.
        stepInput_.emplace(forceProviderInput);
        return;
    }
#endif
    // setup: these instructions in the original patch are BEFORE do_force()
    // now this is called within do_force(), but this does not impact the results
    long int lstep = prepareStepIsCalled_ ? currentStep_ : forceProviderInput.step_;
    // With shell particles, relax_shell_flexcon() calls do_force() again on the same step (passing
    // step 1): PLUMED then calculates its forces again, without updating (depositing, printing...).
    // "Again" means since the last prepareStep(), which MD and reruns call before each step
    repeatedCall_   = prepareStepIsCalled_ && (calculatedStep_ == lstep);
    calculatedStep_ = lstep;
    bindInput(forceProviderInput);

    if (preparedStep_ == lstep)
    {
        // prepareStep() already ran the first half of prepareCalc
        plumed_->cmd("shareData", nullptr);
    }
    else
    {
        plumed_->cmd("prepareCalc", nullptr);
        plumedNeedsEnergy_ = 0;
        plumed_->cmd("isEnergyNeeded", &plumedNeedsEnergy_);
    }
    preparedStep_.reset();

    if (!repeatedCall_)
    {
        plumedWantsToStop_ = 0;
    }
    plumed_->cmd("setStopFlag", &plumedWantsToStop_);

    if (!prepareStepIsCalled_)
    {
        // writeCheckpointData() runs only on the main rank: pass its request to all ranks
        if (forceProviderInput.mpiComm_.isParallel())
        {
            forceProviderInput.mpiComm_.sumReduce(1, &checkpointRequested_);
        }
        if (plumedAPIversion_ > 3)
        {
            int checkp = checkpointRequested_ != 0 ? 1 : 0;
            plumed_->cmd("doCheckPoint", &checkp);
        }
        checkpointRequested_ = 0;
    }

    // Biasing the energy rescales the total MD force, which is only complete once the
    // potential energy has been accumulated: finish in applyAfterPotentialEnergy()
    if (plumedNeedsEnergy_)
    {
        return;
    }

    real* fOut = &(forceProviderOutput->forceWithVirial_.force_.data()->as_vec()[0]);
    plumed_->cmd("setForces", fOut);

    clear_mat(plumedVir_);
    plumed_->cmd("setVirial", &plumedVir_[0][0]);

    // end setup: these instructions in the original patch are BEFORE do_force()
    // in the original patch do_force() was called HERE

    // Do the work
    plumed_->cmd(repeatedCall_ ? "performCalcNoUpdate" : "performCalc", nullptr);

    msmul(plumedVir_, 0.5, plumedVir_);
    forceProviderOutput->forceWithVirial_.addVirialContribution(plumedVir_);
}
catch (const GromacsException&)
{
    // Keep the type and message of the errors raised on purpose, e.g. NotImplementedError
    throw;
}
catch (const std::exception& ex)
{
    GMX_THROW(InternalError(
            std::string("An error occurred while PLUMED was calculating the forces\n:") + ex.what()));
}

bool PlumedForceProvider::requestsStop() const
{
    return plumedWantsToStop_ != 0;
}

bool PlumedForceProvider::usesVirtualSites()
{
    if (isVirtualSite_.empty())
    {
        return false;
    }
    /* PLUMED applies forces only on the atoms it uses. Their list is the same on all ranks,
       so that all take the same decision. It is empty when PLUMED only biases the energy. */
    int numAtoms = 0;
    plumed_->cmd("createFullList", &numAtoms);
    const int* atoms = nullptr;
    plumed_->cmd("getFullList", &atoms);
    const bool found =
            std::any_of(atoms, atoms + numAtoms, [this](int atom) { return isVirtualSite_[atom]; });
    plumed_->cmd("clearFullList", nullptr);
    return found;
}

bool PlumedForceProvider::prepareStep(int64_t step, bool isCheckpointingStep)
try
{
    /* Whether PLUMED needs the energy depends on which actions are active on this step,
       which is only known after prepareDependencies(). That is the first half of
       prepareCalc(); calculateForces() runs the second half (shareData()) once the
       positions are available. */
    long int lstep = step;
    plumed_->cmd("setStepLong", &lstep);
    if (expanded_)
    {
        context_.lambdaForce = 0;
        plumed_->cmd("setExtraCV lambda", &context_.lambdaState);
        plumed_->cmd("setExtraCVForce lambda", &context_.lambdaForce);
    }
    context_.foreignReady = false;
    context_.exchangeReady = false;
    prepareStepIsCalled_ = true;
    currentStep_         = step;
    calculatedStep_.reset();
    if (plumedAPIversion_ > 3)
    {
        // Set on every step: PLUMED resets it only on the steps when some action is active
        int checkp = isCheckpointingStep ? 1 : 0;
        plumed_->cmd("doCheckPoint", &checkp);
    }
    plumed_->cmd("prepareDependencies", nullptr);
    plumedNeedsEnergy_ = 0;
    plumed_->cmd("isEnergyNeeded", &plumedNeedsEnergy_);
    preparedStep_ = step;
    return plumedNeedsEnergy_ != 0;
}
catch (const GromacsException&)
{
    // Keep the type and message of the errors raised on purpose, e.g. NotImplementedError
    throw;
}
catch (const std::exception& ex)
{
    GMX_THROW(InternalError(
            std::string("An error occurred while PLUMED was preparing the step\n:") + ex.what()));
}

#if GMX_LIB_MPI
void PlumedForceProvider::prepareExchange(const gmx_enerdata_t& energy)
{
    if (!usesPreUpdateExchange() || context_.hrex || context_.foreignReady)
        return;
    GMX_RELEASE_ASSERT(stepInput_.has_value(), "Exchange needs this force call's physical inputs");
    const auto& input = *stepInput_;
    // Native lambda terms are still domain-local here. Reduce scratch copies,
    // leaving GROMACS's ordinary energy arrays for their normal later reduction.
    double potentials[2] = { energy.term[InteractionFunction::PotentialEnergy],
                             energy.term[InteractionFunction::PotentialEnergy] };
    if (context_.partnerLambda >= 0)
        potentials[1] += energy.foreignLambdaTerms.deltaPotentialEnergy(context_.partnerLambda);
    mpiComm_->sumReduce(2, potentials);
    double foreignPotential = potentials[1];

    // Assemble coordinates by atom identity, without changing GROMACS's DD,
    // neighbor lists or force buffers. The trial uses the evaluating rank's
    // local Hamiltonian (masses/charges), at the partner's positions and box.
    const std::size_t coordinateCount = std::size_t(DIM) * numAtoms_;
    GMX_RELEASE_ASSERT(coordinateCount <= std::size_t(std::numeric_limits<int>::max()) / sizeof(real),
                       "PLUMED trial coordinate buffer exceeds the MPI count limit");
    std::vector<real> positions(coordinateCount, 0);
    for (int i = 0; i < input.homenr_; ++i)
    {
        const int global = globalAtomIndices_ ? (*globalAtomIndices_)[i] : i;
        for (int d = 0; d < DIM; ++d)
            positions[DIM * global + d] = input.x_[i][d];
    }
    mpiComm_->sumReduce(static_cast<int>(positions.size()), positions.data());
    matrix foreignBox;
    copy_mat(input.box_, foreignBox);
    if (context_.mainRank && context_.partner != multiSim_->simulationIndex_)
    {
        // Send U_partner(x_local); receive U_local(x_partner).
        MPI_Sendrecv_replace(&foreignPotential,
                             1,
                             MPI_DOUBLE,
                             context_.partner,
                             4110,
                             context_.partner,
                             4110,
                             multiSim_->mainRanksComm_,
                             MPI_STATUS_IGNORE);
        MPI_Sendrecv_replace(positions.data(),
                             static_cast<int>(positions.size() * sizeof(real)),
                             MPI_BYTE,
                             context_.partner,
                             4111,
                             context_.partner,
                             4111,
                             multiSim_->mainRanksComm_,
                             MPI_STATUS_IGNORE);
        MPI_Sendrecv_replace(foreignBox,
                             sizeof(matrix),
                             MPI_BYTE,
                             context_.partner,
                             4112,
                             context_.partner,
                             4112,
                             multiSim_->mainRanksComm_,
                             MPI_STATUS_IGNORE);
    }
    MPI_Bcast(&foreignPotential, 1, MPI_DOUBLE, 0, mpiComm_->comm());
    MPI_Bcast(positions.data(), static_cast<int>(positions.size() * sizeof(real)), MPI_BYTE, 0, mpiComm_->comm());
    MPI_Bcast(foreignBox, sizeof(matrix), MPI_BYTE, 0, mpiComm_->comm());
    context_.foreignEnergy = foreignPotential;
    std::vector<RVec> trialPositions(input.homenr_);
    for (int i = 0; i < input.homenr_; ++i)
    {
        const int global = globalAtomIndices_ ? (*globalAtomIndices_)[i] : i;
        for (int d = 0; d < DIM; ++d)
            trialPositions[i][d] = positions[DIM * global + d];
    }
    preparedStep_.reset();
    bindInput(input);
    plumed_->cmd("setPositions", &(trialPositions.data()->as_vec()[0]));
    copy_mat(foreignBox, plumedBox_);
    plumed_->cmd("setBox", plumedBox_);
    context_.trialEnergy = context_.mainRank ? context_.foreignEnergy : 0;
    plumed_->cmd("setEnergy", &context_.trialEnergy);
    int exchange = 1;
    plumed_->cmd("setExchangeStep", &exchange);
    plumed_->cmd("prepareCalc");
    plumed_->cmd("performCalcNoForces");
    plumed_->cmd("getBias", &context_.foreignBias);
    context_.foreignReady = true;
    exchange              = 0;
    plumed_->cmd("setExchangeStep", &exchange);
    // The caller immediately rebinds the local input before the normal force
    // calculation/update; the temporary trialPositions binding is not reused.
    preparedStep_.reset();
}

real PlumedForceProvider::localExchangeBias(real physicalEnergy)
{
    // The strict bias STRIDE policy makes the ordinary local bias complete.
    // savePositions validates that policy and shares inputs, without a calculation.
    long int step = currentStep_;
    plumed_->cmd("setStepLong", &step);
    context_.trialEnergy = physicalEnergy;
    plumed_->cmd("setEnergy", &context_.trialEnergy);
    plumed_->cmd("GREX savePositions");
    real bias = 0;
    plumed_->cmd("getBias", &bias);
    return bias;
}
#endif

bool PlumedForceProvider::applyAfterPotentialEnergy(bool           energyWasComputed,
                                                    real*          potentialEnergy,
                                                    ArrayRef<RVec> force,
                                                    tensor         virial)
try
{
#if GMX_LIB_MPI
    if (usesPreUpdateExchange() && context_.hrexTrial)
    {
        GMX_RELEASE_ASSERT(stepInput_.has_value() && energyWasComputed,
                           "HREX trial must supply coordinates and physical energy");
        preparedStep_.reset();
        bindInput(*stepInput_);
        int exchange = 1;
        plumed_->cmd("setExchangeStep", &exchange);
        plumed_->cmd("prepareCalc");
        plumed_->cmd("setEnergy", potentialEnergy);
        plumed_->cmd("performCalcNoForces");
        plumed_->cmd("getBias", &context_.foreignBias);
        double energy = *potentialEnergy;
        mpiComm_->sumReduce(1, &energy);
        context_.foreignEnergy = energy;
        context_.foreignReady  = true;
        exchange               = 0;
        plumed_->cmd("setExchangeStep", &exchange);
        preparedStep_.reset();
        return false;
    }
    if (usesPreUpdateExchange() && context_.exchangeStep)
    {
        GMX_RELEASE_ASSERT(stepInput_.has_value() && context_.foreignReady,
                           "Exchange requires its foreign bias before the local calculation");
        const long int step = currentStep_;
        repeatedCall_       = calculatedStep_ == step;
        calculatedStep_     = step;
        if (!repeatedCall_)
            plumedWantsToStop_ = 0;
        preparedStep_.reset();
        bindInput(*stepInput_);
        plumed_->cmd("prepareCalc");
        plumedNeedsEnergy_ = 0;
        plumed_->cmd("isEnergyNeeded", &plumedNeedsEnergy_);
        // Run the ordinary local force calculation once, before its bias update.
        plumed_->cmd("setEnergy", potentialEnergy);
        if (!force.empty())
            plumed_->cmd("setForces", &(force.data()->as_vec()[0]));
        msmul(virial, 2.0, plumedVir_);
        plumed_->cmd("setVirial", plumedVir_);
        plumed_->cmd("setStopFlag", &plumedWantsToStop_);
        plumed_->cmd("performCalcNoUpdate");
        const real localBias   = localExchangeBias(*potentialEnergy);
        context_.exchangeDelta = context_.foreignBias - localBias;
        if (context_.hrex)
        {
            double energy = *potentialEnergy;
            mpiComm_->sumReduce(1, &energy);
            context_.exchangeDelta += context_.foreignEnergy - energy;
        }
        context_.exchangeReady = true;
        if (!repeatedCall_)
            plumed_->cmd("update");
        msmul(plumedVir_, 0.5, virial);
        preparedStep_.reset();
        return usesVirtualSites();
    }
#endif
#if GMX_LIB_MPI
    if (context_.hrexTrial)
        return false;
    if (context_.hrex && context_.exchangeStep && energyWasComputed)
        plumed_->cmd("GREX cacheLocalUNow", potentialEnergy);
#endif

    auto saveExchangeInputs = [&]()
    {
#if GMX_LIB_MPI
        if (replex_ && context_.exchangeStep)
        {
            GMX_RELEASE_ASSERT(energyWasComputed && potentialEnergy != nullptr,
                               "Replica exchange requires the physical potential energy");
            // Save only after the normal physical energy has been computed.
            // Older kernels retain those bindings: resetting the step would
            // clear PLUMED 2.9's positions/masses without registering them again.
            if (plumedAPIversion_ >= 13)
            {
                long int savedStep = calculatedStep_.value();
                plumed_->cmd("setStepLong", &savedStep);
                plumed_->cmd("setEnergy", potentialEnergy);
            }
            plumed_->cmd("GREX savePositions");
        }
#endif
    };

    if (!plumedNeedsEnergy_)
    {
        saveExchangeInputs();
        return false;
    }
    if (!energyWasComputed || potentialEnergy == nullptr)
    {
        GMX_THROW(NotImplementedError(
                "The PLUMED input needs the potential energy, but GROMACS did not compute it "
                "on this step. Check the simulator's force/energy scheduling."));
    }

    plumed_->cmd("setEnergy", potentialEnergy);
    if (!force.empty())
    {
        plumed_->cmd("setForces", &(force.data()->as_vec()[0]));
    }
    // Same virial convention as the legacy patch: hand PLUMED 2*virial, then replace it
    msmul(virial, 2.0, plumedVir_);
    plumed_->cmd("setVirial", &plumedVir_[0][0]);
    plumed_->cmd(repeatedCall_ ? "performCalcNoUpdate" : "performCalc", nullptr);
    saveExchangeInputs();


    msmul(plumedVir_, 0.5, virial);
    return usesVirtualSites();
}
catch (const GromacsException&)
{
    // Keep the type and message of the errors raised on purpose, e.g. NotImplementedError
    throw;
}
catch (const std::exception& ex)
{
    GMX_THROW(InternalError(
            std::string("An error occurred while PLUMED was calculating the forces\n:") + ex.what()));
}
} // namespace gmx
