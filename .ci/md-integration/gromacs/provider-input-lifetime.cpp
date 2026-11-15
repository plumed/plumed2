// Validate retained box bindings for GREX snapshots and lambda MCMC trials.
#include <atomic>
#include <cmath>
#include <thread>

#include <fstream>
#include <iostream>
#include <stdexcept>

#include "gromacs/applied_forces/plumed/plumedOptions.h"
#include "gromacs/applied_forces/plumed/plumedforceprovider.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdrun/plumedexchange.h"
#include "gromacs/mdrunutility/handlerestart.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/utility/init.h"
#include "gromacs/utility/mpicomm.h"
#define __PLUMED_WRAPPER_CXX 1
#define __PLUMED_WRAPPER_LIBCXX11 1
#define __PLUMED_WRAPPER_EXTERN 0
#define __PLUMED_WRAPPER_FORTRAN 0
#define __PLUMED_WRAPPER_LINK_RUNTIME 1
#define __PLUMED_WRAPPER_IMPLEMENTATION 1
#define __PLUMED_HAS_DLOPEN
#include "Plumed.h"

int main(int argc, char** argv)
{
    gmx::init(&argc, &argv);
    gmx_omp_nthreads_set(ModuleMultiThread::Default, 1);
    try
    {
        std::ofstream("plumed.dat") << "d: DISTANCE ATOMS=1,2\n"
                                       "r: RESTRAINT ARG=d AT=0 KAPPA=2\n";
        gmx::MpiComm       comm(gmx::MpiComm::SingleRank{});
        gmx::PlumedOptions options;
        options.plumedFile_         = "plumed.dat";
        options.natoms_             = 2;
        options.replex_             = true;
        options.mpiComm_            = &comm;
        options.simulationTimeStep_ = 0.001;
        options.startingBehavior_   = gmx::StartingBehavior::NewSimulation;
        {
            gmx::PlumedForceProvider provider(options);
            gmx::ForceProviders providers;
            providers.addForceProvider(&provider);
            if (gmx::findPlumedForceProvider(&providers) != &provider)
                throw std::runtime_error("PLUMED provider lookup does not match its collection");
            std::vector<gmx::RVec>   x{ { 0, 0, 0 }, { real(0.45), 0, 0 } }, f(2);
            std::vector<real>        mass(2, 1), charge(2, 0);
            matrix                   box = { { 4, 0, 0 }, { 0, 4, 0 }, { 0, 0, 4 } };
            gmx::ForceProviderInput  input(x, 2, charge, mass, 0, 0, box, comm, nullptr);
            gmx::ForceWithVirial     force(f, true);
            gmx_enerdata_t           energy(1, nullptr);
            gmx::ForceProviderOutput output(&force, &energy);
            provider.calculateForces(input, &output);
            PLMD::Plumed plumed(provider.plumedHandle());
#if GMX_LIB_MPI
            plumed.cmd("GREX init");
#endif
            real before = 0, after = 0;
            plumed.cmd("getBias", &before);
            // Model reuse of the short-lived input storage without invoking undefined behavior.
            input.box_[0][0] = input.box_[1][1] = input.box_[2][2] = real(0.5);
            real physicalEnergy = 0;
            long int snapshotStep = 0;
            plumed.cmd("setStepLong", &snapshotStep);
            plumed.cmd("setEnergy", &physicalEnergy);
            real snapshot = before;
#if GMX_LIB_MPI
            plumed.cmd("GREX savePositions");
            int partner = 0;
            plumed.cmd("GREX setPartner", &partner);
            plumed.cmd("GREX calculate");
            plumed.cmd("getBias", &snapshot);
            std::cout << "afterGREXSnapshot=" << snapshot << " expected=0.2025\n";
#endif
            long int step     = 0;
            int      exchange = 1;
            plumed.cmd("setExchangeStep", &exchange);
            plumed.cmd("setStepLong", &step);
            plumed.cmd("prepareCalc");
            plumed.cmd("performCalcNoForces");
            plumed.cmd("getBias", &after);
            std::cout << "before=" << before << " after=" << after << " expected=0.2025\n";
            if (std::abs(before - real(0.2025)) > 1e-6 || std::abs(after - before) > 1e-6
                || std::abs(snapshot - before) > 1e-6)
                throw std::runtime_error("PLUMED retains the temporary force input box");
        }
        // Two live providers/collections must coexist and remain independent.
        gmx::PlumedForceProvider firstProvider(options), secondProvider(options);
        gmx::ForceProviders first, second, empty;
        first.addForceProvider(&firstProvider);
        second.addForceProvider(&secondProvider);
        if (gmx::findPlumedForceProvider(nullptr) || gmx::findPlumedForceProvider(&empty))
            throw std::runtime_error("Empty collection found an unrelated provider");
        struct OtherProvider : gmx::IForceProvider
        {
            void calculateForces(const gmx::ForceProviderInput&, gmx::ForceProviderOutput*) override {}
        } other;
        empty.addForceProvider(&other);
        if (gmx::findPlumedForceProvider(&empty))
            throw std::runtime_error("Non-PLUMED provider was mistaken for PLUMED");
        std::atomic<bool> isolated{true};
        auto probe = [&](gmx::ForceProviders* collection, gmx::PlumedForceProvider* expected, real value)
        {
            for (int i = 0; i < 100000; ++i)
            {
                auto* provider = gmx::findPlumedForceProvider(collection);
                if (provider != expected)
                {
                    isolated = false;
                    return;
                }
                provider->context().lambdaState = value;
                if (provider->context().lambdaState != value)
                    isolated = false;
            }
        };
        std::thread a(probe, &first, &firstProvider, real(0.25));
        std::thread b(probe, &second, &secondProvider, real(0.75));
        a.join();
        b.join();
        if (!isolated || firstProvider.context().lambdaState != real(0.25)
            || secondProvider.context().lambdaState != real(0.75))
            throw std::runtime_error("Provider lookup or buffers are shared between ranks");
        std::cout << "provider ownership and 200000 concurrent collection lookups: pass\n";
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << "\n";
        gmx::finalize();
        return 1;
    }
    gmx::finalize();
}
