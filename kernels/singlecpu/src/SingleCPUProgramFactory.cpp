//
// Created by clonker on 08.04.16.
//

#include <readdy/common/make_unique.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUProgramFactory.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUTestProgram.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUAddParticleProgram.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUEulerBDIntegrator.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUCalculateForces.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUReactionImpls.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUUpdateNeighborList.h>

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            namespace programs {
                SingleCPUProgramFactory::SingleCPUProgramFactory(SingleCPUKernel *kernel) : kernel(kernel) {
                    namespace core_p = readdy::model::programs;
                    factory[core_p::getProgramName<core_p::Test>()] = [] { return new SingleCPUTestProgram(); };
                    factory[core_p::getProgramName<core_p::AddParticle>()] = [kernel] {
                        return new SingleCPUAddParticleProgram(kernel);
                    };
                    factory[core_p::getProgramName<core_p::EulerBDIntegrator>()] = [kernel] {
                        return new SingleCPUEulerBDIntegrator(kernel);
                    };
                    factory[core_p::getProgramName<core_p::UpdateNeighborList>()] = [kernel] {
                        return new SingleCPUUpdateNeighborList(kernel);
                    };
                    factory[core_p::getProgramName<core_p::CalculateForces>()] = [kernel] {
                        return new SingleCPUCalculateForces(kernel);
                    };
                    factory[core_p::getProgramName<core_p::reactions::UncontrolledApproximation>()] = [kernel] {
                        return new reactions::UncontrolledApproximation(kernel);
                    };
                }
            }
        }
    }
}

