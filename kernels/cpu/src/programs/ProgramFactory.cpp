/**
 * << detailed description >>
 *
 * @file CPUProgramFactory.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 23.06.16
 */

#include <readdy/kernel/cpu/programs/ProgramFactory.h>
#include <readdy/kernel/cpu/programs/EulerBDIntegrator.h>
#include <readdy/kernel/cpu/programs/UpdateNeighborList.h>
#include <readdy/kernel/cpu/programs/CalculateForces.h>
#include <readdy/kernel/cpu/programs/Compartments.h>
#include <readdy/kernel/cpu/programs/reactions/Gillespie.h>
#include <readdy/kernel/cpu/programs/reactions/UncontrolledApproximation.h>
#include <readdy/kernel/cpu/programs/reactions/GillespieParallel.h>
#include <readdy/kernel/cpu/programs/reactions/NextSubvolumesReactionScheduler.h>

namespace core_p = readdy::model::programs;

namespace readdy {
namespace kernel {
namespace cpu {
namespace programs {
ProgramFactory::ProgramFactory(Kernel *kernel) {
    factory[core_p::getProgramName<core_p::reactions::UncontrolledApproximation>()] = [kernel] {
        return new reactions::UncontrolledApproximation(kernel);
    };
    factory[core_p::getProgramName<core_p::EulerBDIntegrator>()] = [kernel] {
        return new EulerBDIntegrator(kernel);
    };
    factory[core_p::getProgramName<core_p::UpdateNeighborList>()] = [kernel] {
        return new UpdateNeighborList(kernel);
    };
    factory[core_p::getProgramName<core_p::CalculateForces>()] = [kernel] {
        return new CalculateForces(kernel);
    };
    factory[core_p::getProgramName<core_p::reactions::Gillespie>()] = [kernel] {
        return new reactions::Gillespie(kernel);
    };
    factory[core_p::getProgramName<core_p::reactions::GillespieParallel>()] = [kernel] {
        return new reactions::GillespieParallel(kernel);
    };
    factory[core_p::getProgramName<core_p::reactions::NextSubvolumes>()] = [kernel] {
        return new reactions::NextSubvolumes(kernel);
    };
    factory[core_p::getProgramName<core_p::Compartments>()] = [kernel] {
        return new Compartments(kernel);
    };
}
}
}
}
}