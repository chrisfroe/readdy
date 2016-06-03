//
// Created by Moritz Hoffmann on 18/02/16.
//

#ifndef READDY_SIMULATION_H
#define READDY_SIMULATION_H

#include <memory>

#include <boost/predef.h>
#include <vector>
#include <readdy/model/Vec3.h>
#include <readdy/model/KernelStateModel.h>
#include <functional>

#if BOOST_OS_MACOS
#include <array>
#endif

namespace readdy {
    /**
     * Simulation is the focus of the high-level C++ API of ReaDDy.
     * This is where the system is set up and run for a certain number of
     * timesteps.
     * Things like temperature, boxsize, reactions and potentials belong to the
     * context and are given to the kernel when run() is called.
     */
    class Simulation {
    public:
        Simulation();

        ~Simulation();

        // move
        Simulation(Simulation &&rhs);

        Simulation &operator=(Simulation &&rhs);

        double getKBT() const;

        void setKBT(double kBT);

        std::array<double, 3> getBoxSize() const;

        void setBoxSize(double dx, double dy, double dz);

        std::array<bool, 3> getPeriodicBoundary() const;

        void setPeriodicBoundary(std::array<bool, 3> periodic);

        void registerParticleType(const std::string name, const double diffusionCoefficient);
        //void registerPotential(const Potential& potential);
        //void registerReaction(const Reaction& reaction);
        //void registerReactionByDescriptor(const std::string descriptor);

        void addParticle(double x, double y, double z, const std::string& type);

        const std::vector<readdy::model::Vec3> getParticlePositions() const;

        void setKernel(const std::string& kernel);

        bool isKernelSelected() const;

        const std::string& getSelectedKernelType() const;

        virtual void run(const readdy::model::time_step_type steps, const double timeStep);

    private:
        struct Impl;
        std::unique_ptr<readdy::Simulation::Impl> pimpl;
    };

    class NoKernelSelectedException : public std::runtime_error {
    public:
        NoKernelSelectedException(const std::string &__arg);
    };

}

#endif //READDY_SIMULATION_H
