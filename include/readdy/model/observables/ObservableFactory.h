/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 ********************************************************************/


/**
 * This header file contains the definition of the ObservableFactory. Its purpose is to create observable of different
 * types in the form of unique_ptrs. The actual implementation of an observable can be changed by specializing the
 * dispatcher for its type and invoking a virtual (and then: overridden) method within the factory.
 *
 * @file ObservableFactory.h
 * @brief This header file contains the definition of the ObservableFactory.
 * @author clonker
 * @date 29.04.16
 */
#pragma once

#include <string>
#include <unordered_map>
#include <readdy/common/Utils.h>
#include <readdy/model/observables/Observables.h>

namespace readdy::model {
class Kernel;
namespace observables {

namespace detail {
template<typename T>
using is_observable_type = std::enable_if_t<std::is_base_of<model::observables::ObservableBase, T>::value>;
}

class ObservableFactory {
public:
    template<typename T>
    using ObsCallBack = typename std::function<void(typename T::result_type)>;

    explicit ObservableFactory(Kernel *const kernel) : kernel(kernel) {};

    virtual std::unique_ptr<Energy> energy(Stride stride, ObsCallBack<Energy> callback = [](Energy::result_type){}) const = 0;

    virtual std::unique_ptr<Virial> virial(Stride stride, ObsCallBack<Virial> callback = [](Virial::result_type){}) const = 0;
    
    virtual std::unique_ptr<HistogramAlongAxis>
    histogramAlongAxis(Stride stride, std::vector<scalar> binBorders, std::vector<std::string> typesToCount,
                       unsigned int axis, ObsCallBack<HistogramAlongAxis> callback = [](HistogramAlongAxis::result_type){}) const = 0;
    
    std::unique_ptr<NParticles>
    nParticles(Stride stride, ObsCallBack<NParticles> callback = [](NParticles::result_type){}) const { return nParticles(stride, {}, callback); }
    
    virtual std::unique_ptr<NParticles>
    nParticles(Stride stride, std::vector<std::string> typesToCount, ObsCallBack<NParticles> callback = [](NParticles::result_type){}) const = 0;

    std::unique_ptr<Forces>
    forces(Stride stride, ObsCallBack<Forces> callback = [](Forces::result_type){}) const { return forces(stride, {}, callback); }
    
    virtual std::unique_ptr<Forces>
    forces(Stride stride, std::vector<std::string> typesToCount, ObsCallBack<Forces> callback = [](Forces::result_type){}) const  = 0;

    std::unique_ptr<Positions>
    positions(Stride stride, ObsCallBack <Positions> callback = [](Positions::result_type){}) const { return positions(stride, {}, callback); }

    virtual std::unique_ptr<Positions>
    positions(Stride stride, std::vector<std::string> typesToCount, ObsCallBack <Positions> callback = [](Positions::result_type){}) const = 0;

    virtual std::unique_ptr<RadialDistribution>
    radialDistribution(Stride stride, std::vector<scalar> binBorders, std::vector<std::string> typeCountFrom,
                       std::vector<std::string> typeCountTo, scalar particleDensity,
                       ObsCallBack<RadialDistribution> callback = [](RadialDistribution::result_type){}) const = 0;

    virtual std::unique_ptr<Particles> particles(Stride stride, ObsCallBack<Particles> callback = [](Particles::result_type){}) const = 0;

    virtual std::unique_ptr<Reactions> reactions(Stride stride, ObsCallBack<Reactions> callback = [](Reactions::result_type){}) const = 0;

    virtual std::unique_ptr<ReactionCounts> reactionCounts(Stride stride, ObsCallBack<ReactionCounts> callback = [](ReactionCounts::result_type){}) const = 0;

    std::unique_ptr<Trajectory> trajectory(Stride stride, ObsCallBack<Trajectory> callback = [](Trajectory::result_type){}) const {
        auto obs = std::make_unique<Trajectory>(kernel, stride);
        obs->setCallback(callback);
        return std::move(obs);
    }

    std::unique_ptr<FlatTrajectory> flatTrajectory(Stride stride, ObsCallBack<FlatTrajectory> callback = [](FlatTrajectory::result_type){}) const {
        auto obs = std::make_unique<FlatTrajectory>(kernel, stride);
        obs->setCallback(callback);
        return std::move(obs);
    }

    std::unique_ptr<Topologies> topologies(Stride stride, ObsCallBack<Topologies> callback = [](Topologies::result_type){}) const {
        auto obs = std::make_unique<Topologies>(kernel, stride);
        obs->setCallback(callback);
        return std::move(obs);
    }

protected:
    Kernel *const kernel;
};

}
}
