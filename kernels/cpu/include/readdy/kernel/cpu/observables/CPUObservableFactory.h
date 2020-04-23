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
 * << detailed description >>
 *
 * @file ObservableFactory.h
 * @brief << brief description >>
 * @author clonker
 * @date 21.07.16
 */

#pragma once

#include <readdy/model/observables/ObservableFactory.h>

namespace readdy::kernel::cpu {
class CPUKernel;
namespace observables {

class CPUObservableFactory : public readdy::model::observables::ObservableFactory {

public:
    explicit CPUObservableFactory(CPUKernel* kernel);

    std::unique_ptr<model::observables::Energy>
    energy(Stride stride, ObsCallBack <model::observables::Energy> callback) const override;

    std::unique_ptr<model::observables::Virial>
    virial(Stride stride, ObsCallBack <model::observables::Virial> callback) const override;

    std::unique_ptr<model::observables::HistogramAlongAxis>
    histogramAlongAxis(Stride stride, std::vector<scalar> binBorders, std::vector<std::string> typesToCount,
                       unsigned int axis, ObsCallBack <model::observables::HistogramAlongAxis> callback) const override;

    std::unique_ptr<model::observables::NParticles>
    nParticles(Stride stride, std::vector<std::string> typesToCount,
               ObsCallBack <model::observables::NParticles> callback) const override;

    std::unique_ptr<model::observables::Forces>
    forces(Stride stride, std::vector<std::string> typesToCount,
           ObsCallBack <model::observables::Forces> callback) const override;

    std::unique_ptr<Positions>
    positions(Stride stride, std::vector<std::string> typesToCount, ObsCallBack <Positions> callback) const override;

    std::unique_ptr<RadialDistribution>
    radialDistribution(Stride stride, std::vector<scalar> binBorders, std::vector<std::string> typeCountFrom,
                       std::vector<std::string> typeCountTo, scalar particleDensity,
                       ObsCallBack <RadialDistribution> callback) const override;

    std::unique_ptr<Particles> particles(Stride stride, ObsCallBack <Particles> callback) const override;

    std::unique_ptr<Reactions> reactions(Stride stride, ObsCallBack <Reactions> callback) const override;

    std::unique_ptr<ReactionCounts> reactionCounts(Stride stride, ObsCallBack <ReactionCounts> callback) const override;

private:
    CPUKernel *const kernel;
};

}
}
