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
 * @file SCPUObservableFactory.h
 * @brief Declaration of observable factory for single cpu kernel
 * @author clonker
 * @date 30.06.16
 */

#pragma once
#include <readdy/model/observables/ObservableFactory.h>

namespace readdy::kernel::scpu {
class SCPUKernel;
namespace observables {

class SCPUObservableFactory : public readdy::model::observables::ObservableFactory {
public:
    explicit SCPUObservableFactory(readdy::kernel::scpu::SCPUKernel* kernel);

    std::unique_ptr<readdy::model::observables::Energy>
    energy(Stride stride, ObsCallBack<readdy::model::observables::Energy> callBack) const override;

    std::unique_ptr<readdy::model::observables::Virial>
    virial(Stride stride, ObsCallBack<readdy::model::observables::Virial> callBack) const override;

    std::unique_ptr<readdy::model::observables::HistogramAlongAxis>
    histogramAlongAxis(Stride stride, std::vector<scalar> binBorders, std::vector<std::string> typesToCount,
                       unsigned int axis, ObsCallBack<readdy::model::observables::HistogramAlongAxis> callBack) const override;

    std::unique_ptr<readdy::model::observables::NParticles>
    nParticles(Stride stride, std::vector<std::string> typesToCount,
               ObsCallBack<readdy::model::observables::NParticles> callback) const override;

    std::unique_ptr<readdy::model::observables::Forces>
    forces(Stride stride, std::vector<std::string> typesToCount,
           ObsCallBack<readdy::model::observables::Forces> callback) const override;

    std::unique_ptr<readdy::model::observables::Positions>
    positions(Stride stride, std::vector<std::string> typesToCount,
              ObsCallBack<readdy::model::observables::Positions> callback) const override;

    std::unique_ptr<readdy::model::observables::RadialDistribution>
    radialDistribution(Stride stride, std::vector<scalar> binBorders, std::vector<std::string> typeCountFrom,
                       std::vector<std::string> typeCountTo, scalar particleDensity,
                       ObsCallBack<readdy::model::observables::RadialDistribution> callback) const override;

    std::unique_ptr<readdy::model::observables::Particles>
    particles(Stride stride, ObsCallBack<readdy::model::observables::Particles> callback) const override;

    std::unique_ptr<readdy::model::observables::Reactions>
    reactions(Stride stride, ObsCallBack<readdy::model::observables::Reactions> callback) const override;

    std::unique_ptr<readdy::model::observables::ReactionCounts>
    reactionCounts(Stride stride, ObsCallBack<readdy::model::observables::ReactionCounts> callback) const override;
private:
    readdy::kernel::scpu::SCPUKernel *const kernel;
};

}
}
