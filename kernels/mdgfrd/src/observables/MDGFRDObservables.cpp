/********************************************************************
 * Copyright © 2016 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * This file is part of ReaDDy.                                     *
 *                                                                  *
 * ReaDDy is free software: you can redistribute it and/or modify   *
 * it under the terms of the GNU Lesser General Public License as   *
 * published by the Free Software Foundation, either version 3 of   *
 * the License, or (at your option) any later version.              *
 *                                                                  *
 * This program is distributed in the hope that it will be useful,  *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of   *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
 * GNU Lesser General Public License for more details.              *
 *                                                                  *
 * You should have received a copy of the GNU Lesser General        *
 * Public License along with this program. If not, see              *
 * <http://www.gnu.org/licenses/>.                                  *
 ********************************************************************/


/**
 * << detailed description >>
 *
 * @file Observables.cpp
 * @brief << brief description >>
 * @author clonker
 * @author luigisbailo
 * @date 15.02.18
 */

#include <future>

#include <readdy/common/thread/scoped_async.h>

#include <readdy/kernel/mdgfrd/observables/MDGFRDObservables.h>
#include <readdy/kernel/mdgfrd/MDGFRDKernel.h>

namespace readdy {
namespace kernel {
namespace mdgfrd {
namespace observables {

MDGFRDPositions::MDGFRDPositions(MDGFRDKernel *const kernel, unsigned int stride,
                           const std::vector<std::string> &typesToCount) :
        readdy::model::observables::Positions(kernel, stride, typesToCount), kernel(kernel) {}

void MDGFRDPositions::evaluate() {
    result.clear();
    auto &stateModel = kernel->getMDGFRDKernelStateModel();
    const auto &pd = stateModel.getParticleData();
    if (typesToCount.empty()) {
        result = stateModel.getParticlePositions();
    } else {
        for (const auto &e : *stateModel.getParticleData()) {
            if (!e.deactivated &&
                std::find(typesToCount.begin(), typesToCount.end(), e.type) != typesToCount.end()) {
                result.push_back(e.pos);
            }
        }
    }
}

MDGFRDHistogramAlongAxis::MDGFRDHistogramAlongAxis(MDGFRDKernel *const kernel, unsigned int stride,
                                             const std::vector<scalar> &binBorders,
                                             const std::vector<std::string> &typesToCount, unsigned int axis)
        : readdy::model::observables::HistogramAlongAxis(kernel, stride, binBorders, typesToCount, axis),
          kernel(kernel) {
    size = result.size();
}

void MDGFRDHistogramAlongAxis::evaluate() {
    using Iter = readdy::kernel::mdgfrd::MDGFRDStateModel::data_type::const_iterator;

    std::fill(result.begin(), result.end(), 0);

    const auto binBorders = this->binBorders;
    const auto typesToCount = this->typesToCount;
    const auto resultSize = result.size();
    const auto axis = this->axis;
    const auto data = kernel->getMDGFRDKernelStateModel().getParticleData();

    std::vector<std::future<result_type>> updates;
    updates.reserve(kernel->getNThreads());
    auto worker = [binBorders, typesToCount, resultSize, data, axis](std::size_t, Iter from, Iter to, std::promise<result_type>& update) {
        result_type resultUpdate;
        resultUpdate.resize(resultSize);

        for (auto it = from; it != to; ++it) {
            if (!it->deactivated && typesToCount.find(it->type) != typesToCount.end()) {
                auto upperBound = std::upper_bound(binBorders.begin(), binBorders.end(), it->pos[axis]);
                if (upperBound != binBorders.end()) {
                    auto binBordersIdx = upperBound - binBorders.begin();
                    if (binBordersIdx >= 1 && binBordersIdx < resultSize) {
                        ++resultUpdate[binBordersIdx - 1];
                    }
                }
            }
        }

        update.set_value(std::move(resultUpdate));
    };

    {
        const std::size_t grainSize = data->size() / kernel->getNThreads();

        std::vector<std::promise<result_type>> promises;
        promises.resize(kernel->getNThreads());

        auto &pool = kernel->pool();

        auto workIter = data->cbegin();
        for (unsigned int i = 0; i < kernel->getNThreads() - 1; ++i) {
            updates.push_back(promises.at(i).get_future());
            pool.push(worker, workIter, workIter + grainSize, std::ref(promises.at(i)));
            workIter += grainSize;
        }
        auto& promise = promises.back();
        updates.push_back(promise.get_future());
        pool.push(worker, workIter, data->cend(), std::ref(promise));


        for (auto &update : updates) {
            auto vec = std::move(update.get());
            auto it1 = vec.begin();
            auto it2 = result.begin();
            for (; it1 != vec.end(); ++it1, ++it2) {
                *it2 += *it1;
            }
        }

    }
}


MDGFRDNParticles::MDGFRDNParticles(MDGFRDKernel *const kernel, unsigned int stride, std::vector<std::string> typesToCount)
        : readdy::model::observables::NParticles(kernel, stride, std::move(typesToCount)),
          kernel(kernel) {}

void MDGFRDNParticles::evaluate() {
    std::vector<unsigned long> resultVec = {};
    if (typesToCount.empty()) {
        const auto data = kernel->getMDGFRDKernelStateModel().getParticleData();
        resultVec.push_back(data->size() - data->getNDeactivated());
    } else {
        resultVec.resize(typesToCount.size());
        const auto &pd = kernel->getMDGFRDKernelStateModel().getParticleData();
        for (const auto &e : *pd) {
            if (!e.deactivated) {
                auto typeIt = std::find(typesToCount.begin(), typesToCount.end(), e.type);
                if (typeIt != typesToCount.end()) {
                    ++resultVec[typeIt - typesToCount.begin()];
                }
            }
        }
    }
    result = std::move(resultVec);
}

MDGFRDForces::MDGFRDForces(MDGFRDKernel *const kernel, unsigned int stride, std::vector<std::string> typesToCount) :
        readdy::model::observables::Forces(kernel, stride, std::move(typesToCount)),
        kernel(kernel) {}

void MDGFRDForces::evaluate() {
    result.clear();
    const auto &pd = kernel->getMDGFRDKernelStateModel().getParticleData();
    if (typesToCount.empty()) {
        result.reserve(pd->size());
    }
    for (const auto &e : *pd) {
        if (!e.deactivated) {
            if (typesToCount.empty()) {
                result.push_back(e.force);
            } else {
                for (auto countedParticleType : typesToCount) {
                    if (e.type == countedParticleType) {
                        result.push_back(e.force);
                        break;
                    }
                }
            }
        }
    }
}


MDGFRDParticles::MDGFRDParticles(MDGFRDKernel *const kernel, unsigned int stride)
        : readdy::model::observables::Particles(kernel, stride), kernel(kernel) {}

void MDGFRDParticles::evaluate() {
    auto &resultTypes = std::get<0>(result);
    auto &resultIds = std::get<1>(result);
    auto &resultPositions = std::get<2>(result);
    resultTypes.clear();
    resultIds.clear();
    resultPositions.clear();
    const auto &particleData = kernel->getMDGFRDKernelStateModel().getParticleData();
    resultTypes.reserve(particleData->size());
    resultIds.reserve(particleData->size());
    resultPositions.reserve(particleData->size());
    for (const auto &entry : *particleData) {
        if (!entry.deactivated) {
            resultTypes.push_back(entry.type);
            resultIds.push_back(entry.id);
            resultPositions.push_back(entry.pos);
        }
    }
}

MDGFRDReactions::MDGFRDReactions(MDGFRDKernel *const kernel, unsigned int stride)
        : Reactions(kernel, stride), kernel(kernel) {}

void MDGFRDReactions::evaluate() {
    const auto& model = kernel->getMDGFRDKernelStateModel();
    const auto& records = model.reactionRecords();
    result = records;
}

MDGFRDReactionCounts::MDGFRDReactionCounts(MDGFRDKernel *const kernel, unsigned int stride)
        : ReactionCounts(kernel, stride), kernel(kernel) {}

void MDGFRDReactionCounts::evaluate() {
    result = kernel->getMDGFRDKernelStateModel().reactionCounts();
}

MDGFRDVirial::MDGFRDVirial(MDGFRDKernel *kernel, stride_type stride) : Virial(kernel, stride), kernel(kernel) {}

void MDGFRDVirial::evaluate() {
    result = kernel->getMDGFRDKernelStateModel().virial();
}


}
}
}
}
