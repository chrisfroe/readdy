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
 * @file PositionsObservable.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 13.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/model/observables/Positions.h>
#include <readdy/model/Kernel.h>
#include <readdy/io/DataSet.h>
#include <readdy/model/observables/io/Types.h>

namespace readdy {
namespace model {
namespace observables {

struct Positions::Impl {
    using writer_t = io::DataSet<Vec3, true>;
    std::unique_ptr <writer_t> writer;
};

Positions::Positions(Kernel *const kernel, unsigned int stride,
                     std::vector <std::string> typesToCount) :
        Positions(kernel, stride,
                  _internal::util::transformTypes2(typesToCount, kernel->getKernelContext())) {}

Positions::Positions(Kernel *const kernel, unsigned int stride,
                     std::vector<unsigned int> typesToCount) :
        Observable(kernel, stride), typesToCount(typesToCount), pimpl(std::make_unique<Impl>()) {}

void Positions::append() {
    std::vector <Vec3> podVec(result.begin(), result.end());
    pimpl->writer->append({1}, &podVec);
}

Positions::Positions(Kernel *const kernel, unsigned int stride) : Observable(kernel, stride) {}

void Positions::initializeDataSet(io::File &file, const std::string &dataSetName, unsigned int flushStride) {
    if (!pimpl->writer) {
        std::vector <readdy::io::h5::dims_t> fs = {flushStride};
        std::vector <readdy::io::h5::dims_t> dims = {readdy::io::h5::UNLIMITED_DIMS};
        auto dataSet = std::make_unique < io::DataSet < Vec3,
        true >> (
                dataSetName, file.createGroup(util::OBSERVABLES_GROUP_PATH), fs, dims,
                        util::Vec3MemoryType(), util::Vec3FileType()
        );
        pimpl->writer = std::move(dataSet);
    }
}

void Positions::flush() {
    pimpl->writer->flush();
}

Positions::~Positions() = default;
}
}
}