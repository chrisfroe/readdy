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
 * @file NeighborList.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 08.09.16
 */

#include <set>
#include <sstream>
#include <future>

#include <readdy/common/numeric.h>
#include <readdy/common/thread/scoped_async.h>
#include <readdy/common/thread/barrier.h>
#include <readdy/common/thread/semaphore.h>
#include <readdy/common/Timer.h>

#include <readdy/kernel/cpu/model/CPUNeighborList.h>
#include <readdy/kernel/cpu/util/hilbert.h>
#include <readdy/common/thread/notification_barrier.h>
#include <readdy/common/Utils.h>
#include <readdy/kernel/cpu/util/config.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace model {

namespace thd = readdy::util::thread;

struct CPUNeighborList::Cell {
    std::vector<Cell *> neighbors{};
    double maximal_displacements[2];
    cell_index contiguous_index;
    bool enoughCells;


    // dirty flag indicating whether the cell and its neighboring cells have to be re-created
    bool dirty{false};

    Cell(cell_index i, cell_index j, cell_index k, const std::array<cell_index, 3> &nCells);

    void addNeighbor(Cell *cell);

    void checkDirty(skin_size_t skin);

    friend bool operator==(const Cell &lhs, const Cell &rhs);

    friend bool operator!=(const Cell &lhs, const Cell &rhs);

    void addParticleIndex(const particle_index idx) {
        particleIndices.push_back(idx);
    }

    const std::vector<particle_index> &particles() const {
        return particleIndices;
    }

    std::vector<particle_index> &particles() {
        return particleIndices;
    }

private:
    std::vector<particle_index> particleIndices{};
};

template<typename T, typename Dims>
T get_contiguous_index(T i, T j, T k, Dims I, Dims J) {
    return k + j * J + i * I * J;
};

class ParticleTravelledTooFarException : public std::runtime_error {
public:
    ParticleTravelledTooFarException() : runtime_error("") {}
};

const static std::vector<CPUNeighborList::neighbor_t> no_neighbors{};

CPUNeighborList::CPUNeighborList(const ctx_t *const context, data_t &data, readdy::util::thread::Config const *const config,
                           skin_size_t skin)
        : ctx(context), config(config), cells(std::vector<Cell>()), reorderSignal(std::make_unique<reorder_signal_t>()),
          simBoxSize(ctx->getBoxSize()), skin_size(skin), data(data), groupParticlesOnCreation(false) {}
// fixme also try with groupparticlesoncreation(true)

void CPUNeighborList::setupCells() {
    if (cells.empty()) {
        double maxCutoff = 0;
        for(const auto& entry : ctx->potentials().potentials_order2()) {
            for(const auto& potential : entry.second) {
                maxCutoff = maxCutoff < potential->getCutoffRadius() ? potential->getCutoffRadius() : maxCutoff;
            }
        }
        for (auto &&e : ctx->reactions().order2_flat()) {
            maxCutoff = maxCutoff < e->getEductDistance() ? e->getEductDistance() : maxCutoff;
        }
        CPUNeighborList::maxCutoff = maxCutoff;
        CPUNeighborList::maxCutoffPlusSkin = maxCutoff + skin_size;

        log::debug("got maxCutoff={}, skin={} => r_c + r_s = {}", maxCutoff, skin_size, maxCutoff + skin_size);

        if (maxCutoff > 0) {
            const auto desiredCoarseCellWidth = maxCutoffPlusSkin;

            for (unsigned short i = 0; i < 3; ++i) {
                nCells[i] = std::max(1u, 2 * static_cast<cell_index>(floor(simBoxSize[i] / desiredCoarseCellWidth)));
                if (nCells[i] == 0) nCells[i] = 1;
                cellSize[i] = simBoxSize[i] / static_cast<double>(nCells[i]);
            }
            log::debug("resulting cell size = {}", cellSize);

            for (cell_index i = 0; i < nCells[0]; ++i) {
                for (cell_index j = 0; j < nCells[1]; ++j) {
                    for (cell_index k = 0; k < nCells[2]; ++k) {
                        cells.push_back({i, j, k, nCells});
                        if(i % 2 == 0 && j % 2 == 0 && k % 2 == 0) {
                            // coarseCells.emplace_back();
                        }
                    }
                }
            }

            for (cell_index i = 0; i < nCells[0]; ++i) {
                for (cell_index j = 0; j < nCells[1]; ++j) {
                    for (cell_index k = 0; k < nCells[2]; ++k) {
                        setupNeighboringCells(i, j, k);
                    }
                }
            }
        } else {
            nCells = {{1, 1, 1}};
            cells.push_back({0, 0, 0, nCells});
            setupNeighboringCells(0, 0, 0);
        }
    }
    log::trace("got {} in total with ({}, {}, {}) in each direction", cells.size(), nCells[0], nCells[1], nCells[2]);
}

void
CPUNeighborList::setupNeighboringCells(const signed_cell_index i, const signed_cell_index j, const signed_cell_index k) {
    auto me = getCell(i, j, k);
    for (signed_cell_index _i = -2; _i < 3; ++_i) {
        for (signed_cell_index _j = -2; _j < 3; ++_j) {
            for (signed_cell_index _k = -2; _k < 3; ++_k) {
                // don't add me as neighbor to myself
                if (!(_i == 0 && _j == 0 && _k == 0)) {
                    me->addNeighbor(getCell(i + _i, j + _j, k + _k));
                }
            }
        }
    }
}

void CPUNeighborList::clear() {
    cells.clear();
    for (auto &list : data.neighbors) {
        list.clear();
    }
    initialSetup = true;
}

/**
 * Potentially marks a cell dirty. If there is a particle that has travelled too far, it will return false, otherwise
 * true. Travelling too far means:
 *  - with boxes of width .5(r_c + r_s): displacement > r_c + r_s (this is our case)
 *  - with boxes of width r_c + r_s: displacement > 2r_c + 2r_s
 *
 * @param cell the cell
 * @param data the underlying data
 * @param skin the skin width
 * @return false if everything is ok, true if a particle travelled too far
 */
bool markCell(CPUNeighborList::Cell &cell, const CPUNeighborList::data_t &data, const double r_c, const double r_s) {
    // fixme hier ist alles in ordnung, keine deactivated [zumindest ohne hilbert sorting]
    bool travelledTooFar = false;
    cell.maximal_displacements[0] = 0;
    cell.maximal_displacements[1] = 0;

    for (const auto idx : cell.particles()) {
        const auto &entry = data.entry_at(idx);
        if (!entry.is_deactivated()) {
            const double disp = entry.displacement;
            // fixme is this correct? probably not
            travelledTooFar = disp > r_c + r_s;
            if (travelledTooFar) break;
            if (disp > cell.maximal_displacements[0]) {
                cell.maximal_displacements[0] = disp;
            } else if (disp > cell.maximal_displacements[1]) {
                cell.maximal_displacements[1] = disp;
            }
        } else {
            log::critical("found deactivated particle in my cell, uh oh.");
        }
    }
    cell.checkDirty(r_s);
    return travelledTooFar;
}

bool cellOrNeighborDirty(const CPUNeighborList::Cell &cell) {
    if (cell.dirty) return true;
    for (const auto n : cell.neighbors) {
        if (n->dirty) return true;
    }
    return false;
}

void CPUNeighborList::fillCells() {
    if (maxCutoff > 0) {
        const auto c2 = maxCutoffPlusSkin * maxCutoffPlusSkin;
        auto d2 = ctx->getDistSquaredFun();

        if (initialSetup) {

            for (auto &cell : cells) {
                cell.particles().clear();
            }
            data.neighbors.clear();
            data.neighbors.resize(data.size());

            if (groupParticlesOnCreation) {
                if(!data.empty()) {
                    std::vector<unsigned int> hilbert_indices;
                    std::vector<std::size_t> indices(data.size());
                    {
                        const auto grainSize = data.size() / config->nThreads();
                        std::iota(indices.begin(), indices.end(), 0);
                        hilbert_indices.resize(data.size());

                        auto worker = [&indices, &hilbert_indices, this, grainSize](data_iter_t begin, data_iter_t end,
                                                                                    std::vector<unsigned int>::iterator hilberts_begin) {
                            {
                                auto it = begin;
                                auto hilbert_it = hilberts_begin;
                                for (; it != end; ++it, ++hilbert_it) {
                                    if (!it->is_deactivated()) {
                                        const auto &pos = it->position();
                                        const cell_index i = static_cast<const cell_index>(floor((pos[0] + .5 * simBoxSize[0]) / cellSize[0]));
                                        const cell_index j = static_cast<const cell_index>(floor((pos[1] + .5 * simBoxSize[1]) / cellSize[1]));
                                        const cell_index k = static_cast<const cell_index>(floor((pos[2] + .5 * simBoxSize[2]) / cellSize[2]));
                                        bitmask_t coords[3]{i,j,k};
                                        *hilbert_it = 1 + static_cast<unsigned int>(hilbert_c2i(3, CHAR_BIT, coords));
                                    } else {
                                        *hilbert_it = 0;
                                    }
                                }
                            }
                        };

                        {
                            std::vector<threading_model> threads;
                            threads.reserve(config->nThreads());
                            auto data_it = data.begin();
                            auto hilberts_it = hilbert_indices.begin();
                            for (std::size_t i = 0; i < config->nThreads() - 1; ++i) {
                                threads.emplace_back(worker, data_it, data_it + grainSize, hilberts_it);
                                data_it += grainSize;
                                hilberts_it += grainSize;
                            }
                            threads.emplace_back(worker, data_it, data.end(), hilberts_it);
                        }
                        {
                            std::sort(indices.begin(), indices.end(),
                                      [&hilbert_indices](std::size_t i, std::size_t j) -> bool {
                                          return hilbert_indices[i] < hilbert_indices[j];
                                      });
                        }
                    }

                    data.blanks_moved_to_front();
                    std::vector<std::size_t> inverseIndices(indices.size());
                    for(std::size_t i = 0; i < indices.size(); ++i) {
                        inverseIndices[indices[i]] = i;
                    }
                    reorderSignal->fire_signal(inverseIndices);
                    readdy::util::collections::reorder_destructive(inverseIndices.begin(), inverseIndices.end(), data.begin());
                }
                {
                    data_t::index_t idx = data.getNDeactivated();
                    for (auto it = data.begin() + data.getNDeactivated(); it != data.end(); ++it) {
                        auto cell = getCell(it->position());
                        if (cell) {
                            cell->addParticleIndex(idx);
                        }
                        ++idx;
                    }
                }
            } else {
                data_t::index_t idx = 0;
                for (const auto &entry : data) {
                    if (!entry.is_deactivated()) {
                        auto cell = getCell(entry.pos);
                        if (cell) {
                            cell->addParticleIndex(idx);
                        }
                    }
                    ++idx;
                }
            }

            using cell_it_t = decltype(cells.begin());
            const auto size = cells.size();
            const std::size_t grainSize = size / config->nThreads();

            auto worker = [this, c2, d2](cell_it_t cellsBegin, cell_it_t cellsEnd) -> void {
                for (cell_it_t _b = cellsBegin; _b != cellsEnd; ++_b) {
                    auto &cell = *_b;
                    setUpCell(cell, c2, d2);
                }
            };

            std::vector<threading_model> threads;
            threads.reserve(config->nThreads());

            auto it_cells = cells.begin();
            for (int i = 0; i < config->nThreads() - 1; ++i) {
                threads.emplace_back(worker, it_cells, it_cells + grainSize);
                it_cells += grainSize;
            }
            threads.emplace_back(worker, it_cells, cells.end());
        } else {
            auto dirtyCells = findDirtyCells();

            handleDirtyCells(dirtyCells);
            {
                // auto dirtyCells2 = findDirtyCells();
                // log::warn("dirty cells after: {}", dirtyCells2.size());
            }
        }
    }
}

void CPUNeighborList::create() {
    simBoxSize = ctx->getBoxSize();
    data.setFixPosFun(ctx->getFixPositionFun());
    {
        setupCells();
    }
    bool redoFillCells = false;
    try {
        fillCells();
    } catch (const ParticleTravelledTooFarException &) {
        initialSetup = true;
        redoFillCells = true;
        log::debug("A particle's displacement has been more than r_c + r_s = {} + {} = {}, which means that "
                                     "it might have left its cell linked-list cell. This should, if at all, only happen "
                                     "very rarely and triggers a complete rebuild of the neighbor list.",
                             maxCutoff, skin_size, maxCutoffPlusSkin);
    }
    if (redoFillCells) fillCells();
    initialSetup = false;
}

CPUNeighborList::Cell *CPUNeighborList::getCell(signed_cell_index i, signed_cell_index j, signed_cell_index k) {
    const auto &periodic = ctx->getPeriodicBoundary();
    if (periodic[0]) i = readdy::util::numeric::positive_modulo(i, nCells[0]);
    else if (i < 0 || i >= nCells[0]) return nullptr;
    if (periodic[1]) j = readdy::util::numeric::positive_modulo(j, nCells[1]);
    else if (j < 0 || j >= nCells[1]) return nullptr;
    if (periodic[2]) k = readdy::util::numeric::positive_modulo(k, nCells[2]);
    else if (k < 0 || k >= nCells[2]) return nullptr;
    const auto cix = get_contiguous_index(i, j, k, nCells[1], nCells[2]);
    if (cix < cells.size()) {
        return &cells.at(static_cast<cell_index>(cix));
    } else {
        log::critical("CPUNeighborList::getCell(nonconst): Requested cell ({},{},{})={}, but there are "
                                         "only {} cells.", i, j, k, cix, cells.size());
        throw std::runtime_error("tried to get cell index that was too large");
    }
}

void CPUNeighborList::remove(const particle_index idx) {
    auto cell = getCell(data.pos(idx));
    if (cell != nullptr) {
        {
            if(data.entry_at(idx).is_deactivated()) {
                log::critical("well, shit");
            }
            // remove from own cell
            bool erased = false;
            auto &particles = cell->particles();
            auto find_it = std::find(particles.begin(), particles.end(), idx);
            if (find_it != particles.end()) {
                particles.erase(find_it);
                erased = true;
            } else {
                for (auto neighbor_cell : cell->neighbors) {
                    auto &p2 = neighbor_cell->particles();
                    auto find_it2 = std::find(p2.begin(), p2.end(), idx);
                    if (find_it2 != p2.end()) {
                        p2.erase(find_it2);
                        erased = true;
                        break;
                    }
                }
            }
            if (!erased) {
                log::critical("hier sollte man aber nicht hinkommen");
            }
        }
        {
            // check up on my neighbors and remove myself from them
            for (auto &neighbor : neighbors(idx)) {
                if (data.entry_at(neighbor).deactivated) {
                    log::critical("!!!!!");
                }
                auto &neighbors_2nd = data.neighbors.at(neighbor);
                auto it = std::find(neighbors_2nd.begin(), neighbors_2nd.end(), idx);
                if (it != neighbors_2nd.end()) {
                    neighbors_2nd.erase(it);
                } else {
                    // todo log::debug("abcabc");
                }
            }
        }
    } else {
        log::critical("couldnt find cell");
    }
}

void CPUNeighborList::insert(const particle_index idx) {
    const auto &d2 = ctx->getDistSquaredFun();
    const auto &pos = data.pos(idx);
    const auto cutoffSquared = maxCutoffPlusSkin * maxCutoffPlusSkin;
    auto cell = getCell(pos);
    if (cell) {
        cell->addParticleIndex(idx);
        auto &myNeighbors = data.neighbors.at(idx);
        for (const auto pJ : cell->particles()) {
            if (idx != pJ) {
                const auto distSquared = d2(pos, data.pos(pJ));
                if (distSquared < cutoffSquared) {
                    myNeighbors.push_back(pJ);
                    data.neighbors.at(pJ).push_back(idx);
                }
            }
        }
        for (auto &neighboringCell : cell->neighbors) {
            for (const auto &pJ : neighboringCell->particles()) {
                const auto distSquared = d2(pos, data.pos(pJ));
                if (distSquared < cutoffSquared) {
                    myNeighbors.push_back(pJ);
                    data.neighbors.at(pJ).push_back(idx);
                }
            }
        }
    } else {
        //log::error("could not assign particle (index={}) to any cell!", data.getEntryIndex(idx));
    }
}

CPUNeighborList::Cell *CPUNeighborList::getCell(const readdy::model::Particle::pos_type &pos) {
    const cell_index i = static_cast<const cell_index>(floor((pos[0] + .5 * simBoxSize[0]) / cellSize[0]));
    const cell_index j = static_cast<const cell_index>(floor((pos[1] + .5 * simBoxSize[1]) / cellSize[1]));
    const cell_index k = static_cast<const cell_index>(floor((pos[2] + .5 * simBoxSize[2]) / cellSize[2]));
    return getCell(i, j, k);
}

void CPUNeighborList::updateData(CPUParticleData::update_t &&update) {
    std::unordered_set<CPUNeighborList::Cell *> dirties;
    if (maxCutoff > 0) {
        for (const auto &decayIndex : std::get<1>(update)) {
            // remove(p);
            if(!data.entry_at(decayIndex).is_deactivated()) {
                auto cell = getCell(data.pos(decayIndex));
                if(cell) {
                    {
                        bool found = false;
                        auto &particles = cell->particles();
                        auto find_it = std::find(particles.begin(), particles.end(), decayIndex);
                        if (find_it != particles.end()) {
                            found = true;
                        } else {
                            for (auto neighbor_cell : cell->neighbors) {
                                auto &p2 = neighbor_cell->particles();
                                auto find_it2 = std::find(p2.begin(), p2.end(), decayIndex);
                                if (find_it2 != p2.end()) {
                                    found = true;
                                    cell = neighbor_cell;
                                    break;
                                }
                            }
                        }
                        if(!found) {
                            log::error("did not find particle in any of the surrounding cells");
                        }
                    }

                    cell->dirty = true;
                    for(auto neigh : cell->neighbors) {
                        neigh->dirty = true;
                    }
                } else {
                    log::critical("did not find cell for particle");
                }
            } else {
                log::critical("tried removing deactivated particle");
            }
        }
        for(auto& cell : cells) {
            if(cellOrNeighborDirty(cell)) {
                dirties.insert(&cell);
            }
        }
    }
    auto newEntries = data.update(std::move(update));
    if(!dirties.empty()) {
        log::error("dirties not empty: {} / {}", dirties.size(), cells.size());
        handleDirtyCells(std::move(dirties));
    }
    if (maxCutoff > 0) {
        for (const auto p : newEntries) {
            insert(p);
        }
    }
}

const std::vector<CPUNeighborList::neighbor_t> &CPUNeighborList::neighbors(CPUNeighborList::particle_index const entry) const {
    if (maxCutoff > 0) {
        return data.neighbors.at(entry);
    }
    return no_neighbors;
}

const CPUNeighborList::Cell *const CPUNeighborList::getCell(const readdy::model::Particle::pos_type &pos) const {
    const cell_index i = static_cast<const cell_index>(floor((pos[0] + .5 * simBoxSize[0]) / cellSize[0]));
    const cell_index j = static_cast<const cell_index>(floor((pos[1] + .5 * simBoxSize[1]) / cellSize[1]));
    const cell_index k = static_cast<const cell_index>(floor((pos[2] + .5 * simBoxSize[2]) / cellSize[2]));
    return getCell(i, j, k);
}

const CPUNeighborList::Cell *const CPUNeighborList::getCell(CPUNeighborList::signed_cell_index i,
                                                      CPUNeighborList::signed_cell_index j,
                                                      CPUNeighborList::signed_cell_index k) const {

    const auto &periodic = ctx->getPeriodicBoundary();
    if (periodic[0]) i = readdy::util::numeric::positive_modulo(i, nCells[0]);
    else if (i < 0 || i >= nCells[0]) return nullptr;
    if (periodic[1]) j = readdy::util::numeric::positive_modulo(j, nCells[1]);
    else if (j < 0 || j >= nCells[1]) return nullptr;
    if (periodic[2]) k = readdy::util::numeric::positive_modulo(k, nCells[2]);
    else if (k < 0 || k >= nCells[2]) return nullptr;
    const auto cix = get_contiguous_index(i, j, k, nCells[1], nCells[2]);
    if (cix < cells.size()) {
        return &cells.at(static_cast<cell_index>(cix));
    } else {
        log::critical("CPUNeighborList::getCell(const): Requested cell ({},{},{})={}, but there are "
                                         "only {} cells.", i, j, k, cix, cells.size());
        throw std::out_of_range("tried to access an invalid cell");
    }
}

const std::vector<CPUNeighborList::neighbor_t> &CPUNeighborList::find_neighbors(particle_index const entry) const {
    if (maxCutoff > 0 && entry < data.neighbors.size()) {
        return data.neighbors.at(entry);
    }
    return no_neighbors;
}

CPUNeighborList::~CPUNeighborList() = default;

void CPUNeighborList::Cell::addNeighbor(CPUNeighborList::Cell *cell) {
    if (cell && cell->contiguous_index != contiguous_index
        && (enoughCells || std::find(neighbors.begin(), neighbors.end(), cell) == neighbors.end())) {
        neighbors.push_back(cell);
    }
}

CPUNeighborList::Cell::Cell(cell_index i, cell_index j, cell_index k,
                         const std::array<CPUNeighborList::cell_index, 3> &nCells)
        : contiguous_index(get_contiguous_index(i, j, k, nCells[1], nCells[2])),
          enoughCells(nCells[0] >= 5 && nCells[1] >= 5 && nCells[2] >= 5), maximal_displacements{0, 0} {
}

bool operator==(const CPUNeighborList::Cell &lhs, const CPUNeighborList::Cell &rhs) {
    return lhs.contiguous_index == rhs.contiguous_index;
}

bool operator!=(const CPUNeighborList::Cell &lhs, const CPUNeighborList::Cell &rhs) {
    return !(lhs == rhs);
}

void CPUNeighborList::Cell::checkDirty(skin_size_t skin) {
    dirty = maximal_displacements[0] + maximal_displacements[1] > .5*skin;
}

CPUNeighborList::hilbert_index_t CPUNeighborList::getHilbertIndex(std::size_t i, std::size_t j, std::size_t k) const {
    bitmask_t coords[3]{i, j, k};
    return static_cast<unsigned int>(hilbert_c2i(3, CHAR_BIT, coords));
}

void CPUNeighborList::displace(data_t::Entry &entry, const data_t::particle_type::pos_type &delta) {
    data.displace(entry, delta);
}

void CPUNeighborList::displace(data_iter_t iter, const readdy::model::Vec3 &vec) {
    auto &entry = *iter;
    displace(entry, vec);
}

void CPUNeighborList::displace(data_t::index_t idx, const data_t::particle_type::pos_type &delta) {
    auto &entry = data.entry_at(idx);
    displace(entry, delta);
}

void CPUNeighborList::setPosition(data_t::index_t idx, data_t::particle_type::pos_type &&newPosition) {
    auto &entry = data.entry_at(idx);
    const auto delta = newPosition - entry.position();
    displace(entry, delta);
}

// todo: replace uset w/ vector, check cell bounds with ineq operators
std::unordered_set<CPUNeighborList::Cell *> CPUNeighborList::findDirtyCells() {
    using it_type = decltype(cells.begin());

    using future_t = std::future<std::vector<Cell *>>;
    using promise_t = std::promise<std::vector<Cell *>>;

    std::unordered_set<CPUNeighborList::Cell *> result;


    auto worker = [this](it_type begin, it_type end, promise_t update, const thd::barrier &barrier,
                         std::atomic<bool> &interrupt) {
        // interrupt if there is a particle that travelled too far (r_c + r_s)
        bool shouldInterrupt = false;
        for (it_type it = begin; it != end; ++it) {
            shouldInterrupt |= markCell(*it, data, maxCutoff, skin_size);
        }
        if (shouldInterrupt) {
            interrupt = shouldInterrupt;
        }
        barrier.wait();
        std::vector<Cell *> foundDirtyCells;
        if (!interrupt.load()) {
            for (it_type it = begin; it != end; ++it) {
                if (cellOrNeighborDirty(*it)) {
                    foundDirtyCells.push_back(&*it);
                }
            }
        }
        update.set_value(std::move(foundDirtyCells));
    };

    std::atomic<bool> interrupt(false);
    thd::barrier barrier(config->nThreads());

    std::vector<future_t> dirtyCells;
    dirtyCells.reserve(config->nThreads());
    {
        auto it_cells = cells.begin();
        std::vector<threading_model> threads;
        const std::size_t grainSize = cells.size() / config->nThreads();
        for (int i = 0; i < config->nThreads() - 1; ++i) {
            promise_t promise;
            dirtyCells.push_back(promise.get_future());
            threads.emplace_back(worker, it_cells, it_cells + grainSize, std::move(promise),
                                 std::cref(barrier), std::ref(interrupt));
            it_cells += grainSize;
        }
        promise_t promise;
        dirtyCells.push_back(promise.get_future());
        threads.emplace_back(worker, it_cells, it_cells + grainSize, std::move(promise),
                             std::cref(barrier), std::ref(interrupt));
    }
    if (interrupt.load()) {
        throw ParticleTravelledTooFarException();
    }
    for (auto &&dirties : dirtyCells) {
        auto v = std::move(dirties.get());
        std::copy(v.begin(), v.end(), std::inserter(result, result.end()));
    }
    return std::move(result);
}

void CPUNeighborList::setUpCell(CPUNeighborList::Cell &cell, const double cutoffSquared, const ctx_t::dist_squared_fun &d2) {
    cell.dirty = false;
    cell.maximal_displacements[0] = 0;
    cell.maximal_displacements[1] = 0;
    for (const auto &pI : cell.particles()) {
        auto &entry_i = data.entry_at(pI);
        entry_i.displacement = 0;
        auto &neighbors_i = data.neighbors_at(pI);
        neighbors_i.clear();
        for (const auto &pJ : cell.particles()) {
            if (pI != pJ) {
                const auto distSquared = d2(entry_i.position(), data.pos(pJ));
                if (distSquared < cutoffSquared) {
                    neighbors_i.push_back(pJ);
                }
            }
        }
        for (const auto &neighboringCell : cell.neighbors) {
            for (const auto pJ : neighboringCell->particles()) {
                const auto distSquared = d2(entry_i.position(), data.pos(pJ));
                if (distSquared < cutoffSquared) {
                    neighbors_i.push_back(pJ);
                }
            }
        }

    }
}

int CPUNeighborList::getCellIndex(const readdy::model::Particle::pos_type &pos) const {
    signed_cell_index i = static_cast<const signed_cell_index>(floor((pos[0] + .5 * simBoxSize[0]) / cellSize[0]));
    signed_cell_index j = static_cast<const signed_cell_index>(floor((pos[1] + .5 * simBoxSize[1]) / cellSize[1]));
    signed_cell_index k = static_cast<const signed_cell_index>(floor((pos[2] + .5 * simBoxSize[2]) / cellSize[2]));
    const auto &periodic = ctx->getPeriodicBoundary();
    if (periodic[0]) i = readdy::util::numeric::positive_modulo(i, nCells[0]);
    else if (i < 0 || i >= nCells[0]) return -1;
    if (periodic[1]) j = readdy::util::numeric::positive_modulo(j, nCells[1]);
    else if (j < 0 || j >= nCells[1]) return -1;
    if (periodic[2]) k = readdy::util::numeric::positive_modulo(k, nCells[2]);
    else if (k < 0 || k >= nCells[2]) return -1;
    return get_contiguous_index(i, j, k, nCells[1], nCells[2]);
}

void CPUNeighborList::setSkinSize(CPUNeighborList::skin_size_t skin_size) {
    initialSetup = true;
    cells.clear();
    CPUNeighborList::skin_size = skin_size;
}

bool CPUNeighborList::isInCell(const CPUNeighborList::Cell &cell, const data_t::particle_type::pos_type &pos) const {
    return hash_pos(pos) == cell.contiguous_index;
}

std::size_t CPUNeighborList::hash_pos(const data_t::particle_type::pos_type &pos) const {
    const cell_index i = static_cast<const cell_index>(floor((pos[0] + .5 * simBoxSize[0]) / cellSize[0]));
    const cell_index j = static_cast<const cell_index>(floor((pos[1] + .5 * simBoxSize[1]) / cellSize[1]));
    const cell_index k = static_cast<const cell_index>(floor((pos[2] + .5 * simBoxSize[2]) / cellSize[2]));
    return get_contiguous_index(i, j, k, nCells[1], nCells[2]);
}

void CPUNeighborList::setGroupParticlesOnCreation(bool groupParticlesOnCreation) {
    CPUNeighborList::groupParticlesOnCreation = groupParticlesOnCreation;
}

std::tuple<CPUNeighborList::cell_index, CPUNeighborList::cell_index, CPUNeighborList::cell_index> CPUNeighborList::mapPositionToCell(
        const readdy::model::Particle::pos_type &pos) const {
    const cell_index i = static_cast<const cell_index>(floor((pos[0] + .5 * simBoxSize[0]) / cellSize[0]));
    const cell_index j = static_cast<const cell_index>(floor((pos[1] + .5 * simBoxSize[1]) / cellSize[1]));
    const cell_index k = static_cast<const cell_index>(floor((pos[2] + .5 * simBoxSize[2]) / cellSize[2]));
    return std::make_tuple(i, j, k);
}

readdy::signals::scoped_connection CPUNeighborList::registerReorderEventListener(const reorder_signal_t::slot_type &slot) {
    return reorderSignal->connect_scoped(slot);
}

void CPUNeighborList::handleDirtyCells(std::unordered_set<Cell *> dirtyCells) {
    const auto c2 = maxCutoffPlusSkin * maxCutoffPlusSkin;
    auto d2 = ctx->getDistSquaredFun();
    const double fraction = .9;
    if (dirtyCells.size() > cells.size() * fraction) {
        initialSetup = true;
        // log::debug("had more than {}% ({} / {}) dirty cells, recreate neighbor list", fraction * 100, dirtyCells.size(), cells.size());
        setupCells();
        fillCells();
        initialSetup = false;
    } else {

        const std::size_t grainSize = dirtyCells.size() / config->nThreads();

        using cell_it_t = decltype(dirtyCells.begin());

        auto worker = [this, c2, d2, grainSize](cell_it_t begin, cell_it_t end, const thd::barrier &barrier) {
            // first gather updates
            std::vector<std::vector<CPUNeighborList::particle_index>> cellUpdates;
            cellUpdates.reserve(grainSize);
            {
                for (cell_it_t _b = begin; _b != end; ++_b) {
                    auto cell = *_b;
                    cellUpdates.push_back({});
                    auto &updatedIndices = cellUpdates.back();
                    updatedIndices.reserve(cell->particles().size());
                    for (const auto idx : cell->particles()) {
                        const auto &entry = data.entry_at(idx);
                        if (!entry.is_deactivated() &&
                            cell->contiguous_index == getCellIndex(entry.position())) {
                            updatedIndices.push_back(idx);
                        }
                    }
                    for (const auto neigh : cell->neighbors) {
                        if (cell->dirty || neigh->dirty) {
                            for (const auto idx : neigh->particles()) {
                                const auto &entry = data.entry_at(idx);
                                if (!entry.is_deactivated()
                                    && cell->contiguous_index == getCellIndex(entry.position())) {
                                    updatedIndices.push_back(idx);
                                }
                            }
                        }
                    }
                }
            }

            // wait until all threads have their updates gathered
            barrier.wait();
            // apply updates
            {
                auto it = begin;
                for (auto &&update : cellUpdates) {
                    auto cell = *it;
                    cell->particles() = std::move(update);
                    ++it;
                }
            }
            // wait until all threads have applied their updates
            barrier.wait();
            // update neighbor list in respective cells
            {
                for (cell_it_t _b = begin; _b != end; ++_b) {
                    auto cell = *_b;
                    setUpCell(*cell, c2, d2);
                }
            }
        };

        // readdy::util::Timer t("        update dirty cells");
        thd::barrier b(config->nThreads());

        std::vector<threading_model> threads;
        threads.reserve(config->nThreads());

        log::warn("got dirty cells {} vs total cells {}", dirtyCells.size(), cells.size());

        {
            auto it_cells = dirtyCells.begin();
            for (int i = 0; i < config->nThreads() - 1; ++i) {
                auto advanced = std::next(it_cells, grainSize);
                threads.emplace_back(worker, it_cells, advanced, std::cref(b));
                it_cells = advanced;
            }
            threads.emplace_back(worker, it_cells, dirtyCells.end(), std::cref(b));
        }
    }
}


}
}
}
}
