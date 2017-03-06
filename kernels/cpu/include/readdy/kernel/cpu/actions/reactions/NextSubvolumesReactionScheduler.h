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
 * @file SingleCPUNextSubvolumesReactionScheduler.h
 * @brief << brief description >>
 * @author clonker
 * @date 09.09.16
 */

#ifndef READDY_KERNEL_CPU_NEXTSUBVOLUMESREACTIONSCHEDULER_H
#define READDY_KERNEL_CPU_NEXTSUBVOLUMESREACTIONSCHEDULER_H

#include <readdy/model/actions/Actions.h>
#include "readdy/kernel/cpu/CPUKernel.h"

namespace readdy {
namespace kernel {
namespace cpu {
namespace actions {
namespace reactions {

class CPUNextSubvolumes : public readdy::model::actions::reactions::NextSubvolumes {
using cell_index_t = unsigned int;
using signed_cell_index_t = typename std::make_signed<cell_index_t>::type;
    using super = readdy::model::actions::reactions::NextSubvolumes;
public:
    CPUNextSubvolumes(const CPUKernel *const kernel, double timeStep);
    ~CPUNextSubvolumes();

    virtual void perform() override;

    double getMaxReactionRadius() const;
private:
    struct ReactionEvent;
    struct GridCell;

    CPUKernel const* const kernel;

    // sets up a grid cell (rate, timestamp, next event)
    void setUpCell(GridCell& cell);
    // sets up the computational grid with spacing >= max(reaction radii)
    void setUpGrid();
    // assigns particles to the computational grid
    void assignParticles();
    // schedules the reactions
    void setUpEventQueue();
    // evaluates the collected events
    void evaluateReactions();
    // sets up the neighbor-linked-list structure
    void setUpNeighbors(GridCell& cell);
    GridCell * getCell(const readdy::model::Vec3& particlePosition);
    // fetches a cell at (i,j,k)
    GridCell * getCell(signed_cell_index_t i, signed_cell_index_t j, signed_cell_index_t k);

    // array holding the number of cells in each spatial direction
    std::array<unsigned int, 3> nCells;
    // size of each box
    readdy::model::Vec3 cellSize;

    std::vector<GridCell> cells;
    std::vector<GridCell*> eventQueue;
};

}
}
}
}
}

#endif //READDY_KERNEL_CPU_NEXTSUBVOLUMESREACTIONSCHEDULER_H