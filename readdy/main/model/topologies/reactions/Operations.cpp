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
 * @file Operations.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 05.04.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/model/topologies/reactions/Operation.h>
#include <readdy/model/topologies/GraphTopology.h>

namespace readdy {
namespace model {
namespace top {
namespace reactions {
namespace op {

Operation::Operation(GraphTopology *const topology) : topology(topology){ }

ChangeParticleType::ChangeParticleType(GraphTopology *const topology, const graph::vertex_ref &v,
                                       const particle_type_type &type_to)
        : Operation(topology), vertex(v), type_to(type_to), previous_type(type_to){}


AddEdge::AddEdge(GraphTopology *const topology, const graph::edge& edge) : Operation(topology), edge(edge) {}


void AddEdge::execute() {
    topology->graph().addEdge(edge);
}

void AddEdge::undo() {
    topology->graph().removeEdge(edge);
}


RemoveEdge::RemoveEdge(GraphTopology *const topology, const graph::edge& edge) : Operation(topology), edge(edge) {}

void RemoveEdge::execute() {
    topology->graph().removeEdge(edge);
}

void RemoveEdge::undo() {
    topology->graph().addEdge(edge);
}



}
}
}
}
}
