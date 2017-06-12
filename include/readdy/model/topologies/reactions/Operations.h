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
 * @file Operations.h
 * @brief << brief description >>
 * @author clonker
 * @date 13.04.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <memory>
#include <readdy/common/macros.h>
#include "TopologyReactionAction.h"
#include "TopologyReactionActionFactory.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(top)
class GraphTopology;
NAMESPACE_BEGIN(reactions)
NAMESPACE_BEGIN(op)

class Operation {
public:
    using Ref = std::shared_ptr<Operation>;
    using factory_ref = const actions::TopologyReactionActionFactory *const;
    using topology_ref = GraphTopology* const;
    using graph_t = actions::TopologyReactionAction::graph_t;
    using action_ptr = std::unique_ptr<actions::TopologyReactionAction>;

    using label_edge = graph_t::label_edge;
    using label_vertex = graph_t::label;
    virtual action_ptr create_action(topology_ref topology, factory_ref factory) const = 0;
};

class ChangeParticleType : public Operation {
public:
    ChangeParticleType(const label_vertex &vertex, particle_type_type type_to);

    virtual action_ptr create_action(topology_ref topology, factory_ref factory) const override;
private:
    label_vertex label_vertex_;
    particle_type_type type_to_;
};

NAMESPACE_END(op)
NAMESPACE_END(reactions)
NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
