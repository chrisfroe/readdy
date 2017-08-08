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
 * Base class for all possible types of reaction. Currently:
 *  - Conversion: A -> B
 *  - Enzymatic: A + C -> B + C where C is a catalyst
 *  - Fission: A -> B + C
 *  - Fusion: A + B -> C
 *
 * @file Reactions.h
 * @brief Reaction base class.
 * @author clonker
 * @date 17.06.16
 */

#pragma once
#include <string>
#include <ostream>
#include <utility>
#include <spdlog/fmt/ostr.h>
#include <readdy/model/Particle.h>
#include <readdy/model/RandomProvider.h>
#include <readdy/common/make_unique.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(reactions)

enum class ReactionType { Conversion, Fusion, Fission, Enzymatic, Decay };

std::ostream& operator<<(std::ostream& os, const ReactionType& reactionType);

template<unsigned int N_EDUCTS>
class Reaction {
protected:
    static short counter;
    using particle_type_type = readdy::model::Particle::type_type;
public:

    using rnd_normal = std::function<Vec3(const scalar, const scalar)>;
    // static constexpr unsigned int n_educts = N_EDUCTS;

    Reaction(std::string name, const scalar rate, const scalar eductDistance,
             const scalar productDistance, const unsigned int n_products) :
            name(std::move(name)),
            id(counter++),
            rate(rate),
            eductDistance(eductDistance),
            eductDistanceSquared(eductDistance * eductDistance),
            productDistance(productDistance),
            _n_products(n_products) {}

    Reaction(const Reaction&) = default;
    Reaction& operator=(const Reaction&) = default;
    Reaction(Reaction&&) = default;
    Reaction& operator=(Reaction&&) = default;

    virtual ~Reaction() = default;

    virtual const ReactionType getType() = 0;

    const std::string &getName() const {
        return name;
    }

    const short getId() const {
        return id;
    }

    const scalar getRate() const {
        return rate;
    }

    const unsigned int getNEducts() const {
        return _n_educts;
    }

    const unsigned int getNProducts() const {
        return _n_products;
    }

    const scalar getEductDistance() const {
        return eductDistance;
    }

    const scalar getEductDistanceSquared() const {
        return eductDistanceSquared;
    }

    const scalar getProductDistance() const {
        return productDistance;
    }

    friend std::ostream &operator<<(std::ostream &os, const Reaction &reaction) {
        os << "Reaction(\"" << reaction.name << "\", N_Educts=" << reaction._n_educts << ", N_Products="
           << reaction._n_products << ", (";
        for (unsigned int i = 0; i < reaction._n_educts; i++) {
            if (i > 0) os << ",";
            os << reaction.educts[i];
        }
        os << ") -> (";
        for (unsigned int i = 0; i < reaction._n_products; i++) {
            if (i > 0) os << ",";
            os << reaction.products[i];
        }
        os << "), rate=" << reaction.rate << ", eductDist=" << reaction.eductDistance << ", prodDist="
           << reaction.productDistance << ")";
        return os;
    }

    const std::array<particle_type_type, N_EDUCTS> &getEducts() const {
        return educts;
    }

    const std::array<particle_type_type, 2> &getProducts() const {
        return products;
    }

    const scalar getWeight1() const {
        return weight1;
    }

    const scalar getWeight2() const {
        return weight2;
    }


protected:
    const unsigned int _n_educts = N_EDUCTS;
    const unsigned int _n_products;
    std::array<particle_type_type, N_EDUCTS> educts;
    std::array<particle_type_type, 2> products {{0, 0}};
    const std::string name;
    const short id;
    const scalar rate;
    const scalar eductDistance, eductDistanceSquared;
    const scalar productDistance;

    scalar weight1 = .5, weight2 = .5;
};

NAMESPACE_END(reactions)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
