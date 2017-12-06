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
 * @file ParticleTypeRegistry.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 29.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/model/ParticleTypeRegistry.h>
#include <readdy/model/Utils.h>

namespace readdy {
namespace model {


ParticleTypeInfo::ParticleTypeInfo(const std::string &name, const scalar diffusionConstant,
                                   const particle_flavor flavor, const Particle::type_type typeId)
        : name(name), diffusionConstant(diffusionConstant), flavor(flavor), typeId(typeId) {}


const ParticleTypeInfo &ParticleTypeRegistry::infoOf(const std::string &name) const {
    return infoOf(_idOf(name));
}

const ParticleTypeInfo &ParticleTypeRegistry::infoOf(Particle::type_type type) const {
    return particle_info_.at(type);
}

const ParticleTypeRegistry::type_map &ParticleTypeRegistry::typeMapping() const {
    return type_mapping_;
}

std::string ParticleTypeRegistry::nameOf(particle_type_type id) const {
    for (const auto &e : type_mapping_) {
        if (e.second == id) return e.first;
    }
    return "";
}

std::vector<particle_type_type> ParticleTypeRegistry::typesFlat() const {
    std::vector<particle_type_type> v;
    for (auto &&entry : type_mapping_) {
        v.push_back(entry.second);
    }
    return v;
}

scalar ParticleTypeRegistry::diffusionConstantOf(particle_type_type particleType) const {
    return particle_info_.at(particleType).diffusionConstant;
}

scalar ParticleTypeRegistry::diffusionConstantOf(const std::string &particleType) const {
    return diffusionConstantOf(idOf(particleType));
}

void ParticleTypeRegistry::add(const std::string &name, const scalar diffusionConst, const particle_flavor flavor) {
    util::validateTypeName(name);
    particle_type_type t_id = type_counter_++;
    type_mapping_.emplace(name, t_id);
    particle_info_.emplace(std::make_pair(t_id, ParticleTypeInfo{name, diffusionConst, flavor, t_id}));
    n_types_++;
}

particle_type_type ParticleTypeRegistry::idOf(const std::string &name) const {
    return _idOf(name);
}

const std::size_t &ParticleTypeRegistry::nTypes() const {
    return n_types_;
}

particle_type_type ParticleTypeRegistry::_idOf(const std::string &name) const {
    auto it = type_mapping_.find(name);
    if (it == type_mapping_.end()) {
        throw std::invalid_argument(
                fmt::format("Could not find type \"{}\", did you forget to register it before accessing it?", name)
        );
    }
    return it->second;
}

std::string ParticleTypeRegistry::describe() const {
    namespace rus = readdy::util::str;
    std::string description;
    description += fmt::format(" - particle types:{}", rus::newline);
    for (const auto &entry : particle_info_) {
        description += fmt::format("     * particle type \"{}\" with D={}, flavor={}, id={}{}", entry.second.name,
                                   entry.second.diffusionConstant,
                                   particleflavor::particle_flavor_to_str(entry.second.flavor), entry.second.typeId,
                                   rus::newline);
    }
    return description;
}

void ParticleTypeRegistry::configure() { /*no op*/ }

particle_type_type ParticleTypeRegistry::operator()(const std::string &name) const {
    return idOf(name);
}

}
}