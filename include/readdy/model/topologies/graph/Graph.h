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
 * @file Graph.h
 * @brief << brief description >>
 * @author clonker
 * @date 16.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <functional>
#include <stdexcept>
#include <list>
#include <unordered_map>
#include <readdy/common/macros.h>
#include "Vertex.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(top)
NAMESPACE_BEGIN(graph)

class Graph {
public:

    using vertex_list = std::list<Vertex>;
    using vertex_ref = vertex_list::iterator;
    using vertex_cref = vertex_list::const_iterator;

    using edge = std::tuple<vertex_ref, vertex_ref>;
    using cedge = std::tuple<vertex_cref, vertex_cref>;

    using path_len_2 = std::tuple<vertex_ref, vertex_ref, vertex_ref>;
    using cpath_len_2 = std::tuple<vertex_cref, vertex_cref, vertex_cref>;

    using path_len_3 = std::tuple<vertex_ref, vertex_ref, vertex_ref, vertex_ref>;
    using cpath_len_3 = std::tuple<vertex_cref, vertex_cref, vertex_cref, vertex_cref>;

    using edge_callback = std::function<void(const edge &)>;
    using path_len_2_callback = std::function<void(const path_len_2 &)>;
    using path_len_3_callback = std::function<void(const path_len_3 &)>;

    Graph() = default;

    explicit Graph(vertex_list vertexList) : _vertices(std::move(vertexList)) {}

    virtual ~Graph() = default;

    Graph(const Graph &) = delete;

    Graph &operator=(const Graph &) = delete;

    Graph(Graph &&) = default;

    Graph &operator=(Graph &&) = default;

    const vertex_list &vertices() const {
        return _vertices;
    }

    vertex_list &vertices() {
        return _vertices;
    }

    vertex_ref firstVertex() {
        return vertices().begin();
    }

    vertex_ref lastVertex() {
        return --vertices().end();
    }

    bool containsEdge(const cedge& edge) const {
        const auto& v1 = std::get<0>(edge);
        const auto& v2 = std::get<1>(edge);
        const auto& v1Neighbors = v1->neighbors();
        const auto& v2Neighbors = v2->neighbors();
        return std::find(v1Neighbors.begin(), v1Neighbors.end(), v2) != v1Neighbors.end()
               && std::find(v2Neighbors.begin(), v2Neighbors.end(), v1) != v2Neighbors.end();
    }

    bool containsEdge(vertex_cref v1, vertex_cref v2) const {
        return containsEdge(std::tie(v1, v2));
    }

    const Vertex &vertexForParticleIndex(std::size_t particleIndex) const {
        auto it = std::find_if(_vertices.begin(), _vertices.end(), [particleIndex](const Vertex &vertex) {
            return vertex.particleIndex == particleIndex;
        });
        if (it != _vertices.end()) {
            return *it;
        }
        throw std::invalid_argument("graph did not contain the particle index " + std::to_string(particleIndex));
    }

    void addVertex(std::size_t particleIndex, particle_type_type particleType) {
        _vertices.emplace_back(particleIndex, particleType);
    }

    void addEdge(vertex_ref v1, vertex_ref v2) {
        v1->addNeighbor(v2);
        v2->addNeighbor(v1);
    }

    void addEdge(const edge& edge) {
        addEdge(std::get<0>(edge), std::get<1>(edge));
    }

    void addEdgeBetweenParticles(std::size_t particleIndex1, std::size_t particleIndex2) {
        auto it1 = vertexItForParticleIndex(particleIndex1);
        auto it2 = vertexItForParticleIndex(particleIndex2);
        if (it1 != _vertices.end() && it2 != _vertices.end()) {
            it1->addNeighbor(it2);
            it2->addNeighbor(it1);
        } else {
            throw std::invalid_argument("the particles indices did not exist...");
        }
    }

    void removeEdge(vertex_ref v1, vertex_ref v2) {
        assert(v1 != v2);
        v1->removeNeighbor(v2);
        v2->removeNeighbor(v1);
    }

    void removeEdge(const edge& edge) {
        removeEdge(std::get<0>(edge), std::get<1>(edge));
    }

    void removeVertex(vertex_ref vertex) {
        removeNeighborsEdges(vertex);
        _vertices.erase(vertex);
    }

    void removeParticle(std::size_t particleIndex) {
        auto v = vertexItForParticleIndex(particleIndex);
        if (v != _vertices.end()) {
            removeNeighborsEdges(v);
            _vertices.erase(v);
        } else {
            throw std::invalid_argument(
                    "the vertex corresponding to the particle with topology index " + std::to_string(particleIndex) +
                    " did not exist in the graph");
        }
    }

    bool isConnected();
    
    std::vector<std::tuple<vertex_ref, vertex_ref>> edges() {
        std::vector<std::tuple<Graph::vertex_ref, Graph::vertex_ref>> result;
        findEdges([&result](const edge& tup) {
            result.push_back(tup);
        });
        return result;
    };

    void findEdges(const edge_callback &edgeCallback);
    
    void findNTuples(const edge_callback &tuple_callback,
                     const path_len_2_callback &triple_callback,
                     const path_len_3_callback &quadruple_callback);

    std::tuple<std::vector<edge>, std::vector<path_len_2>, std::vector<path_len_3>>
    findNTuples() {
        auto tuple = std::make_tuple(std::vector<edge>(), std::vector<path_len_2>(), std::vector<path_len_3>());
        findNTuples([&](const edge& tup) {
            std::get<0>(tuple).push_back(tup);
        }, [&](const path_len_2& path2) {
            std::get<1>(tuple).push_back(path2);
        }, [&](const path_len_3& path3) {
            std::get<2>(tuple).push_back(path3);
        });
        return tuple;
    };

    /**
     * Returns the connected components, invalidates this graph
     * @return connected components
     */
    std::vector<Graph> connectedComponentsDestructive();

private:
    vertex_list _vertices {};

    void removeNeighborsEdges(vertex_ref vertex) {
        std::for_each(std::begin(vertex->neighbors()), std::end(vertex->neighbors()), [vertex](const auto neighbor) {
            neighbor->removeNeighbor(vertex);
        });
    }

    auto vertexItForParticleIndex(std::size_t particleIndex) -> decltype(_vertices.begin()) {
        return std::find_if(_vertices.begin(), _vertices.end(), [particleIndex](const Vertex &vertex) {
            return vertex.particleIndex == particleIndex;
        });
    }
};

NAMESPACE_END(graph)
NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
