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
 * @file TestObservables.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 02.05.16
 */

#include <readdy/plugin/KernelProvider.h>
#include <readdy/api/Simulation.h>
#include <readdy/testing/KernelTest.h>
#include <readdy/testing/Utils.h>

namespace m = readdy::model;

namespace {
class TestObservables : public KernelTest {

};

TEST_P(TestObservables, TestParticlePositions) {
    const unsigned int n_particles = 100;
    kernel->context().particle_types().add("type", 1.);
    const readdy::scalar  timeStep = 1.0;
    const auto particleTypeId = kernel->context().particle_types().idOf("type");
    const auto particles = std::vector<m::Particle>(n_particles, m::Particle(0, 0, 0, particleTypeId));
    kernel->stateModel().addParticles(particles);
    auto &&obs = kernel->observe().positions(3);
    auto &&connection = kernel->connectObservable(obs.get());

    auto &&integrator = kernel->actions().createIntegrator("EulerBDIntegrator", timeStep);
    using update_nl = readdy::model::actions::UpdateNeighborList;
    auto &&neighborListInit = kernel->actions().updateNeighborList(update_nl::Operation::init, 0);
    auto &&neighborList = kernel->actions().updateNeighborList(update_nl::Operation::update, -1);
    neighborListInit->perform();
    for (readdy::time_step_type t = 0; t < 100; t++) {
        integrator->perform();
        neighborList->perform();
        kernel->evaluateObservables(t);
    }

    const auto &result = obs->getResult();
    const auto &&positions = kernel->stateModel().getParticlePositions();
    auto it_pos = positions.begin();
    int j = 0;
    for (auto it = result.begin(); it != result.end(); it = std::next(it)) {
        EXPECT_EQ(*it, *it_pos);
        it_pos++;
        ++j;
    }
    EXPECT_TRUE(j == 100);
    connection.disconnect();
}

TEST_P(TestObservables, Topologies) {
    using namespace readdy;
    auto &ctx = kernel->context();
    ctx.particle_types().add("Topology A", 1.0, readdy::model::particleflavor::TOPOLOGY);
    ctx.particle_types().add("Topology B", 1.0, readdy::model::particleflavor::TOPOLOGY);
    ctx.particle_types().add("Topology Invalid Type", 1.0, readdy::model::particleflavor::TOPOLOGY);
    ctx.particle_types().add("A", 1.0, readdy::model::particleflavor::NORMAL);

    ctx.topology_registry().configureBondPotential("Topology A", "Topology A", {10, 10});
    ctx.topology_registry().configureBondPotential("Topology A", "Topology B", {10, 10});
    ctx.topology_registry().configureBondPotential("Topology B", "Topology B", {10, 10});

    ctx.boxSize() = {{10, 10, 10}};

    std::size_t n_chain_elements = 50;
    auto &toptypes = ctx.topology_registry();
    toptypes.addType("TA");

    ctx.boxSize() = {{10, 10, 10}};
    std::vector<readdy::model::TopologyParticle> topologyParticles;
    {
        topologyParticles.reserve(n_chain_elements);
        for (std::size_t i = 0; i < n_chain_elements; ++i) {
            const auto id = ctx.particle_types().idOf("Topology A");
            topologyParticles.emplace_back(-5 + i * 10. / static_cast<readdy::scalar>(n_chain_elements), 0, 0, id);
        }
    }
    auto topology = kernel->stateModel().addTopology(toptypes.idOf("TA"), topologyParticles);
    {
        auto it = topology->graph().vertices().begin();
        auto it2 = ++topology->graph().vertices().begin();
        while(it2 != topology->graph().vertices().end()) {
            topology->graph().addEdge(it, it2);
            std::advance(it, 1);
            std::advance(it2, 1);
        }
    }

    {
        // split reaction
        auto reactionFunction = [&](model::top::GraphTopology &top) {
            model::top::reactions::Recipe recipe (top);
            auto& vertices = top.graph().vertices();
            auto current_n_vertices = vertices.size();
            if(current_n_vertices > 1) {
                auto edge = readdy::model::rnd::uniform_int<>(0, static_cast<int>(current_n_vertices - 2));
                auto it1 = vertices.begin();
                auto it2 = ++vertices.begin();
                for(int i = 0; i < edge; ++i) {
                    ++it1;
                    ++it2;
                }
                recipe.removeEdge(it1, it2);
            }

            return recipe;
        };
        auto rateFunction = [](const model::top::GraphTopology &top) {
            return top.getNParticles() > 1 ? top.getNParticles()/50. : 0;
        };
        model::top::reactions::StructuralTopologyReaction reaction {reactionFunction, rateFunction};
        reaction.create_child_topologies_after_reaction();
        reaction.roll_back_if_invalid();

        toptypes.addStructuralReaction("TA", reaction);
    }
    {
        // decay reaction
        auto reactionFunction = [&](model::top::GraphTopology &top) {
            model::top::reactions::Recipe recipe (top);
            if(top.graph().vertices().size() == 1) {
                recipe.changeParticleType(top.graph().vertices().begin(),
                                          kernel->context().particle_types().idOf("A"));
            } else {
                throw std::logic_error("this reaction should only be executed when there is exactly "
                                               "one particle in the topology");
            }
            return recipe;
        };
        auto rateFunction = [](const model::top::GraphTopology &top) {
            return top.getNParticles() > 1 ? 0 : 1;
        };
        model::top::reactions::StructuralTopologyReaction reaction {reactionFunction, rateFunction};
        reaction.create_child_topologies_after_reaction();
        reaction.roll_back_if_invalid();
        toptypes.addStructuralReaction("TA", reaction);
    }

    {
        auto integrator = kernel->actions().createIntegrator("EulerBDIntegrator", 1.0);
        auto forces = kernel->actions().calculateForces();
        auto topReactions = kernel->actions().evaluateTopologyReactions(1.0);

        std::size_t time = 0;
        std::size_t n_time_steps = 500;

        auto obs = kernel->observe().topologies(1);
        obs->setCallback([&](const model::observables::Topologies::result_type &value) {
            auto tops = kernel->stateModel().getTopologies();
            ASSERT_EQ(value.size(), tops.size());
            for(auto its = std::make_pair(tops.begin(), value.begin()); its.first != tops.end(); ++its.first, ++its.second) {
                auto topPtr = *its.first;
                const auto& record = *its.second;
                const auto &topParticles = topPtr->getParticles();
                const auto &recordParticles = record.particleIndices;
                auto contains1 = std::all_of(recordParticles.begin(), recordParticles.end(), [&](auto idx) {
                    return std::find(topParticles.begin(), topParticles.end(), idx) != topParticles.end();
                });
                auto contains2 = std::all_of(topParticles.begin(), topParticles.end(), [&](auto idx) {
                    return std::find(recordParticles.begin(), recordParticles.end(), idx) != recordParticles.end();
                });
                ASSERT_TRUE(contains1) << "record.particleIndices was not contained in topology particles";
                ASSERT_TRUE(contains2) << "topology particles were not contained in record.particleIndices";

                auto topEdges = topPtr->graph().edges();
                ASSERT_EQ(topEdges.size(), record.edges.size());

                contains1 = std::all_of(record.edges.begin(), record.edges.end(), [&](const auto &edge) {
                    std::size_t ix1 = std::get<0>(edge);
                    std::size_t ix2 = std::get<1>(edge);
                    for(const auto &topEdge : topEdges) {
                        const auto &v1 = std::get<0>(topEdge);
                        const auto &v2 = std::get<1>(topEdge);
                        if(v1->particleIndex == ix1 && v2->particleIndex == ix2) {
                            return true;
                        }
                        if(v1->particleIndex == ix2 && v2->particleIndex == ix1) {
                            return true;
                        }
                    }
                    return false;
                });

                contains2 = std::all_of(topEdges.begin(), topEdges.end(), [&](const auto &e) {
                    auto vtup1 = std::make_tuple(std::get<0>(e)->particleIndex, std::get<1>(e)->particleIndex);
                    auto vtup2 = std::make_tuple(std::get<1>(e)->particleIndex, std::get<0>(e)->particleIndex);

                    auto find1 = std::find(record.edges.begin(), record.edges.end(), vtup1);
                    auto find2 = std::find(record.edges.begin(), record.edges.end(), vtup2);

                    return find1 != record.edges.end() || find2 != record.edges.end();
                });

                ASSERT_TRUE(contains1) << "there were non existent edges in the record";
                ASSERT_TRUE(contains2) << "some edges were not in the record";
            }
        });
        auto connection = kernel->connectObservable(obs.get());

        kernel->initialize();

        forces->perform();
        kernel->evaluateObservables(time);
        for(time = 1; time < n_time_steps; ++time) {
            integrator->perform();
            topReactions->perform();
            forces->perform();
            kernel->evaluateObservables(time);

        }
        kernel->finalize();
    }
}

TEST_P(TestObservables, TestForcesObservable) {
    // Setup particles
    kernel->context().particle_types().add("A", 42.);
    kernel->context().particle_types().add("B", 1337.);
    const auto typeIdA = kernel->context().particle_types().idOf("A");
    const auto typeIdB = kernel->context().particle_types().idOf("B");
    const unsigned int n_particles = 2; // There will be 55 Bs
    const auto particlesA = std::vector<m::Particle>(n_particles, m::Particle(0, 0, 0, typeIdA));
    const auto particlesB = std::vector<m::Particle>(n_particles + 5, m::Particle(0, 0, 0, typeIdB));
    kernel->stateModel().addParticles(particlesA);
    kernel->stateModel().addParticles(particlesB);
    {
        // Check if result has correct size
        // Check that empty particleType argument gives correct object, namely all forces
        auto &&obsA = kernel->observe().forces(1, std::vector<std::string>{"A"});
        auto &&obsB = kernel->observe().forces(1, std::vector<std::string>{"B"});
        auto &&obsBoth = kernel->observe().forces(1);
        auto &&connectionA = kernel->connectObservable(obsA.get());
        auto &&connectionB = kernel->connectObservable(obsB.get());
        auto &&connectionBoth = kernel->connectObservable(obsBoth.get());
        // Evaluate twice to ensure that results do not accumulate
        kernel->evaluateObservables(0);
        kernel->evaluateObservables(1);
        const auto &resA = obsA->getResult();
        const auto &resB = obsB->getResult();
        const auto &resBoth = obsBoth->getResult();
        EXPECT_EQ(resA.size(), n_particles);
        EXPECT_EQ(resB.size(), n_particles+5);
        EXPECT_EQ(resBoth.size(), n_particles + n_particles + 5);
        readdy::Vec3 zero = readdy::Vec3(0, 0, 0);
        for (auto force : resBoth) {
            EXPECT_TRUE(force == zero);
        }
    }
    // Two particles C and C with radius 1 and harmonic repulsion at distance 1.5 -> force = kappa * (radiiSum - 1.5)
    kernel->context().periodicBoundaryConditions() = {{false, false, false}};
    kernel->context().boxSize() = {{5, 5, 5}};
    kernel->context().particle_types().add("C", 1.);
    kernel->context().potentials().addBox("A", .001, {-2.4, -2.4, -2.4}, {4.8, 4.8, 4.8});
    kernel->context().potentials().addBox("B", .001, {-2.4, -2.4, -2.4}, {4.8, 4.8, 4.8});
    kernel->context().potentials().addBox("C", .001, {-2.4, -2.4, -2.4}, {4.8, 4.8, 4.8});
    const auto typeIdC = kernel->context().particle_types().idOf("C");
    const auto particlesC = std::vector<m::Particle>{m::Particle(0, 0, 0, typeIdC), m::Particle(0, -1.5, 0, typeIdC)};
    kernel->stateModel().addParticles(particlesC);

    kernel->context().potentials().addHarmonicRepulsion("C", "C", 2.0, 2.0);

    using update_nl = readdy::model::actions::UpdateNeighborList;
    auto &&nl = kernel->actions().updateNeighborList();
    auto &&forces = kernel->actions().calculateForces();
    kernel->context().configure();
    kernel->initialize();
    {
        auto obsC = kernel->observe().forces(1, std::vector<std::string>{"C"});
        auto connectionC = kernel->connectObservable(obsC.get());
        nl->perform();
        forces->perform();
        kernel->evaluateObservables(2);
        const auto &resC = obsC->getResult();
        readdy::Vec3 force0 = readdy::Vec3(0., 1., 0.);
        readdy::Vec3 force1 = readdy::Vec3(0., -1., 0.);
        EXPECT_EQ(resC.size(), 2);
        EXPECT_TRUE(resC[0] == force0 || resC[1] == force0);
        EXPECT_TRUE(resC[1] == force1 || resC[0] == force1);
    }
}

INSTANTIATE_TEST_CASE_P(TestObservablesKernel, TestObservables,
                        ::testing::ValuesIn(readdy::testing::getKernelsToTest()));
}

