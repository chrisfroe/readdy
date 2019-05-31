/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 ********************************************************************/


/**
 * << detailed description >>
 *
 * @file TestReactions.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 21.06.16
 */

#include <catch2/catch.hpp>

#include <readdy/plugin/KernelProvider.h>
#include <readdy/api/SimulationLoop.h>
#include <readdy/testing/KernelTest.h>
#include <readdy/testing/Utils.h>

using namespace readdytesting::kernel;

static constexpr const char *const REACTION_HANDLERS[] = {
        "UncontrolledApproximation", "Gillespie", "DetailedBalance"
};

// Compare two vectors via their distance to the (1,1,1) plane in 3d space. This is quite arbitrary and only used
// to construct an ordered set.
struct Vec3ProjectedLess {
    bool operator()(const readdy::Vec3& lhs, const readdy::Vec3& rhs) const {
        return (lhs.x + lhs.y + lhs.z) < (rhs.x + rhs.y + rhs.z);
    }
};

TEMPLATE_TEST_CASE("Test reaction handlers", "[reactions]", SingleCPU, CPU) {
    auto kernel = create<TestType>();
    auto &ctx = kernel->context();
    for (const auto &handler : REACTION_HANDLERS) {
        SECTION(handler) {
            if (kernel->name() == "CPU" && std::string(handler) == "DetailedBalance") {
                INFO("Skipping tests for the combination (CPU, DetailedBalance)");
                continue;
            }
            SECTION("Constant number of particles") {
                // scenario: two particle types A and B, which can form a complex AB which after a time is going
                // to dissolve back into A and B. Therefore,
                // the numbers #(A) + #(AB) and #(B) + #(AB) must remain constant.

                using n_particles_obs = readdy::model::observables::NParticles;

                std::random_device rd;
                std::uniform_real_distribution<readdy::scalar> dist(-2.5, 2.5);
                auto stdRand = [&dist, &rd]() -> readdy::scalar {
                    return dist(rd);
                };

                ctx.particleTypes().add("A", 1.0);
                ctx.particleTypes().add("B", 1.0);
                ctx.particleTypes().add("AB", 0.0);
                ctx.periodicBoundaryConditions() = {{true, true, true}};
                ctx.boxSize() = {{5, 5, 5}};
                ctx.reactions().addFusion("Form complex", "A", "B", "AB", .5, 1.0);
                ctx.reactions().addFission("Dissolve", "AB", "A", "B", .5, 1.0);

                unsigned long n_A = 50;
                unsigned long n_B = n_A;

                for (unsigned long i = 0; i < n_A; ++i) {
                    kernel->addParticle("A", {stdRand(), stdRand(), stdRand()});
                    kernel->addParticle("B", {stdRand(), stdRand(), stdRand()});
                }


                auto obs = kernel->observe().nParticles(1, std::vector<std::string>({"A", "B", "AB"}));
                auto conn = kernel->connectObservable(obs.get());
                obs->callback() = [&n_A, &n_B](const n_particles_obs::result_type &result) {
                    if (result.size() == 2) {
                        REQUIRE(n_A == result[0] + result[2]);
                        REQUIRE(n_B == result[1] + result[2]);
                    }
                };

                readdy::api::SimulationLoop loop(kernel.get(), 1);
                loop.useReactionScheduler(handler);
                loop.run(10);
            }

            SECTION("Fusion and Fission weights") {
                // A diffuses, F is fixed
                // Idea: position F particles and remember their positions (ordered set),
                // do ONE time-step and check if current positions are still the same.
                ctx.particleTypes().add("A", 0.5);
                ctx.particleTypes().add("F", 0.0);
                ctx.periodicBoundaryConditions() = {{true, true, true}};
                ctx.boxSize() = {{20, 20, 20}};

                const readdy::scalar weightF {static_cast<readdy::scalar>(0)};
                const readdy::scalar weightA  {static_cast<readdy::scalar>(1.)};
                ctx.reactions().addFusion("F+A->F", "F", "A", "F", 1.0, 2.0, weightF, weightA);
                ctx.reactions().addFission("F->F+A", "F", "F", "A", 1.0, 2.0, weightF, weightA);

                std::set<readdy::Vec3, Vec3ProjectedLess> fPositions;
                auto n3 = readdy::model::rnd::normal3<readdy::scalar>;
                for (std::size_t i = 0; i < 15; ++i) {
                    auto fPos = n3(static_cast<readdy::scalar>(0.), static_cast<readdy::scalar>(0.8));
                    kernel->addParticle("F", fPos);
                    // fPos
                    fPositions.emplace(kernel->stateModel().getParticles().back().pos());
                    kernel->addParticle("A", n3(static_cast<readdy::scalar>(0.), static_cast<readdy::scalar>(1.)));
                }

                auto obs = kernel->observe().positions(1, std::vector<std::string>({"F"}));
                obs->callback() =
                        [&](const readdy::model::observables::Positions::result_type &result) {
                            std::set<readdy::Vec3, Vec3ProjectedLess> checklist;
                            for (const auto &pos : result) {
                                checklist.emplace(pos);
                            }
                            REQUIRE(fPositions.size() == checklist.size());
                            auto itPos = fPositions.begin();
                            auto itCheck = checklist.begin();
                            while (itPos != fPositions.end()) {
                                readdy::testing::vec3eq(*itPos, *itCheck, kernel->doublePrecision() ? 1e-8 : 1e-6);
                                ++itPos;
                                ++itCheck;
                            }
                        };
                auto connection = kernel->connectObservable(obs.get());
                {
                    readdy::api::SimulationLoop loop (kernel.get(), .1);
                    loop.useReactionScheduler(handler);
                    loop.run(1);
                }
            }

            SECTION("Fusion through periodic boundary") {
                ctx.boxSize() = {{10, 10, 10}};
                ctx.periodicBoundaryConditions() = {true, true, true};
                ctx.particleTypes().add("A", 0.);
                ctx.reactions().add("fus: A +(2) A -> A", 1e16);

                auto &&reactions = kernel->actions().createReactionScheduler(handler, 1);
                auto &&initNeighborList = kernel->actions().createNeighborList(ctx.calculateMaxCutoff());
                auto &&neighborList = kernel->actions().updateNeighborList();
                // resulting particle should be at 4.9
                kernel->stateModel().addParticle({4.5, 4.5, 4.5, ctx.particleTypes().idOf("A")});
                kernel->stateModel().addParticle({-4.7, -4.7, -4.7, ctx.particleTypes().idOf("A")});

                initNeighborList->perform();
                neighborList->perform();
                reactions->perform();

                const auto particles = kernel->stateModel().getParticles();
                REQUIRE(particles.size() == 1);
                const auto &pos = particles[0].pos();
                REQUIRE(pos.x == Approx(4.9));
                REQUIRE(pos.y == Approx(4.9));
                REQUIRE(pos.z == Approx(4.9));
            }

            SECTION("Fission near boundary") {
                ctx.boxSize() = {{10, 10, 10}};
                ctx.periodicBoundaryConditions() = {true, true, true};
                ctx.particleTypes().add("A", 0.);
                ctx.reactions().add("fis: A -> A +(2) A", 1e16);

                auto &&reactions = kernel->actions().createReactionScheduler(handler, 1);
                // one product will be in negative-x halfspace, the other in positive-x halfspace
                kernel->stateModel().addParticle({-4.9999999, 0, 0, ctx.particleTypes().idOf("A")});

                reactions->perform();

                const auto particles = kernel->stateModel().getParticles();
                REQUIRE(particles.size() == 2);
                const auto &pos1 = particles[0].pos();
                const auto &pos2 = particles[1].pos();
                if (pos1.x <= 0.) {
                    REQUIRE(pos2.x > 0.);
                } else if (pos1.x >= 0.) {
                    REQUIRE(pos2.x< 0.);
                }
            }

            SECTION("Michaelis Menten") {
                //if (std::string(handler) == "UncontrolledApproximation") {
                    /**
                    * Since comparing the value of a stochastic process (number of particles over time) is not well
                    * suited for testing, here integrate the stochastic process over the simulated time,
                    * giving a quite robust observable. This can be compared to an analytic value from ODE kinetics.
                    */
                    namespace rnd = readdy::model::rnd;

                    readdy::scalar length {10.};
                    ctx.boxSize() = {{length, length, length}};
                    ctx.periodicBoundaryConditions() = {true, true, true};

                    ctx.particleTypes().add("E", 0.5);
                    ctx.particleTypes().add("S", 0.5);
                    ctx.particleTypes().add("ES", 0.5);
                    ctx.particleTypes().add("P", 0.5);

                    ctx.reactions().add("fwd: E +(1.0) S -> ES", 0.0023417691399750494);
                    ctx.reactions().add("back: ES -> E +(1.0) S", 1.);
                    ctx.reactions().add("prod: ES -> E +(1.0) P", 1.);

                    readdy::scalar dt {8e-3};
                    auto &&reactions = kernel->actions().createReactionScheduler(handler, dt);
                    auto &&bd = kernel->actions().eulerBDIntegrator(dt);
                    auto &&initNeighborList = kernel->actions().createNeighborList(1.0);
                    auto &&neighborList = kernel->actions().updateNeighborList();

                    std::size_t nParticlesE {70};
                    std::size_t nParticlesS {700};

                    for (int i = 0; i < nParticlesE; ++i) {
                        readdy::Vec3 pos{rnd::uniform_real() * length - 0.5 * length,
                                         rnd::uniform_real() * length - 0.5 * length,
                                         rnd::uniform_real() * length - 0.5 * length};
                        kernel->addParticle("E", pos);
                    }

                    for (int i = 0; i < nParticlesS; ++i) {
                        readdy::Vec3 pos{rnd::uniform_real() * length - 0.5 * length,
                                         rnd::uniform_real() * length - 0.5 * length,
                                         rnd::uniform_real() * length - 0.5 * length};
                        kernel->addParticle("S", pos);
                    }

                    std::size_t nSteps {37500};

                    std::vector<std::string> typesToCount {{"S"}};
                    auto &&obs = kernel->observe().nParticles(1, typesToCount);
                    readdy::scalar meanS {0.};

                    obs->callback() = [&](const readdy::model::observables::NParticles::result_type &result) {
                        meanS += static_cast<readdy::scalar>(result[0]);
                    };

                    auto &&connection = kernel->connectObservable(obs.get());

                    initNeighborList->perform();
                    neighborList->perform();
                    kernel->evaluateObservables(0);
                    auto t1 = std::chrono::high_resolution_clock::now();
                    for (readdy::time_step_type t = 1; t < nSteps+1; t++) {
                        if (t == (nSteps+1)/100) {
                            auto t2 = std::chrono::high_resolution_clock::now();
                            auto dur = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
                            readdy::log::warn("1% done in {} ms, ETA is {} seconds", dur, dur / 10.);
                        }
                        bd->perform();
                        neighborList->perform();
                        reactions->perform();
                        neighborList->perform();
                        kernel->evaluateObservables(t);
                    }

                    meanS /= static_cast<readdy::scalar>(nSteps+1);

                    CHECK(meanS < 665.2 + 0.03 * 665.2);
                    CHECK(meanS > 665.2 - 0.03 * 665.2);
                //}
            }
        }
    }
}
