#!/usr/bin/env python
# coding=utf-8

# Copyright © 2017 Computational Molecular Biology Group,
#                  Freie Universität Berlin (GER)
#
# This file is part of ReaDDy.
#
# ReaDDy is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General
# Public License along with this program. If not, see
# <http://www.gnu.org/licenses/>.

"""
Execute performance measurement and generate some plots
"""

import numpy as np
import matplotlib.pyplot as plt

import readdy.examples.performance_scenarios as ps

__author__ = "chrisfroe"
__license__ = "LGPL"


def sample_and_plot_scenario(scenario_type=ps.PerformanceScenario, kernel="SingleCPU"):
    n_samples = 10
    n_time_steps = 20
    number_factors = np.logspace(0, 3.5, n_samples)
    times, counts, system_vars = ps.sample_n_particles_const_density(number_factors, n_time_steps, scenario_type=ps.Collisive, kernel=kernel)
    result = (times, counts, system_vars)
    print("density", system_vars["density"])
    print("reactivity", system_vars["reactivity"])
    print("displacement", system_vars["displacement"])
    ps.plot_times(times, counts, system_vars["n_particles"])
    factor = 1e-6
    plt.plot(system_vars["n_particles"], np.fromiter(map(lambda x: factor * x, system_vars["n_particles"]), dtype=np.float), "--", label="linear")
    plt.legend(loc="best")
    plt.title(scenario_type.describe() + " at constant density, " + kernel)
    plt.xlabel("number of particles")
    plt.xscale("log")
    plt.yscale("log")
    plt.savefig(scenario_type.describe() + kernel + ".pdf")
    plt.clf()
    return result, kernel + scenario_type.describe()


if __name__ == '__main__':
    scpu_coll = sample_and_plot_scenario(ps.Collisive, "SingleCPU")
    cpu_coll = sample_and_plot_scenario(ps.Collisive, "CPU")
    scpu_react = sample_and_plot_scenario(ps.Reactive, "SingleCPU")
    cpu_react = sample_and_plot_scenario(ps.Reactive, "CPU")

    for result, label in [scpu_coll, cpu_coll, scpu_react, cpu_react]:
        ps.plot_times(result[0], result[1], result[2]["n_particles"], total_label=label)

    plt.title("SCPU vs CPU")
    plt.xlabel("number of particles")
    plt.xscale("log")
    plt.yscale("log")
    plt.legend(loc="best")
    plt.savefig("SCPU-vs-CPU.pdf")
    plt.clf()
