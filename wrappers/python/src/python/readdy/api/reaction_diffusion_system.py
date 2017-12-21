# coding=utf-8

# Copyright © 2016 Computational Molecular Biology Group,
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
Created on 08.09.17

@author: clonker
@author: chrisfroe
"""

import numpy as _np

from readdy._internal.readdybinding.api import Context as _Context
from readdy._internal.readdybinding.api import ParticleTypeFlavor as _ParticleTypeFlavor

from readdy.api.utils import vec3_of as _v3_of
from readdy.api.conf.UnitConfiguration import UnitConfiguration as _UnitConfiguration
from readdy.api.conf.UnitConfiguration import NoUnitConfiguration as _NoUnitConfiguration

from readdy.api.registry.compartments import CompartmentRegistry as _CompartmentRegistry
from readdy.api.registry.topologies import TopologyRegistry as _TopologyRegistry
from readdy.api.registry.potentials import PotentialRegistry as _PotentialRegistry
from readdy.api.registry.reactions import ReactionRegistry as _ReactionRegistry
from readdy.api.simulation import Simulation as _Simulation

__all__ = ['ReactionDiffusionSystem']


class ReactionDiffusionSystem(object):
    def __init__(self, box_size, temperature=None, periodic_boundary_conditions=None, unit_system='default'):
        """
        Constructs a new reaction diffusion system, starting point for all reaction diffusion simulations.

        If the unit_system is set to 'default', the units
          - nanometer for lengths,
          - nanosecond for time,
          - kilojoule/mole for energy,
          - kelvin for temperature are used.
        If the unit_system is set to `None`, all values will be interpreted in units of 1. In this case
        setting or reading the temperature is not supported, as this would require setting a specific value
        for the Boltzmann constant. Use kbt directly in this case. By default, if unit_system is `None`, then kbt=1.

        :param box_size: size of the simulation box [length]
        :param temperature: temperature, per default room temperature in the selected units
        :param periodic_boundary_conditions: periodic boundary conditions, per default all periodic
        :param unit_system: can be set to 'default', None, or dictionary
        """
        self._context = _Context()
        if unit_system is None or not unit_system:
            self._unit_conf = _NoUnitConfiguration()
        else:
            args = {}
            if isinstance(unit_system, dict):
                args = unit_system
            self._unit_conf = _UnitConfiguration(**args)
        self._compartment_registry = _CompartmentRegistry(self._context.compartments, self._unit_conf)
        self._topology_registry = _TopologyRegistry(self._context.topologies, self._unit_conf)
        self._potential_registry = _PotentialRegistry(self._context.potentials, self._unit_conf)
        self._reaction_registry = _ReactionRegistry(self._context.reactions, self._unit_conf)

        if (temperature is not None) and (unit_system is None):
            raise ValueError(
                "Setting the temperature without a unit system is not supported. "
                "If working without units, set kbt instead."
            )
        if (temperature is None) and (unit_system is None):
            self.kbt = 1.
        if (temperature is None) and (unit_system is not None):
            temperature = 293. * self.units.kelvin
            self.temperature = temperature
        if (temperature is not None) and (unit_system is not None):
            self.temperature = temperature

        if periodic_boundary_conditions is None:
            periodic_boundary_conditions = [True, True, True]

        self.box_size = box_size
        self.periodic_boundary_conditions = periodic_boundary_conditions

    @property
    def units(self):
        """
        Returns access to the unit registry, so that one can produce quantities like `5. * system.units.meter`.

        For more detail, see the pint documentation: pint.readthedocs.io
        :return: the unit registry
        """
        return self._unit_conf.reg

    @property
    def length_unit(self):
        """
        Returns the system's length unit. This is a property that can only be set in the constructor.
        :return: the length unit
        """
        return self._unit_conf.length_unit

    @property
    def energy_unit(self):
        """
        Returns the system's energy unit. This is a property that can only be set in the constructor.
        :return: the energy unit
        """
        return self._unit_conf.energy_unit

    @property
    def time_unit(self):
        """
        Returns the system's time unit. This is a property that can only be set in the constructor.
        :return: the time unit
        """
        return self._unit_conf.time_unit

    @property
    def temperature_unit(self):
        """
        Returns the system's temperature unit. This is a property that can only be set in the constructor.
        :return: the temperature unit
        """
        return self._unit_conf.temperature_unit

    @property
    def kbt(self):
        """
        Returns the thermal energy of the system.
        :return: the thermal energy [energy]
        """
        return self._context.kbt * self.energy_unit

    @kbt.setter
    def kbt(self, value):
        """
        Sets the thermal energy of the system.
        :param value: the thermal energy [energy]
        """
        if not isinstance(self._unit_conf, _NoUnitConfiguration):
            raise ValueError("Setting kbt is only supported in a unitless system. Set the temperature instead.")
        value = self._unit_conf.convert(value, self.energy_unit)
        self._context.kbt = value

    @property
    def temperature(self):
        """
        Retrieves the set temperature.
        :return: the temperature [temperature]
        """
        kbt = self.kbt
        if self.temperature_unit != 1.:
            return (kbt / (self._unit_conf.boltzmann * self._unit_conf.avogadro)).to(self.temperature_unit)
        else:
            raise ValueError("No temperature unit was set. In a unitless system, refer to kbt instead.")

    @temperature.setter
    def temperature(self, value):
        """
        Sets the temperature of the system.
        :param value: the new temperature [temperature]
        """
        value = self._unit_conf.convert(value, self.temperature_unit)
        if self.temperature_unit != 1:
            kbt = self._unit_conf.convert(value * self.temperature_unit * self._unit_conf.boltzmann
                                          * self._unit_conf.avogadro, self.energy_unit)
        else:
            raise ValueError("No temperature unit was set. In a unitless system, refer to kbt instead.")
        self._context.kbt = kbt

    @property
    def box_size(self):
        """
        Retrieves the system's current box size
        :return: the box size [length]
        """
        return self.length_unit * _np.array(self._context.box_size)

    @box_size.setter
    def box_size(self, value):
        """
        Sets the system's box size
        :param value: the box size (list or tuple of length 3, numpy scalar array of squeezed shape (3,)) [length]
        """
        value = self._unit_conf.convert(value, self.length_unit)
        v3 = _v3_of(value)
        self._context.box_size = [v3[0], v3[1], v3[2]]

    @property
    def periodic_boundary_conditions(self):
        """
        The periodic boundary conditions of the system.

        If boundaries are set to non-periodic, one has to take care that the particle cannot escape the simulation
        box, e.g., with a box potential that is inside of the simulation box
        at least within the directions of non-periodicity.

        :return: a list of length 3 containing boolean values indicating whether x-, y-, or z-boundaries are periodic
                 or not
        """
        return self._context.pbc

    @periodic_boundary_conditions.setter
    def periodic_boundary_conditions(self, value):
        """
        Sets the periodic boundary conditions of the system.
        :param value: List or Tuple or NDarray of length 3 containing boolean values.
        """
        if isinstance(value, _np.ndarray):
            if value.squeeze().ndim != 1:
                raise ValueError("Invalid shape for box size!")
            value = value.astype(bool).squeeze().tolist()

        value = list(value)
        if len(value) != 3:
            raise ValueError("Invalid length for box size list! Length can only be 3 but was {}.".format(len(value)))
        self._context.pbc = value

    def calculate_max_cutoff(self):
        """
        Calculates the maximal cutoff distance w.r.t. pair potentials, binary reactions and spatially
        dependent topology reactions.

        This quantity is used for constructing neighbor lists and can have a great impact on the performance.

        :return: The maximal cutoff distance [length]
        """
        self._context.configure(False)
        return self._context.calculate_max_cutoff() * self.length_unit

    @property
    def box_volume(self):
        """
        Calculates the box volume.
        :return: the box volume [length**3]
        """
        return self._context.box_volume() * (self.length_unit ** 3)

    def add_species(self, name, diffusion_constant=1.):
        """
        Adds a species to the system.
        :param name: The species' name.
        :param diffusion_constant: The species' diffusion constant [length**2/time]
        """
        diffusion_constant = self._unit_conf.convert(diffusion_constant, self._unit_conf.diffusion_constant_unit)
        self._context.particle_types.add(name, diffusion_constant, _ParticleTypeFlavor.NORMAL)

    def add_topology_species(self, name, diffusion_constant=1.):
        """
        Adds a species to the system that can be used with topologies.
        :param name: The species' name
        :param diffusion_constant: The species' diffusion constant [length**2/time]
        """
        diffusion_constant = self._unit_conf.convert(diffusion_constant, self._unit_conf.diffusion_constant_unit)
        self._context.particle_types.add(name, diffusion_constant, _ParticleTypeFlavor.TOPOLOGY)

    def registered_species(self):
        """
        Returns a list of thus far registered species.
        """
        return [k for k in self._context.particle_types.type_mapping.keys()]

    @property
    def compartments(self):
        return self._compartment_registry

    @property
    def reactions(self):
        return self._reaction_registry

    @property
    def topologies(self):
        return self._topology_registry

    @property
    def potentials(self):
        return self._potential_registry

    def simulation(self, kernel="SingleCPU", output_file="", integrator="EulerBDIntegrator",
                   reaction_handler="Gillespie", evaluate_topology_reactions=True, evaluate_forces=True,
                   evaluate_observables=True, skin=0):
        """
        Generates a simulation object for this reaction diffusion system configuration. The configuration is copied
        into the simulation object, so subsequent changes to the reaction diffusion system will not propagate into
        the simulation.

        :param kernel: The kernel that is used to perform the simulation
        :param output_file: the output file
        :param integrator: the integrator
        :param reaction_handler: the reaction handler
        :param evaluate_topology_reactions: whether to evaluate topology reactions
        :param evaluate_forces: whether to evaluate forces
        :param evaluate_observables: whether to evaluate observables
        :param skin: the skin size [length]
        :return: the simulation object
        """
        return _Simulation(kernel, self._context, self._unit_conf, output_file=output_file, integrator=integrator,
                           reaction_handler=reaction_handler, evaluate_topology_reactions=evaluate_topology_reactions,
                           evaluate_forces=evaluate_forces, evaluate_observables=evaluate_observables, skin=skin)

    def __str__(self):
        return self._context.describe()
