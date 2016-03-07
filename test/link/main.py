#!/usr/bin/env python

import numpy as np
import plumed


if plumed.installed():
    print 'PLUMED API version {:d}'.format(plumed.get_API_version())
else:
    raise StandardError('PLUMED not installed.')


n_atoms = 10
positions = np.zeros([n_atoms, 3])
forces = np.zeros_like(positions)
masses = np.ones(n_atoms)
charges = np.zeros(n_atoms)
box = np.diag(np.ones(3))
virial = np.zeros((3, 3))
positions[:,0] = np.arange(n_atoms)
positions[:,1] = 1.0
positions[:,2] = 2.0


#
# use the global PLUMED instance
#

plumed.gcreate()

print 'global instance initialized:', plumed.ginitialized()

plumed.gcmd('setPlumedDat', 'plumed.dat')

plumed.gfinalize()

print 'global instance initialized:', plumed.ginitialized()


#
# use the object interface
#

plmd = plumed.Plumed('plumed.dat', 'plumed.log', 'FancyMD', 1.0, 1.0, 1.0, 0.1, n_atoms, masses, charges)

print plmd.check_action('DISTANCE')
print plmd.check_action('XYZ')

plmd.set_step(42)
plmd.set_masses(masses)
plmd.set_charges(charges)
plmd.set_box(box.reshape(-1))
plmd.set_positions(positions.reshape(-1))
plmd.set_forces(forces.reshape(-1))
plmd.set_virial(virial.reshape(-1))

plmd.calc()
