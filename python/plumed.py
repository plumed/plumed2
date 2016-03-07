import ctypes as ct

"""
Python interface to PLUMED 2

author: Ondrej Marsalek
        ondrej.marsalek@gmail.com

Features that are currently not supported:
- MPI, should not be a problem if needed
- callback functions, should be possible:
  https://docs.python.org/2/library/ctypes.html#callback-functions

Links to useful documentation:

https://docs.python.org/2/library/ctypes.html
http://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.ctypes.html
http://plumed.github.io/doc-v2.1/developer-doc/html/_how_to_plumed_your_m_d.html

"""

#
# load the PLUMED shared library
#

# TODO: might need `ctypes.util.find_library` for more robustness?

_libplumed = ct.CDLL('libplumed.so')


#
# thinly wrap global instance methods
#

def gcreate():
    _libplumed.plumed_gcreate()


def ginitialized():
    return bool(_libplumed.plumed_ginitialized())


def gfinalize():
    _libplumed.plumed_gfinalize()


def gcmd(cmd, ptr=None):
    _libplumed.plumed_gcmd(cmd, ptr)


def installed():
    return bool(_libplumed.plumed_installed())

def get_API_version():
    version = ct.c_int()
    do_finalize = False
    if not ginitialized():
        do_finalize = True
        gcreate()
    _libplumed.plumed_gcmd('getApiVersion', ct.byref(version))
    if do_finalize:
        gfinalize()
    return version.value


#
# create Pythonic interface
#

class Plumed(object):

    # TODO: allow construction with access to the global instance?

    def __init__(self, fn_in, fn_out, engine,
                 energy_units, length_units, time_units, dt,
                 n_atoms, masses, charges, no_virial=False):

        self.n_atoms = n_atoms

        self.main = _libplumed.plumed_create()

        _libplumed.plumed_cmd(self.main, 'setPlumedDat', fn_in)
        _libplumed.plumed_cmd(self.main, 'setLogFile', fn_out)
        _libplumed.plumed_cmd(self.main, 'setMDEngine', engine)
        _libplumed.plumed_cmd(self.main, 'setRealPrecision', ct.byref(ct.c_int(8)))
        _libplumed.plumed_cmd(self.main, 'setMDEnergyUnits', ct.byref(ct.c_double(energy_units)))
        _libplumed.plumed_cmd(self.main, 'setMDLengthUnits', ct.byref(ct.c_double(length_units)))
        _libplumed.plumed_cmd(self.main, 'setMDTimeUnits', ct.byref(ct.c_double(time_units)))
        _libplumed.plumed_cmd(self.main, 'setTimestep', ct.byref(ct.c_double(dt)))
        _libplumed.plumed_cmd(self.main, 'setNatoms', ct.byref(ct.c_int(n_atoms)))

        if no_virial:
            _libplumed.plumed_cmd(self.main, 'setNoVirial', None)

        _libplumed.plumed_cmd(self.main, 'init', None)

    def __del__(self):
        _libplumed.plumed_finalize(self.main)

    def cmd(self, cmd, ptr=None):
        _libplumed.plumed_cmd(self.main, cmd, ptr)

    def check_action(self, action):
        check = ct.c_int()
        _libplumed.plumed_cmd(self.main, 'checkAction %s'.format(action), ct.byref(check))
        return check.value

    def set_masses(self, masses):
        if masses.flags.c_contiguous and (masses.shape == (self.n_atoms,)):
            _libplumed.plumed_cmd(self.main, 'setMasses', masses.ctypes.data)
        else:
            raise ValueError('Incompatible array passed to `set_masses`.')

    def set_charges(self, charges):
        if charges.flags.c_contiguous and (charges.shape == (self.n_atoms,)):
            _libplumed.plumed_cmd(self.main, 'setCharges', charges.ctypes.data)
        else:
            print charges.flags
            print charges.shape
            raise ValueError('Incompatible array passed to `set_charges`.')

    def set_step(self, step):
        # TODO
        # Does this memory location need to be persistent until calc is called?
        # If so, store the int in self.
        _libplumed.plumed_cmd(self.main, 'setStep', ct.byref(ct.c_int(step)))

    def set_box(self, box):
        if box.flags.c_contiguous and (box.shape == (9,)):
            _libplumed.plumed_cmd(self.main, 'setBox', box.ctypes.data)
        else:
            print box.flags
            print box.shape
            raise ValueError('Incompatible array passed to `set_box`.')

    def set_positions(self, positions):
        if positions.flags.c_contiguous and (positions.shape == (3 * self.n_atoms,)):
            _libplumed.plumed_cmd(self.main, 'setPositions', positions.ctypes.data)
        else:
            print positions.flags
            print positions.shape
            raise ValueError('Incompatible array passed to `set_positions`.')

    def set_position_xyz(self, x, y, z):
        raise NotImplementedError

    def set_forces(self, forces):
        if forces.flags.c_contiguous and forces.flags.writeable and (forces.shape == (3 * self.n_atoms,)):
            _libplumed.plumed_cmd(self.main, 'setForces', forces.ctypes.data)
        else:
            print forces.flags
            print forces.shape
            raise ValueError('Incompatible array passed to `set_forces`.')

    def set_forces_xyz(self, f_x, f_y, f_z):
        raise NotImplementedError

    def set_virial(self, virial):
        if virial.flags.c_contiguous and virial.flags.writeable and (virial.shape == (9,)):
            _libplumed.plumed_cmd(self.main, 'setVirial', virial.ctypes.data)
        else:
            print virial.flags
            print virial.shape
            raise ValueError('Incompatible array passed to `set_virial`.')

    def calc(self):
        _libplumed.plumed_cmd(self.main, 'calc', None)
