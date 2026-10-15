"""Assert that pycv and metatomic execute in the same Python interpreter."""

import builtins
import sys

import plumedCommunications as PLMD

init = {"Value": PLMD.defaults.COMPONENT_NODEV}


def function(action):
    assert builtins._metatomic_test_interpreter == id(sys.modules)
    return 0.0
