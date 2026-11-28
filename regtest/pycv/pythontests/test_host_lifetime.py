import os
import subprocess
import sys
import textwrap

import pytest

# PLUMED runs within Python, which keeps plumedCommunications after the Plumed
# object has been finalized (and the pycv library released by PLUMED): using it,
# creating a new Plumed object and exiting must work. In a subprocess, so that a
# crash at exit is a failure of this test.

USER_MODULE = """
import plumedCommunications

plumedInit = {"Value": plumedCommunications.defaults.COMPONENT_NODEV}


def plumedCalculate(action):
    return 2.0
"""

SCRIPT = """
from plumed import Plumed


def run():
    plmd = Plumed()
    plmd.cmd("setNatoms", 1)
    plmd.cmd("setLogFile", "test.log")
    plmd.cmd("init")
    plmd.cmd("readInputLine", "LOAD MODULE=pycv")
    plmd.cmd("readInputLine", "x: CONSTANT VALUE=2")
    plmd.cmd("readInputLine", "f: PYFUNCTION ARG=x IMPORT=lifetime_user")
    plmd.finalize()


run()
import plumedCommunications

kept = plumedCommunications
run()
print("kept:", kept.defaults.COMPONENT_NODEV["derivative"])
"""


# with autoload, and with the explicit LOAD MODULE only
@pytest.mark.parametrize("autoload", ["1", "0"])
def test_module_kept_after_finalize(tmp_path, autoload):
    (tmp_path / "lifetime_user.py").write_text(textwrap.dedent(USER_MODULE))
    (tmp_path / "lifetime_script.py").write_text(textwrap.dedent(SCRIPT))
    result = subprocess.run(
        [sys.executable, "lifetime_script.py"],
        cwd=tmp_path,
        env=dict(os.environ, PLUMED_AUTOLOAD=autoload),
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "kept: False" in result.stdout
