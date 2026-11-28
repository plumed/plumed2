import os
import subprocess
import sys
import textwrap

# PLUMED runs within Python: actions that cannot be created, Plumed objects used
# by several threads at the same time and Plumed objects destroyed in another
# thread. In a subprocess, so that a crash or a deadlock is a failure of the test.

USER_MODULE = """
import numpy as np
import plumedCommunications

plumedInit = {"Value": plumedCommunications.defaults.COMPONENT_NODEV}


def plumedCalculate(action):
    at = action.getPositions()
    return np.linalg.norm(at[0] - at[1])


def square(action):
    return action.argument(0) ** 2
"""

HELPERS = """
import numpy as np
from plumed import Plumed


def make(log):
    plmd = Plumed()
    plmd.cmd("setNatoms", 2)
    plmd.cmd("setLogFile", log)
    plmd.cmd("init")
    plmd.cmd("readInputLine", "LOAD MODULE=pycv")
    return plmd


def value(plmd, line):
    name = line.split(":")[0]
    plmd.cmd("readInputLine", line)
    shape = np.zeros(1, dtype=np.int_)
    plmd.cmd("getDataRank " + name, shape)
    data = np.zeros(1)
    plmd.cmd("setMemoryForData " + name, data)
    return data


def step(plmd, i):
    positions = np.array([[0.0, 0.0, 0.0], [1.0 + 0.01 * i, 0.0, 0.0]])
    plmd.cmd("setStep", i)
    plmd.cmd("setBox", np.diag(10.0 * np.ones(3)))
    plmd.cmd("setMasses", np.ones(2))
    plmd.cmd("setCharges", np.zeros(2))
    plmd.cmd("setPositions", positions)
    plmd.cmd("setForces", np.zeros((2, 3)))
    plmd.cmd("setVirial", np.zeros((3, 3)))
    plmd.cmd("calc")
    return positions[1, 0]
"""

FAILURES = """
failing = [
    "c: PYCVINTERFACE ATOMS=1,2 IMPORT=no_such_module",
    "c: PYCVINTERFACE ATOMS=1,2 IMPORT=threads_user INIT=no_such_function",
    "c: PYCVINTERFACE ATOMS=1,2 IMPORT=threads_user CALCULATE=no_such_function",
    "f: PYFUNCTION ARG=x IMPORT=no_such_module",
    "f: PYFUNCTION ARG=x IMPORT=threads_user CALCULATE=no_such_function",
]
for line in failing:
    plmd = make("failure.log")
    plmd.cmd("readInputLine", "x: CONSTANT VALUE=2")
    try:
        plmd.cmd("readInputLine", line)
    except Exception:
        print("refused:", line)
    else:
        raise SystemExit("accepted: " + line)
    plmd.finalize()

plmd = make("working.log")
c = value(plmd, "c: PYCVINTERFACE ATOMS=1,2 IMPORT=threads_user")
f = value(plmd, "f: PYFUNCTION ARG=c IMPORT=threads_user CALCULATE=square")
d = step(plmd, 3)
assert abs(c[0] - d) < 1e-12 and abs(f[0] - d * d) < 1e-12, (c, f, d)
plmd.finalize()
print("done")
"""

THREADS = """
import threading

errors = []


def work(name, objects):
    try:
        plmd = make(name + ".log")
        c = value(plmd, "c: PYCVINTERFACE ATOMS=1,2 IMPORT=threads_user")
        f = value(plmd, "f: PYFUNCTION ARG=c IMPORT=threads_user CALCULATE=square")
        for i in range(200):
            d = step(plmd, i)
            if abs(c[0] - d) > 1e-12 or abs(f[0] - d * d) > 1e-12:
                errors.append(f"{name}: step {i} {c[0]} {f[0]} {d}")
        # destroyed by the main thread
        objects.append(plmd)
    except Exception as e:
        errors.append(f"{name}: {e}")


objects = []
threads = [threading.Thread(target=work, args=(f"t{i}", objects)) for i in range(4)]
for t in threads:
    t.start()
for t in threads:
    t.join()
assert not errors, errors
assert len(objects) == 4
for plmd in objects:
    plmd.finalize()

# a Plumed object created in the main thread and destroyed in another one
plmd = make("main.log")
c = value(plmd, "c: PYCVINTERFACE ATOMS=1,2 IMPORT=threads_user")
step(plmd, 0)
t = threading.Thread(target=plmd.finalize)
t.start()
t.join()
print("done")
"""


def run(tmp_path, script):
    (tmp_path / "threads_user.py").write_text(textwrap.dedent(USER_MODULE))
    (tmp_path / "script.py").write_text(
        textwrap.dedent(HELPERS) + textwrap.dedent(script)
    )
    result = subprocess.run(
        [sys.executable, "script.py"],
        cwd=tmp_path,
        env=dict(os.environ, PLUMED_MAXBACKUP="0"),
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "done" in result.stdout, result.stdout + result.stderr
    return result


def test_failing_constructors(tmp_path):
    result = run(tmp_path, FAILURES)
    assert result.stdout.count("refused:") == 5, result.stdout


def test_threads(tmp_path):
    run(tmp_path, THREADS)
