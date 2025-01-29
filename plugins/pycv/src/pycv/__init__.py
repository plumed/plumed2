# checked and formatted with ruff 0.6.9


def getPythonCVInterface():
    """returns the location of the pycv shared object"""
    # older version
    # import inspect
    # import pycv

    # path_of_this = inspect.getfile(pycv)
    import importlib.util

    path_of_this = importlib.util.find_spec("pycv").origin
    path_of_lib = path_of_this[: path_of_this.rfind("/")]
    return path_of_lib + "/PythonCVInterface.so"


def main():
    from argparse import ArgumentParser

    parser = ArgumentParser(
        prog="pyCV",
        description="""shows the path for the pycv shared object.

        Just run this with no arguments to see the path.
        """,
    )
    parser.add_argument(
        "-p",
        "--plumed",
        help="print on screen an example plumed.dat",
        action="store_true",
    )
    args = parser.parse_args()

    if args.plumed:
        print(
            f"""LOAD FILE={getPythonCVInterface()}
cvPy:PYCVINTERFACE IMPORT=mypycv
fPy: PYFUNCTION IMPORT=mypycvfunc ARG=cvPy
PRINT FILE=colvar.out ARG=*"""
        )
        return 0

    print(getPythonCVInterface())
    return 0
