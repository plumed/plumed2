#!/usr/bin/env python3
"""Compare restart output while allowing OpenMP reduction roundoff."""

from itertools import zip_longest
import math
from pathlib import Path
import sys


def compare(reference, restarted):
    data_rows = 0
    with Path(reference).open() as left, Path(restarted).open() as right:
        for lineno, (a, b) in enumerate(zip_longest(left, right), 1):
            if a is None or b is None:
                raise ValueError(f"line {lineno}: different number of lines")
            if a.startswith("#") or b.startswith("#"):
                if a != b:
                    raise ValueError(f"line {lineno}: different COLVAR headers")
                continue
            fields_a, fields_b = a.split(), b.split()
            if len(fields_a) != len(fields_b):
                raise ValueError(f"line {lineno}: different number of columns")
            for column, (x, y) in enumerate(zip(fields_a, fields_b), 1):
                x, y = float(x), float(y)
                if not (math.isfinite(x) and math.isfinite(y)):
                    raise ValueError(f"line {lineno}, column {column}: nonfinite value")
                if not math.isclose(x, y, rel_tol=1e-12, abs_tol=1e-12):
                    raise ValueError(f"line {lineno}, column {column}: {x} != {y}")
            data_rows += bool(fields_a)
    if not data_rows:
        raise ValueError("no COLVAR data rows")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("usage: compare-colvar.py REFERENCE RESTARTED")
    try:
        compare(*sys.argv[1:])
    except (OSError, ValueError) as error:
        sys.exit(str(error))
