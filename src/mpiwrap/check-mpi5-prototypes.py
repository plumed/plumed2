#!/usr/bin/env python3
"""Compile exact manifest prototypes, retained constants and types against selected local MPI headers."""
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / "support-matrix.tsv"

def split_arguments(text):
    result, start, depth = [], 0, 0
    for index, char in enumerate(text):
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
        elif char == "," and depth == 0:
            result.append(text[start:index].strip())
            start = index + 1
    tail = text[start:].strip()
    if tail and tail != "void":
        result.append(tail)
    return result

def function_pointer_type(signature):
    match = re.match(r"^(.+?)\s+(MPI_[A-Za-z0-9_]+)\s*\((.*)\)$", signature.strip())
    if not match:
        raise ValueError("cannot parse function declaration: " + signature)
    return_type, name, raw_args = match.groups()
    params = []
    for argument in split_arguments(raw_args):
        argument = argument.strip()
        if argument == "...":
            params.append("...")
            continue
        array = re.search(r"([A-Za-z_][A-Za-z0-9_]*)\s*((?:\[[^]]*\])+)$", argument)
        if array:
            declarator = array.group(2)
            base = argument[:array.start()].strip()
            dimensions = re.findall(r"\[([^]]*)\]", declarator)
            # Function-parameter arrays decay by one dimension.
            if len(dimensions) == 1:
                params.append(base + " *")
            else:
                remaining = "".join("[" + d + "]" for d in dimensions[1:])
                params.append(base + " (*)" + remaining)
            continue
        # Remove the final parameter identifier, keeping pointer stars and qualifiers.
        ident = re.search(r"([A-Za-z_][A-Za-z0-9_]*)\s*$", argument)
        if not ident:
            raise ValueError("cannot parse parameter: " + argument)
        base = argument[:ident.start()].rstrip()
        # A lone typedef-like parameter is only valid for abstract declarations; all
        # manifest prototypes carry names, so treating that as a failure is safer.
        if not base:
            raise ValueError("parameter lacks a type: " + argument)
        params.append(base)
    return name, f"{return_type} (*)({', '.join(params) if params else 'void'})"

def manifest_functions():
    lines = MANIFEST.read_text().splitlines()
    functions = []
    for line in lines[1:]:
        row = line.split("\t")
        if row[1] == "function" and row[0].startswith("MPI_") and row[4] not in {
            "excluded-mpi5-removed", "excluded-deprecated"
        } and row[0] != "MPI_Comm_from_native":
            functions.append(row)
    return functions

def source_text(rows):
    body = ["#include <mpi.h>", "#include <type_traits>"]
    count = 0
    for row in rows:
        try:
            name, expected = function_pointer_type(row[2])
        except ValueError as error:
            raise SystemExit(str(error))
        if name == "MPI_Comm_f2c":
            # MPI permits this conversion routine to be a macro; exercise the call
            # expression instead of taking the macro's address.
            body.append("static_assert(std::is_same<decltype(MPI_Comm_f2c(MPI_Fint{})), MPI_Comm>::value, \"MPI_Comm_f2c result type\");")
            count += 1
            continue
        body.append(f"using expected_{name} = {expected};")
        body.append(f"static_assert(std::is_same<decltype(&{name}), expected_{name}>::value, \"MPI-5 signature: {name}\");")
        count += 1
    body.extend(["namespace PLMD { int run_probe() { return 0; } }", "int main() { return PLMD::run_probe(); }"])
    return "\n".join(body) + "\n", count

def source_types_and_constants():
    body = ["#include <mpi.h>", "#include <type_traits>"]
    lines = MANIFEST.read_text().splitlines()
    count = 0
    for line in lines[1:]:
        row = line.split("\t")
        if row[4] in {"excluded-mpi5-removed", "excluded-deprecated", "excluded-non-C-binding"}:
            continue
        name, kind = row[0], row[1]
        if kind == "type" and name in {
            "MPI_User_function", "MPI_Comm_errhandler_function"
        }:
            body.append(f"static_assert(std::is_function<{name}>::value, \"callback type: {name}\");")
            expected = {
                "MPI_User_function": "void(void *, void *, int *, MPI_Datatype *)",
                "MPI_Comm_errhandler_function": "void(MPI_Comm *, int *, ...)",
            }[name]
            body.append(f"static_assert(std::is_same<{name}, {expected}>::value, \"MPI-5 callback signature: {name}\");")
            count += 1
        elif kind == "type" and name.startswith("MPI_"):
            body.append(f"static_assert(sizeof({name}) > 0, \"MPI type: {name}\");")
            count += 1
        elif kind == "constant" and name.startswith("MPI_"):
            body.append(f"static_assert(sizeof({name}) > 0, \"MPI constant: {name}\");")
            count += 1
    body.extend(["namespace PLMD { int run_values_probe() { return 0; } }", "int main() { return PLMD::run_values_probe(); }"])
    return "\n".join(body) + "\n", count

def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: check-mpi5-prototypes.py MPICH_MPICXX OPENMPI_MPICXX")
    source, count = source_text(manifest_functions())
    value_source, value_count = source_types_and_constants()
    for label, wrapper, variable in (("mpich", sys.argv[1], "MPICH_CXX"),
                                     ("openmpi", sys.argv[2], "OMPI_CXX")):
        env = os.environ.copy()
        env.setdefault(variable, "/usr/bin/clang++")
        with tempfile.TemporaryDirectory(prefix="mpiw-signature-") as temp:
            cpp = Path(temp) / "signatures.cpp"
            values_cpp = Path(temp) / "values.cpp"
            exe = Path(temp) / "signatures"
            values_exe = Path(temp) / "values"
            cpp.write_text(source)
            values_cpp.write_text(value_source)
            cmd = [wrapper, "-std=c++17", "-Wall", "-Wextra", str(cpp), "-o", str(exe)]
            print(f"[{label}] {wrapper}: compiling {count} MPI declarations")
            subprocess.run(cmd, check=True, env=env)
            subprocess.run([str(exe)], check=True)
            value_cmd = [wrapper, "-std=c++17", "-Wall", "-Wextra", str(values_cpp), "-o", str(values_exe)]
            print(f"[{label}] compiling {value_count} retained MPI types/constants")
            subprocess.run(value_cmd, check=True, env=env)
            subprocess.run([str(values_exe)], check=True)
    print(f"verified {count} active function signatures and {value_count} active type/constant names on both headers")

if __name__ == "__main__":
    main()
