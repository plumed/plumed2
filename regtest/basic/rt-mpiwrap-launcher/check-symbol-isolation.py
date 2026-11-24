#!/usr/bin/env python3
"""Build and run upstream host tMPI and namespaced PLUMED tMPI together."""
import argparse
import pathlib
import subprocess
import tarfile
import tempfile

PINNED = "6c69f720aa7378fddfee67f99d51aff21ccfb0a2"
SOURCE_CPP = """alltoall.cpp atomic.cpp barrier.cpp bcast.cpp collective.cpp comm.cpp
errhandler.cpp event.cpp gather.cpp group.cpp list.cpp lock.cpp numa_malloc.cpp
once.cpp p2p_protocol.cpp p2p_send_recv.cpp p2p_wait.cpp profile.cpp
pthreads.cpp reduce.cpp reduce_fast.cpp scan.cpp scatter.cpp tmpi_init.cpp
tmpi_malloc.cpp topology.cpp type.cpp""".split()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gromacs", type=pathlib.Path)
    ap.add_argument("--source", type=pathlib.Path,
                    default=pathlib.Path(__file__).resolve().parent)
    args = ap.parse_args()
    root = args.gromacs.resolve()
    rev = subprocess.check_output(["git", "-C", str(root), "rev-parse", "v2026.4^{commit}"], text=True).strip()
    if rev != PINNED:
        raise SystemExit(f"expected v2026.4 at {PINNED}, found {rev}")
    archive = subprocess.check_output(["git", "-C", str(root), "archive", PINNED,
                                       "src/external/thread_mpi"])
    source = args.source.resolve()
    with tempfile.TemporaryDirectory(prefix="plumed-tmpi-isolation-") as tmp:
        work = pathlib.Path(tmp)
        with tarfile.open(fileobj=__import__("io").BytesIO(archive)) as tf:
            tf.extractall(work)
        host = work / "src/external/thread_mpi"
        cpp = "\n".join(f"  {host / 'src' / f}" for f in SOURCE_CPP)
        cmake = f'''cmake_minimum_required(VERSION 3.16)
project(thread_mpi_symbol_isolation LANGUAGES CXX)
find_package(Threads REQUIRED)
add_subdirectory("{source.as_posix()}" imported)
add_library(host_thread_mpi STATIC
{cpp}
)
target_include_directories(host_thread_mpi PUBLIC "{host / 'include'}" PRIVATE "{host / 'src'}")
target_compile_features(host_thread_mpi PUBLIC cxx_std_17)
target_compile_definitions(host_thread_mpi PRIVATE THREAD_PTHREADS HAVE_UNISTD_H HAVE_SYS_TIME_H HAVE_PTHREAD_H)
target_link_libraries(host_thread_mpi PUBLIC Threads::Threads)
configure_file("{(source / 'host-runner.cpp.in').as_posix()}" "${{CMAKE_CURRENT_BINARY_DIR}}/host-runner.cpp" COPYONLY)
configure_file("{(source / 'internal-runner.cpp.in').as_posix()}" "${{CMAKE_CURRENT_BINARY_DIR}}/internal-runner.cpp" COPYONLY)
configure_file("{(source / 'isolation-main.cpp.in').as_posix()}" "${{CMAKE_CURRENT_BINARY_DIR}}/isolation-main.cpp" COPYONLY)
add_library(host_runner OBJECT "${{CMAKE_CURRENT_BINARY_DIR}}/host-runner.cpp")
target_include_directories(host_runner PRIVATE "{host / 'include'}")
add_library(internal_runner OBJECT "${{CMAKE_CURRENT_BINARY_DIR}}/internal-runner.cpp")
target_link_libraries(internal_runner PRIVATE plumed_thread_mpi)
add_executable(symbol_isolation
  $<TARGET_OBJECTS:host_runner> $<TARGET_OBJECTS:internal_runner>
  "${{CMAKE_CURRENT_BINARY_DIR}}/isolation-main.cpp")
target_link_libraries(symbol_isolation PRIVATE host_thread_mpi plumed_thread_mpi)
'''
        (work / "CMakeLists.txt").write_text(cmake)
        build = work / "build"
        subprocess.run(["cmake", "-S", str(work), "-B", str(build)], check=True)
        subprocess.run(["cmake", "--build", str(build), "-j4"], check=True)
        subprocess.run([str(build / "symbol_isolation")], check=True, timeout=20)
        print(f"host/internal tMPI isolation passed against {PINNED}")


if __name__ == "__main__":
    main()
