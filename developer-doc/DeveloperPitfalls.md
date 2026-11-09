\page DeveloperPitfalls Pitfalls and checks for developers

This page collects problems that were found while revising the code for PLUMED 2.11, and the checks that found them.
Read it before changing an action or a tool that other actions use.

\section DeveloperPitfallsRegtests Regression tests

- **A test for a fix must fail without the fix.** Check it by reverting the change of the source file only
  (e.g. `git stash push src/module/File.cpp`), rebuilding and running the test, then restoring the fix.
- **Wrong inputs** are tested in `plumed_regtest_after` of the `config` file: each bad input is written to a file, run with the driver,
  and the exit code and the first `ERROR in input` line are appended to a file `errors`, which is compared with `errors.reference`
  (see e.g. `regtest/multicolvar/rt-multicolvar-errors` or `regtest/basic/rt-molinfo-errors`).
  Check the text of the message, not whether it comes from `error()` or from an assertion.
  Set `export PLUMED_STACK_TRACE=no` in such a `config`: the CI runs with the stack trace on, and for an error raised
  outside an action the stack dump is printed where the test expects the message.
- **Do not print full precision** in the reference files (use `FMT=%8.4f` or similar): the last digits change with the compiler,
  the optimization flags and the number of threads. Results that depend on `std::shuffle`, on the distributions of `<random>`
  or on cancellations that change with fused multiply-add are not portable either.
- **Forces.** Check new or changed forces with `driver --debug-forces`, and on inputs that are not trivial: sparse matrices,
  matrices that are not square, rows with no elements, two OpenMP threads, `AVOID_THREAD_GATHER`.
  For models computed in single precision (e.g. pytorch), compare with an equivalent input in double precision instead.
- **Actions that collect data** (`COLLECT`, `COLLECT_FRAMES`, `ACCUMULATE` and the actions built on them) must also be run
  on a real trajectory through to the end (the jobs at the end of the calculation run there) and with an output at step 0.
- **Memory errors.** `make valgrind` inside a regtest directory runs that test under valgrind (the summary is in `tmp/err`).
  Tests with `mpiprocs>1` are skipped without an MPI build.
- **Unit tests that compile against the headers** (`type=make`, e.g. `regtest/basic/rt65-rmsd2` and `regtest/tools/rt-make-grid`)
  use functions of the tools that no action calls. Run them when changing the interface of a tool.

\section DeveloperPitfallsDriver The driver

- `plumed driver --noatoms` without a trajectory runs forever: give a small trajectory with `--ixyz`.
  Inside a regtest use `--noatoms` only with an input that reads its data with READ.
- A `--trajectory-stride` larger than the number of frames is an error.
- Test a new trajectory reader with valid and corrupted files under valgrind before trusting its regtest
  (see regtest/basic/rt-trajectory-malformed).

\section DeveloperPitfallsInput Input errors and messages

- A wrong input must stop with `error()`, naming the action and the keyword, and never at a `plumed_assert`
  or with an exception from the standard library (e.g. `map::at`): the user cannot tell what to change.
- Existing messages are part of the interface (regtests, tutorials and scripts look for them): do not reword them without a reason.
- When a keyword is renamed, register the old name with `Keywords::addDeprecatedKeyword`. When a keyword is removed,
  keep it as a hidden keyword that stops with a message giving the replacement (see `CUTOFF` in RDF).
- After changing an example in the manual, check every keyword it uses with `plumed manual --action NAME`.
- `Action::keywords` refers to the keywords of the options the action was created with: use it only in the constructor
  (e.g. not in `checkRead()` or later), where those options may not exist any more. The same holds for the functions
  that use it, e.g. `getkBT()`, which stops with an assertion when it is called after the constructor: store its result.
- `.ci/audit-keywords.py` adds each registered keyword to an example of the manual that uses the action, and reports
  the keywords that are refused as not understood, and the crashes; `--static-only` lists the keywords read with `parse`
  that no action of the file registers. `.ci/audit-inputs.py` changes the values of the keywords of the examples.

\section DeveloperPitfallsCompat Back compatibility

- An input that ran with the previous release must still run, or stop with a message that says what to write instead.
  The inputs of PLUMED-NEST and of the PLUMED tutorials are a good test: the build of the NEST site stores the output of
  every input for the last release and for master.
- Code that users compile with LOAD uses the headers of `src/tools` and `src/core`. When a function is removed or changed,
  add a line to the list of changes on the manual page of LOAD.

\section DeveloperPitfallsMisc Other limits

- Restarting with RESTART is exact only for the actions that store their whole state in their output files.
  CHECKPOINT gives an exact restart.
- The random number generator (`tools/Random.h`) has a period of about 2^31: the noise of `simplemd` or `pesmd`
  repeats after about 5.6e5 steps with 1000 atoms.
