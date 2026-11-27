# Helpers shared by the REPLEX_BIASVALUE regtests (copied with extra_files)

# make_traj NREP [NFRAMES]: writes traj.<i>.xyz for each replica i, two atoms
# whose distance oscillates differently on each replica, so that some exchanges
# are accepted and some rejected
make_traj() {
  local nrep=$1 nframes=${2:-40} i
  for ((i = 0; i < nrep; i++)); do
    awk -v r="$i" -v n="$nframes" 'BEGIN{
      for (f = 0; f < n; f++) {
        d = 0.25 + 0.12 * sin(1.7 * f + 2.1 * r) + 0.06 * cos(3.1 * f + 0.7 * r)
        print 2; print "10 10 10"; print "X 0 0 0"; printf "X %f 0 0\n", d
      }
    }' > "traj.$i.xyz"
  done
}

# split_traj NREP NFIRST: first.<i>.xyz with the first NFIRST frames of
# traj.<i>.xyz, and second.<i>.xyz with the others (2 atoms: 4 lines per frame)
split_traj() {
  local nrep=$1 nfirst=$2 i
  for ((i = 0; i < nrep; i++)); do
    head -n $((4 * nfirst)) "traj.$i.xyz" > "first.$i.xyz"
    tail -n +$((4 * nfirst + 1)) "traj.$i.xyz" > "second.$i.xyz"
  done
}

# compare DIR FILE...: whether each FILE is bitwise identical to DIR/FILE
compare() {
  local dir=$1 f
  shift
  for f in "$@"; do
    if cmp -s "$f" "$dir/$f"; then echo "$dir/$f IDENTICAL"; else echo "$dir/$f DIFFERENT"; fi
  done
}

# expect_error NAME TEXT DRIVER_ARGS...: runs the driver on NREP replicas in
# directory NAME (which must contain plumed.dat) and reports whether it stopped
# with an error containing TEXT. A timeout means that some replicas did not
# stop (deadlock), which is what these tests are meant to catch.
expect_error() {
  local name=$1 text=$2 code limit=
  shift 2
  if command -v timeout > /dev/null; then limit="timeout 120"; fi
  (
    cd "$name" || exit
    $limit $mpi $plumed driver --plumed plumed.dat --multi "$NREP" "$@" \
      --log-file log > out 2> err
    code=$?
    if ((code == 124)); then
      echo "$name: TIMEOUT (some replicas did not stop)"
    elif ((code != 0)) && cat out err log.* 2> /dev/null | grep -q -F "$text"; then
      echo "$name: rejected"
    else
      echo "$name: NOT REJECTED PROPERLY (exit code $code)"
    fi
  )
}

# extract_log LOGFILE: the lines written by REPLEX_BIASVALUE (setup, exchanges,
# statistics), without the timings of the rest of the log
extract_log() {
  sed 's/^PLUMED: //' "$1" | awk '
    /Cycles +Total +Average/ { exit }
    /^(Action )?REPLEX_BIASVALUE/ { keep = 1; print; next }
    keep && /^  / { print; next }
    { keep = 0 }'
}
